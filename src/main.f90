! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> dopamine: fractional-step FV LES solver. Parses args, reads mesh, runs time loop.
program dopamine_main
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gmsh_reader
  use mod_foam_reader
  use mod_partitioner
  use mod_mesh_io
  use mod_mesh_validate, only: validate_processor_faces, validate_periodic_ghost_centroids, report_mesh_quality, &
                              report_mesh_integrity
  use mod_nonconforming, only: enforce_face_conformity
  use mod_bc, only: BCConfig
  use mod_bc_setup, only: build_bc_configs
  use mod_fv_solver
  use mod_input_parser
  use mod_runtime_env
  use mod_warn_log
  use mod_output_writer
  use mod_checkpoint_io
  use mod_stats
  use mod_sampling
  use mod_initconditions
  use mod_device, only: device_bind, device_mpi_gpu_aware
  use mpi_f08
  implicit none

  ! -- MPI
  integer(ip)    :: rank, nranks, ierr

  ! -- Objects
  type(Mesh)           :: mesh_global, mesh_local
  type(FVSolver)       :: solver
  type(InputConfig)    :: config
  type(OutputWriter)   :: writer
  type(StatsCollector) :: stats
  type(SamplingCollector) :: sampler
  type(BCConfig), allocatable :: bc_configs(:)
  type(BCConfig), allocatable :: bc_configs_global(:)

  ! -- Book-keeping
  character(256) :: input_file, mesh_file_local
  integer(ip)    :: step
  real(rp)       :: dt
  logical        :: do_write
  logical        :: do_stats
  logical        :: do_checkpoint
  logical        :: write_time_hit
  logical        :: use_write_time
  integer(ip)    :: restart_step
  character(256) :: checkpoint_latest_file
  integer(ip)    :: iunit_tmp, ios_tmp
  character(16)  :: rank_str

  ! Initialise MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

  ! Parse command line
  if (command_argument_count() < 1) then
    if (rank == 0) write(*, '(A)') 'Usage: dopamine <input_file>'
    call MPI_Finalize(ierr)
    stop
  end if
  call get_command_argument(1, input_file)

  ! Parse input file
  call parse_input_file(trim(input_file), config)
  use_write_time = (config%params%write_interval < 0)
  if (use_write_time .and. config%params%write_time <= 0.0_rp) &
    error stop 'write_interval < 0 requires write_time > 0'
  call configure_hdf5_file_locking(config%hdf5_use_file_locking, rank)
  call warn_log_init(MPI_COMM_WORLD, rank, nranks, config%warn_log_file, &
                     config%warn_flush_interval, 'input='//trim(input_file)//'  mesh='//trim(config%mesh_file))

  if (rank == 0) then
    write(*, '(A)') '================================================='
    write(*, '(A)') ' DOPAMINE -- Finite Volume LES solver'
    write(*, '(A)') '================================================='
    write(*, '(A,A)')  ' Mesh file    : ', trim(config%mesh_file)
    write(*, '(A,A)')  ' Mesh type    : ', trim(config%mesh_type)
    write(*, '(A,ES10.3)') ' nu           :', config%params%nu
    write(*, '(A,ES10.3)') ' dt           :', config%params%dt
    write(*, '(A,ES10.3)') ' t_end        :', config%params%t_end
    write(*, '(A,A)')  ' Time scheme  : ', 'RK3W'
    write(*, '(A,I6)')  ' n_scalars    :', config%params%n_scalars
    write(*, '(A,A)')  ' Initial cond.: ', trim(config%ic_type)
    write(*, '(A,I6)')  ' n_ranks      :', nranks
    write(*, '(A)') '-------------------------------------------------'
  end if

  ! One rank per GPU; must precede every !$acc region and PETSc initialisation.
  call device_bind(MPI_COMM_WORLD)

  ! Read mesh (skip for preprocessed type; will be read per-rank later)
  if (trim(config%mesh_type) /= 'preprocessed') then
    mesh_global%comm   = MPI_COMM_WORLD
    mesh_global%rank   = rank
    mesh_global%nranks = nranks

    if (rank == 0) then
      write(*, '(A)') ' Reading mesh on rank 0...'
      select case (trim(config%mesh_type))
      case ('gmsh')
        call gmsh_read(trim(config%mesh_file), mesh_global)
      case ('openfoam', 'foam')
        call foam_read(trim(config%mesh_file), mesh_global)
      case default
        write(*, '(A)') 'ERROR: unknown mesh_type '//trim(config%mesh_type)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end select
      write(*, '(A,I0,A,I0,A)') ' Read mesh: ', mesh_global%n_cells, &
                                  ' cells, ', mesh_global%n_faces, ' faces'
    end if

    ! Build connectivity and geometry (rank 0 only; then partition)
    if (rank == 0) then
      write(*, '(A)') ' Building global mesh connectivity and geometry...'
      ! Detect and repair non-conforming (hanging-node) T-junction faces before
      ! build_connectivity constructs the cell->face CSR from the face array.
      call enforce_face_conformity(mesh_global)
      call mesh_global%build_connectivity()
      call mesh_global%compute_geometry()

      ! Build periodic face pairing on the global mesh before partitioning.
      ! Local (post-partition) geometry cannot always infer cross-rank periodic pairs.
      ! Only rank 0 gets here, so the global mesh must not run its collectives on
      ! MPI_COMM_WORLD (the other ranks are already waiting in partition_mesh).
      if (nranks > 1) then
        mesh_global%comm   = MPI_COMM_SELF
        mesh_global%nranks = 1
        call build_bc_configs(config, mesh_global, bc_configs_global)
        deallocate(bc_configs_global)
        mesh_global%comm   = MPI_COMM_WORLD
        mesh_global%nranks = nranks
      end if
    end if
  end if

  ! Partition and distribute mesh (or read preprocessed mesh)
  if (trim(config%mesh_type) == 'preprocessed') then
    ! Each rank reads its own preprocessed mesh file
    write(rank_str, '(I0)') rank
    mesh_file_local = trim(config%mesh_file)//'/processor'//trim(rank_str)//'/mesh.bin'
    if (rank == 0) write(*, '(A)') ' Reading preprocessed mesh (parallel)...'
    call read_local_mesh_file(mesh_local, trim(mesh_file_local))
    mesh_local%comm   = MPI_COMM_WORLD
    mesh_local%rank   = rank
    mesh_local%nranks = nranks
    ! Decomposer already computed geometry; build connectivity for cell_faces and halo operators.
    call mesh_local%build_connectivity()
    ! skew/is_refine are derived, not serialised, so they must be rebuilt after a preprocessed read.
    call mesh_local%tag_refinement_faces()
  else if (nranks > 1) then
    call partition_mesh(mesh_global, mesh_local, MPI_COMM_WORLD)
    ! Build local connectivity and geometry
    call mesh_local%build_connectivity()
    call mesh_local%compute_geometry()
    ! mesh_global is destroyed inside partition_mesh_morton on rank 0
  else
    if (rank == 0) write(*, '(A)') ' Single-rank: no partitioning needed.'
    mesh_local = mesh_global
    mesh_local%n_cells_global = int(mesh_global%n_cells, lp)  ! For single-rank: global = local
    ! For single-rank: destroy the original reference to free memory
    ! (mesh_local and mesh_global are now separate copies)
    call mesh_global%destroy()
  end if

  if (rank == 0) &
    write(*, '(A,I0,A,I0,A)') ' Local mesh (rank 0): ', mesh_local%n_cells, &
                                ' cells, ', mesh_local%n_ghost, ' ghost cells'

  ! Map input BC configs to mesh patches
  if (rank == 0) write(*, '(A)') ' Building BC configs...'
  call build_bc_configs(config, mesh_local, bc_configs)
  if (rank == 0) write(*, '(A)') ' BC configs complete.'

  ! Build LS stencils after periodic pairing is known.
  if (rank == 0) write(*, '(A)') ' Building LS stencils...'
  call mesh_local%build_ls_stencil()
  if (rank == 0) write(*, '(A)') ' LS stencils complete...'

  ! Report cells whose LS gradient stencil needed regularisation (sliver or near-coplanar cells).
  block
    integer(ip)     :: n_deg_local, n_deg_global, ierr_deg
    character(64)   :: deg_text
    n_deg_local = mesh_ls_degenerate_count()
    call MPI_Allreduce(n_deg_local, n_deg_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr_deg)
    if (n_deg_global > 0) then
      write(deg_text, '(I0,A)') n_deg_global, ' cells: LS stencil regularised (sliver/coplanar)'
      call warn_log_message('[mesh]      ', trim(deg_text))
    end if
  end block

  ! Coloring only pays off when the face loops are threaded; serial loops lose cache locality to it.
  if (THREADED_FACE_LOOPS) call mesh_local%build_face_coloring()

  ! Every colored !$acc/!$omp face loop falls back to silent single-thread host execution when
  ! use_color is false; catch that here instead of it quietly showing up as a huge slowdown.
  if (mesh_local%n_faces > 0 .and. mesh_local%n_colors <= 0) then
    call warn_log_message('[mesh]      ', 'face coloring unavailable: colored kernels ' // &
      '(gradients, advection, viscous) will run single-threaded on host, not on GPU/OpenMP')
  end if

  ! One rank losing coloring (e.g. a tiny/degenerate local partition) desyncs GPU execution across
  ! the same collective even though each rank's own warning above looks harmless in isolation.
  block
    integer(ip) :: n_uncolored_local, n_uncolored_global, ierr_col
    character(64) :: col_text
    n_uncolored_local = merge(1_ip, 0_ip, mesh_local%n_faces > 0 .and. mesh_local%n_colors <= 0)
    call MPI_Allreduce(n_uncolored_local, n_uncolored_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr_col)
    if (n_uncolored_global > 0 .and. n_uncolored_global < nranks) then
      write(col_text, '(I0,A,I0,A)') n_uncolored_global, ' of ', nranks, ' ranks lack face coloring'
      call warn_log_message('[mesh]      ', trim(col_text) // &
        ': those ranks run colored kernels on host while others use GPU/OpenMP')
    end if
  end block

  ! SoA mirror of face/cell geometry: needed by !$acc kernels that can't map the AoS mesh type.
  call mesh_local%build_soa_mirror()

  ! Static mesh/stencil geometry never changes; staging it once means every !$acc kernel's own copyin (checked "if not already present") finds it already resident and skips the transfer.
#ifdef _OPENACC
  !$acc enter data copyin(mesh_local%face_owner_soa, mesh_local%face_neighbour_soa, &
  !$acc      mesh_local%face_sf_soa, mesh_local%face_gc_soa, mesh_local%face_eps_soa, &
  !$acc      mesh_local%cell_volume_soa, mesh_local%ls_stencil_offset, mesh_local%ls_stencil, &
  !$acc      mesh_local%ls_weights, mesh_local%ls_stencil_face, &
  !$acc      mesh_local%cell_centroid_soa, mesh_local%face_normal_soa, mesh_local%face_area_soa, &
  !$acc      mesh_local%face_centre_soa, mesh_local%face_dcf_soa, mesh_local%face_patch_id_soa, &
  !$acc      mesh_local%patch_bc_id_soa, mesh_local%periodic_neighbor_soa, &
  !$acc      mesh_local%periodic_pair_soa, mesh_local%periodic_pair_dcf_soa, mesh_local%periodic_gc_soa, &
  !$acc      mesh_local%face_skew_soa, mesh_local%cell_face_offset_soa, mesh_local%cell_faces)
  if (allocated(mesh_local%face_color_perm)) then
    !$acc enter data copyin(mesh_local%face_color_perm)
  end if
#endif

  ! Cross-rank processor-face geometry self-test
  call validate_processor_faces(mesh_local)
  call validate_periodic_ghost_centroids(mesh_local)
  block
    real(rp) :: nonortho_max_deg
    integer(lp) :: n_refine_faces
    logical :: corr_is_default
    corr_is_default = (config%params%n_nonortho_corr <= 0)  ! 0 = auto (parser default)
    call report_mesh_quality(mesh_local, nonortho_max_deg_out=nonortho_max_deg, n_refine_out=n_refine_faces)
    call report_mesh_integrity(mesh_local)
    ! Auto-upgrade n_nonortho_corr based on mesh non-orthogonality (only when the user left it at auto, 0).
    if (corr_is_default) then
      if (nonortho_max_deg > 60.0_rp) then
        config%params%n_nonortho_corr = 3
        if (rank == 0) write(*, '(A)') &
          ' [nonortho] auto-set n_nonortho_corr=3 (max non-ortho > 60 deg)'
      else if (nonortho_max_deg > 45.0_rp) then
        config%params%n_nonortho_corr = 2
        if (rank == 0) write(*, '(A)') &
          ' [nonortho] auto-set n_nonortho_corr=2 (max non-ortho > 45 deg)'
      end if
      ! A 2:1 face leaves ~47% of the Poisson operator explicit; two sweeps leave ~20% of the divergence.
      if (n_refine_faces > 0_lp .and. config%params%n_nonortho_corr < 4) then
        config%params%n_nonortho_corr = 4
        if (rank == 0) write(*, '(A)') &
          ' [nonortho] auto-set n_nonortho_corr=4 (refinement interfaces present)'
      end if
    end if
  end block

  ! Initialise solver
  if (rank == 0) write(*, '(A)') ' Initializing solver...'
  call solver%init(mesh_local, config%params, bc_configs, config%abl_config)
  if (rank == 0) write(*, '(A)') ' Solver initialization complete.'

  ! Handle restart or initial conditions
  restart_step = 0
  if (config%restart_enabled) then
    ! Restart from checkpoint
    restart_step = config%restart_step
    if (restart_step < 0) then
      ! Auto-read checkpoint_latest
      checkpoint_latest_file = trim(config%restart_dir)//'/checkpoint_latest'
      if (rank == 0) then
        open(newunit=iunit_tmp, file=trim(checkpoint_latest_file), form='formatted', &
             action='read', status='old', iostat=ios_tmp)
        if (ios_tmp /= 0) then
          write(*, '(A)') 'ERROR: cannot read checkpoint_latest from '//trim(checkpoint_latest_file)
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        read(iunit_tmp, *, iostat=ios_tmp) restart_step
        close(iunit_tmp)
        if (ios_tmp /= 0 .or. restart_step < 0) then
          write(*, '(A)') 'ERROR: invalid checkpoint_latest file'
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
      end if
      call MPI_Bcast(restart_step, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end if

    if (rank == 0) write(*, '(A,I0)') ' Reading checkpoint from step ', restart_step
    call read_checkpoint(solver, mesh_local, trim(config%restart_dir), restart_step, rank)
    if (rank == 0) write(*, '(A)') ' Checkpoint restored. Continuing from previous state...'
  else
    ! Normal initial conditions
    ! Temperature initial condition is uniform reference temperature.
    solver%T = config%params%T_ref

    if (rank == 0) then
      write(*, '(A,ES12.4)') ' IC: ic_uref = ', config%ic_uref
      write(*, '(A,ES12.4)') ' IC: ic_lref = ', config%ic_lref
      write(*, '(A,ES12.4)') ' IC: noise_level = ', config%ic_noise_level
    end if

    call apply_initial_conditions(solver, trim(config%ic_type), &
                                  config%ic_uref, config%ic_lref, config%ic_noise_level, &
                                  config%wall_normal_direction)
  end if

  ! A consistent projection is idempotent; a large change here is a discretisation defect.
  if (config%params%debug_diagnostics) call solver%projection_idempotency()

  ! Initialise output writer
  call output_writer_init(writer, trim(config%output_prefix), mesh_local, solver, &
                          restart=config%restart_enabled)
  if (config%stats_global_frequency > 0) then
    call stats_init(stats, trim(config%output_prefix), &
                    config%stats_global_frequency, rank, &
                    restart=config%restart_enabled)
  end if

  call sampling_init(sampler, trim(config%output_prefix), mesh_local, rank, &
                     config%n_planes, config%planes, &
                     config%n_lines, config%lines, &
                     config%n_curves, config%curves, &
                     config%n_surfaces, config%surfaces)

  ! Time loop
  call warn_log_set_step_start(restart_step + 1_ip)

  ! From here on the solver fields live on the device (OpenACC builds); the host copy is
  ! refreshed with sync_host only when output, statistics, sampling or checkpoints need it.
  call solver%to_device(gpu_aware_mpi=device_mpi_gpu_aware)

  if (rank == 0) write(*, '(A)') ' Starting time loop...'
  call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! Synchronize all ranks before time loop

  dt = config%params%dt

  do step = restart_step + 1, huge(1_ip) - 1_ip
    if (config%params%nstep_max > 0 .and. step > config%params%nstep_max) then
      if (rank == 0) write(*, '(A,I0,A)') &
        ' Reached nstep_max=', config%params%nstep_max, '; stopping time loop.'
      exit
    end if
    ! Natural dt for this step: adaptive (check CFL from step 1 to catch overshoots
    ! early) or fixed. Recomputed fresh every iteration -- write_time snapping below
    ! only ever shrinks dt, so a fixed dt must not be left shrunk on later steps.
    if (config%params%adapt_dt) then
      dt = solver%compute_dt()
    else
      dt = config%params%dt
    end if

    ! Snap dt down (never up, so this cannot relax the CFL/diffusive stability
    ! limit already applied above) to land exactly on the next write_time boundary.
    write_time_hit = .false.
    if (use_write_time) &
      call snap_dt_to_write_time(solver%time, config%params%write_time, dt, write_time_hit)

    ! Advance one step
    call solver%advance(dt)

    do_stats = .false.
    if (config%stats_global_frequency > 0) &
      do_stats = (mod(step, config%stats_global_frequency) == 0)
    do_checkpoint = .false.
    if (config%checkpoint_interval > 0) &
      do_checkpoint = (mod(step, config%checkpoint_interval) == 0)

    ! Output
    if (use_write_time) then
      do_write = write_time_hit
    else
      do_write = (mod(step, config%params%write_interval) == 0)
    end if
    if (do_write .or. do_stats .or. do_checkpoint .or. sampling_is_due(sampler, step)) &
      call solver%sync_host()
    if (do_write) then
      call output_writer_write(writer, solver, solver%time, step)
    end if

    if (do_stats) then
      call stats_sample(stats, solver, solver%time, step, rank, dt=dt, output_written=do_write)
    end if

    call sampling_maybe_sample(sampler, solver, solver%time, step)

    call warn_log_flush(step, solver%time)

    ! Checkpoint
    if (do_checkpoint) then
      call write_checkpoint(solver, mesh_local, trim(config%output_prefix), step, rank)
      if (rank == 0) write(*, '(A,I0)') ' Checkpoint written at step ', step
    end if

    if (solver%time >= config%params%t_end) exit
  end do

  if (rank == 0) write(*, '(A)') ' Time loop complete.....'

  call warn_log_finalize(step, solver%time)
  if (rank == 0 .and. warn_log_total() > 0) &
    write(*, '(A,I0,A,A)') ' Solver warnings: ', warn_log_total(), ' recorded in ', warn_log_filename()

  ! Finalise
  call output_writer_finalize(writer)
  if (config%stats_global_frequency > 0) then
    call stats_finalize(stats, rank)
  end if
  call sampling_finalize(sampler)
  call solver%destroy()
  call mesh_local%destroy()
  deallocate(bc_configs)

  call MPI_Finalize(ierr)
end program dopamine_main
