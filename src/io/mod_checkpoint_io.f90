! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Per-rank binary checkpoint serialization/deserialization of FVSolver state.
module mod_checkpoint_io
  use mod_kinds
  use iso_fortran_env, only: int64
  implicit none
  private
  public :: write_checkpoint, read_checkpoint

contains

  !> Write checkpoint file for a single rank
  !> File contains: header, time, global_id, fields (U,p,T,C)
  subroutine write_checkpoint(solver, msh, checkpoint_dir, step, rank)
    ! Import types locally to avoid module-level dependency issues
    use mod_mesh,    only: Mesh
    use mod_wall_model, only: wall_model_get_filter_state
    use mod_sem_inflow, only: sem_get_rng_state
    use mod_fv_solver, only: FVSolver
    use mpi_f08
    type(FVSolver),       intent(in) :: solver
    type(Mesh),           intent(in) :: msh
    character(len=*),     intent(in) :: checkpoint_dir
    integer(ip),          intent(in) :: step
    integer(ip),          intent(in) :: rank

    integer(ip) :: iunit, ios, nc, ns, i
    character(len=512) :: filename, dir_cmd
    character(len=16)  :: rank_str, step_str
    integer(ip) :: header(8)
    integer(ip) :: has_C
    real(rp), allocatable :: wm_filt(:,:)
    logical :: wm_ready, sem_ready
    logical :: dir_exists

    nc = msh%n_cells
    ns = solver%params%n_scalars

    ! Determine which optional fields to write
    has_C = 0;  if (ns > 0) has_C = 1

    ! Pack header
    header(1) = nc
    header(2) = ns
    ! Layout shared with dopamine-checkpoint-partition: header(3) is the history-record count
    ! (RK3W keeps none, so 0), header(4) the step. header(7) flags a trailing U_face record.
    header(3) = 0
    header(4) = solver%step
    header(5) = has_C
    header(6) = 0  ! reserved
    header(7) = 0
    if (solver%u_face_ready .and. allocated(solver%U_face)) header(7) = 1
    ! header(8) masks trailing wall-model records after U_face: 1 tau_w_bnd, 2 q_w_bnd, 4 L_bnd, 8 stress filter.
    header(8) = 0
    if (allocated(solver%tau_w_bnd)) header(8) = header(8) + 1
    if (allocated(solver%q_w_bnd))   header(8) = header(8) + 2
    if (allocated(solver%L_bnd))     header(8) = header(8) + 4
    call wall_model_get_filter_state(wm_filt, wm_ready)
    if (wm_ready) header(8) = header(8) + 8  ! 8 = wall-stress time-filter state
    sem_ready = solver%abl_config%use_sem .and. solver%sem_state%initialized
    if (sem_ready) sem_ready = allocated(solver%sem_state%pos) .and. allocated(solver%sem_state%amp)
    if (sem_ready) header(8) = header(8) + 16  ! 16 = synthetic-eddy inflow state

    ! Create processorR directory if needed (each rank creates its own to avoid race condition)
    write(rank_str, '(I0)') rank
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)
    inquire(file=trim(filename), exist=dir_exists)
    if (.not. dir_exists) then
      write(dir_cmd, '(A)') 'mkdir -p '//trim(filename)
      call execute_command_line(trim(dir_cmd), wait=.true., exitstat=ios)
      if (ios /= 0) error stop 'write_checkpoint: mkdir failed'
    end if

    ! Construct checkpoint filename: checkpoint_XXXXXXXX.bin (8-digit zero-padded step)
    write(step_str, '(I8.8)') step
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)//'/checkpoint_'//trim(step_str)//'.bin'

    ! Open file for writing
    open(newunit=iunit, file=trim(filename), form='unformatted', access='sequential', &
         action='write', status='replace', iostat=ios)
    if (ios /= 0) error stop 'write_checkpoint: cannot open '//trim(filename)

    ! Record 1: header
    write(iunit) header

    ! Record 2: time
    write(iunit) solver%time

    ! Record 3: global_id (owned cells only, 1..nc)
    if (nc > 0) then
      if (allocated(msh%global_id)) then
        write(iunit) msh%global_id(1:nc)
      else
        ! Fallback: write 1-based indices if global_id not allocated
        write(iunit) [(int(i, lp), i=1,nc)]
      end if
    end if

    ! Record 4-7: Primary fields (owned cells only)
    if (nc > 0) then
      write(iunit) solver%U(:, 1:nc)
      write(iunit) solver%p(1:nc)
      write(iunit) solver%T(1:nc)
      write(iunit) solver%nu_sgs(1:nc)
    end if

    ! Record 8: Passive scalars (if present)
    if (has_C == 1 .and. nc > 0) then
      write(iunit) solver%C(1:ns, 1:nc)
    end if

    ! Divergence-free face flux: advection reads it on the first stage after a restart.
    if (header(7) == 1) write(iunit) solver%U_face
    if (iand(header(8), 1_ip) /= 0) write(iunit) solver%tau_w_bnd
    if (iand(header(8), 2_ip) /= 0) write(iunit) solver%q_w_bnd
    if (iand(header(8), 4_ip) /= 0) write(iunit) solver%L_bnd
    if (iand(header(8), 8_ip) /= 0) write(iunit) wm_filt
    if (iand(header(8), 16_ip) /= 0) then
      write(iunit) solver%sem_state%n_eddies, sem_get_rng_state()
      write(iunit) solver%sem_state%pos
      write(iunit) solver%sem_state%amp
    end if

    close(iunit)

    ! All ranks must finish their files before the manifest advertises this step.
    call MPI_Barrier(msh%comm, ios)

    ! Rank 0 writes checkpoint_latest manifest
    if (rank == 0) then
      filename = trim(checkpoint_dir)//'/checkpoint_latest'
      open(newunit=iunit, file=trim(filename), form='formatted', &
           action='write', status='replace', iostat=ios)
      if (ios /= 0) error stop 'write_checkpoint: cannot open checkpoint_latest'
      write(iunit, '(I0)') step
      close(iunit)
    end if

  end subroutine write_checkpoint


  !> Read checkpoint file for a single rank
  !> Populates FVSolver fields from checkpoint_XXXXXXXX.bin
  subroutine read_checkpoint(solver, msh, checkpoint_dir, step, rank)
    ! Import types locally to avoid module-level dependency issues
    use mod_mesh,    only: Mesh
    use mod_fv_solver, only: FVSolver
    use mod_wall_model, only: wall_model_set_filter_state
    use mod_sem_inflow, only: sem_set_rng_state
    use mod_fv_solver, only: update_sem_inflow
    type(FVSolver),       intent(inout) :: solver
    type(Mesh),           intent(in) :: msh
    character(len=*),     intent(in) :: checkpoint_dir
    integer(ip),          intent(in) :: step
    integer(ip),          intent(in) :: rank

    real(rp), allocatable :: wm_filt(:,:)
    integer(ip) :: n_eddies_ckpt
    integer(int64) :: rng_ckpt
    real(rp), allocatable :: eddy_pos(:,:), eddy_amp(:,:)
    integer(ip) :: iunit, ios, nc, ns, nc_ckpt, ns_ckpt, step_ckpt
    character(len=512) :: filename
    character(len=16)  :: rank_str, step_str
    integer(ip) :: header(8)
    integer(lp), allocatable :: global_id_ckpt(:)

    nc = msh%n_cells
    ns = solver%params%n_scalars

    ! Construct checkpoint filename
    write(rank_str, '(I0)') rank
    write(step_str, '(I8.8)') step
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)//'/checkpoint_'//trim(step_str)//'.bin'

    ! Open file for reading
    open(newunit=iunit, file=trim(filename), form='unformatted', access='sequential', &
         action='read', status='old', iostat=ios)
    if (ios /= 0) error stop 'read_checkpoint: cannot open '//trim(filename)

    ! Record 1: header
    read(iunit) header
    nc_ckpt = header(1)
    ns_ckpt = header(2)
    step_ckpt = header(4)

    ! Validate dimensions
    if (nc_ckpt /= nc) then
      write(*, '(A,I0,A,I0)') 'ERROR: checkpoint n_cells=', nc_ckpt, &
                               ' does not match mesh n_cells=', nc
      error stop 'read_checkpoint: dimension mismatch (n_cells)'
    end if
    if (ns_ckpt /= ns) then
      write(*, '(A,I0,A,I0)') 'ERROR: checkpoint n_scalars=', ns_ckpt, &
                               ' does not match solver n_scalars=', ns
      error stop 'read_checkpoint: dimension mismatch (n_scalars)'
    end if

    ! Record 2: time
    read(iunit) solver%time

    ! Record 3: global_id (just for validation/tracking; not used in read)
    if (nc > 0) then
      allocate(global_id_ckpt(nc))
      read(iunit) global_id_ckpt
      if (allocated(msh%global_id)) then
        if (any(global_id_ckpt /= msh%global_id(1:nc))) &
          error stop 'read_checkpoint: global_id mismatch (different partition)'
      end if
      deallocate(global_id_ckpt)
    end if

    ! Record 4-7: Primary fields
    if (nc > 0) then
      read(iunit) solver%U(1:3, 1:nc)
      read(iunit) solver%p(1:nc)
      read(iunit) solver%T(1:nc)
      read(iunit) solver%nu_sgs(1:nc)
    end if

    ! Record 8: Passive scalars
    if (ns_ckpt > 0 .and. nc > 0) then
      read(iunit) solver%C(1:ns, 1:nc)
    end if

    if (header(7) == 1 .and. allocated(solver%U_face)) then
      read(iunit) solver%U_face
      solver%u_face_ready = .true.
      !$acc update device(solver%U_face) if_present
    else
      write(*, '(A)') 'WARNING: checkpoint has no face flux; it is rebuilt from the restored velocity'
    end if

    ! Wall-model state (absent in older checkpoints and in repartitioned ones: header(8) = 0 there).
    if (iand(header(8), 1_ip) /= 0) then
      if (allocated(solver%tau_w_bnd)) then
        read(iunit) solver%tau_w_bnd
      else
        read(iunit)
      end if
    end if
    if (iand(header(8), 2_ip) /= 0) then
      if (allocated(solver%q_w_bnd)) then
        read(iunit) solver%q_w_bnd
      else
        read(iunit)
      end if
    end if
    if (iand(header(8), 4_ip) /= 0) then
      if (allocated(solver%L_bnd)) then
        read(iunit) solver%L_bnd
      else
        read(iunit)
      end if
    end if
    if (iand(header(8), 8_ip) /= 0) then
      allocate(wm_filt(3, size(solver%tau_w_bnd, 2)))
      read(iunit) wm_filt
      call wall_model_set_filter_state(wm_filt)
    end if
    if (iand(header(8), 16_ip) /= 0) then
      read(iunit) n_eddies_ckpt, rng_ckpt
      allocate(eddy_pos(3, n_eddies_ckpt), eddy_amp(3, n_eddies_ckpt))
      read(iunit) eddy_pos
      read(iunit) eddy_amp
      if (solver%abl_config%use_sem .and. solver%sem_state%initialized .and. &
          solver%sem_state%n_eddies == n_eddies_ckpt) then
        solver%sem_state%pos = eddy_pos
        solver%sem_state%amp = eddy_amp
        call sem_set_rng_state(rng_ckpt)
        call update_sem_inflow(solver, solver%time)  ! inlet fluctuations for the restored eddies
      else
        write(*, '(A)') 'WARNING: checkpoint synthetic-eddy state does not match this run; eddies restart fresh'
      end if
    end if

    ! Restore step counter
    solver%step = step_ckpt

    close(iunit)

  end subroutine read_checkpoint

end module mod_checkpoint_io
