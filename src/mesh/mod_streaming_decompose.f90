! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Streaming multi-pass Morton-ordered mesh decomposition for large meshes.
module mod_streaming_decompose
  use mod_kinds
  use mod_mesh
  use mod_gmsh_reader
  use mod_foam_reader
  use mod_mesh_io
  use mod_bc,           only: BCConfig
  use mod_bc_setup,     only: build_bc_configs
  use mod_input_parser
  use mod_partitioner,  only: build_morton_partition, pack_global_mesh, sort_and_offset_partition, &
                              build_local_partition_mesh
  use mod_nonconforming, only: enforce_face_conformity
  use mpi_f08
  !$ use omp_lib
  implicit none
  private

  public :: decompose_streaming

contains

  !> Main entry point for streaming mesh decomposition.
  subroutine decompose_streaming(input_file, config, nranks_target, output_dir, chunk_size, max_memory_gb, &
                                 weight_by_volume, n_threads)
    character(len=*), intent(in) :: input_file, output_dir
    type(InputConfig), intent(in) :: config
    integer(ip), intent(in) :: nranks_target, chunk_size
    real(rp), intent(in), optional :: max_memory_gb
    logical, intent(in), optional :: weight_by_volume  ! default: equal cell count
    integer(ip), intent(in), optional :: n_threads  ! writer threads; default: all available
    ! input_file is kept for CLI compatibility; mesh path comes from config%mesh_file.
    ! chunk_size is used as the baseline face scan chunk for large meshes.

    type(Mesh) :: mesh_global, tmp_mesh
    type(BCConfig), allocatable :: bc_configs_global(:)
    integer(ip), allocatable :: part(:), cell_order(:), rank_offsets(:)
    integer(ip), allocatable :: face_owner(:), face_nei(:), face_patch(:)
    integer(ip), allocatable :: patch_start(:), patch_nfaces(:), patch_bcid(:)
    integer(ip), allocatable :: periodic_pair_global(:)
    character(len=128), allocatable :: patch_name(:)
    character(len=64),  allocatable :: patch_type(:)
    real(rp), allocatable :: cell_centroid(:,:), cell_volume(:)
    real(rp), allocatable :: face_area(:), face_centre(:,:), face_normal(:,:)
    real(rp), allocatable :: face_sf(:,:), face_gc(:)
    real(rp), allocatable :: face_dcf(:)
    real(rp), allocatable :: vertices(:,:)
    integer(ip), allocatable :: cell_offsets(:), cell_vertices(:), cell_types(:)
    integer(ip), allocatable :: face_offsets(:), face_vertices(:)
    integer(ip) :: nc, nf, np, nv, nintf, nbndf, topo_flag
    integer(ip) :: face_topo_flag
    integer(ip) :: r, cmd_stat, exit_status, i, face_chunk_size
    integer(ip) :: staged_geom_write_unit, staged_geom_read_unit
    integer(ip) :: cell_chunk_size, last_pct
    real(rp) :: mem_budget_gb
    logical :: use_staged_geom, vol_weight
    integer(ip) :: writer_threads
    ! Progress-bar temporaries hoisted out of a `block` (nvfortran rejects BLOCK inside a parallel directive); scoped via !$omp private() instead.
    integer(ip) :: pct_now
    character(len=32) :: bar_now
    character(len=16) :: rank_str
    character(len=512) :: cmd, outdir, meshbin, mesh_path, staged_geom_file
    integer(lp) :: t_wall_start, t_step_start, t_now, clock_rate
    integer(ip), allocatable :: rf_count(:), rf_ptr(:), rf_list(:)
    integer(ip) :: total_rf, rf_pos, r1, r2

    call system_clock(t_wall_start, clock_rate)
    t_step_start = t_wall_start

    mem_budget_gb = 0.0_rp
    if (present(max_memory_gb)) mem_budget_gb = max(0.0_rp, max_memory_gb)
    vol_weight = .false.
    if (present(weight_by_volume)) vol_weight = weight_by_volume
    writer_threads = default_writer_threads()
    if (present(n_threads)) then
      if (n_threads > 0) writer_threads = n_threads
    end if
    face_chunk_size = max(100000_ip, chunk_size)
    if (mem_budget_gb > 0.0_rp) then
      face_chunk_size = int(max(100000.0_rp, min(real(face_chunk_size, rp), mem_budget_gb * 5.0e5_rp)), ip)
    end if
    cell_chunk_size = max(100000_ip, face_chunk_size)
    use_staged_geom = (mem_budget_gb > 0.0_rp)
    staged_geom_write_unit = -1
    staged_geom_read_unit = -1

    mesh_path = trim(config%mesh_file)
    if (len_trim(mesh_path) == 0 .and. len_trim(input_file) > 0) mesh_path = trim(input_file)
    write(cmd, '(A)') 'mkdir -p ' // trim(output_dir)
    call execute_command_line(trim(cmd), wait=.true., exitstat=exit_status, cmdstat=cmd_stat)

    ! Initialize mesh
    mesh_global%comm   = MPI_COMM_WORLD
    mesh_global%rank   = 0
    mesh_global%nranks = 1  ! Single process for decomposition

    write(*, '(A)') ' '
    write(*, '(A)') ' [1/6]  Reading mesh...'
    select case (trim(config%mesh_type))
    case ('gmsh')
      call gmsh_read(trim(mesh_path), mesh_global)
    case ('openfoam', 'foam')
      call foam_read(trim(mesh_path), mesh_global)
    case default
      write(*, '(A)') 'ERROR: unknown mesh_type '//trim(config%mesh_type)
      error stop 'decompose_streaming: unknown mesh_type'
    end select
    write(*, '(A,I0,A,I0,A)') ' Read mesh: ', mesh_global%n_cells, &
                              ' cells, ', mesh_global%n_faces, ' faces'
    call check_index_range(mesh_global)
    call report_decompose_memory_estimate(mesh_global%n_cells, mesh_global%n_faces, mem_budget_gb, use_staged_geom)
    call system_clock(t_now)
    write(*, '(2X,A,F8.2,A,F8.2,A)') 'Step time: ', real(t_now - t_step_start, rp)/real(clock_rate, rp), &
        ' s  |  Total elapsed: ', real(t_now - t_wall_start, rp)/real(clock_rate, rp), ' s'
    t_step_start = t_now

    ! Pass 1: Build mesh connectivity and compute geometry
    write(*, '(A)') ' '
    write(*, '(A)') ' [2/6]  Building connectivity and geometry...'
    ! Repair non-conforming (hanging-node) T-junction boundary faces before
    ! build_connectivity, which constructs cell->face CSR from the face array.
    call enforce_face_conformity(mesh_global)
    call mesh_global%build_connectivity()

    call mesh_global%compute_geometry()
    call system_clock(t_now)
    write(*, '(2X,A,F8.2,A,F8.2,A)') 'Step time: ', real(t_now - t_step_start, rp)/real(clock_rate, rp), &
        ' s  |  Total elapsed: ', real(t_now - t_wall_start, rp)/real(clock_rate, rp), ' s'
    t_step_start = t_now

    ! Pass 2: Identify periodic boundary pairs
    write(*, '(A)') ' '
    write(*, '(A)') ' [3/6]  Resolving boundary conditions...'
    call build_bc_configs(config, mesh_global, bc_configs_global)
    deallocate(bc_configs_global)

    ! Extract mesh connectivity (topology-first: skip expensive geometry packing)
    nc  = mesh_global%n_cells
    nf  = mesh_global%n_faces
    np  = mesh_global%n_patches
    nv  = mesh_global%n_vertices
    nintf = mesh_global%n_int_faces
    nbndf = mesh_global%n_bnd_faces

    ! Allocate and extract topology only (no geometry - computed per-rank)
    allocate(cell_centroid(3, nc), cell_volume(nc))
    allocate(face_owner(nf), face_nei(nf), face_patch(nf))
    allocate(patch_start(np), patch_nfaces(np), patch_bcid(np))
    allocate(patch_name(np), patch_type(np))

    ! Extract cell data
    do i = 1, nc
      cell_centroid(:, i) = mesh_global%cells(i)%centroid
      cell_volume(i) = mesh_global%cells(i)%volume
    end do

    ! Extract face connectivity
    face_owner = [(mesh_global%faces(i)%owner, i=1,nf)]
    face_nei = [(mesh_global%faces(i)%neighbour, i=1,nf)]
    face_patch = [(mesh_global%faces(i)%patch_id, i=1,nf)]

    ! Extract patch info
    do i = 1, np
      patch_start(i) = mesh_global%patches(i)%start_face
      patch_nfaces(i) = mesh_global%patches(i)%n_faces
      patch_bcid(i) = mesh_global%patches(i)%bc_id
      patch_name(i) = mesh_global%patches(i)%name
      patch_type(i) = mesh_global%patches(i)%bc_type
    end do

    ! Extract all face geometry (all faces already computed by mesh_global%compute_geometry())
    allocate(face_area(nf))
    if (.not. use_staged_geom) then
      allocate(face_centre(3, nf), face_normal(3, nf), face_sf(3, nf), face_gc(nf))
    else
      staged_geom_file = trim(output_dir)//'/.dopamine_face_geom.bin'
      open(newunit=staged_geom_write_unit, file=trim(staged_geom_file), form='unformatted', access='stream', &
           action='write', status='replace', iostat=cmd_stat)
      if (cmd_stat /= 0) then
        write(*, '(A)') ' WARNING: Could not open staged geometry file; using in-memory geometry arrays instead.'
        use_staged_geom = .false.
        allocate(face_centre(3, nf), face_normal(3, nf), face_sf(3, nf), face_gc(nf))
      end if
    end if
    allocate(face_dcf(nf))
    do i = 1, nf
      face_area(i) = mesh_global%faces(i)%area
      ! Persist sf and gc (S7 / F1): without these the partition
      ! mesh would silently fall back to sf=area*normal and gc=1.0, which
      if (use_staged_geom) then
        write(staged_geom_write_unit) mesh_global%faces(i)%centre, &
                                      mesh_global%faces(i)%normal, &
                                      mesh_global%faces(i)%sf, &
                                      mesh_global%faces(i)%gc
      else
        face_centre(:, i) = mesh_global%faces(i)%centre
        face_normal(:, i) = mesh_global%faces(i)%normal
        face_sf(:, i)     = mesh_global%faces(i)%sf
        face_gc(i)        = mesh_global%faces(i)%gc
      end if
      face_dcf(i) = mesh_global%faces(i)%d_cf
    end do

    if (use_staged_geom) then
      close(staged_geom_write_unit)
      open(newunit=staged_geom_read_unit, file=trim(staged_geom_file), form='unformatted', access='stream', &
           action='read', status='old', iostat=cmd_stat)
      if (cmd_stat /= 0) then
        write(*, '(A)') 'ERROR: failed to reopen staged geometry for reading.'
        error stop 'decompose_streaming: staged geometry reopen failed'
      end if
    end if

    ! Get optional topology
    allocate(periodic_pair_global(nbndf))
    if (allocated(mesh_global%periodic_pair)) then
      periodic_pair_global = mesh_global%periodic_pair
    else
      periodic_pair_global = 0
    end if

    if (allocated(mesh_global%vertices)) then
      topo_flag = 1
      nv = size(mesh_global%vertices, 2)
      allocate(vertices(3, nv))
      vertices = mesh_global%vertices
      allocate(cell_offsets(nc+1))
      cell_offsets = mesh_global%cell_offsets
      allocate(cell_vertices(size(mesh_global%cell_vertices)))
      cell_vertices = mesh_global%cell_vertices
      allocate(cell_types(nc))
      cell_types = mesh_global%cell_types
    else
      topo_flag = 0
    end if

    if (allocated(mesh_global%face_offsets)) then
      face_topo_flag = 1
      allocate(face_offsets(nf+1))
      face_offsets = mesh_global%face_offsets
      allocate(face_vertices(size(mesh_global%face_vertices)))
      face_vertices = mesh_global%face_vertices
    else
      face_topo_flag = 0
    end if

    ! Release heavyweight mesh object now that extracted arrays are ready.
    ! This lowers peak memory during partition sorting for very large meshes.
    call mesh_global%destroy()
    call system_clock(t_now)
    write(*, '(2X,A,F8.2,A,F8.2,A)') 'Step time: ', real(t_now - t_step_start, rp)/real(clock_rate, rp), &
        ' s  |  Total elapsed: ', real(t_now - t_wall_start, rp)/real(clock_rate, rp), ' s'
    t_step_start = t_now

    ! Pass 3: Compute Morton partition with volume-weighted load balancing
    write(*, '(A)') ' '
    write(*, '(A)') ' [4/6]  Partitioning mesh (Morton SFC)...'
    allocate(part(nc))
    call build_morton_partition(nc, nranks_target, 0, cell_centroid, part, cell_volume, &
                                 weight_by_volume=vol_weight)

    ! Sort cells by partition and compute offsets
    allocate(cell_order(nc), rank_offsets(0:nranks_target))
    call sort_and_offset_partition(nc, nranks_target, part, cell_order, rank_offsets)
    call system_clock(t_now)
    write(*, '(2X,A,F8.2,A,F8.2,A)') 'Step time: ', real(t_now - t_step_start, rp)/real(clock_rate, rp), &
        ' s  |  Total elapsed: ', real(t_now - t_wall_start, rp)/real(clock_rate, rp), ' s'
    t_step_start = t_now

    ! Pre-create all output directories in a single shell command
    write(cmd, '(A,I0,A,A,A)') &
        'for r in $(seq 0 ', nranks_target - 1, '); do mkdir -p "', trim(output_dir), '/processor$r"; done'
    call execute_command_line(trim(cmd), wait=.true., exitstat=exit_status, cmdstat=cmd_stat)
    if (cmd_stat /= 0 .or. exit_status /= 0) then
      write(*, '(A)') 'WARNING: pre-mkdir for output directories may have failed.'
    end if

    ! Build per-rank interior face lists: one O(nintf) pass instead of O(N*nintf) scans
    allocate(rf_count(0:nranks_target-1), source=0_ip)
    do i = 1, nintf
      r1 = part(face_owner(i))
      r2 = part(face_nei(i))
      rf_count(r1) = rf_count(r1) + 1
      if (r2 /= r1) rf_count(r2) = rf_count(r2) + 1
    end do
    allocate(rf_ptr(0:nranks_target))
    rf_ptr(0) = 1
    do i = 0, nranks_target - 1
      rf_ptr(i+1) = rf_ptr(i) + rf_count(i)
    end do
    total_rf = rf_ptr(nranks_target) - 1
    allocate(rf_list(max(1_ip, total_rf)))
    rf_count = 0
    do i = 1, nintf
      r1 = part(face_owner(i))
      rf_pos = rf_ptr(r1) + rf_count(r1)
      rf_list(rf_pos) = i
      rf_count(r1) = rf_count(r1) + 1
      r2 = part(face_nei(i))
      if (r2 /= r1) then
        rf_pos = rf_ptr(r2) + rf_count(r2)
        rf_list(rf_pos) = i
        rf_count(r2) = rf_count(r2) + 1
      end if
    end do
    deallocate(rf_count)

    ! Pass 4: Build and write local meshes
    write(*, '(A)') ' '
    write(*, '(A)') ' [5/6]  Writing local meshes...'

    last_pct = -1
    write(*, '(2X,A,I0,A)') 'Writer threads: ', writer_threads, ' (override with --threads)'
    if (.not. use_staged_geom) then
      !$omp parallel do num_threads(writer_threads) schedule(dynamic) default(shared) &
      !$omp& private(tmp_mesh, rank_str, outdir, meshbin, pct_now, bar_now)
      do r = 0, nranks_target - 1
        if (rank_offsets(r) > rank_offsets(r + 1) - 1) then
          !$omp critical (decompose_io)
          write(*,'(A,I0,A)') '   Rank ', r, ' has no owned cells (skipped).'
          !$omp end critical (decompose_io)
          cycle
        end if

        tmp_mesh%comm   = MPI_COMM_WORLD
        tmp_mesh%rank   = r
        tmp_mesh%nranks = nranks_target

        call build_local_partition_mesh(tmp_mesh, MPI_COMM_WORLD, r, nranks_target, &
                        nc, nf, np, part, &
                        cell_centroid=cell_centroid, cell_volume=cell_volume, &
                        face_owner=face_owner, face_nei=face_nei, face_patch=face_patch, &
                        face_area=face_area, face_dcf=face_dcf, &
                        patch_start=patch_start, patch_nfaces=patch_nfaces, patch_bcid=patch_bcid, &
                        patch_name=patch_name, patch_type=patch_type, &
                        periodic_pair_global=periodic_pair_global, &
                        topo_flag=topo_flag, n_vertices_global=nv, vertices_global=vertices, &
                        cell_offsets_global=cell_offsets, cell_vertices_global=cell_vertices, &
                        cell_types_global=cell_types, &
                        face_topo_flag=face_topo_flag, face_offsets_global=face_offsets, &
                        face_vertices_global=face_vertices, &
                        cell_order=cell_order, cell_start=rank_offsets(r), &
                        cell_end=rank_offsets(r + 1) - 1, &
                        face_chunk_size=face_chunk_size, face_centre=face_centre, face_normal=face_normal, &
                        face_sf=face_sf, face_gc=face_gc, &
                        int_face_list=rf_list(rf_ptr(r):rf_ptr(r+1)-1), n_int_faces_global=nintf)

        write(rank_str, '(I0)') r
        outdir = trim(output_dir)//'/processor'//trim(rank_str)
        meshbin = trim(outdir)//'/mesh.bin'
        call write_local_mesh_file(tmp_mesh, trim(meshbin))
        call tmp_mesh%destroy()
        pct_now = int(100.0_rp * real(r + 1_ip, rp) / real(max(1_ip, nranks_target), rp), ip)
        if (pct_now > 100_ip) pct_now = 100_ip
        !$omp critical (decompose_io)
        if (pct_now >= last_pct + 5_ip .or. r + 1_ip == nranks_target) then
          bar_now = make_bar(pct_now, 20_ip)
          write(*, '(2X,A,2X,I0,A,I0,2X,A,I3,A)') trim(bar_now), r+1, '/', nranks_target, '(', pct_now, '%)'
          last_pct = pct_now
        end if
        !$omp end critical (decompose_io)
      end do
      !$omp end parallel do
      ! Early deallocation: face geometry no longer needed after partition mesh building
      if (allocated(face_centre)) deallocate(face_centre)
      if (allocated(face_normal)) deallocate(face_normal)
      if (allocated(face_sf))     deallocate(face_sf)
      if (allocated(face_gc))     deallocate(face_gc)
    else
      ! Staged geometry uses positional reads, so every thread needs its own unit on the file.
      !$omp parallel default(shared) num_threads(writer_threads) &
      !$omp& private(tmp_mesh, rank_str, outdir, meshbin, r, staged_geom_read_unit, cmd_stat, pct_now, bar_now)
      open(newunit=staged_geom_read_unit, file=trim(staged_geom_file), form='unformatted', &
           access='stream', action='read', status='old', iostat=cmd_stat)
      if (cmd_stat /= 0) error stop 'decompose: per-thread staged geometry open failed'
      !$omp do schedule(dynamic)
      do r = 0, nranks_target - 1
        if (rank_offsets(r) > rank_offsets(r + 1) - 1) then
          !$omp critical (decompose_io)
          write(*,'(A,I0,A)') '   Rank ', r, ' has no owned cells (skipped).'
          !$omp end critical (decompose_io)
          cycle
        end if

        tmp_mesh%comm   = MPI_COMM_WORLD
        tmp_mesh%rank   = r
        tmp_mesh%nranks = nranks_target

        call build_local_partition_mesh(tmp_mesh, MPI_COMM_WORLD, r, nranks_target, &
                        nc, nf, np, part, &
                        cell_centroid=cell_centroid, cell_volume=cell_volume, &
                        face_owner=face_owner, face_nei=face_nei, face_patch=face_patch, &
                        face_area=face_area, face_dcf=face_dcf, &
                        patch_start=patch_start, patch_nfaces=patch_nfaces, patch_bcid=patch_bcid, &
                        patch_name=patch_name, patch_type=patch_type, &
                        periodic_pair_global=periodic_pair_global, &
                        topo_flag=topo_flag, n_vertices_global=nv, vertices_global=vertices, &
                        cell_offsets_global=cell_offsets, cell_vertices_global=cell_vertices, &
                        cell_types_global=cell_types, &
                        face_topo_flag=face_topo_flag, face_offsets_global=face_offsets, &
                        face_vertices_global=face_vertices, &
                        cell_order=cell_order, cell_start=rank_offsets(r), &
                        cell_end=rank_offsets(r + 1) - 1, &
                        face_chunk_size=face_chunk_size, staged_geom_unit=staged_geom_read_unit, &
                        int_face_list=rf_list(rf_ptr(r):rf_ptr(r+1)-1), n_int_faces_global=nintf)

        write(rank_str, '(I0)') r
        outdir = trim(output_dir)//'/processor'//trim(rank_str)
        meshbin = trim(outdir)//'/mesh.bin'
        call write_local_mesh_file(tmp_mesh, trim(meshbin))
        call tmp_mesh%destroy()
        pct_now = int(100.0_rp * real(r + 1_ip, rp) / real(max(1_ip, nranks_target), rp), ip)
        if (pct_now > 100_ip) pct_now = 100_ip
        !$omp critical (decompose_io)
        if (pct_now >= last_pct + 5_ip .or. r + 1_ip == nranks_target) then
          bar_now = make_bar(pct_now, 20_ip)
          write(*, '(2X,A,2X,I0,A,I0,2X,A,I3,A)') trim(bar_now), r+1, '/', nranks_target, '(', pct_now, '%)'
          last_pct = pct_now
        end if
        !$omp end critical (decompose_io)
      end do
      !$omp end do
      close(staged_geom_read_unit)
      !$omp end parallel
    end if
    write(*, '(2X,A,I0,A,A,A)') 'Wrote ', nranks_target, ' meshes to ', trim(output_dir), '/'
    call system_clock(t_now)
    write(*, '(2X,A,F8.2,A,F8.2,A)') 'Step time: ', real(t_now - t_step_start, rp)/real(clock_rate, rp), &
        ' s  |  Total elapsed: ', real(t_now - t_wall_start, rp)/real(clock_rate, rp), ' s'

    ! Cleanup (face geometry already deallocated in partition mesh building loop)
    deallocate(rf_ptr, rf_list)
    deallocate(part, cell_order, rank_offsets)
    deallocate(face_owner, face_nei, face_patch)
    deallocate(face_area, face_dcf)
    deallocate(patch_start, patch_nfaces, patch_bcid, patch_name, patch_type)
    deallocate(periodic_pair_global)
    if (allocated(cell_centroid)) deallocate(cell_centroid, cell_volume)
    if (allocated(vertices)) deallocate(vertices)
    if (allocated(cell_offsets)) deallocate(cell_offsets, cell_vertices, cell_types)
    if (allocated(face_offsets)) deallocate(face_offsets, face_vertices)

    if (use_staged_geom .and. staged_geom_read_unit > 0) then
      close(staged_geom_read_unit)
      write(cmd, '(A)') 'rm -f ' // trim(staged_geom_file)
      call execute_command_line(trim(cmd), wait=.true., exitstat=exit_status, cmdstat=cmd_stat)
    end if

  end subroutine decompose_streaming

  !> Fail loudly before ip overflows; connectivity CSR is the first array to wrap.
  subroutine check_index_range(msh)
    type(Mesh), intent(in) :: msh

    integer(lp) :: n_conn, limit

    limit = int(huge(1_ip), lp)
    n_conn = 0_lp
    if (allocated(msh%face_vertices)) n_conn = max(n_conn, int(size(msh%face_vertices), lp))
    if (allocated(msh%cell_vertices)) n_conn = max(n_conn, int(size(msh%cell_vertices), lp))
    n_conn = max(n_conn, 4_lp * int(msh%n_faces, lp))

    if (n_conn > limit) then
      write(*,'(A,I0,A,I0)') ' ERROR: connectivity size ', n_conn, ' exceeds integer limit ', limit
      error stop 'mesh too large for ip; rebuild with 64-bit ip'
    else if (real(n_conn, rp) > 0.8_rp * real(limit, rp)) then
      write(*,'(A,F5.1,A)') ' WARNING: connectivity index at ', &
        100.0_rp * real(n_conn, rp) / real(limit, rp), '% of the integer limit'
    end if
  end subroutine check_index_range

  !> Writer threads available; 1 when built without OpenMP.
  integer(ip) function default_writer_threads() result(n)
    n = 1_ip
    !$ n = int(omp_get_max_threads(), ip)
    n = max(1_ip, n)
  end function default_writer_threads

  subroutine report_decompose_memory_estimate(nc, nf, mem_budget_gb, use_staged_geom)
    integer(ip), intent(in) :: nc, nf
    real(rp), intent(in) :: mem_budget_gb
    logical, intent(in) :: use_staged_geom

    real(rp) :: bytes_cell, bytes_face_geom, bytes_partition, total_gb

    bytes_cell = real(nc, rp) * (4.0_rp * 8.0_rp)
    if (use_staged_geom) then
      bytes_face_geom = real(nf, rp) * (2.0_rp * 8.0_rp)
    else
      bytes_face_geom = real(nf, rp) * (8.0_rp * 8.0_rp)
    end if
    bytes_partition = real(nc, rp) * (4.0_rp + 4.0_rp + 8.0_rp)
    total_gb = (bytes_cell + bytes_face_geom + bytes_partition) / (1024.0_rp**3)

    write(*,'(2X,A,F10.2,A)') 'Working-set estimate: ', total_gb, ' GB'
    if (mem_budget_gb > 0.0_rp .and. total_gb > mem_budget_gb) then
      write(*,'(A,F8.2,A)') ' WARNING: estimate exceeds requested memory budget of ', mem_budget_gb, ' GB'
    end if
    if (use_staged_geom) then
      write(*,'(2X,A)') 'Geometry staging: face data streamed from disk'
    end if
  end subroutine report_decompose_memory_estimate

  function make_bar(pct, width) result(bar)
    integer(ip), intent(in) :: pct, width
    character(len=32) :: bar
    integer(ip) :: i, nfill, use_width

    bar = '[..............................]'
    use_width = max(1_ip, min(30_ip, width))
    nfill = int(real(use_width, rp) * real(max(0_ip, min(100_ip, pct)), rp) / 100.0_rp, ip)

    bar(1:1) = '['
    do i = 1, use_width
      if (i <= nfill) then
        bar(i+1:i+1) = '|'
      else
        bar(i+1:i+1) = '.'
      end if
    end do
    bar(use_width+2:use_width+2) = ']'
    bar(use_width+3:) = ' '
  end function make_bar

end module mod_streaming_decompose
