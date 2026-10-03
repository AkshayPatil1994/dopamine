! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> NetCDF output writer: appends one time record per call to "{prefix}.nc".
module mod_output_writer
  use mod_kinds
  use mod_mesh
  use mod_fv_solver, only: FVSolver
  use mpi_f08
  use netcdf
  implicit none
  private

  public :: OutputWriter
  public :: output_writer_init
  public :: output_writer_write
  public :: output_writer_finalize

  ! VTK element type codes written to cell_type variable.
  integer(ip), parameter :: VTK_VERTEX     =  1_ip
  integer(ip), parameter :: VTK_POLYGON    =  7_ip
  integer(ip), parameter :: VTK_TETRA      = 10_ip
  integer(ip), parameter :: VTK_HEXAHEDRON = 12_ip
  integer(ip), parameter :: VTK_WEDGE      = 13_ip
  integer(ip), parameter :: VTK_PYRAMID    = 14_ip
  integer(ip), parameter :: VTK_POLYHEDRON = 42_ip

  type, public :: OutputWriter
    character(256)      :: prefix    = 'dopamine'
    integer(ip)         :: snap_count = 0
    type(Mesh), pointer :: mesh      => null()
    integer(ip)         :: comm_rank = 0
    integer(ip)         :: comm_size = 1
    logical             :: has_unstructured_topology = .false.
    logical             :: boussinesq = .false.
    logical             :: use_sgs    = .false.  ! nu_sgs is only defined/written when true
    integer(ip), allocatable :: recvcounts(:)  ! < cells per rank
    integer(ip), allocatable :: recvdispls(:)  ! < cell offset per rank (0-based)
    integer(ip), allocatable :: recv_nverts(:)  ! < used vertices per rank
    integer(ip), allocatable :: recv_connsizes(:)  ! < connectivity entries per rank
    integer(ip), allocatable :: recv_toposizes(:)  ! < XDMF topology entries per rank
    real(rp), allocatable    :: times(:)
    integer(ip), allocatable :: steps(:)
    integer(ip)              :: snapshot_capacity = 0
    integer(ip)              :: n_verts_dedup = 0  ! < deduplicated unique vertex count after merge
    ! NetCDF IDs (-1 = not open / not applicable)
    integer(ip) :: ncid        = -1
    integer(ip) :: dim_time    = -1
    integer(ip) :: dim_cells   = -1
    integer(ip) :: dim_verts   = -1
    integer(ip) :: dim_conn    = -1
    integer(ip) :: dim_cp1     = -1  ! n_cells + 1
    integer(ip) :: dim_three   = -1
    integer(ip) :: dim_xdmf_topo = -1
    integer(ip) :: xdmf_topo_size = 0
    integer(ip) :: var_time    = -1
    integer(ip) :: var_step    = -1
    integer(ip) :: var_coords  = -1
    integer(ip) :: var_coords_xdmf = -1
    integer(ip) :: var_conn    = -1
    integer(ip) :: var_offsets = -1
    integer(ip) :: var_ctype   = -1
    integer(ip) :: var_xdmf_topology = -1
    integer(ip) :: var_p       = -1
    integer(ip) :: var_p_xdmf  = -1
    integer(ip) :: var_T       = -1
    integer(ip) :: var_T_xdmf  = -1
    integer(ip) :: var_U       = -1
    integer(ip) :: var_U_xdmf  = -1
    integer(ip) :: var_nu_sgs  = -1
    integer(ip) :: n_scalars = 0
    integer(ip), allocatable :: var_scalar(:)
    integer(ip), allocatable :: var_scalar_xdmf(:)
    character(len=64), allocatable :: scalar_var_names(:)
  end type OutputWriter

contains

  ! Public interface

  subroutine output_writer_init(writer, prefix, msh, solver, restart, restart_step)
    type(OutputWriter), intent(inout) :: writer
    character(len=*),   intent(in)    :: prefix
    type(Mesh), target, intent(in)    :: msh
    type(FVSolver),     intent(in), optional :: solver
    logical,            intent(in), optional :: restart
    integer(ip),        intent(in), optional :: restart_step  ! snapshots written after this step are discarded

    integer(ip) :: ierr, nc_local, nv_local, conn_local
    integer(ip) :: n_cells_total, n_verts_total, n_conn_total, n_topo_total
    integer(ip) :: r
    integer(ip), allocatable :: vgeom_counts(:), vgeom_displs(:)
    integer(ip), allocatable :: conn_counts(:),  conn_displs(:)
    integer(ip), allocatable :: topo_counts(:),  topo_displs(:)
    real(rp),    allocatable :: geom_local(:),   geom_global(:)
    integer(ip), allocatable :: conn_local_buf(:), conn_global(:)
    integer(ip), allocatable :: topo_local_buf(:), topo_global(:)
    integer(ip), allocatable :: vc_local(:), vc_global(:)
    integer(ip), allocatable :: ct_local(:), ct_global(:)
    integer(ip), allocatable :: offsets_global(:)
    integer(ip), allocatable :: vtx_map(:)
    real(rp),    allocatable :: geom_dedup(:)
    integer(ip) :: n_verts_dedup
    integer(ip), allocatable :: vtx_remap(:)
    character(len=512) :: fname
    logical :: do_restart, nc_exists

    writer%prefix     = trim(prefix)
    writer%mesh       => msh
    writer%snap_count = 0
    writer%n_scalars = 0
    do_restart = .false.
    if (present(restart)) do_restart = restart
    if (present(solver)) then
      writer%use_sgs = solver%params%use_sgs
      writer%n_scalars = max(0_ip, solver%params%n_scalars)
      allocate(writer%var_scalar(writer%n_scalars), writer%var_scalar_xdmf(writer%n_scalars))
      allocate(writer%scalar_var_names(writer%n_scalars))
      if (writer%n_scalars > 0) then
        writer%var_scalar = -1
        writer%var_scalar_xdmf = -1
        writer%scalar_var_names = ''
        do r = 1, writer%n_scalars
          ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
          ! must never be evaluated in the same expression as the allocated() check.
          writer%scalar_var_names(r) = sanitize_netcdf_name('scalar', r)
          if (allocated(solver%params%scalar_names)) then
            if (r <= size(solver%params%scalar_names)) &
              writer%scalar_var_names(r) = sanitize_netcdf_name(solver%params%scalar_names(r), r)
          end if
        end do
      end if
    end if

    call MPI_Comm_rank(msh%comm, writer%comm_rank, ierr)
    call MPI_Comm_size(msh%comm, writer%comm_size, ierr)

    ! ---- 1. Gather per-rank cell counts ----
    nc_local = msh%n_cells
    allocate(writer%recvcounts(writer%comm_size), writer%recvdispls(writer%comm_size))
    call MPI_Allgather(nc_local, 1, MPI_INTEGER, &
                       writer%recvcounts, 1, MPI_INTEGER, msh%comm, ierr)
    writer%recvdispls(1) = 0
    do r = 2, writer%comm_size
      writer%recvdispls(r) = writer%recvdispls(r-1) + writer%recvcounts(r-1)
    end do
    n_cells_total = sum(writer%recvcounts)

    allocate(writer%recv_nverts(writer%comm_size), writer%recv_connsizes(writer%comm_size), &
         writer%recv_toposizes(writer%comm_size))

    writer%has_unstructured_topology = &
        allocated(msh%vertices)     .and. allocated(msh%cell_offsets) .and. &
        allocated(msh%cell_vertices).and. allocated(msh%cell_types)

    ! ---- 2. Build local geometry / topology ----
    if (writer%has_unstructured_topology) then
      call build_local_geometry(msh, geom_local, nv_local, vtx_map)
      call build_local_verts_per_cell(msh, vc_local)
      call build_local_conn(msh, vtx_map, conn_local_buf)
      call build_local_cell_types(msh, ct_local)
      call build_local_xdmf_topology(msh, vtx_map, conn_local_buf, ct_local, topo_local_buf)
      conn_local = size(conn_local_buf)
      deallocate(vtx_map)
    else
      nv_local   = msh%n_cells
      conn_local = msh%n_cells
      allocate(geom_local(3 * nv_local))
      call pack_centroids(msh, geom_local)
      allocate(vc_local(msh%n_cells), source=1_ip)
      allocate(conn_local_buf(msh%n_cells))
      call iota(conn_local_buf, 0_ip)
      allocate(ct_local(msh%n_cells), source=VTK_VERTEX)
      allocate(topo_local_buf(3 * msh%n_cells))
      do r = 1, msh%n_cells
        topo_local_buf(3*(r-1)+1) = 1_ip
        topo_local_buf(3*(r-1)+2) = 1_ip
        topo_local_buf(3*(r-1)+3) = conn_local_buf(r)
      end do
    end if

    ! ---- 3. Exchange sizes ----
    call MPI_Allgather(nv_local,   1, MPI_INTEGER, writer%recv_nverts,    1, MPI_INTEGER, msh%comm, ierr)
    call MPI_Allgather(conn_local, 1, MPI_INTEGER, writer%recv_connsizes, 1, MPI_INTEGER, msh%comm, ierr)
    call MPI_Allgather(size(topo_local_buf), 1, MPI_INTEGER, writer%recv_toposizes, 1, MPI_INTEGER, msh%comm, ierr)
    n_verts_total = sum(writer%recv_nverts)
    n_conn_total  = sum(writer%recv_connsizes)
    n_topo_total  = sum(writer%recv_toposizes)

    ! ---- 4. Gather geometry to rank 0 (3 doubles per vertex) ----
    allocate(vgeom_counts(writer%comm_size), vgeom_displs(writer%comm_size))
    vgeom_counts = 3 * writer%recv_nverts
    vgeom_displs(1) = 0
    do r = 2, writer%comm_size
      vgeom_displs(r) = vgeom_displs(r-1) + vgeom_counts(r-1)
    end do
    if (writer%comm_rank == 0) then
      allocate(geom_global(3 * n_verts_total))
    else
      allocate(geom_global(1))
    end if
    call MPI_Gatherv(geom_local, 3*nv_local, MPI_DOUBLE_PRECISION, &
                     geom_global, vgeom_counts, vgeom_displs, MPI_DOUBLE_PRECISION, &
                     0, msh%comm, ierr)
    deallocate(geom_local, vgeom_counts, vgeom_displs)

    ! ---- 5. Gather connectivity ----
    allocate(conn_counts(writer%comm_size), conn_displs(writer%comm_size))
    conn_counts = writer%recv_connsizes
    conn_displs(1) = 0
    do r = 2, writer%comm_size
      conn_displs(r) = conn_displs(r-1) + conn_counts(r-1)
    end do
    if (writer%comm_rank == 0) then
      allocate(conn_global(n_conn_total))
    else
      allocate(conn_global(1))
    end if
    call MPI_Gatherv(conn_local_buf, conn_local, MPI_INTEGER, &
                     conn_global, conn_counts, conn_displs, MPI_INTEGER, &
                     0, msh%comm, ierr)
    deallocate(conn_local_buf, conn_counts, conn_displs)

    ! ---- 6. Gather XDMF mixed-topology stream ----
    allocate(topo_counts(writer%comm_size), topo_displs(writer%comm_size))
    topo_counts = writer%recv_toposizes
    topo_displs(1) = 0
    do r = 2, writer%comm_size
      topo_displs(r) = topo_displs(r-1) + topo_counts(r-1)
    end do
    if (writer%comm_rank == 0) then
      allocate(topo_global(n_topo_total))
    else
      allocate(topo_global(1))
    end if
    call MPI_Gatherv(topo_local_buf, size(topo_local_buf), MPI_INTEGER, &
                     topo_global, topo_counts, topo_displs, MPI_INTEGER, &
                     0, msh%comm, ierr)
    deallocate(topo_local_buf, topo_counts, topo_displs)

    ! ---- 7. Gather verts-per-cell and cell types ----
    if (writer%comm_rank == 0) then
      allocate(vc_global(n_cells_total), ct_global(n_cells_total))
    else
      allocate(vc_global(1), ct_global(1))
    end if
    call MPI_Gatherv(vc_local, nc_local, MPI_INTEGER, &
                     vc_global, writer%recvcounts, writer%recvdispls, MPI_INTEGER, &
                     0, msh%comm, ierr)
    call MPI_Gatherv(ct_local, nc_local, MPI_INTEGER, &
                     ct_global, writer%recvcounts, writer%recvdispls, MPI_INTEGER, &
                     0, msh%comm, ierr)
    deallocate(vc_local, ct_local)

    ! ---- 8. Rank 0: remap connectivity, build offsets, create NetCDF ----
    if (writer%comm_rank == 0) then
      ! First globalise per-rank-local vertex indices by adding cumulative rank offsets.
      ! After this step every vertex index in conn/topo is in [0, n_verts_total).
      call remap_connectivity_vertex_offsets(conn_global, &
                                             writer%recv_connsizes, writer%recv_nverts)
      call remap_xdmf_topology_vertex_offsets(topo_global, writer%recv_toposizes, writer%recv_nverts)

      ! Deduplicate geometry and remap the now-global vertex indices to deduped indices.
      call deduplicate_geometry(geom_global, n_verts_total, geom_dedup, &
                               n_verts_dedup, vtx_remap)
      writer%n_verts_dedup = n_verts_dedup
      call apply_vertex_dedup_to_connectivity(conn_global, vtx_remap)
      call apply_vertex_dedup_to_xdmf_topology(topo_global, vtx_remap)

      deallocate(geom_global)
      allocate(geom_global(3 * n_verts_dedup))
      geom_global = geom_dedup
      deallocate(geom_dedup, vtx_remap)
      n_verts_total = n_verts_dedup
      allocate(offsets_global(n_cells_total + 1))
      offsets_global(1) = 0
      do r = 1, n_cells_total
        offsets_global(r + 1) = offsets_global(r) + vc_global(r)
      end do

      fname = trim(prefix)//'.nc'
      if (do_restart) then
        ! A restart into a fresh output directory (e.g. after repartitioning) has no file to continue.
        inquire(file=trim(fname), exist=nc_exists)
        if (.not. nc_exists) then
          write(*, '(A)') ' Output: '//trim(fname)//' not found; starting a new output file for the restart.'
          do_restart = .false.
        end if
      end if
      if (do_restart) then
        if (present(restart_step)) then
          call reopen_netcdf_file(writer, trim(fname), restart_step)
        else
          call reopen_netcdf_file(writer, trim(fname))
        end if
      else
        call create_netcdf_file(writer, trim(fname), n_cells_total, n_verts_total, n_conn_total, &
                                geom_global, conn_global, offsets_global, ct_global, topo_global)
      end if
      deallocate(offsets_global)
    end if

    deallocate(geom_global, conn_global, topo_global, vc_global, ct_global)

  end subroutine output_writer_init


  subroutine output_writer_write(writer, solver, time, step)
    type(OutputWriter), intent(inout) :: writer
    type(FVSolver),     intent(in)    :: solver
    real(rp),           intent(in)    :: time
    integer(ip),        intent(in)    :: step

    integer(ip) :: ierr, nc_local, n_cells_total, r, c, tidx
    integer(ip), allocatable :: u_counts(:), u_displs(:)
    real(rp),    allocatable :: p_global(:), T_global(:), U_global(:), u_local(:)
    real(rp),    allocatable :: S_global(:,:), S_local(:), nu_sgs_global(:)
    integer(ip) :: s

    writer%boussinesq = solver%params%boussinesq
    call append_snapshot_metadata(writer, time, step)
    nc_local      = solver%mesh%n_cells
    n_cells_total = sum(writer%recvcounts)

    ! Pack U into flat contiguous buffer (Ux,Uy,Uz per cell).
    allocate(u_local(3 * nc_local))
    do c = 1, nc_local
      u_local(3*(c-1)+1:3*c) = solver%U(1:3, c)
    end do

    ! ---- Gather p and T (one scalar per cell) ----
    if (writer%comm_rank == 0) then
      allocate(p_global(n_cells_total))
    else
      allocate(p_global(1))
    end if
    call MPI_Gatherv(solver%p(1:nc_local), nc_local, MPI_DOUBLE_PRECISION, &
                     p_global, writer%recvcounts, writer%recvdispls, MPI_DOUBLE_PRECISION, &
                     0, solver%mesh%comm, ierr)
    if (writer%boussinesq) then
      if (writer%comm_rank == 0) then
        allocate(T_global(n_cells_total))
      else
        allocate(T_global(1))
      end if
      call MPI_Gatherv(solver%T(1:nc_local), nc_local, MPI_DOUBLE_PRECISION, &
                       T_global, writer%recvcounts, writer%recvdispls, MPI_DOUBLE_PRECISION, &
                       0, solver%mesh%comm, ierr)
    else
      allocate(T_global(1), source=0.0_rp)
    end if

    if (writer%use_sgs) then
      if (writer%comm_rank == 0) then
        allocate(nu_sgs_global(n_cells_total))
      else
        allocate(nu_sgs_global(1))
      end if
      call MPI_Gatherv(solver%nu_sgs(1:nc_local), nc_local, MPI_DOUBLE_PRECISION, &
                       nu_sgs_global, writer%recvcounts, writer%recvdispls, MPI_DOUBLE_PRECISION, &
                       0, solver%mesh%comm, ierr)
    else
      allocate(nu_sgs_global(1), source=0.0_rp)
    end if

    if (writer%n_scalars > 0) then
      if (writer%comm_rank == 0) then
        allocate(S_global(writer%n_scalars, n_cells_total))
      else
        allocate(S_global(1, 1))
      end if
      allocate(S_local(nc_local))
      do s = 1, writer%n_scalars
        S_local = solver%C(s, 1:nc_local)
        if (writer%comm_rank == 0) then
          call MPI_Gatherv(S_local, nc_local, MPI_DOUBLE_PRECISION, &
                           S_global(s, :), writer%recvcounts, writer%recvdispls, MPI_DOUBLE_PRECISION, &
                           0, solver%mesh%comm, ierr)
        else
          call MPI_Gatherv(S_local, nc_local, MPI_DOUBLE_PRECISION, &
                           S_global(1, :), writer%recvcounts, writer%recvdispls, MPI_DOUBLE_PRECISION, &
                           0, solver%mesh%comm, ierr)
        end if
      end do
      deallocate(S_local)
    end if

    ! ---- Gather U (3 values per cell) ----
    allocate(u_counts(writer%comm_size), u_displs(writer%comm_size))
    u_counts = 3 * writer%recvcounts
    u_displs(1) = 0
    do r = 2, writer%comm_size
      u_displs(r) = u_displs(r-1) + u_counts(r-1)
    end do
    if (writer%comm_rank == 0) then
      allocate(U_global(3 * n_cells_total))
    else
      allocate(U_global(1))
    end if
    call MPI_Gatherv(u_local, 3*nc_local, MPI_DOUBLE_PRECISION, &
                     U_global, u_counts, u_displs, MPI_DOUBLE_PRECISION, &
                     0, solver%mesh%comm, ierr)
    deallocate(u_local, u_counts, u_displs)

    ! ---- Rank 0: append time record ----
    tidx = writer%snap_count + 1
    writer%snap_count = tidx
    if (writer%comm_rank == 0) then
      call append_netcdf_record(writer, tidx, time, step, &
                                p_global, T_global, U_global, n_cells_total, S_global, &
                                nu_sgs_global)
      ! Force check for nc data file
      call nc_check(nf90_sync(writer%ncid), 'nf90_sync')
      call write_temporal_xmf(writer)
    end if
    deallocate(p_global, T_global, U_global, nu_sgs_global)
    if (allocated(S_global)) deallocate(S_global)
  end subroutine output_writer_write


  subroutine output_writer_finalize(writer)
    type(OutputWriter), intent(inout) :: writer

    if (writer%comm_rank == 0 .and. writer%ncid >= 0) then
      call nc_check(nf90_close(writer%ncid), 'nf90_close')
      writer%ncid = -1
    end if
    if (allocated(writer%recvcounts))     deallocate(writer%recvcounts)
    if (allocated(writer%recvdispls))     deallocate(writer%recvdispls)
    if (allocated(writer%recv_nverts))    deallocate(writer%recv_nverts)
    if (allocated(writer%recv_connsizes)) deallocate(writer%recv_connsizes)
    if (allocated(writer%recv_toposizes)) deallocate(writer%recv_toposizes)
    if (allocated(writer%times))          deallocate(writer%times)
    if (allocated(writer%steps))          deallocate(writer%steps)
    if (allocated(writer%var_scalar))     deallocate(writer%var_scalar)
    if (allocated(writer%var_scalar_xdmf)) deallocate(writer%var_scalar_xdmf)
    if (allocated(writer%scalar_var_names)) deallocate(writer%scalar_var_names)
    if (associated(writer%mesh)) writer%mesh => null()
  end subroutine output_writer_finalize

  ! Private: local topology builders

  !> Build a contiguous geometry buffer of used vertices and a 0-based vertex
  !> remapping table for the local partition.
  subroutine build_local_geometry(msh, geom_buf, n_used, vtx_map)
    type(Mesh),              intent(in)  :: msh
    real(rp), allocatable,   intent(out) :: geom_buf(:)  ! (3 * n_used)
    integer(ip),             intent(out) :: n_used
    integer(ip), allocatable,intent(out) :: vtx_map(:)  ! (n_vertices) -> 0-based output ID

    integer(ip) :: c, n0, n1, v, vid, new_id

    allocate(vtx_map(msh%n_vertices), source=-1_ip)
    new_id = -1
    do c = 1, msh%n_cells
      n0 = msh%cell_offsets(c)
      n1 = msh%cell_offsets(c+1) - 1
      do v = n0, n1
        vid = msh%cell_vertices(v)
        if (vtx_map(vid) < 0) then
          new_id = new_id + 1
          vtx_map(vid) = new_id
        end if
      end do
    end do

    n_used = new_id + 1
    allocate(geom_buf(3 * n_used))
    do v = 1, msh%n_vertices
      if (vtx_map(v) >= 0) then
        geom_buf(3*vtx_map(v)+1:3*vtx_map(v)+3) = msh%vertices(1:3, v)
      end if
    end do
  end subroutine build_local_geometry

  !> Build per-cell vertex count array (needed to reconstruct global CSR offsets).
  subroutine build_local_verts_per_cell(msh, vc)
    type(Mesh),              intent(in)  :: msh
    integer(ip), allocatable,intent(out) :: vc(:)
    integer(ip) :: c

    allocate(vc(msh%n_cells))
    do c = 1, msh%n_cells
      vc(c) = msh%cell_offsets(c+1) - msh%cell_offsets(c)
    end do
  end subroutine build_local_verts_per_cell

  !> Build the local 0-based connectivity array using the vtx_map.
  !> Input vtx_map maps 1-based vertex IDs to 0-based output IDs.
  subroutine build_local_conn(msh, vtx_map, conn)
    type(Mesh),              intent(in)  :: msh
    integer(ip),             intent(in)  :: vtx_map(:)
    integer(ip), allocatable,intent(out) :: conn(:)

    integer(ip) :: c, n0, n1, k, n_conn

    n_conn = msh%cell_offsets(msh%n_cells+1) - msh%cell_offsets(1)
    allocate(conn(n_conn))
    k = 0
    do c = 1, msh%n_cells
      n0 = msh%cell_offsets(c)
      n1 = msh%cell_offsets(c+1) - 1
      conn(k+1:k+n1-n0+1) = vtx_map(msh%cell_vertices(n0:n1))
      k = k + n1 - n0 + 1
    end do
  end subroutine build_local_conn

  !> Map vertex count per cell to VTK element type codes.
  subroutine build_local_cell_types(msh, ct)
    type(Mesh),              intent(in)  :: msh
    integer(ip), allocatable,intent(out) :: ct(:)
    integer(ip) :: c, nverts

    allocate(ct(msh%n_cells))
    ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
    ! must never be evaluated in the same expression as the allocated() check.
    if (allocated(msh%cell_types)) then
      if (size(msh%cell_types) >= msh%n_cells) then
        ct = msh%cell_types(1:msh%n_cells)
        return
      end if
    end if

    do c = 1, msh%n_cells
      nverts = msh%cell_offsets(c+1) - msh%cell_offsets(c)
      select case (nverts)
      case (4);      ct(c) = VTK_TETRA
      case (5);      ct(c) = VTK_PYRAMID
      case (6);      ct(c) = VTK_WEDGE
      case (8);      ct(c) = VTK_HEXAHEDRON
      case default;  ct(c) = VTK_POLYHEDRON  ! Arbitrary polyhedra with nverts vertices
      end select
    end do
  end subroutine build_local_cell_types

  !> Fallback: pack cell centroids as a flat (3*n_cells) geometry buffer.
  subroutine pack_centroids(msh, buf)
    type(Mesh), intent(in)  :: msh
    real(rp),   intent(out) :: buf(:)
    integer(ip) :: c
    do c = 1, msh%n_cells
      buf(3*(c-1)+1:3*c) = msh%cells(c)%centroid
    end do
  end subroutine pack_centroids

  !> Fill arr with start_val, start_val+1, ..., start_val+n-1.
  subroutine iota(arr, start_val)
    integer(ip), intent(out) :: arr(:)
    integer(ip), intent(in)  :: start_val
    integer(ip) :: i
    do i = 1, size(arr)
      arr(i) = start_val + i - 1
    end do
  end subroutine iota

  ! Private: connectivity remapping and NetCDF I/O

  !> Apply per-rank vertex offsets to the concatenated connectivity array so
  !> that every index becomes a global 0-based vertex ID.
  subroutine remap_connectivity_vertex_offsets(conn, recv_connsizes, recv_nverts)
    integer(ip), intent(inout) :: conn(:)
    integer(ip), intent(in)    :: recv_connsizes(:), recv_nverts(:)

    integer(ip) :: r, conn_off, vert_off, k

    conn_off = 0
    vert_off = 0
    do r = 1, size(recv_connsizes)
      if (vert_off > 0) then
        do k = conn_off + 1, conn_off + recv_connsizes(r)
          conn(k) = conn(k) + vert_off
        end do
      end if
      conn_off = conn_off + recv_connsizes(r)
      vert_off = vert_off + recv_nverts(r)
    end do
  end subroutine remap_connectivity_vertex_offsets

  !> Re-open an existing NetCDF output file for appending after a restart.
  !> Queries all dimension/variable IDs and restores the times/steps arrays.
  subroutine reopen_netcdf_file(writer, fname, restart_step)
    type(OutputWriter), intent(inout) :: writer
    character(len=*),   intent(in)    :: fname
    integer(ip),        intent(in), optional :: restart_step

    integer(ip) :: existing_snaps, s, keep
    character(len=64) :: vname

    ! Open existing file for read/write access
    call nc_check(nf90_open(trim(fname), nf90_write, writer%ncid), &
                  'nf90_open('//trim(fname)//')')

    ! Re-acquire dimension IDs
    call nc_check(nf90_inq_dimid(writer%ncid, 'time',               writer%dim_time),      'inq_dimid:time')
    call nc_check(nf90_inq_dimid(writer%ncid, 'n_cells',            writer%dim_cells),     'inq_dimid:n_cells')
    call nc_check(nf90_inq_dimid(writer%ncid, 'n_verts',            writer%dim_verts),     'inq_dimid:n_verts')
    call nc_check(nf90_inq_dimid(writer%ncid, 'connectivity_size',  writer%dim_conn),      'inq_dimid:connectivity_size')
    call nc_check(nf90_inq_dimid(writer%ncid, 'n_cells_p1',         writer%dim_cp1),       'inq_dimid:n_cells_p1')
    call nc_check(nf90_inq_dimid(writer%ncid, 'three',              writer%dim_three),     'inq_dimid:three')
    call nc_check(nf90_inq_dimid(writer%ncid, 'xdmf_topology_size', writer%dim_xdmf_topo),'inq_dimid:xdmf_topology_size')

    ! Query sizes needed at runtime
    call nc_check(nf90_inquire_dimension(writer%ncid, writer%dim_xdmf_topo, len=writer%xdmf_topo_size), &
                  'inquire_dim:xdmf_topo_size')

    ! How many time records already exist?
    call nc_check(nf90_inquire_dimension(writer%ncid, writer%dim_time, len=existing_snaps), &
                  'inquire_dim:time')
    writer%snap_count = existing_snaps

    ! Re-acquire variable IDs
    call nc_check(nf90_inq_varid(writer%ncid, 'time',          writer%var_time),           'inq_varid:time')
    call nc_check(nf90_inq_varid(writer%ncid, 'step',          writer%var_step),           'inq_varid:step')
    call nc_check(nf90_inq_varid(writer%ncid, 'coords',        writer%var_coords),         'inq_varid:coords')
    call nc_check(nf90_inq_varid(writer%ncid, 'coords_xdmf',   writer%var_coords_xdmf),   'inq_varid:coords_xdmf')
    call nc_check(nf90_inq_varid(writer%ncid, 'connectivity',  writer%var_conn),           'inq_varid:connectivity')
    call nc_check(nf90_inq_varid(writer%ncid, 'offsets',       writer%var_offsets),        'inq_varid:offsets')
    call nc_check(nf90_inq_varid(writer%ncid, 'cell_type',     writer%var_ctype),          'inq_varid:cell_type')
    call nc_check(nf90_inq_varid(writer%ncid, 'xdmf_topology', writer%var_xdmf_topology), 'inq_varid:xdmf_topology')
    call nc_check(nf90_inq_varid(writer%ncid, 'p',             writer%var_p),              'inq_varid:p')
    call nc_check(nf90_inq_varid(writer%ncid, 'p_xdmf',        writer%var_p_xdmf),        'inq_varid:p_xdmf')
    call nc_check(nf90_inq_varid(writer%ncid, 'T',             writer%var_T),              'inq_varid:T')
    call nc_check(nf90_inq_varid(writer%ncid, 'T_xdmf',        writer%var_T_xdmf),        'inq_varid:T_xdmf')
    call nc_check(nf90_inq_varid(writer%ncid, 'U',             writer%var_U),              'inq_varid:U')
    call nc_check(nf90_inq_varid(writer%ncid, 'U_xdmf',        writer%var_U_xdmf),        'inq_varid:U_xdmf')

    ! Absent when the original run had no SGS model: fall back to not writing it.
    if (writer%use_sgs) then
      if (nf90_inq_varid(writer%ncid, 'nu_sgs', writer%var_nu_sgs) /= nf90_noerr) then
        writer%var_nu_sgs = -1
        writer%use_sgs    = .false.
      end if
    end if

    if (writer%n_scalars > 0) then
      do s = 1, writer%n_scalars
        vname = writer%scalar_var_names(s)
        call nc_check(nf90_inq_varid(writer%ncid, trim(vname), writer%var_scalar(s)), &
                      'inq_varid:'//trim(vname))
        call nc_check(nf90_inq_varid(writer%ncid, trim(vname)//'_xdmf', writer%var_scalar_xdmf(s)), &
                      'inq_varid:'//trim(vname)//'_xdmf')
      end do
    end if

    ! Restore snapshot metadata so write_temporal_xmf regenerates the full XMF.
    if (existing_snaps > 0) then
      writer%snapshot_capacity = max(64_ip, existing_snaps + 64_ip)
      allocate(writer%times(writer%snapshot_capacity))
      allocate(writer%steps(writer%snapshot_capacity))
      writer%times = 0.0_rp
      writer%steps = 0_ip
      call nc_check(nf90_get_var(writer%ncid, writer%var_time, writer%times(1:existing_snaps), &
                    start=(/1/), count=(/existing_snaps/)), 'get:time_history')
      call nc_check(nf90_get_var(writer%ncid, writer%var_step, writer%steps(1:existing_snaps), &
                    start=(/1/), count=(/existing_snaps/)), 'get:step_history')

      ! Snapshots past the checkpoint belong to the abandoned continuation; the next write overwrites
      ! them in place (the unlimited dimension cannot shrink, so a shorter rerun leaves a stale tail
      ! that the XDMF index, which only lists snap_count entries, does not reference).
      if (present(restart_step)) then
        keep = existing_snaps
        do while (keep > 0)
          if (writer%steps(keep) <= restart_step) exit
          keep = keep - 1
        end do
        if (keep < existing_snaps) then
          write(*, '(A,I0,A,I0,A)') ' Output: dropping ', existing_snaps - keep, &
            ' snapshot(s) written after restart step ', restart_step, '.'
          writer%snap_count = keep
        end if
      end if
    end if
  end subroutine reopen_netcdf_file

  !> Create the NetCDF file, define all dimensions and variables, end define
  !> mode, and write the time-invariant geometry in a single pass.
  subroutine create_netcdf_file(writer, fname, n_cells, n_verts, n_conn, &
                                geom, connectivity, offsets, cell_type, xdmf_topology)
    type(OutputWriter), intent(inout) :: writer
    character(len=*),   intent(in)    :: fname
    integer(ip),        intent(in)    :: n_cells, n_verts, n_conn
    real(rp),           intent(in)    :: geom(:)  ! (3 * n_verts) flat
    integer(ip),        intent(in)    :: connectivity(:)  ! (n_conn)
    integer(ip),        intent(in)    :: offsets(:)  ! (n_cells + 1)
    integer(ip),        intent(in)    :: cell_type(:)  ! (n_cells)
    integer(ip),        intent(in)    :: xdmf_topology(:)

    real(rp), allocatable :: coords_2d(:,:), coords_xdmf(:,:)
    integer(ip) :: n_xdmf_topo, s
    character(len=64) :: vname

    ! ---- Open (NetCDF-4 / HDF5 backend) ----
    call nc_check(nf90_create(trim(fname), ior(nf90_clobber, nf90_netcdf4), writer%ncid), &
                  'nf90_create('//trim(fname)//')')

    ! ---- Global metadata ----
    call nc_check(nf90_put_att(writer%ncid, nf90_global, 'Conventions', 'VTK-CSR, XDMF'), &
                  'att:Conventions')
    call nc_check(nf90_put_att(writer%ncid, nf90_global, 'description', &
                  'DOPAMINE LES solver output - merged unstructured mesh'), 'att:description')
    call nc_check(nf90_put_att(writer%ncid, nf90_global, 'cell_type_convention', &
                  'VTK codes: 1=vertex,7=polygon,10=tet,12=hex,13=wedge,14=pyramid,42=polyhedron'), &
                  'att:cell_type_convention')
    call nc_check(nf90_put_att(writer%ncid, nf90_global, 'connectivity_format', &
                  '0-based VTK CSR: offsets(n_cells+1), connectivity(total_conn)'), &
                  'att:connectivity_format')
    call nc_check(nf90_put_att(writer%ncid, nf90_global, 'n_mpi_ranks', &
                  writer%comm_size), 'att:n_mpi_ranks')

    ! ---- Dimensions ----
    call nc_check(nf90_def_dim(writer%ncid, 'time',               nf90_unlimited, writer%dim_time),  'dim:time')
    call nc_check(nf90_def_dim(writer%ncid, 'n_cells',            n_cells,        writer%dim_cells), 'dim:n_cells')
    call nc_check(nf90_def_dim(writer%ncid, 'n_verts',            n_verts,        writer%dim_verts), 'dim:n_verts')
    call nc_check(nf90_def_dim(writer%ncid, 'connectivity_size',  n_conn,         writer%dim_conn),  'dim:conn')
    call nc_check(nf90_def_dim(writer%ncid, 'n_cells_p1',         n_cells + 1,    writer%dim_cp1),   'dim:n_cells_p1')
    call nc_check(nf90_def_dim(writer%ncid, 'three',              3,              writer%dim_three), 'dim:three')
    n_xdmf_topo = size(xdmf_topology)
    writer%xdmf_topo_size = n_xdmf_topo
    call nc_check(nf90_def_dim(writer%ncid, 'xdmf_topology_size', n_xdmf_topo,    writer%dim_xdmf_topo), 'dim:xdmf_topology')

    ! ---- Geometry variables (time-independent) ----
    ! (three, n_verts) -> Python: coords.shape = (n_verts, 3)
    call nc_check(nf90_def_var(writer%ncid, 'coords', nf90_double, &
                  (/writer%dim_three, writer%dim_verts/), writer%var_coords), 'def:coords')
      call nc_check(nf90_def_var(writer%ncid, 'coords_xdmf', nf90_double, &
              (/writer%dim_verts, writer%dim_three/), writer%var_coords_xdmf), 'def:coords_xdmf')
      call nc_check(nf90_put_att(writer%ncid, writer%var_coords_xdmf, 'long_name', &
              'vertex coordinates (XDMF order n_verts x 3)'), 'att:coords_xdmf:long_name')

    call nc_check(nf90_put_att(writer%ncid, writer%var_coords, 'long_name', 'vertex coordinates'), &
                  'att:coords:long_name')
    call nc_check(nf90_put_att(writer%ncid, writer%var_coords, 'units', 'm'), &
                  'att:coords:units')

    call nc_check(nf90_def_var(writer%ncid, 'connectivity', nf90_int, &
                  (/writer%dim_conn/), writer%var_conn), 'def:connectivity')
    call nc_check(nf90_put_att(writer%ncid, writer%var_conn, 'long_name', &
                  '0-based global vertex indices (CSR, concatenated per cell)'), &
                  'att:connectivity:long_name')

    call nc_check(nf90_def_var(writer%ncid, 'offsets', nf90_int, &
                  (/writer%dim_cp1/), writer%var_offsets), 'def:offsets')
    call nc_check(nf90_put_att(writer%ncid, writer%var_offsets, 'long_name', &
                  'CSR row offsets (size n_cells+1): conn(offsets(c):offsets(c+1)-1) are cell c vertices'), &
                  'att:offsets:long_name')

    call nc_check(nf90_def_var(writer%ncid, 'cell_type', nf90_int, &
                  (/writer%dim_cells/), writer%var_ctype), 'def:cell_type')
      call nc_check(nf90_def_var(writer%ncid, 'xdmf_topology', nf90_int, &
              (/writer%dim_xdmf_topo/), writer%var_xdmf_topology), 'def:xdmf_topology')
      call nc_check(nf90_put_att(writer%ncid, writer%var_xdmf_topology, 'long_name', &
              'XDMF Mixed topology stream'), 'att:xdmf_topology:long_name')

    call nc_check(nf90_put_att(writer%ncid, writer%var_ctype, 'long_name', &
                  'VTK element type per cell'), 'att:cell_type:long_name')

    ! ---- Time / step variables ----
    call nc_check(nf90_def_var(writer%ncid, 'time', nf90_double, &
                  (/writer%dim_time/), writer%var_time), 'def:time')
    call nc_check(nf90_put_att(writer%ncid, writer%var_time, 'long_name', 'simulation time'), &
                  'att:time:long_name')
    call nc_check(nf90_put_att(writer%ncid, writer%var_time, 'units', 's'), &
                  'att:time:units')

    call nc_check(nf90_def_var(writer%ncid, 'step', nf90_int, &
                  (/writer%dim_time/), writer%var_step), 'def:step')
    call nc_check(nf90_put_att(writer%ncid, writer%var_step, 'long_name', 'time-step index'), &
                  'att:step:long_name')

    ! ---- Field variables (time- and cell-varying) ----
    ! (n_cells, time) -> Python: p.shape = (n_time, n_cells)
    call nc_check(nf90_def_var(writer%ncid, 'p', nf90_double, &
                  (/writer%dim_cells, writer%dim_time/), writer%var_p), 'def:p')
      call nc_check(nf90_def_var(writer%ncid, 'p_xdmf', nf90_double, &
              (/writer%dim_time, writer%dim_cells/), writer%var_p_xdmf), 'def:p_xdmf')

    call nc_check(nf90_put_att(writer%ncid, writer%var_p, 'long_name', 'pressure'), 'att:p:long_name')
    call nc_check(nf90_put_att(writer%ncid, writer%var_p, 'units', 'Pa'), 'att:p:units')
    call nc_check(nf90_put_att(writer%ncid, writer%var_p, 'cell_data', 'true'), 'att:p:cell_data')

    call nc_check(nf90_def_var(writer%ncid, 'T', nf90_double, &
                  (/writer%dim_cells, writer%dim_time/), writer%var_T), 'def:T')
      call nc_check(nf90_def_var(writer%ncid, 'T_xdmf', nf90_double, &
              (/writer%dim_time, writer%dim_cells/), writer%var_T_xdmf), 'def:T_xdmf')

    call nc_check(nf90_put_att(writer%ncid, writer%var_T, 'long_name', 'temperature'), 'att:T:long_name')
    call nc_check(nf90_put_att(writer%ncid, writer%var_T, 'units', 'K'), 'att:T:units')
    call nc_check(nf90_put_att(writer%ncid, writer%var_T, 'cell_data', 'true'), 'att:T:cell_data')

    ! (three, n_cells, time) -> Python: U.shape = (n_time, n_cells, 3)
    call nc_check(nf90_def_var(writer%ncid, 'U', nf90_double, &
                  (/writer%dim_three, writer%dim_cells, writer%dim_time/), writer%var_U), 'def:U')
      call nc_check(nf90_def_var(writer%ncid, 'U_xdmf', nf90_double, &
              (/writer%dim_time, writer%dim_cells, writer%dim_three/), writer%var_U_xdmf), 'def:U_xdmf')

    call nc_check(nf90_put_att(writer%ncid, writer%var_U, 'long_name', &
                  'velocity vector (Ux, Uy, Uz)'), 'att:U:long_name')
    call nc_check(nf90_put_att(writer%ncid, writer%var_U, 'units', 'm/s'), 'att:U:units')
    call nc_check(nf90_put_att(writer%ncid, writer%var_U, 'cell_data', 'true'), 'att:U:cell_data')

    ! Defined only when a model is active, so readers can test for its presence.
    if (writer%use_sgs) then
      call nc_check(nf90_def_var(writer%ncid, 'nu_sgs', nf90_double, &
                    (/writer%dim_cells, writer%dim_time/), writer%var_nu_sgs), 'def:nu_sgs')
      call nc_check(nf90_put_att(writer%ncid, writer%var_nu_sgs, 'long_name', &
                    'subgrid-scale eddy viscosity'), 'att:nu_sgs:long_name')
      call nc_check(nf90_put_att(writer%ncid, writer%var_nu_sgs, 'units', 'm2/s'), 'att:nu_sgs:units')
      call nc_check(nf90_put_att(writer%ncid, writer%var_nu_sgs, 'cell_data', 'true'), 'att:nu_sgs:cell_data')
    end if

    if (writer%n_scalars > 0) then
      do s = 1, writer%n_scalars
        ! Split into nested ifs: .and. does not short-circuit in Fortran, so len_trim()
        ! must never be evaluated in the same expression as the allocated() check.
        vname = ''
        if (allocated(writer%scalar_var_names)) then
          if (len_trim(writer%scalar_var_names(s)) > 0) vname = trim(writer%scalar_var_names(s))
        end if
        if (len_trim(vname) == 0) then
          write(vname, '(A,I0)') 'scalar_', s
          writer%scalar_var_names(s) = trim(vname)
        end if
        call nc_check(nf90_def_var(writer%ncid, trim(vname), nf90_double, &
                      (/writer%dim_cells, writer%dim_time/), writer%var_scalar(s)), 'def:'//trim(vname))
        call nc_check(nf90_put_att(writer%ncid, writer%var_scalar(s), 'long_name', trim(vname)), &
                      'att:'//trim(vname)//':long_name')
        call nc_check(nf90_put_att(writer%ncid, writer%var_scalar(s), 'cell_data', 'true'), &
                      'att:'//trim(vname)//':cell_data')

        call nc_check(nf90_def_var(writer%ncid, trim(vname)//'_xdmf', nf90_double, &
                      (/writer%dim_time, writer%dim_cells/), writer%var_scalar_xdmf(s)), &
                      'def:'//trim(vname)//'_xdmf')
      end do
    end if

    call nc_check(nf90_enddef(writer%ncid), 'nf90_enddef')

    ! ---- Write time-invariant geometry ----
    allocate(coords_2d(3, n_verts))
    coords_2d = reshape(geom, (/3, n_verts/))
    call nc_check(nf90_put_var(writer%ncid, writer%var_coords, coords_2d), 'put:coords')
    allocate(coords_xdmf(n_verts, 3))
    coords_xdmf = transpose(coords_2d)
    call nc_check(nf90_put_var(writer%ncid, writer%var_coords_xdmf, coords_xdmf), 'put:coords_xdmf')
    deallocate(coords_2d)
    deallocate(coords_xdmf)

    call nc_check(nf90_put_var(writer%ncid, writer%var_conn,    connectivity), 'put:connectivity')
    call nc_check(nf90_put_var(writer%ncid, writer%var_offsets, offsets),      'put:offsets')
    call nc_check(nf90_put_var(writer%ncid, writer%var_ctype,   cell_type),    'put:cell_type')
    call nc_check(nf90_put_var(writer%ncid, writer%var_xdmf_topology, xdmf_topology), 'put:xdmf_topology')
  end subroutine create_netcdf_file

  !> Append one time record.
  !> p, T : (n_cells) flat.
  subroutine append_netcdf_record(writer, tidx, time, step, p, T, U_flat, n_cells, S_all, nu_sgs)
    type(OutputWriter), intent(inout) :: writer
    integer(ip),        intent(in)    :: tidx, step, n_cells
    real(rp),           intent(in)    :: time
    real(rp),           intent(in)    :: p(:), T(:), U_flat(:)
    real(rp), optional, intent(in)    :: S_all(:,:)
    real(rp), optional, intent(in)    :: nu_sgs(:)

    real(rp)     :: time_buf(1)
    integer(ip)  :: step_buf(1)
    integer(ip)  :: s
    real(rp), allocatable :: u_2d(:,:), u_xdmf(:,:)

    time_buf(1) = time
    step_buf(1) = step

    call nc_check(nf90_put_var(writer%ncid, writer%var_time, time_buf, &
                  start=(/tidx/), count=(/1/)), 'put:time')
    call nc_check(nf90_put_var(writer%ncid, writer%var_step, step_buf, &
                  start=(/tidx/), count=(/1/)), 'put:step')

    ! (n_cells, time): write slice at time index tidx.
    call nc_check(nf90_put_var(writer%ncid, writer%var_p, p, &
                  start=(/1, tidx/), count=(/n_cells, 1/)), 'put:p')
    call nc_check(nf90_put_var(writer%ncid, writer%var_p_xdmf, p, &
                  start=(/tidx, 1/), count=(/1, n_cells/)), 'put:p_xdmf')
    if (writer%boussinesq) then
      call nc_check(nf90_put_var(writer%ncid, writer%var_T, T, &
                    start=(/1, tidx/), count=(/n_cells, 1/)), 'put:T')
      call nc_check(nf90_put_var(writer%ncid, writer%var_T_xdmf, T, &
                    start=(/tidx, 1/), count=(/1, n_cells/)), 'put:T_xdmf')
    end if

    ! (three, n_cells, time): reshape flat buffer to (3, n_cells).
    allocate(u_2d(3, n_cells))
    u_2d = reshape(U_flat, (/3, n_cells/))
    call nc_check(nf90_put_var(writer%ncid, writer%var_U, u_2d, &
                  start=(/1, 1, tidx/), count=(/3, n_cells, 1/)), 'put:U')
    allocate(u_xdmf(n_cells, 3))
    u_xdmf = transpose(u_2d)
    call nc_check(nf90_put_var(writer%ncid, writer%var_U_xdmf, u_xdmf, &
                  start=(/tidx, 1, 1/), count=(/1, n_cells, 3/)), 'put:U_xdmf')
    deallocate(u_xdmf)

    if (writer%use_sgs .and. present(nu_sgs)) then
      call nc_check(nf90_put_var(writer%ncid, writer%var_nu_sgs, nu_sgs, &
                    start=(/1, tidx/), count=(/n_cells, 1/)), 'put:nu_sgs')
    end if

    if (writer%n_scalars > 0 .and. present(S_all)) then
      do s = 1, writer%n_scalars
        call nc_check(nf90_put_var(writer%ncid, writer%var_scalar(s), S_all(s, :), &
                      start=(/1, tidx/), count=(/n_cells, 1/)), 'put:scalar')
        call nc_check(nf90_put_var(writer%ncid, writer%var_scalar_xdmf(s), S_all(s, :), &
                      start=(/tidx, 1/), count=(/1, n_cells/)), 'put:scalar_xdmf')
      end do
    end if
    deallocate(u_2d)
  end subroutine append_netcdf_record


  subroutine append_snapshot_metadata(writer, time, step)
    type(OutputWriter), intent(inout) :: writer
    real(rp),           intent(in)    :: time
    integer(ip),        intent(in)    :: step

    real(rp), allocatable :: times_tmp(:)
    integer(ip), allocatable :: steps_tmp(:)
    integer(ip) :: next_idx, new_cap

    next_idx = writer%snap_count + 1

    if (.not. allocated(writer%times)) then
      writer%snapshot_capacity = 64_ip
      allocate(writer%times(writer%snapshot_capacity), writer%steps(writer%snapshot_capacity))
      writer%times = 0.0_rp
      writer%steps = 0_ip
    end if

    if (writer%snapshot_capacity <= 0) writer%snapshot_capacity = size(writer%times)

    if (next_idx > writer%snapshot_capacity) then
      new_cap = max(2_ip * writer%snapshot_capacity, next_idx)
      allocate(times_tmp(new_cap), steps_tmp(new_cap))
      times_tmp = 0.0_rp
      steps_tmp = 0_ip
      times_tmp(1:writer%snapshot_capacity) = writer%times(1:writer%snapshot_capacity)
      steps_tmp(1:writer%snapshot_capacity) = writer%steps(1:writer%snapshot_capacity)
      call move_alloc(times_tmp, writer%times)
      call move_alloc(steps_tmp, writer%steps)
      writer%snapshot_capacity = new_cap
    end if

    writer%times(next_idx) = time
    writer%steps(next_idx) = step
  end subroutine append_snapshot_metadata

  subroutine write_temporal_xmf(writer)
    type(OutputWriter), intent(in) :: writer

    integer(ip) :: unit_id, s, n_cells, n_verts, n_topo
    character(len=512) :: xmf_file, nc_rel

    if (writer%snap_count <= 0) return

    n_cells = sum(writer%recvcounts)
    n_verts = writer%n_verts_dedup
    n_topo  = writer%xdmf_topo_size
    xmf_file = trim(writer%prefix)//'.xmf'
    nc_rel = trim(writer%prefix)//'.nc'

    open(newunit=unit_id, file=trim(xmf_file), status='replace', action='write')
    write(unit_id,'(A)') '<?xml version="1.0" ?>'
    write(unit_id,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(unit_id,'(A)') '<Xdmf Version="3.0">'
    write(unit_id,'(A)') '  <Domain>'
    write(unit_id,'(A)') '    <Grid Name="dopamine" GridType="Collection" CollectionType="Temporal">'

    do s = 1, writer%snap_count
      call write_xmf_grid(unit_id, nc_rel, n_cells, n_verts, n_topo, writer%snap_count, s, writer%times(s), &
                          writer%boussinesq, writer%n_scalars, writer%scalar_var_names, writer%use_sgs)
    end do

    write(unit_id,'(A)') '    </Grid>'
    write(unit_id,'(A)') '  </Domain>'
    write(unit_id,'(A)') '</Xdmf>'
    close(unit_id)
  end subroutine write_temporal_xmf

  subroutine write_xmf_grid(unit_id, nc_file, n_cells, n_verts, n_topo, ntime, tidx, time_val, has_temp, n_scalars, scalar_names, has_nu_sgs)
    integer(ip),      intent(in) :: unit_id, n_cells, n_verts, n_topo, ntime, tidx
    character(len=*), intent(in) :: nc_file
    real(rp),         intent(in) :: time_val
    logical,          intent(in) :: has_temp
    integer(ip),      intent(in) :: n_scalars
    character(len=*), intent(in) :: scalar_names(:)
    logical,          intent(in) :: has_nu_sgs
    integer(ip) :: s

    write(unit_id,'(A)') '      <Grid Name="step'//trim(i2s(tidx))//'" GridType="Uniform">'
    write(unit_id,'(A,ES24.16,A)') '        <Time Value="', time_val, '"/>'
    write(unit_id,'(A)') '        <Topology TopologyType="Mixed" NumberOfElements="'//trim(i2s(n_cells))//'">'
    write(unit_id,'(A)') '          <DataItem Format="HDF" NumberType="Int" Precision="4" Dimensions="'// &
      trim(i2s(n_topo))//'">'//trim(nc_file)//':/xdmf_topology</DataItem>'
    write(unit_id,'(A)') '        </Topology>'
    write(unit_id,'(A)') '        <Geometry GeometryType="XYZ">'
    ! NetCDF/XDMF expect coords as (n_verts, 3), so use that dimension order.
    write(unit_id,'(A)') '          <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="'// &
      trim(i2s(n_verts))//' '//trim(i2s(3))//'">'// &
      trim(nc_file)//':/coords</DataItem>'
    write(unit_id,'(A)') '        </Geometry>'

    call write_xmf_hyperslab_scalar(unit_id, 'pressure', trim(nc_file)//':/p', ntime, tidx, n_cells)
    if (has_temp) call write_xmf_hyperslab_scalar(unit_id, 'temperature', trim(nc_file)//':/T', ntime, tidx, n_cells)
    call write_xmf_hyperslab_vector(unit_id, 'velocity', trim(nc_file)//':/U', ntime, tidx, n_cells)
    if (has_nu_sgs) call write_xmf_hyperslab_scalar(unit_id, 'nu_sgs', trim(nc_file)//':/nu_sgs', ntime, tidx, n_cells)
    do s = 1, n_scalars
      call write_xmf_hyperslab_scalar(unit_id, trim(scalar_names(s)), trim(nc_file)//':/'//trim(scalar_names(s)), ntime, tidx, n_cells)
    end do

    write(unit_id,'(A)') '      </Grid>'
  end subroutine write_xmf_grid

  subroutine write_xmf_hyperslab_scalar(unit_id, name, hdf_path, ntime, tidx, n_cells)
    integer(ip),      intent(in) :: unit_id, ntime, tidx, n_cells
    character(len=*), intent(in) :: name, hdf_path

    write(unit_id,'(A)') '        <Attribute Name="'//trim(name)//'" AttributeType="Scalar" Center="Cell">'
    write(unit_id,'(A)') '          <DataItem ItemType="HyperSlab" Dimensions="'//trim(i2s(n_cells))//'" Type="HyperSlab">'
    write(unit_id,'(A)') '            <DataItem Dimensions="3 2" Format="XML">'
    write(unit_id,'(A)') '              '//trim(i2s(tidx-1))//' 0  1 1  1 '//trim(i2s(n_cells))
    write(unit_id,'(A)') '            </DataItem>'
    write(unit_id,'(A)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="'// &
      trim(i2s(ntime))//' '//trim(i2s(n_cells))//'">'//trim(hdf_path)//'</DataItem>'
    write(unit_id,'(A)') '          </DataItem>'
    write(unit_id,'(A)') '        </Attribute>'
  end subroutine write_xmf_hyperslab_scalar

  subroutine write_xmf_hyperslab_vector(unit_id, name, hdf_path, ntime, tidx, n_cells)
    integer(ip),      intent(in) :: unit_id, ntime, tidx, n_cells
    character(len=*), intent(in) :: name, hdf_path

    write(unit_id,'(A)') '        <Attribute Name="'//trim(name)//'" AttributeType="Vector" Center="Cell">'
    write(unit_id,'(A)') '          <DataItem ItemType="HyperSlab" Dimensions="'//trim(i2s(n_cells))//' 3" Type="HyperSlab">'
    write(unit_id,'(A)') '            <DataItem Dimensions="3 3" Format="XML">'
    write(unit_id,'(A)') '              '//trim(i2s(tidx-1))//' 0 0  1 1 1  1 '//trim(i2s(n_cells))//' 3'
    write(unit_id,'(A)') '            </DataItem>'
    write(unit_id,'(A)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="'// &
      trim(i2s(ntime))//' '//trim(i2s(n_cells))//' 3">'//trim(hdf_path)//'</DataItem>'
    write(unit_id,'(A)') '          </DataItem>'
    write(unit_id,'(A)') '        </Attribute>'
  end subroutine write_xmf_hyperslab_vector

  subroutine build_local_xdmf_topology(msh, vtx_map, connectivity, cell_type, topo)
    type(Mesh),              intent(in)  :: msh
    integer(ip),             intent(in)  :: vtx_map(:)
    integer(ip),             intent(in)  :: connectivity(:)
    integer(ip),             intent(in)  :: cell_type(:)
    integer(ip), allocatable,intent(out) :: topo(:)

    integer(ip) :: c, p, idx, lf, nfv, nfaces, xcode
    integer(ip) :: c0, c1, nverts, total_size

    total_size = 0
    do c = 1, msh%n_cells
      if (allocated(msh%face_offsets) .and. allocated(msh%face_vertices) .and. cell_type(c) == VTK_POLYHEDRON) then
        nfaces = msh%cells(c)%n_faces
        total_size = total_size + 2
        do idx = msh%cells(c)%face_offset, msh%cells(c)%face_offset + nfaces - 1
          lf = msh%cell_faces(idx)
          nfv = msh%face_offsets(lf + 1) - msh%face_offsets(lf)
          total_size = total_size + 1 + nfv
        end do
      else
        c0 = msh%cell_offsets(c)
        c1 = msh%cell_offsets(c + 1) - 1
        nverts = c1 - c0 + 1
        xcode = xdmf_code_from_vtk(cell_type(c))
        if (xcode == 1_ip) then
          total_size = total_size + 2 + nverts
        else
          total_size = total_size + 1 + nverts
        end if
      end if
    end do

    allocate(topo(total_size))
    p = 1
    do c = 1, msh%n_cells
      if (allocated(msh%face_offsets) .and. allocated(msh%face_vertices) .and. cell_type(c) == VTK_POLYHEDRON) then
        nfaces = msh%cells(c)%n_faces
        topo(p) = 16_ip
        p = p + 1
        topo(p) = nfaces
        p = p + 1
        do idx = msh%cells(c)%face_offset, msh%cells(c)%face_offset + nfaces - 1
          lf = msh%cell_faces(idx)
          c0 = msh%face_offsets(lf)
          c1 = msh%face_offsets(lf + 1) - 1
          nfv = c1 - c0 + 1
          topo(p) = nfv
          p = p + 1
          if (msh%faces(lf)%owner == c) then
            ! Face orientation is stored owner->neighbour in OpenFOAM.
            ! For owner cell keep order; for neighbour cell reverse order so
            topo(p:p+nfv-1) = vtx_map(msh%face_vertices(c0:c1))
          else
            topo(p:p+nfv-1) = vtx_map(msh%face_vertices(c1:c0:-1))
          end if
          p = p + nfv
        end do
      else
        c0 = msh%cell_offsets(c)
        c1 = msh%cell_offsets(c + 1) - 1
        nverts = c1 - c0 + 1
        xcode = xdmf_code_from_vtk(cell_type(c))
        topo(p) = xcode
        p = p + 1
        if (xcode == 1_ip) then
          topo(p) = nverts
          p = p + 1
        end if
        topo(p:p+nverts-1) = connectivity(c0:c1)
        p = p + nverts
      end if
    end do
  end subroutine build_local_xdmf_topology

  subroutine remap_xdmf_topology_vertex_offsets(topo, recv_toposizes, recv_nverts)
    integer(ip), intent(inout) :: topo(:)
    integer(ip), intent(in)    :: recv_toposizes(:), recv_nverts(:)

    integer(ip) :: r, p0, p1, p, xcode, n, nfaces, f
    integer(ip) :: vert_off

    p0 = 1
    vert_off = 0
    do r = 1, size(recv_toposizes)
      p1 = p0 + recv_toposizes(r) - 1
      p = p0
      do while (p <= p1)
        xcode = topo(p)
        p = p + 1
        select case (xcode)
        case (1_ip, 2_ip, 3_ip)
          n = topo(p)
          p = p + 1
          if (n > 0) topo(p:p+n-1) = topo(p:p+n-1) + vert_off
          p = p + n
        case (4_ip)
          topo(p:p+2) = topo(p:p+2) + vert_off
          p = p + 3
        case (5_ip)
          topo(p:p+3) = topo(p:p+3) + vert_off
          p = p + 4
        case (6_ip)
          topo(p:p+3) = topo(p:p+3) + vert_off
          p = p + 4
        case (7_ip)
          topo(p:p+4) = topo(p:p+4) + vert_off
          p = p + 5
        case (8_ip)
          topo(p:p+5) = topo(p:p+5) + vert_off
          p = p + 6
        case (9_ip)
          topo(p:p+7) = topo(p:p+7) + vert_off
          p = p + 8
        case (16_ip)
          nfaces = topo(p)
          p = p + 1
          do f = 1, nfaces
            n = topo(p)
            p = p + 1
            if (n > 0) topo(p:p+n-1) = topo(p:p+n-1) + vert_off
            p = p + n
          end do
        case default
          error stop 'output_writer_init: unsupported XDMF mixed topology code'
        end select
      end do
      p0 = p1 + 1
      vert_off = vert_off + recv_nverts(r)
    end do
  end subroutine remap_xdmf_topology_vertex_offsets

  integer(ip) function xdmf_code_from_vtk(vtk_code) result(xdmf_code)
    integer(ip), intent(in) :: vtk_code
    select case (vtk_code)
    case (VTK_TETRA)
      xdmf_code = 6_ip
    case (VTK_PYRAMID)
      xdmf_code = 7_ip
    case (VTK_WEDGE)
      xdmf_code = 8_ip
    case (VTK_HEXAHEDRON)
      xdmf_code = 9_ip
    case (VTK_POLYHEDRON)
      xdmf_code = 16_ip
    case default
      ! Fallback for unsupported/variable-arity cells.
      xdmf_code = 1_ip
    end select
  end function xdmf_code_from_vtk

  character(len=64) function i2s(v) result(s)
    integer(ip), intent(in) :: v
    write(s, '(i0)') v
  end function i2s

  character(len=64) function sanitize_netcdf_name(name_in, idx) result(name_out)
    character(len=*), intent(in) :: name_in
    integer(ip), intent(in) :: idx
    integer(ip) :: i, n

    name_out = adjustl(trim(name_in))
    n = len_trim(name_out)
    if (n <= 0) then
      write(name_out, '(A,I0)') 'scalar_', idx
      return
    end if

    do i = 1, n
      select case (name_out(i:i))
      case ('a':'z', 'A':'Z', '0':'9', '_')
      case default
        name_out(i:i) = '_'
      end select
    end do

    if (.not. ((name_out(1:1) >= 'A' .and. name_out(1:1) <= 'Z') .or. &
               (name_out(1:1) >= 'a' .and. name_out(1:1) <= 'z') .or. name_out(1:1) == '_')) then
      name_out = 's_'//trim(name_out)
    end if
  end function sanitize_netcdf_name

  ! Vertex deduplication for decomposed meshes:
  ! In partitioned OpenFOAM grids, shared boundary vertices appear twice with

  subroutine deduplicate_geometry(geom_in, n_in, geom_out, n_out, vertex_remap_out)
    real(rp),       intent(in)  :: geom_in(:)  ! (3 * n_in) concatenated coords
    integer(ip),    intent(in)  :: n_in  ! total vertices before dedup
    real(rp), allocatable, intent(out) :: geom_out(:)  ! (3 * n_out) deduplicated
    integer(ip),    intent(out) :: n_out
    integer(ip), allocatable, intent(out) :: vertex_remap_out(:)  ! (n_in) old_idx -> new_idx (0-based)

    integer(ip) :: i, j, k, new_id
    real(rp) :: dx, dy, dz, dist_sq
    real(rp), parameter :: COORD_TOL_SQ = 1.0e-20_rp  ! (1e-10)^2
    integer(ip), allocatable :: sort_idx(:), remap(:)

    ! 1. Build index array and sort by (x, y, z) - O(n log n).
    ! After lexicographic sort, spatial duplicates are guaranteed adjacent.
    allocate(sort_idx(0:n_in-1))
    do i = 0, n_in - 1
      sort_idx(i) = i
    end do
    call sort_vert_idx(sort_idx, geom_in, 0, n_in - 1)

    ! 2. Linear scan: each vertex compared only to its sorted predecessor - O(n).
    allocate(remap(0:n_in-1))
    new_id = 0
    remap(sort_idx(0)) = 0
    do k = 1, n_in - 1
      i = sort_idx(k)
      j = sort_idx(k - 1)
      dx = geom_in(3*i+1) - geom_in(3*j+1)
      dy = geom_in(3*i+2) - geom_in(3*j+2)
      dz = geom_in(3*i+3) - geom_in(3*j+3)
      dist_sq = dx*dx + dy*dy + dz*dz
      if (dist_sq < COORD_TOL_SQ) then
        remap(i) = new_id  ! duplicate of predecessor
      else
        new_id = new_id + 1
        remap(i) = new_id
      end if
    end do
    n_out = new_id + 1
    vertex_remap_out = remap

    ! 3. Build deduplicated geometry - O(n).
    ! Reverse traversal so the first sorted occurrence of each unique vertex wins.
    allocate(geom_out(3 * n_out))
    do k = n_in - 1, 0, -1
      i = sort_idx(k)
      geom_out(3*remap(i)+1:3*remap(i)+3) = geom_in(3*i+1:3*i+3)
    end do
    deallocate(sort_idx, remap)
  end subroutine deduplicate_geometry

  subroutine apply_vertex_dedup_to_connectivity(conn, vertex_remap)
    integer(ip), intent(inout) :: conn(:)
    integer(ip), intent(in)    :: vertex_remap(0:)  ! 0-based like the vertex indices
    integer(ip) :: i

    do i = 1, size(conn)
      conn(i) = vertex_remap(conn(i))
    end do
  end subroutine apply_vertex_dedup_to_connectivity

  ! Structure-aware XDMF mixed-topology remap: walks type codes, counts and
  ! vertex indices correctly, only remapping actual vertex slots.
  subroutine apply_vertex_dedup_to_xdmf_topology(topo, vertex_remap)
    integer(ip), intent(inout) :: topo(:)
    integer(ip), intent(in)    :: vertex_remap(0:)  ! 0-based
    integer(ip) :: p, xcode, n, nfaces, f, i, tlen

    tlen = size(topo)
    p = 1
    do while (p <= tlen)
      xcode = topo(p)
      p = p + 1
      select case (xcode)
      case (1_ip, 2_ip, 3_ip)  ! variable-length: next int is vertex count
        n = topo(p)
        p = p + 1
        do i = 0, n - 1
          topo(p + i) = vertex_remap(topo(p + i))
        end do
        p = p + n
      case (4_ip)  ! triangle: 3 verts
        do i = 0, 2; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 3
      case (5_ip)  ! quad: 4 verts
        do i = 0, 3; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 4
      case (6_ip)  ! tet: 4 verts
        do i = 0, 3; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 4
      case (7_ip)  ! pyramid: 5 verts
        do i = 0, 4; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 5
      case (8_ip)  ! wedge: 6 verts
        do i = 0, 5; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 6
      case (9_ip)  ! hex: 8 verts
        do i = 0, 7; topo(p+i) = vertex_remap(topo(p+i)); end do
        p = p + 8
      case (16_ip)  ! polyhedron
        nfaces = topo(p)
        p = p + 1
        do f = 1, nfaces
          n = topo(p)
          p = p + 1
          do i = 0, n - 1
            topo(p + i) = vertex_remap(topo(p + i))
          end do
          p = p + n
        end do
      case default
        error stop 'apply_vertex_dedup_to_xdmf_topology: unsupported XDMF code'
      end select
    end do
  end subroutine apply_vertex_dedup_to_xdmf_topology

  subroutine nc_check(status, where)
    integer(ip),      intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(*, '(A,1X,A,1X,A)') 'ERROR(NetCDF):', trim(where), trim(nf90_strerror(status))
      error stop 'NetCDF failure in output writer'
    end if
  end subroutine nc_check

  ! Lexicographic (x, y, z) quicksort of a 0-based vertex index array.
  ! Median-of-three pivot + insertion sort below threshold -> O(n log n) average.
  recursive subroutine sort_vert_idx(idx, geom, lo, hi)
    integer(ip), intent(inout) :: idx(0:)
    real(rp),    intent(in)    :: geom(:)
    integer(ip), intent(in)    :: lo, hi
    integer(ip) :: i, j, k, tmp, pivot_val

    if (hi - lo < 16) then
      ! Insertion sort for small subarrays
      do i = lo + 1, hi
        tmp = idx(i)
        j = i - 1
        do while (j >= lo .and. vert_lt_coord(tmp, idx(j), geom))
          idx(j + 1) = idx(j)
          j = j - 1
        end do
        idx(j + 1) = tmp
      end do
      return
    end if

    ! Median-of-three: sort lo, mid, hi to choose a good pivot
    k = (lo + hi) / 2
    if (vert_lt_coord(idx(hi), idx(lo), geom)) then; tmp = idx(lo); idx(lo) = idx(hi); idx(hi) = tmp; end if
    if (vert_lt_coord(idx(k),  idx(lo), geom)) then; tmp = idx(lo); idx(lo) = idx(k);  idx(k)  = tmp; end if
    if (vert_lt_coord(idx(hi), idx(k),  geom)) then; tmp = idx(k);  idx(k)  = idx(hi); idx(hi) = tmp; end if
    pivot_val = idx(k)  ! capture before it may be swapped

    ! Hoare partition
    i = lo; j = hi
    do while (i <= j)
      do while (vert_lt_coord(idx(i), pivot_val, geom)); i = i + 1; end do
      do while (vert_lt_coord(pivot_val, idx(j), geom)); j = j - 1; end do
      if (i <= j) then
        tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
        i = i + 1; j = j - 1
      end if
    end do

    if (lo < j) call sort_vert_idx(idx, geom, lo, j)
    if (i < hi) call sort_vert_idx(idx, geom, i, hi)
  end subroutine sort_vert_idx

  pure logical function vert_lt_coord(a, b, geom)
    integer(ip), intent(in) :: a, b
    real(rp),    intent(in) :: geom(:)
    real(rp) :: ax, ay, az, bx, by, bz
    logical :: ex, ey
    ax = geom(3*a+1); ay = geom(3*a+2); az = geom(3*a+3)
    bx = geom(3*b+1); by = geom(3*b+2); bz = geom(3*b+3)
    ! Use IEEE-safe "not less and not greater" instead of == to avoid -Wcompare-reals
    ex = .not. (ax < bx) .and. .not. (bx < ax)
    ey = .not. (ay < by) .and. .not. (by < ay)
    vert_lt_coord = (ax < bx) .or. (ex .and. ay < by) .or. (ex .and. ey .and. az < bz)
  end function vert_lt_coord

end module mod_output_writer
