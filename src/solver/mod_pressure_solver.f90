! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Pressure Poisson solver with PETSc backend
module mod_pressure_solver
  use mod_kinds
  use mod_parameters, only: BC_OUTLET, BC_PERIODIC, SMALL, THREADED_FACE_LOOPS
  use mod_mesh
  use mod_bc, only: BCConfig
  use mod_warn_log, only: warn_log_bump, warn_log_message, WARN_PSOLVER_RESIDUAL
  use iso_c_binding
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: pressure_solver_handle_t
  public :: pressure_solver_init
  public :: pressure_solver_solve
  public :: pressure_solver_free

  logical, save :: petsc_runtime_active = .false.

  interface
    integer(c_int) function fvm_petsc_backend_available() bind(C)
      use iso_c_binding
    end function

    integer(c_int) function fvm_petsc_backend_initialize() bind(C)
      use iso_c_binding
    end function

    integer(c_int) function fvm_petsc_backend_finalize() bind(C)
      use iso_c_binding
    end function

    !> Set the preconditioner type used for all subsequent pressure solves.
    !> Must be called after fvm_petsc_backend_initialize and before the first solve.
    subroutine fvm_petsc_set_pressure_solver(solver_name) bind(C)
      use iso_c_binding
      character(kind=c_char), intent(in) :: solver_name(*)
    end subroutine

    integer(c_int) function fvm_petsc_device_io_supported() bind(C)
      use iso_c_binding
    end function

    !> Solve with rhs/solution given as device pointers (0 = host pointers again).
    subroutine fvm_petsc_set_io_on_device(on) bind(C)
      use iso_c_binding
      integer(c_int), value :: on
    end subroutine

    !> Insert a single option string into PETSc's global options database.
    !> opt_str must be null-terminated, e.g. "-pc_hypre_boomeramg_strong_threshold 0.5"
    subroutine fvm_petsc_insert_option(opt_str) bind(C)
      use iso_c_binding
      character(kind=c_char), intent(in) :: opt_str(*)
    end subroutine

    integer(c_int) function fvm_petsc_solve_pressure(n_cells_local, n_cells_global, n_faces, local_to_global, face_owner, &
      face_neighbor_local, face_neighbor_global, skip_face, diagonal, rhs, face_coeff, reference_cell, max_iterations, &
      rtol, atol, dirichlet_rows, dirichlet_values, n_dirichlet_rows, solution, iteration_count, residual_norm) bind(C)
      use iso_c_binding
      integer(c_int), value :: n_cells_local, n_cells_global, n_faces, reference_cell, max_iterations, n_dirichlet_rows
      integer(c_int), intent(in) :: local_to_global(*), face_owner(*), face_neighbor_local(*), face_neighbor_global(*), skip_face(*)
      real(c_double), intent(in) :: diagonal(*), rhs(*), face_coeff(*)
      real(c_double), value :: rtol, atol
      integer(c_int), intent(in) :: dirichlet_rows(*)
      real(c_double), intent(in) :: dirichlet_values(*)
      real(c_double), intent(out) :: solution(*)
      integer(c_int), intent(out) :: iteration_count
      real(c_double), intent(out) :: residual_norm
    end function
  end interface

  type, public :: pressure_solver_handle_t
    integer(ip) :: n_cells = 0
    integer(ip) :: n_cells_global = 0
    integer(ip) :: n_faces = 0
    integer(ip) :: reference_cell = 1
    integer(ip) :: n_dirichlet_rows = 0
    integer(c_int), allocatable :: local_to_global(:)
    integer(c_int), allocatable :: face_owner(:), face_neighbor_local(:), face_neighbor_global(:), skip_face(:)
    real(c_double), allocatable :: diagonal(:), face_coeff(:)
    real(c_double), allocatable :: rhs_work(:), p_work(:)
    ! Outlet (Dirichlet) BC: stiffness c_b = area/d_perp added to diagonal
    ! at solver assembly; rhs_work(o) += c_b * p_val applied per outlet face at solve time.
    integer(c_int), allocatable :: dirichlet_rows(:)    ! local owner index per outlet face
    real(c_double), allocatable :: dirichlet_values(:)  ! c_b * p_val per outlet face
    integer(ip), allocatable :: face_mesh_idx(:)
    integer(ip), allocatable :: face_neighbor_ext(:)
    ! Own idx-space coloring for fv_solver's scatter-add loops; idx differs from msh%face_color_perm (periodic edges compacted).
    integer(ip)              :: n_colors = 0
    integer(ip), allocatable :: color_offsets(:)     ! size n_colors+1
    integer(ip), allocatable :: idx_color_perm(:)    ! permutation of 1..n_faces
  end type pressure_solver_handle_t

contains

  subroutine pressure_solver_init(handle, msh, bc_configs, solver_type, petsc_options)
    type(pressure_solver_handle_t), intent(inout) :: handle
    type(Mesh),                     intent(in)    :: msh
    type(BCConfig),       optional, intent(in)    :: bc_configs(:)
    character(len=*),     optional, intent(in)    :: solver_type
    character(len=*),     optional, intent(in)    :: petsc_options(:)

    integer(ip) :: nc, nf_int, nf, q, f_loc, f_bnd, f_glob, o, idx, n_outlet_faces, nb
    integer(ip) :: n_periodic_edges, pair_bnd, pair_face
    integer(ip) :: n_skip_int, n_skip_per, n_used_per, n_dortho_fallback
    integer(ip) :: ref_local, ref_global, ierr_local, i_opt
    real(rp) :: coeff, area_f, d_per, d_perp_b, c_b
    real(rp) :: pair_dcf
    real(rp) :: d_on_vec(3), d_ortho
    integer(ip) :: rank, nranks

    if (fvm_petsc_backend_available() == 0_c_int) then
      error stop 'pressure_solver_init: PETSc backend unavailable'
    end if
    if (.not. petsc_runtime_active) then
      if (fvm_petsc_backend_initialize() /= 0_c_int) then
        error stop 'pressure_solver_init: PETSc initialization failed'
      end if
      petsc_runtime_active = .true.
    end if

    ! Set preconditioner before any KSP is created. Supported: 'boomeramg' | 'gamg' | 'gmres'.
    if (present(solver_type)) then
      if (trim(solver_type) == 'gmres') then
        ! GMRES matrix-free: use ASM preconditioner
        call fvm_petsc_set_pressure_solver('gmres'//c_null_char)
      else
        call fvm_petsc_set_pressure_solver(trim(solver_type)//c_null_char)
      end if
    else
      call fvm_petsc_set_pressure_solver('boomeramg'//c_null_char)
    end if

    ! Inject user-supplied PETSc options so KSPSetFromOptions picks them up.
    if (present(petsc_options)) then
      do i_opt = 1, size(petsc_options)
        if (len_trim(petsc_options(i_opt)) > 0) then
          call fvm_petsc_insert_option(trim(petsc_options(i_opt))//c_null_char)
        end if
      end do
    end if

    nc = msh%n_cells
    nf_int = msh%n_int_faces

    ! Audit periodic mapping (single-rank and MPI modes); prefers msh%periodic_neighbor.
    if (present(bc_configs)) then
      call check_periodic_consistency(msh, bc_configs)
    end if

    n_periodic_edges = 0
    if (present(bc_configs) .and. (allocated(msh%periodic_neighbor) .or. allocated(msh%periodic_pair))) then
      do q = 1, msh%n_patches
        if (bc_configs(q)%bc_id /= BC_PERIODIC) cycle
        do f_loc = 1, msh%patches(q)%n_faces
          f_bnd = msh%patches(q)%start_face + f_loc - 1
          nb = 0
          pair_bnd = 0
          if (allocated(msh%periodic_neighbor)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) nb = msh%periodic_neighbor(f_bnd)
          end if
          if (nb <= 0 .and. allocated(msh%periodic_pair)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) pair_bnd = msh%periodic_pair(f_bnd)
            if (pair_bnd > 0) then
              pair_face = msh%n_int_faces + pair_bnd
              if (pair_face >= 1 .and. pair_face <= size(msh%faces)) nb = msh%faces(pair_face)%owner
            end if
          end if
          if (nb <= 0) cycle
          f_glob = msh%n_int_faces + f_bnd
          o = msh%faces(f_glob)%owner
          if (o < 1 .or. nb < 1 .or. o == nb) cycle
          ! Canonicalize: each same-rank periodic pair contributes one matrix edge.
          if (nb <= nc) then
            if (msh%global_id(o) >= msh%global_id(nb)) cycle
          end if
          n_periodic_edges = n_periodic_edges + 1
        end do
      end do
    end if
    nf = nf_int + n_periodic_edges

    handle%n_cells = nc
    handle%n_cells_global = int(msh%n_cells_global, ip)
    handle%n_faces = nf
    if (nc > 0 .and. allocated(msh%compact_id)) then
      ref_local = minval(int(msh%compact_id(1:nc), ip))
    else
      ref_local = huge(1_ip)
    end if
    call MPI_Allreduce(ref_local, ref_global, 1, MPI_INTEGER, MPI_MIN, msh%comm, ierr_local)
    if (ref_global < 1 .or. ref_global > handle%n_cells_global) ref_global = 1
    handle%reference_cell = ref_global

    allocate(handle%local_to_global(nc))
    allocate(handle%face_owner(nf), handle%face_neighbor_local(nf), handle%face_neighbor_global(nf), handle%skip_face(nf))
    allocate(handle%diagonal(nc), handle%face_coeff(nf))
    allocate(handle%rhs_work(nc), handle%p_work(nc))
    allocate(handle%face_mesh_idx(nf), handle%face_neighbor_ext(nf))

    if (allocated(msh%compact_id)) then
      handle%local_to_global = int(msh%compact_id(1:nc), c_int)
    else
      do idx = 1, nc
        handle%local_to_global(idx) = int(idx, c_int)
      end do
    end if
    handle%face_neighbor_local = 0_c_int
    handle%face_neighbor_global = 0_c_int
    handle%skip_face = 0_c_int
    handle%diagonal = 0.0_c_double
    handle%face_coeff = 0.0_c_double
    handle%face_mesh_idx = 0_ip
    handle%face_neighbor_ext = 0_ip

    do idx = 1, nf_int
      handle%face_owner(idx) = int(msh%faces(idx)%owner, c_int)
      handle%face_mesh_idx(idx) = idx
      if (msh%faces(idx)%owner < 1 .or. msh%faces(idx)%owner > nc) then
        handle%skip_face(idx) = 1_c_int
      else if (msh%faces(idx)%neighbour < 1 .or. msh%faces(idx)%neighbour > nc + msh%n_ghost) then
        handle%skip_face(idx) = 1_c_int
      else
        nb = msh%faces(idx)%neighbour
        if (nb <= nc) handle%face_neighbor_local(idx) = int(nb, c_int)
        ! Same identity fallback as local_to_global above: without compact_id (e.g. a
        ! hand-built single-rank mesh), the local index already is the global id.
        if (allocated(msh%compact_id)) then
          handle%face_neighbor_global(idx) = int(msh%compact_id(nb), c_int)
        else
          handle%face_neighbor_global(idx) = int(nb, c_int)
        end if
        handle%face_neighbor_ext(idx) = nb
      end if
    end do

    idx = nf_int
    if (n_periodic_edges > 0) then
      do q = 1, msh%n_patches
        if (bc_configs(q)%bc_id /= BC_PERIODIC) cycle
        do f_loc = 1, msh%patches(q)%n_faces
          f_bnd = msh%patches(q)%start_face + f_loc - 1
          pair_bnd = 0
          pair_face = 0
          nb = 0
          if (allocated(msh%periodic_neighbor)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) nb = msh%periodic_neighbor(f_bnd)
          end if
          if (nb <= 0 .and. allocated(msh%periodic_pair)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) pair_bnd = msh%periodic_pair(f_bnd)
            if (pair_bnd > 0) then
              pair_face = msh%n_int_faces + pair_bnd
              if (pair_face >= 1 .and. pair_face <= size(msh%faces)) nb = msh%faces(pair_face)%owner
            end if
          end if
          if (nb <= 0) cycle
          f_glob = msh%n_int_faces + f_bnd
          o = msh%faces(f_glob)%owner
          if (o < 1 .or. nb < 1 .or. o == nb) cycle
          ! Canonicalize: each same-rank periodic pair contributes one matrix edge.
          if (nb <= nc) then
            if (msh%global_id(o) >= msh%global_id(nb)) cycle
          end if

          if (idx >= nf) then
            error stop 'pressure_solver_init: periodic edge count mismatch'
          end if
          idx = idx + 1
          handle%face_owner(idx) = int(o, c_int)
          handle%face_mesh_idx(idx) = f_glob
          handle%face_neighbor_ext(idx) = nb
          ! Range-check nb before indexing compact_id (sized n_cells + n_ghost).
          if (o > nc .or. nb > nc + msh%n_ghost) then
            handle%skip_face(idx) = 1_c_int
            cycle
          end if
          if (nb >= 1 .and. nb <= nc) handle%face_neighbor_local(idx) = int(nb, c_int)
          if (allocated(msh%compact_id)) then
            handle%face_neighbor_global(idx) = int(msh%compact_id(nb), c_int)
          else
            handle%face_neighbor_global(idx) = int(nb, c_int)
          end if
          ! One coefficient for both sides of a periodic edge: un-projected d_cf(own) + d_cf(partner), like the periodic operators.
          pair_dcf = 0.0_rp
          if (allocated(msh%periodic_pair_dcf)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair_dcf)) pair_dcf = msh%periodic_pair_dcf(f_bnd)
          end if
          if (allocated(msh%periodic_pair)) then
            if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) pair_bnd = msh%periodic_pair(f_bnd)
          end if
          pair_face = 0
          if (pair_bnd > 0) pair_face = msh%n_int_faces + pair_bnd
          if (pair_dcf <= SMALL .and. pair_face >= 1 .and. pair_face <= size(msh%faces)) &
            pair_dcf = msh%faces(pair_face)%d_cf
          if (pair_dcf <= SMALL) pair_dcf = msh%faces(f_glob)%d_cf
          d_per = msh%faces(f_glob)%d_cf + pair_dcf
          area_f = msh%faces(f_glob)%area
          if (pair_face >= 1 .and. pair_face <= size(msh%faces)) &
            area_f = 0.5_rp * (area_f + msh%faces(pair_face)%area)
          if (.not. ieee_is_finite(area_f) .or. .not. ieee_is_finite(d_per)) then
            handle%skip_face(idx) = 1_c_int
            cycle
          end if
          if (area_f <= SMALL .or. d_per <= SMALL) then
            handle%skip_face(idx) = 1_c_int
            cycle
          end if
          coeff = area_f / (d_per + SMALL)
          if (.not. ieee_is_finite(coeff) .or. coeff <= SMALL) then
            handle%skip_face(idx) = 1_c_int
            cycle
          end if
          handle%face_coeff(idx) = real(coeff, c_double)
          handle%diagonal(o) = handle%diagonal(o) + real(coeff, c_double)
          if (nb >= 1 .and. nb <= nc) handle%diagonal(nb) = handle%diagonal(nb) + real(coeff, c_double)
        end do
      end do
    end if

    n_dortho_fallback = 0
    do idx = 1, nf_int
      o = msh%faces(idx)%owner
      nb = msh%faces(idx)%neighbour
      if (o < 1 .or. nb < 1) cycle
      if (o > nc .or. nb > nc + msh%n_ghost) cycle
      if (.not. ieee_is_finite(msh%faces(idx)%area) .or. .not. ieee_is_finite(msh%faces(idx)%d_cf)) then
        handle%skip_face(idx) = 1_c_int
        cycle
      end if
      if (msh%faces(idx)%area <= SMALL .or. msh%faces(idx)%d_cf <= SMALL) then
        handle%skip_face(idx) = 1_c_int
        cycle
      end if
      ! FPS Sect. 9.7.1: face Laplacian coefficient |S|^2 / (S . d_oc), positive on a well-formed mesh.
      d_on_vec(1) = msh%cells(nb)%centroid(1) - msh%cells(o)%centroid(1)
      d_on_vec(2) = msh%cells(nb)%centroid(2) - msh%cells(o)%centroid(2)
      d_on_vec(3) = msh%cells(nb)%centroid(3) - msh%cells(o)%centroid(3)
      d_ortho = dot_product(d_on_vec, msh%faces(idx)%normal)
      if (d_ortho > SMALL) then
        coeff = msh%faces(idx)%area / d_ortho
      else
        ! Non-convex or inverted cell: this substitutes a different operator, so do not hide it.
        coeff = msh%faces(idx)%area / (msh%faces(idx)%d_cf + SMALL)
        n_dortho_fallback = n_dortho_fallback + 1_ip
      end if
      if (.not. ieee_is_finite(coeff) .or. coeff <= SMALL) then
        handle%skip_face(idx) = 1_c_int
        cycle
      end if
      handle%face_coeff(idx) = real(coeff, c_double)
      handle%diagonal(o) = handle%diagonal(o) + real(coeff, c_double)
      if (nb >= 1 .and. nb <= nc) handle%diagonal(nb) = handle%diagonal(nb) + real(coeff, c_double)
    end do

    ! Outlet Dirichlet: c_b = area/d_perp added to diagonal; stored for rhs accumulation at solve.
    n_outlet_faces = 0
    if (present(bc_configs)) then
      do q = 1, msh%n_patches
        if (bc_configs(q)%bc_id /= BC_OUTLET) cycle
        do f_loc = 1, msh%patches(q)%n_faces
          f_bnd = msh%patches(q)%start_face + f_loc - 1
          f_glob = msh%n_int_faces + f_bnd
          o = msh%faces(f_glob)%owner
          if (o >= 1 .and. o <= nc) n_outlet_faces = n_outlet_faces + 1
        end do
      end do
    end if

    handle%n_dirichlet_rows = n_outlet_faces

    if (n_outlet_faces > 0) then
      allocate(handle%dirichlet_rows(n_outlet_faces))
      allocate(handle%dirichlet_values(n_outlet_faces))
      idx = 0
      do q = 1, msh%n_patches
        if (bc_configs(q)%bc_id /= BC_OUTLET) cycle
        do f_loc = 1, msh%patches(q)%n_faces
          f_bnd = msh%patches(q)%start_face + f_loc - 1
          f_glob = msh%n_int_faces + f_bnd
          o = msh%faces(f_glob)%owner
          if (o < 1 .or. o > nc) cycle
          d_perp_b = abs(dot_product(msh%faces(f_glob)%centre - msh%cells(o)%centroid, &
                                     msh%faces(f_glob)%normal))
          if (d_perp_b <= SMALL) d_perp_b = msh%faces(f_glob)%d_cf
          if (d_perp_b <= SMALL) cycle
          c_b = msh%faces(f_glob)%area / d_perp_b
          if (.not. ieee_is_finite(c_b) .or. c_b <= SMALL) cycle
          handle%diagonal(o) = handle%diagonal(o) + real(c_b, c_double)
          idx = idx + 1
          handle%dirichlet_rows(idx)   = int(o, c_int)
          handle%dirichlet_values(idx) = real(c_b * bc_configs(q)%p_val, c_double)
        end do
      end do
      handle%n_dirichlet_rows = idx
    end if

    rank = msh%rank
    nranks = msh%nranks

    do idx = 1, nc
      if (.not. ieee_is_finite(handle%diagonal(idx)) .or. handle%diagonal(idx) <= real(SMALL, c_double)) then
        handle%diagonal(idx) = 1.0_c_double
      end if
    end do

    n_skip_int = 0
    n_skip_per = 0
    n_used_per = 0
    if (nf_int > 0) n_skip_int = count(handle%skip_face(1:nf_int) /= 0_c_int)
    if (n_periodic_edges > 0) then
      n_skip_per = count(handle%skip_face(nf_int+1:nf) /= 0_c_int)
      n_used_per = n_periodic_edges - n_skip_per
    end if
    if (msh%rank == 0) then
      write(*, '(A,I0,A,I0,A,I0,A,I0,A,I0)') ' Pressure assembly summary: internal_faces=', nf_int, &
        ', internal_skipped=', n_skip_int, ', periodic_edges=', n_periodic_edges, &
        ', periodic_used=', n_used_per, ', periodic_skipped=', n_skip_per
    end if
    if (n_dortho_fallback > 0) then
      block
        character(80) :: msg_text
        write(msg_text, '(I0,A)') n_dortho_fallback, ' faces: Sf.d_on<=0, Poisson used the d_cf fallback'
        call warn_log_message('[psolver]   ', trim(msg_text))
      end block
    end if

    if (THREADED_FACE_LOOPS) call build_psolver_face_coloring(handle, nc + msh%n_ghost)
  end subroutine pressure_solver_init

  !> Greedy face coloring over the handle's idx-space (mirrors mesh_build_face_coloring).
  subroutine build_psolver_face_coloring(handle, ntot)
    type(pressure_solver_handle_t), intent(inout) :: handle
    integer(ip),                    intent(in)    :: ntot  ! n_cells + n_ghost

    integer(ip), allocatable :: face_color(:), counts(:)
    logical, allocatable :: cell_color_used(:,:)
    integer(ip) :: idx, o, nb, col, max_colors, c

    if (handle%n_faces <= 0 .or. ntot <= 0) return

    max_colors = 32_ip  ! safety slack; grown below if ever exceeded
    allocate(face_color(handle%n_faces), source=0_ip)
    allocate(cell_color_used(max_colors, ntot), source=.false.)

    do idx = 1, handle%n_faces
      o  = int(handle%face_owner(idx), ip)
      nb = handle%face_neighbor_ext(idx)
      do col = 1, max_colors
        if (o >= 1 .and. o <= ntot) then
          if (cell_color_used(col, o)) cycle
        end if
        if (nb >= 1 .and. nb <= ntot) then
          if (cell_color_used(col, nb)) cycle
        end if
        exit
      end do
      if (col > max_colors) then
        block
          logical, allocatable :: tmp(:,:)
          integer(ip) :: new_size
          new_size = max_colors * 2_ip
          allocate(tmp(new_size, ntot), source=.false.)
          tmp(1:max_colors, :) = cell_color_used
          call move_alloc(tmp, cell_color_used)
          col = max_colors + 1_ip
          max_colors = new_size
        end block
      end if
      face_color(idx) = col
      if (o  >= 1 .and. o  <= ntot) cell_color_used(col, o)  = .true.
      if (nb >= 1 .and. nb <= ntot) cell_color_used(col, nb) = .true.
    end do

    handle%n_colors = maxval(face_color)
    if (handle%n_colors <= 0) then
      deallocate(face_color, cell_color_used)
      return
    end if

    allocate(counts(handle%n_colors), source=0_ip)
    do idx = 1, handle%n_faces
      counts(face_color(idx)) = counts(face_color(idx)) + 1
    end do

    allocate(handle%color_offsets(handle%n_colors + 1))
    handle%color_offsets(1) = 1
    do c = 1, handle%n_colors
      handle%color_offsets(c+1) = handle%color_offsets(c) + counts(c)
    end do

    allocate(handle%idx_color_perm(handle%n_faces))
    counts = 0_ip
    do idx = 1, handle%n_faces
      col = face_color(idx)
      handle%idx_color_perm(handle%color_offsets(col) + counts(col)) = idx
      counts(col) = counts(col) + 1_ip
    end do

    deallocate(face_color, cell_color_used, counts)
  end subroutine build_psolver_face_coloring

  subroutine check_periodic_consistency(msh, bc_configs)
    type(Mesh),      intent(in) :: msh
    type(BCConfig),  intent(in) :: bc_configs(:)

    integer(ip) :: q, f_loc, f_bnd, f_glob, o, nb, pair_bnd, pair_face
    integer(ip) :: local_faces, local_resolved, local_unresolved
    integer(ip) :: local_neighbor_oob, local_pair_oob, local_self_pair
    integer(ip) :: local_recip_mismatch, local_nb_cell_oob
    integer(ip) :: g_faces, g_resolved, g_unresolved
    integer(ip) :: g_neighbor_oob, g_pair_oob, g_self_pair
    integer(ip) :: g_recip_mismatch, g_nb_cell_oob
    integer(ip) :: ierr

    local_faces = 0
    local_resolved = 0
    local_unresolved = 0
    local_neighbor_oob = 0
    local_pair_oob = 0
    local_self_pair = 0
    local_recip_mismatch = 0
    local_nb_cell_oob = 0

    do q = 1, msh%n_patches
      if (bc_configs(q)%bc_id /= BC_PERIODIC) cycle
      do f_loc = 1, msh%patches(q)%n_faces
        f_bnd = msh%patches(q)%start_face + f_loc - 1
        f_glob = msh%n_int_faces + f_bnd
        if (f_glob < 1 .or. f_glob > size(msh%faces)) cycle
        local_faces = local_faces + 1

        o = msh%faces(f_glob)%owner
        nb = 0

        if (allocated(msh%periodic_neighbor)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) then
            nb = msh%periodic_neighbor(f_bnd)
          else
            local_neighbor_oob = local_neighbor_oob + 1
          end if
        end if

        if (nb <= 0 .and. allocated(msh%periodic_pair)) then
          pair_bnd = 0
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) then
            pair_bnd = msh%periodic_pair(f_bnd)
          else
            local_pair_oob = local_pair_oob + 1
          end if

          if (pair_bnd > 0) then
            if (pair_bnd < 1 .or. pair_bnd > msh%n_bnd_faces) then
              local_pair_oob = local_pair_oob + 1
            else
              pair_face = msh%n_int_faces + pair_bnd
              if (pair_face < 1 .or. pair_face > size(msh%faces)) then
                local_pair_oob = local_pair_oob + 1
              else
                nb = msh%faces(pair_face)%owner
                if (allocated(msh%periodic_pair)) then
                  if (pair_bnd >= 1 .and. pair_bnd <= size(msh%periodic_pair)) then
                    if (msh%periodic_pair(pair_bnd) /= f_bnd) local_recip_mismatch = local_recip_mismatch + 1
                  else
                    local_recip_mismatch = local_recip_mismatch + 1
                  end if
                end if
              end if
            end if
          end if
        end if

        if (nb <= 0) then
          local_unresolved = local_unresolved + 1
          cycle
        end if

        if (nb == o) then
          local_self_pair = local_self_pair + 1
          cycle
        end if

        if (nb > msh%n_cells + msh%n_ghost) then
          local_nb_cell_oob = local_nb_cell_oob + 1
          cycle
        end if

        local_resolved = local_resolved + 1
      end do
    end do

    call MPI_Allreduce(local_faces, g_faces, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_resolved, g_resolved, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_unresolved, g_unresolved, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_neighbor_oob, g_neighbor_oob, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_pair_oob, g_pair_oob, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_self_pair, g_self_pair, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_recip_mismatch, g_recip_mismatch, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_nb_cell_oob, g_nb_cell_oob, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)

    if (msh%rank == 0) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
        ' Periodic consistency check: faces=', g_faces, ', resolved=', g_resolved, &
        ', unresolved=', g_unresolved, ', neighbor_oob=', g_neighbor_oob, &
        ', pair_oob=', g_pair_oob, ', self_pair=', g_self_pair, &
        ', reciprocity_mismatch=', g_recip_mismatch, ', nb_cell_oob=', g_nb_cell_oob
    end if

    if (g_unresolved > 0 .or. g_self_pair > 0 .or. g_neighbor_oob > 0 .or. &
        g_pair_oob > 0 .or. g_recip_mismatch > 0 .or. g_nb_cell_oob > 0) then
      error stop 'check_periodic_consistency: invalid periodic mapping detected'
    end if
  end subroutine check_periodic_consistency

  subroutine pressure_solver_solve(handle, rhs, p, tol, max_iter, n_iter, residual, atol_in, on_device, p_is_total)
    type(pressure_solver_handle_t), intent(inout), target :: handle
    real(rp),                       intent(in)    :: rhs(:)
    real(rp),                       intent(inout) :: p(:)
    !> rhs and p are device-resident (OpenACC): the solve then runs on device pointers with no host copy.
    logical,              optional, intent(in)    :: on_device
    real(rp),             optional, intent(in)    :: tol
    integer(ip),          optional, intent(in)    :: max_iter
    integer(ip),          optional, intent(out)   :: n_iter
    real(rp),             optional, intent(out)   :: residual
    real(rp),             optional, intent(in)    :: atol_in  ! absolute floor, used alongside rtol when present
    !> p is the whole pressure (initial projection): outlet faces then carry p_val. Default: p is the
    !> increment phi, which is zero on the outlet because the boundary pressure already equals p_val.
    logical,              optional, intent(in)    :: p_is_total

    real(c_double) :: c_tol, c_atol, c_residual
    integer(c_int) :: c_max_iter, c_n_iter, status
    integer(ip) :: nc, o, idx, nd
    integer(c_int) :: empty_rows(1)
    real(c_double) :: empty_vals(1)
    logical :: dev, host_sync, dev_io_ok, use_outlet_val
    real(c_double), pointer :: rw(:), pw(:), dvals(:)
    integer(c_int), pointer :: drows(:)

    nc = handle%n_cells
    use_outlet_val = .false.
    if (present(p_is_total)) use_outlet_val = p_is_total
    ! rtol: primary convergence criterion; atol_in: dimensional floor for near-zero RHS only.
    c_tol = real(1.0e-8_rp, c_double)
    if (present(tol)) c_tol = real(tol, c_double)
    if (present(atol_in)) then
      c_atol = real(atol_in, c_double)
    else
      c_atol = real(1.0e-50_rp, c_double)  ! effectively disabled; rtol governs
    end if
    c_max_iter = int(500, c_int)
    if (present(max_iter)) c_max_iter = int(max_iter, c_int)

    ! Split into two ifs: .or. does not short-circuit in Fortran, so size() must never
    ! be evaluated in the same expression as the allocated() check.
    block
      logical :: needs_alloc
      needs_alloc = .not. allocated(handle%rhs_work)
      if (.not. needs_alloc) needs_alloc = (size(handle%rhs_work) < nc)
      if (needs_alloc) then
        if (allocated(handle%rhs_work)) deallocate(handle%rhs_work)
        allocate(handle%rhs_work(nc))
      end if
      needs_alloc = .not. allocated(handle%p_work)
      if (.not. needs_alloc) needs_alloc = (size(handle%p_work) < nc)
      if (needs_alloc) then
        if (allocated(handle%p_work)) deallocate(handle%p_work)
        allocate(handle%p_work(nc))
      end if
    end block

    dev = .false.
    if (present(on_device)) dev = on_device
    ! Backend without CUDA-runtime copies: fall back to host arrays with explicit updates.
    host_sync = .false.
    dev_io_ok = .true.
    if (dev) dev_io_ok = (fvm_petsc_device_io_supported() /= 0_c_int)
    if (dev .and. .not. dev_io_ok) then
      dev = .false.
      host_sync = .true.
      !$acc update host(rhs, p)
    end if
    rw => handle%rhs_work
    pw => handle%p_work

    if (dev) then
      ! Device path: build the work vectors on the device (with the NaN scrub the host path does
      ! in the backend) and hand device pointers to PETSc below.
      nd = 0
      if (use_outlet_val .and. handle%n_dirichlet_rows > 0 .and. allocated(handle%dirichlet_rows)) nd = handle%n_dirichlet_rows
      if (nd > 0) then
        drows => handle%dirichlet_rows
        dvals => handle%dirichlet_values
      else
        drows => null()
        dvals => null()
      end if
      !$acc enter data create(rw(1:nc), pw(1:nc))
      !$acc parallel loop present(rhs, p, rw, pw)
      do idx = 1, nc
        rw(idx) = real(rhs(idx), c_double)
        pw(idx) = real(p(idx), c_double)
        if (.not. (abs(rw(idx)) <= huge(1.0_c_double))) rw(idx) = 0.0_c_double
        if (.not. (abs(pw(idx)) <= huge(1.0_c_double))) pw(idx) = 0.0_c_double
      end do
      if (nd > 0) then
        ! Outlet Dirichlet rhs increments; outlet faces can share an owner, hence atomics.
        !$acc parallel loop copyin(drows(1:nd), dvals(1:nd)) present(rw) private(o)
        do idx = 1, nd
          o = int(drows(idx), ip)
          if (o >= 1 .and. o <= nc) then
            !$acc atomic update
            rw(o) = rw(o) + dvals(idx)
          end if
        end do
      end if
    else
      handle%rhs_work(1:nc) = real(rhs(1:nc), c_double)
      handle%p_work(1:nc) = real(p(1:nc), c_double)

      ! Outlet Dirichlet rhs increments: rhs(o) += c_b * p_val per outlet face.
      if (use_outlet_val .and. handle%n_dirichlet_rows > 0 .and. allocated(handle%dirichlet_rows)) then
        do idx = 1, handle%n_dirichlet_rows
          o = int(handle%dirichlet_rows(idx), ip)
          if (o >= 1 .and. o <= nc) then
            handle%rhs_work(o) = handle%rhs_work(o) + handle%dirichlet_values(idx)
          end if
        end do
      end if
    end if

    empty_rows(1) = int(0, c_int)
    empty_vals(1) = 0.0_c_double

    if (dev) call fvm_petsc_set_io_on_device(1_c_int)
    !$acc host_data use_device(rw, pw) if(dev)
    if (handle%n_dirichlet_rows > 0) then
      status = fvm_petsc_solve_pressure(int(nc, c_int), int(handle%n_cells_global, c_int), int(handle%n_faces, c_int), &
                        handle%local_to_global, handle%face_owner, handle%face_neighbor_local, &
                        handle%face_neighbor_global, handle%skip_face, &
                                        handle%diagonal, rw, handle%face_coeff, &
                        int(handle%reference_cell, c_int), c_max_iter, c_tol, c_atol, &
                        handle%dirichlet_rows, handle%dirichlet_values, int(handle%n_dirichlet_rows, c_int), &
                        pw, c_n_iter, c_residual)
    else
      status = fvm_petsc_solve_pressure(int(nc, c_int), int(handle%n_cells_global, c_int), int(handle%n_faces, c_int), &
                        handle%local_to_global, handle%face_owner, handle%face_neighbor_local, &
                        handle%face_neighbor_global, handle%skip_face, &
                                        handle%diagonal, rw, handle%face_coeff, &
                        int(handle%reference_cell, c_int), c_max_iter, c_tol, c_atol, &
                        empty_rows, empty_vals, int(0, c_int), &
                        pw, c_n_iter, c_residual)
    end if
    !$acc end host_data
    if (dev) call fvm_petsc_set_io_on_device(0_c_int)
    if (status /= 0_c_int) then
      write(*,'(A,I0,A,I0,A,ES12.4)') 'pressure_solver_solve: PETSc solve failed status=', &
           int(status, ip), ' iter=', int(c_n_iter, ip), ' residual=', real(c_residual, rp)
      error stop 'pressure_solver_solve: PETSc solve failed'
    end if

    if (dev) then
      !$acc parallel loop present(p, pw)
      do idx = 1, nc
        p(idx) = real(pw(idx), rp)
      end do
      !$acc exit data delete(rw(1:nc), pw(1:nc))
    else
      p(1:nc) = real(handle%p_work(1:nc), rp)
      if (host_sync) then
        !$acc update device(p)
      end if
    end if

    if (present(n_iter)) n_iter = int(c_n_iter, ip)
    if (present(residual)) residual = real(c_residual, rp)

    ! Warn if residual unexpectedly exceeds atol_in floor (true rtol*||r0|| not directly available here).
    if (present(atol_in)) then
      if (real(c_residual, rp) > 1.0e4_rp * atol_in) &
        call warn_log_bump(WARN_PSOLVER_RESIDUAL, 1_ip, real(c_residual, rp))
    end if

  end subroutine pressure_solver_solve

  subroutine pressure_solver_free(handle)
    type(pressure_solver_handle_t), intent(inout) :: handle

    if (allocated(handle%local_to_global)) deallocate(handle%local_to_global)
    if (allocated(handle%face_owner)) deallocate(handle%face_owner)
    if (allocated(handle%face_neighbor_local)) deallocate(handle%face_neighbor_local)
    if (allocated(handle%face_neighbor_global)) deallocate(handle%face_neighbor_global)
    if (allocated(handle%skip_face)) deallocate(handle%skip_face)
    if (allocated(handle%diagonal)) deallocate(handle%diagonal)
    if (allocated(handle%face_coeff)) deallocate(handle%face_coeff)
    if (allocated(handle%rhs_work)) deallocate(handle%rhs_work)
    if (allocated(handle%p_work)) deallocate(handle%p_work)
    if (allocated(handle%face_mesh_idx)) deallocate(handle%face_mesh_idx)
    if (allocated(handle%face_neighbor_ext)) deallocate(handle%face_neighbor_ext)
    if (allocated(handle%dirichlet_rows)) deallocate(handle%dirichlet_rows)
    if (allocated(handle%dirichlet_values)) deallocate(handle%dirichlet_values)

    if (petsc_runtime_active) then
      if (fvm_petsc_backend_finalize() /= 0_c_int) then
        error stop 'pressure_solver_free: PETSc finalize failed'
      end if
      petsc_runtime_active = .false.
    end if

    handle%n_cells = 0
    handle%n_cells_global = 0
    handle%n_faces = 0
    handle%reference_cell = 1
  end subroutine pressure_solver_free

end module mod_pressure_solver
