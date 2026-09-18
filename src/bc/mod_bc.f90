! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Boundary condition application for velocity, pressure, temperature, and scalars.
module mod_bc
  use mod_kinds
  use mod_parameters
  use mod_scalar
  use mod_mesh
  use mod_abl_profiles
  implicit none
  private

  public :: BCConfig
  public :: apply_velocity_bc
  public :: apply_pressure_bc
  public :: apply_temperature_bc
  public :: apply_passive_scalar_bc

  type, public :: BCConfig
    integer(ip) :: bc_id    = BC_ZEROGRADIENT
    real(rp)    :: U_val(3) = 0.0_rp  ! inlet velocity
    real(rp)    :: T_val    = 0.0_rp  ! inlet/fixed temperature
    real(rp)    :: p_val    = 0.0_rp  ! fixed pressure (outlet)
    real(rp)    :: a_wall(3) = 0.0_rp  ! wall acceleration for Kim-Moin p BC
    integer(ip), allocatable :: scalar_bc_mode(:)
    real(rp), allocatable :: scalar_value(:)
    real(rp), allocatable :: scalar_flux(:)
    real(rp), allocatable :: scalar_flux_omega(:)
    real(rp), allocatable :: scalar_flux_phase(:)
    ! ABL inflow fields
    real(rp)    :: abl_z0_m       = 0.1_rp  ! momentum roughness [m]
    real(rp)    :: abl_z0_h       = 0.01_rp  ! thermal roughness [m]
    real(rp)    :: abl_u_ref      = 8.0_rp  ! reference speed [m/s]
    real(rp)    :: abl_z_ref      = 100.0_rp  ! reference height [m]
    real(rp)    :: abl_T_surface  = 300.0_rp  ! surface temperature [K]
    real(rp)    :: abl_T_ref      = 300.0_rp  ! reference temperature at z_ref [K]
    real(rp)    :: abl_wind_dir   = 270.0_rp  ! wind direction [degrees]
    logical     :: abl_wall       = .false.  ! MOST wall model flag
  end type BCConfig

  ! Periodic-partner cell cache, keyed by boundary-face index; resolved once instead of once per field per RK stage. Rebuilt if n_bnd_faces changes.
  integer(ip), allocatable, save :: periodic_nb_cache(:)

  ! Flattened, allocatable-component-free mirror of bc_configs(:); rebuilt only when n_patches changes, so bc_configs must never be mutated post-init without also changing patch count.
  real(rp),    allocatable, save :: bc_U_val_soa(:,:), bc_a_wall_soa(:,:)  ! (3, n_patches)
  real(rp),    allocatable, save :: bc_T_val_soa(:), bc_p_val_soa(:)
  real(rp),    allocatable, save :: bc_abl_z0_m_soa(:), bc_abl_z0_h_soa(:), bc_abl_u_ref_soa(:), &
                                     bc_abl_z_ref_soa(:), bc_abl_T_surface_soa(:), bc_abl_T_ref_soa(:), &
                                     bc_abl_wind_dir_soa(:)
  integer(ip), allocatable, save :: bc_scalar_mode_soa(:,:)  ! (n_scalars, n_patches)
  real(rp),    allocatable, save :: bc_scalar_value_soa(:,:), bc_scalar_flux_soa(:,:), &
                                     bc_scalar_flux_omega_soa(:,:), bc_scalar_flux_phase_soa(:,:)
  integer(ip), save :: bc_soa_n_patches = -1
  integer(ip), save :: bc_soa_n_scalars = 0

contains

  !> Populate periodic_nb_cache on first use or after a topology change.
  subroutine ensure_periodic_nb_cache(msh)
    type(Mesh), intent(in) :: msh
    integer(ip) :: f_bnd, f_global, o

    if (allocated(periodic_nb_cache)) then
      if (size(periodic_nb_cache) == msh%n_bnd_faces) return
      deallocate(periodic_nb_cache)
    end if

    allocate(periodic_nb_cache(msh%n_bnd_faces))
    do f_bnd = 1, msh%n_bnd_faces
      f_global = msh%n_int_faces + f_bnd
      o = msh%faces(f_global)%owner
      periodic_nb_cache(f_bnd) = resolve_periodic_nb(msh, f_bnd, o)
    end do
  end subroutine ensure_periodic_nb_cache

  !> Flatten bc_configs(:) into plain SoA arrays on first use or after a patch-count change.
  subroutine ensure_bc_config_soa(bc_configs, msh)
    type(BCConfig), intent(in) :: bc_configs(:)
    type(Mesh),     intent(in) :: msh

    integer(ip) :: p, n_patches, n_scalars

    n_patches = msh%n_patches
    if (bc_soa_n_patches == n_patches) return

    if (allocated(bc_U_val_soa)) deallocate(bc_U_val_soa, bc_a_wall_soa, bc_T_val_soa, bc_p_val_soa, &
      bc_abl_z0_m_soa, bc_abl_z0_h_soa, bc_abl_u_ref_soa, bc_abl_z_ref_soa, &
      bc_abl_T_surface_soa, bc_abl_T_ref_soa, bc_abl_wind_dir_soa)
    if (allocated(bc_scalar_mode_soa)) deallocate(bc_scalar_mode_soa, bc_scalar_value_soa, &
      bc_scalar_flux_soa, bc_scalar_flux_omega_soa, bc_scalar_flux_phase_soa)

    allocate(bc_U_val_soa(3, n_patches), bc_a_wall_soa(3, n_patches))
    allocate(bc_T_val_soa(n_patches), bc_p_val_soa(n_patches))
    allocate(bc_abl_z0_m_soa(n_patches), bc_abl_z0_h_soa(n_patches), bc_abl_u_ref_soa(n_patches), &
             bc_abl_z_ref_soa(n_patches), bc_abl_T_surface_soa(n_patches), bc_abl_T_ref_soa(n_patches), &
             bc_abl_wind_dir_soa(n_patches))

    do p = 1, n_patches
      bc_U_val_soa(:, p)      = bc_configs(p)%U_val
      bc_a_wall_soa(:, p)     = bc_configs(p)%a_wall
      bc_T_val_soa(p)         = bc_configs(p)%T_val
      bc_p_val_soa(p)         = bc_configs(p)%p_val
      bc_abl_z0_m_soa(p)      = bc_configs(p)%abl_z0_m
      bc_abl_z0_h_soa(p)      = bc_configs(p)%abl_z0_h
      bc_abl_u_ref_soa(p)     = bc_configs(p)%abl_u_ref
      bc_abl_z_ref_soa(p)     = bc_configs(p)%abl_z_ref
      bc_abl_T_surface_soa(p) = bc_configs(p)%abl_T_surface
      bc_abl_T_ref_soa(p)     = bc_configs(p)%abl_T_ref
      bc_abl_wind_dir_soa(p)  = bc_configs(p)%abl_wind_dir
    end do

    ! build_bc_configs allocates scalar_* uniformly (same length) across every patch, or not at all.
    n_scalars = 0
    if (n_patches >= 1) then
      if (allocated(bc_configs(1)%scalar_bc_mode)) n_scalars = size(bc_configs(1)%scalar_bc_mode)
    end if
    bc_soa_n_scalars = n_scalars
    if (n_scalars > 0) then
      allocate(bc_scalar_mode_soa(n_scalars, n_patches), source=SCALAR_BC_ZERO_GRADIENT)
      allocate(bc_scalar_value_soa(n_scalars, n_patches), source=0.0_rp)
      allocate(bc_scalar_flux_soa(n_scalars, n_patches), source=0.0_rp)
      allocate(bc_scalar_flux_omega_soa(n_scalars, n_patches), source=0.0_rp)
      allocate(bc_scalar_flux_phase_soa(n_scalars, n_patches), source=0.0_rp)
      do p = 1, n_patches
        if (allocated(bc_configs(p)%scalar_bc_mode) .and. size(bc_configs(p)%scalar_bc_mode) /= n_scalars) &
          error stop 'ensure_bc_config_soa: scalar_bc_mode size mismatch across patches'
        if (allocated(bc_configs(p)%scalar_value) .and. size(bc_configs(p)%scalar_value) /= n_scalars) &
          error stop 'ensure_bc_config_soa: scalar_value size mismatch across patches'
        if (allocated(bc_configs(p)%scalar_flux) .and. size(bc_configs(p)%scalar_flux) /= n_scalars) &
          error stop 'ensure_bc_config_soa: scalar_flux size mismatch across patches'
        if (allocated(bc_configs(p)%scalar_flux_omega) .and. size(bc_configs(p)%scalar_flux_omega) /= n_scalars) &
          error stop 'ensure_bc_config_soa: scalar_flux_omega size mismatch across patches'
        if (allocated(bc_configs(p)%scalar_flux_phase) .and. size(bc_configs(p)%scalar_flux_phase) /= n_scalars) &
          error stop 'ensure_bc_config_soa: scalar_flux_phase size mismatch across patches'

        if (allocated(bc_configs(p)%scalar_bc_mode))    bc_scalar_mode_soa(:, p)       = bc_configs(p)%scalar_bc_mode
        if (allocated(bc_configs(p)%scalar_value))      bc_scalar_value_soa(:, p)      = bc_configs(p)%scalar_value
        if (allocated(bc_configs(p)%scalar_flux))       bc_scalar_flux_soa(:, p)       = bc_configs(p)%scalar_flux
        if (allocated(bc_configs(p)%scalar_flux_omega)) bc_scalar_flux_omega_soa(:, p) = bc_configs(p)%scalar_flux_omega
        if (allocated(bc_configs(p)%scalar_flux_phase)) bc_scalar_flux_phase_soa(:, p) = bc_configs(p)%scalar_flux_phase
      end do
    end if

    bc_soa_n_patches = n_patches
  end subroutine ensure_bc_config_soa

  !> Apply velocity BC: sets U_bnd, uses grad_U for O(h^2) extrapolation on zero-gradient faces
  subroutine apply_velocity_bc(U, U_bnd, bc_configs, msh, periodic_only, grad_U, u_fluc_bnd, u_mean_bnd)
    real(rp),      intent(inout) :: U(:,:)  ! (3, n_cells + n_ghost)
    real(rp),      intent(inout) :: U_bnd(:,:)  ! (3, n_bnd_faces)
    type(BCConfig), intent(in)   :: bc_configs(:)  ! per-patch configs
    type(Mesh),    intent(in)    :: msh
    logical, optional, intent(in) :: periodic_only
    real(rp), optional, intent(in) :: grad_U(:,:,:)  ! (3, 3, n_cells+n_ghost)
    real(rp), optional, intent(in) :: u_fluc_bnd(:,:)  ! (3, n_bnd_faces) SEM inflow fluctuation
    real(rp), optional, allocatable, intent(in) :: u_mean_bnd(:,:)  ! (3, n_bnd_faces) prescribed inlet mean, overrides the log law

    integer(ip) :: p, f_bnd, f_global, o, nb, patch_id, bc_id
    real(rp)    :: n_hat(3), U_o(3), U_n(3), d_cf(3)
    real(rp)    :: u_star(msh%n_patches)
    logical     :: have_mean, have_grad, have_fluc, filter_periodic, filter_nonperiodic

    have_mean = .false.
    if (present(u_mean_bnd)) then
      if (allocated(u_mean_bnd)) have_mean = .true.
    end if
    have_grad = present(grad_U)
    have_fluc = present(u_fluc_bnd)

    call ensure_periodic_nb_cache(msh)
    call ensure_bc_config_soa(bc_configs, msh)

    filter_periodic    = .false.
    filter_nonperiodic = .false.
    if (present(periodic_only)) then
      if (periodic_only) then
        filter_periodic = .true.
      else
        filter_nonperiodic = .true.
      end if
    end if

    ! u_star is identical for every face of a patch (same abl_u_ref/z_ref/z0_m); compute once per patch here.
    do p = 1, msh%n_patches
      u_star(p) = 0.0_rp
      if (msh%patch_bc_id_soa(p) == BC_ABL_INLET .and. .not. have_mean) &
        u_star(p) = abl_ustar_from_ref(bc_abl_u_ref_soa(p), bc_abl_z_ref_soa(p), bc_abl_z0_m_soa(p))
    end do

    do f_bnd = 1, msh%n_bnd_faces
      f_global = msh%n_int_faces + f_bnd
      patch_id = msh%face_patch_id_soa(f_global)
      if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
        bc_id = msh%patch_bc_id_soa(patch_id)
      else
        bc_id = BC_ZEROGRADIENT  ! unassigned patch: fall back rather than read out of bounds
      end if
      if (filter_periodic    .and. bc_id /= BC_PERIODIC) cycle
      if (filter_nonperiodic .and. bc_id == BC_PERIODIC) cycle

      o     = msh%face_owner_soa(f_global)
      n_hat = msh%face_normal_soa(:, f_global)

      select case (bc_id)

      case (BC_INLET)
        U_bnd(:, f_bnd) = bc_U_val_soa(:, patch_id)

      case (BC_OUTLET, BC_ZEROGRADIENT)
        ! Linear extrapolation to face centre (O(h^2)); fall back to zeroth-order if
        ! grad_U not provided (e.g. during init before first gradient sweep).
        if (have_grad .and. size(grad_U, 3) >= o) then
          d_cf = msh%face_centre_soa(:, f_global) - msh%cell_centroid_soa(:, o)
          U_bnd(1, f_bnd) = U(1, o) + dot_product(grad_U(:, 1, o), d_cf)
          U_bnd(2, f_bnd) = U(2, o) + dot_product(grad_U(:, 2, o), d_cf)
          U_bnd(3, f_bnd) = U(3, o) + dot_product(grad_U(:, 3, o), d_cf)
        else
          U_bnd(:, f_bnd) = U(:, o)
        end if

      case (BC_WALL_NOSLIP, BC_WALL_MODEL)
        U_bnd(:, f_bnd) = 0.0_rp

      case (BC_WALL_VELOCITY)
        ! Prescribed moving wall velocity with impermeability enforced.
        U_o = bc_U_val_soa(:, patch_id)
        U_n = dot_product(U_o, n_hat) * n_hat
        U_bnd(:, f_bnd) = U_o - U_n

      case (BC_PERIODIC)
        nb = periodic_nb_cache(f_bnd)
        if (nb >= 1 .and. nb <= size(U, 2)) then
          ! Face-centre value: gc-weighted interpolation along the owner-to-neighbour
          ! centroid vector (FPS Sect. 7.2).  Reduces to 0.5 on uniform grids.
          U_bnd(:, f_bnd) = msh%periodic_gc_soa(f_bnd) * U(:, o) &
                          + (1.0_rp - msh%periodic_gc_soa(f_bnd)) * U(:, nb)
        else
          U_bnd(:, f_bnd) = U(:, o)
        end if

      case (BC_SLIP, BC_SYMMETRY)
        ! Remove normal component
        U_o = U(:, o)
        U_n = dot_product(U_o, n_hat) * n_hat
        U_bnd(:, f_bnd) = U_o - U_n

      case (BC_ABL_INLET)
        ! A prescribed profile supersedes the analytic log law when one was loaded.
        if (have_mean) then
          U_bnd(:, f_bnd) = u_mean_bnd(:, f_bnd)
        else
          call abl_inlet_velocity(bc_abl_z0_m_soa(patch_id), bc_abl_wind_dir_soa(patch_id), &
                                  u_star(patch_id), msh%face_centre_soa(:, f_global), U_bnd(:, f_bnd))
        end if
        ! SEM fluctuation is precomputed by the solver and already flux-corrected.
        if (have_fluc) U_bnd(:, f_bnd) = U_bnd(:, f_bnd) + u_fluc_bnd(:, f_bnd)

      end select
    end do
  end subroutine apply_velocity_bc

  !> Apply pressure BC: OUTLET Dirichlet, others Neumann with O(h^2) extrapolation if grad_p provided
  subroutine apply_pressure_bc(p, p_bnd, bc_configs, msh, periodic_only, grad_p)
    real(rp),      intent(inout) :: p(:)  ! (n_cells + n_ghost)
    real(rp),      intent(inout) :: p_bnd(:)  ! (n_bnd_faces)
    type(BCConfig), intent(in)   :: bc_configs(:)
    type(Mesh),    intent(in)    :: msh
    logical, optional, intent(in) :: periodic_only
    real(rp), optional, intent(in) :: grad_p(:,:)  ! (3, n_cells+n_ghost)

    integer(ip) :: f_bnd, f_global, o, nb, patch_id, bc_id
    real(rp) :: an, d_ortho, d_cf(3)
    logical :: have_grad, filter_periodic, filter_nonperiodic

    call ensure_periodic_nb_cache(msh)
    call ensure_bc_config_soa(bc_configs, msh)

    have_grad = present(grad_p)
    filter_periodic    = .false.
    filter_nonperiodic = .false.
    if (present(periodic_only)) then
      if (periodic_only) then
        filter_periodic = .true.
      else
        filter_nonperiodic = .true.
      end if
    end if

    do f_bnd = 1, msh%n_bnd_faces
      f_global = msh%n_int_faces + f_bnd
      patch_id = msh%face_patch_id_soa(f_global)
      if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
        bc_id = msh%patch_bc_id_soa(patch_id)
      else
        bc_id = BC_ZEROGRADIENT  ! unassigned patch: fall back rather than read out of bounds
      end if
      if (filter_periodic    .and. bc_id /= BC_PERIODIC) cycle
      if (filter_nonperiodic .and. bc_id == BC_PERIODIC) cycle

      o = msh%face_owner_soa(f_global)

      select case (bc_id)
      case (BC_OUTLET)
        p_bnd(f_bnd) = bc_p_val_soa(patch_id)  ! fixed (usually 0)
      case (BC_PERIODIC)
        nb = periodic_nb_cache(f_bnd)
        if (nb >= 1 .and. nb <= size(p)) then
          p_bnd(f_bnd) = msh%periodic_gc_soa(f_bnd) * p(o) &
                       + (1.0_rp - msh%periodic_gc_soa(f_bnd)) * p(nb)
        else
          p_bnd(f_bnd) = p(o)
        end if
      case (BC_WALL_NOSLIP, BC_WALL_VELOCITY, BC_WALL_MODEL, BC_SLIP, BC_SYMMETRY)
        ! Neumann wall p BC (Kim & Moin 1985): dp/dn = -rho*(a_wall.n); viscous term omitted (LES).
        an = dot_product(bc_a_wall_soa(:, patch_id), msh%face_normal_soa(:, f_global))
        d_ortho = abs(dot_product(msh%face_centre_soa(:, f_global) - msh%cell_centroid_soa(:, o), &
                                  msh%face_normal_soa(:, f_global)))
        p_bnd(f_bnd) = p(o) - an * d_ortho
      case default
        if (have_grad .and. size(grad_p, 2) >= o) then
          d_cf = msh%face_centre_soa(:, f_global) - msh%cell_centroid_soa(:, o)
          p_bnd(f_bnd) = p(o) + dot_product(grad_p(:, o), d_cf)
        else
          p_bnd(f_bnd) = p(o)  ! zero-gradient fallback
        end if
      end select
    end do
  end subroutine apply_pressure_bc

  !> Apply temperature boundary conditions.
  !> If grad_T is provided, the zero-gradient (default) case uses linear extrapolation.
  subroutine apply_temperature_bc(T, T_bnd, bc_configs, msh, periodic_only, grad_T)
    real(rp),      intent(inout) :: T(:)
    real(rp),      intent(inout) :: T_bnd(:)
    type(BCConfig), intent(in)   :: bc_configs(:)
    type(Mesh),    intent(in)    :: msh
    logical, optional, intent(in) :: periodic_only
    real(rp), optional, intent(in) :: grad_T(:,:)  ! (3, n_cells+n_ghost)

    integer(ip) :: p, f_bnd, f_global, o, nb, patch_id, bc_id
    real(rp)    :: d_cf(3)
    real(rp)    :: T_star(msh%n_patches)
    logical     :: have_grad, filter_periodic, filter_nonperiodic

    call ensure_periodic_nb_cache(msh)
    call ensure_bc_config_soa(bc_configs, msh)

    have_grad = present(grad_T)
    filter_periodic    = .false.
    filter_nonperiodic = .false.
    if (present(periodic_only)) then
      if (periodic_only) then
        filter_periodic = .true.
      else
        filter_nonperiodic = .true.
      end if
    end if

    ! T_star is identical for every face of a patch; compute once per patch here.
    do p = 1, msh%n_patches
      T_star(p) = 0.0_rp
      if (msh%patch_bc_id_soa(p) == BC_ABL_INLET) &
        T_star(p) = abl_tstar_from_ref(bc_abl_T_ref_soa(p), bc_abl_T_surface_soa(p), &
                                       bc_abl_z_ref_soa(p), bc_abl_z0_h_soa(p))
    end do

    do f_bnd = 1, msh%n_bnd_faces
      f_global = msh%n_int_faces + f_bnd
      patch_id = msh%face_patch_id_soa(f_global)
      if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
        bc_id = msh%patch_bc_id_soa(patch_id)
      else
        bc_id = BC_ZEROGRADIENT  ! unassigned patch: fall back rather than read out of bounds
      end if
      if (filter_periodic    .and. bc_id /= BC_PERIODIC) cycle
      if (filter_nonperiodic .and. bc_id == BC_PERIODIC) cycle

      o = msh%face_owner_soa(f_global)

      select case (bc_id)
      case (BC_INLET, BC_FIXEDTEMP)
        T_bnd(f_bnd) = bc_T_val_soa(patch_id)
      case (BC_ABL_INLET)
        ! MOST temperature profile (neutral form for inlet)
        call abl_inlet_temperature(bc_abl_T_surface_soa(patch_id), bc_abl_z0_h_soa(patch_id), &
                                   T_star(patch_id), msh%face_centre_soa(:, f_global), T_bnd(f_bnd))
      case (BC_PERIODIC)
        nb = periodic_nb_cache(f_bnd)
        if (nb >= 1 .and. nb <= size(T)) then
          T_bnd(f_bnd) = msh%periodic_gc_soa(f_bnd) * T(o) &
                       + (1.0_rp - msh%periodic_gc_soa(f_bnd)) * T(nb)
        else
          T_bnd(f_bnd) = T(o)
        end if
      case default
        if (have_grad .and. size(grad_T, 2) >= o) then
          d_cf = msh%face_centre_soa(:, f_global) - msh%cell_centroid_soa(:, o)
          T_bnd(f_bnd) = T(o) + dot_product(grad_T(:, o), d_cf)
        else
          T_bnd(f_bnd) = T(o)  ! zero-gradient fallback
        end if
      end select
    end do
  end subroutine apply_temperature_bc

  !> Apply boundary conditions for one passive scalar index.
  !> Neumann mode expects scalar_flux in units of scalar-flux per area.
  subroutine apply_passive_scalar_bc(C, C_bnd, bc_configs, msh, scalar_id, alpha_mol, time, periodic_only)
    real(rp), intent(inout) :: C(:)
    real(rp), intent(inout) :: C_bnd(:)
    type(BCConfig), intent(in) :: bc_configs(:)
    type(Mesh), intent(in) :: msh
    integer(ip), intent(in) :: scalar_id
    real(rp), intent(in) :: alpha_mol
    real(rp), intent(in) :: time
    logical, optional, intent(in) :: periodic_only

    integer(ip) :: f_bnd, f_global, o, nb, patch_id, bc_id, mode
    real(rp) :: flux_val, alpha_face, d_cf
    logical :: have_scalar, filter_periodic, filter_nonperiodic

    call ensure_periodic_nb_cache(msh)
    call ensure_bc_config_soa(bc_configs, msh)

    have_scalar = (scalar_id >= 1 .and. scalar_id <= bc_soa_n_scalars)
    filter_periodic    = .false.
    filter_nonperiodic = .false.
    if (present(periodic_only)) then
      if (periodic_only) then
        filter_periodic = .true.
      else
        filter_nonperiodic = .true.
      end if
    end if

    do f_bnd = 1, msh%n_bnd_faces
      f_global = msh%n_int_faces + f_bnd
      patch_id = msh%face_patch_id_soa(f_global)
      if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
        bc_id = msh%patch_bc_id_soa(patch_id)
      else
        bc_id = BC_ZEROGRADIENT  ! unassigned patch: fall back rather than read out of bounds
      end if
      if (filter_periodic    .and. bc_id /= BC_PERIODIC) cycle
      if (filter_nonperiodic .and. bc_id == BC_PERIODIC) cycle

      o = msh%face_owner_soa(f_global)
      if (o < 1 .or. o > size(C)) cycle

      mode = SCALAR_BC_ZERO_GRADIENT
      if (have_scalar) mode = bc_scalar_mode_soa(scalar_id, patch_id)

      select case (bc_id)
      case (BC_PERIODIC)
        nb = periodic_nb_cache(f_bnd)
        if (nb >= 1 .and. nb <= size(C)) then
          C_bnd(f_bnd) = msh%periodic_gc_soa(f_bnd) * C(o) &
                       + (1.0_rp - msh%periodic_gc_soa(f_bnd)) * C(nb)
        else
          C_bnd(f_bnd) = C(o)
        end if
      case default
        select case (mode)
        case (SCALAR_BC_DIRICHLET)
          C_bnd(f_bnd) = C(o)
          if (have_scalar) C_bnd(f_bnd) = bc_scalar_value_soa(scalar_id, patch_id)
        case (SCALAR_BC_NEUMANN)
          flux_val = 0.0_rp
          if (have_scalar) then
            flux_val = bc_scalar_flux_soa(scalar_id, patch_id)
            if (abs(bc_scalar_flux_omega_soa(scalar_id, patch_id)) > SMALL) &
              flux_val = flux_val * cos(bc_scalar_flux_omega_soa(scalar_id, patch_id) * time &
                                        + bc_scalar_flux_phase_soa(scalar_id, patch_id))
          end if
          alpha_face = max(alpha_mol, SMALL)
          d_cf = max(msh%face_dcf_soa(f_global), SMALL)
          C_bnd(f_bnd) = C(o) + flux_val * d_cf / alpha_face
        case default
          C_bnd(f_bnd) = C(o)
        end select
      end select
    end do
  end subroutine apply_passive_scalar_bc

  ! Resolve periodic partner owner cell robustly for BC sampling.
  ! Prefer periodic_neighbor in MPI runs; fall back to same-rank periodic_pair.
  integer(ip) function resolve_periodic_nb(msh, f_bnd, owner) result(nb)
    type(Mesh),   intent(in) :: msh
    integer(ip),  intent(in) :: f_bnd, owner

    integer(ip) :: pair_bnd, pair_face

    nb = 0

    if (allocated(msh%periodic_neighbor)) then
      if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) then
        nb = msh%periodic_neighbor(f_bnd)
        if (nb >= 1 .and. nb /= owner) return
      end if
    end if

    if (allocated(msh%periodic_pair)) then
      if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) then
        pair_bnd = msh%periodic_pair(f_bnd)
        if (pair_bnd >= 1 .and. pair_bnd <= msh%n_bnd_faces) then
          pair_face = msh%n_int_faces + pair_bnd
          if (pair_face >= 1 .and. pair_face <= size(msh%faces)) then
            nb = msh%faces(pair_face)%owner
            if (nb >= 1 .and. nb /= owner) return
          end if
        end if
      end if
    end if

    ! For preprocessed distributed meshes the decomposer stores the periodic partner
    ! as a ghost cell directly in faces%neighbour; use it as the authoritative fallback.
    pair_face = msh%n_int_faces + f_bnd
    if (pair_face >= 1 .and. pair_face <= size(msh%faces)) then
      nb = msh%faces(pair_face)%neighbour
      if (nb >= 1 .and. nb /= owner) return
    end if

    nb = 0
  end function resolve_periodic_nb

  !> Neutral MOST velocity profile at an ABL inlet face; SEM fluctuations are added by the solver.
  !> u_star is identical for every face of the patch; the caller computes it once, not per face.
  subroutine abl_inlet_velocity(z0_m, wind_dir, u_star, face_centre, U_bnd)
    real(rp), intent(in)  :: z0_m, wind_dir, u_star
    real(rp), intent(in)  :: face_centre(3)
    real(rp), intent(out) :: U_bnd(3)

    real(rp) :: u_mag

    u_mag  = abl_u_profile(face_centre(3), u_star, z0_m, 0.0_rp)
    U_bnd  = abl_wind_vector(u_mag, wind_dir)
  end subroutine abl_inlet_velocity

  !> Neutral MOST temperature profile at an ABL inlet face.
  !> T_star is identical for every face of the patch; the caller computes it once, not per face.
  subroutine abl_inlet_temperature(T_surface, z0_h, T_star, face_centre, T_bnd)
    real(rp), intent(in)  :: T_surface, z0_h, T_star
    real(rp), intent(in)  :: face_centre(3)
    real(rp), intent(out) :: T_bnd

    T_bnd  = abl_T_profile(face_centre(3), T_surface, T_star, z0_h, 0.0_rp)
  end subroutine abl_inlet_temperature

end module mod_bc
