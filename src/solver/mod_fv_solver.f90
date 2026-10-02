! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Fractional-step finite volume LES solver (Boussinesq incompressible NS).
module mod_fv_solver
  use mod_kinds
  use mod_parameters
  use mod_scalar
  use mod_mesh
  use mod_gradient
  use mod_nonortho
  use mod_advection
  use mod_viscous
  use mod_vreman
  use mod_dynamic_smagorinsky
  use mod_wall_model
  use mod_abl_types
  use mod_abl_profiles
  use mod_sem_inflow
  use mod_inflow_planes
  use mod_inflow_gen, only: inflow_write_grid_file
  use mod_abl_weather
  use mod_abl_wall
  use mod_pressure_solver
  use mod_bc
  use mod_warn_log
  use mod_profiling, only: prof_push, prof_pop
  use iso_c_binding, only: c_int, c_double
  use mpi_f08
  implicit none
  private

  integer(ip), parameter :: PERF_REPORT_INTERVAL = 1000_ip

  public :: FVSolver
  public :: solver_params_t
  public :: snap_dt_to_write_time

  type, public :: solver_params_t
    real(rp)    :: nu           = 1.0e-5_rp
    real(rp)    :: alpha_T      = 1.4e-5_rp
    real(rp)    :: beta         = 3.4e-3_rp
    real(rp)    :: T_ref        = 300.0_rp
    real(rp)    :: g(3)         = [0.0_rp, 0.0_rp, 0.0_rp]
    real(rp)    :: f_const(3)   = [0.0_rp, 0.0_rp, 0.0_rp]
    !> Dynamic body force recomputed each step to hold the volume-averaged velocity at U_bulk_target;
    !> mutually exclusive with f_const (enforced in the input parser).
    logical     :: constant_mass_flux_forcing = .false.
    real(rp)    :: U_bulk_target(3) = [0.0_rp, 0.0_rp, 0.0_rp]
    logical     :: use_rotation = .false.
    real(rp)    :: omega_rot(3) = [0.0_rp, 0.0_rp, 0.0_rp]
    logical     :: use_harmonic_forcing = .false.
    real(rp)    :: forcing_amp   = 0.0_rp
    real(rp)    :: forcing_omega = 0.0_rp
    real(rp)    :: forcing_phase = 0.0_rp
    real(rp)    :: forcing_dir(3) = [1.0_rp, 0.0_rp, 0.0_rp]
    real(rp)    :: Pr_t         = 0.85_rp
    real(rp)    :: dt           = 1.0e-3_rp
    real(rp)    :: t_end        = 100.0_rp
    !> hard cap on number of integrator steps; 0 disables the cap.
    integer(ip) :: nstep_max    = 0_ip
    real(rp)    :: cfl_max      = 0.5_rp
    integer(ip) :: scheme       = SCHEME_KEP
    character(16) :: time_integrator = 'RK3W'  ! 'RK3W' is the only integrator
    logical     :: use_skew_symmetric_advection = .false.  ! Deprecated no-op; use advection_scheme = KEP
    integer(ip) :: write_interval = 100
    !> Physical-time output cadence; active when write_interval < 0 and write_time > 0.
    real(rp)    :: write_time     = 0.0_rp
    real(rp)    :: p_tol        = 1.0e-10_rp
    integer(ip) :: p_max_iter   = 1000
    real(rp)    :: p_target_divU = 1.0e-8_rp  ! parsed but currently unused; atol is computed from a fixed machine-precision floor (see p_atol_factor)
    character(32) :: pressure_solver = 'boomeramg'  ! boomeramg | gamg
    logical     :: boussinesq   = .false.
    logical     :: adapt_dt     = .false.
    integer(ip) :: n_scalars = 0
    integer(ip) :: scalar_scheme = SCHEME_VANLEER_TVD
    integer(ip) :: temperature_scheme = SCHEME_VANLEER_TVD
    character(64), allocatable :: scalar_names(:)
    real(rp), allocatable :: scalar_diffusivity(:)
    real(rp), allocatable :: scalar_pr_t(:)
    integer(ip) :: n_scalar_sources = 0
    type(ScalarSourceDef), allocatable :: scalar_sources(:)
    character(16) :: gradient_scheme = 'LEAST_SQUARES'
    integer(ip)   :: gg_iterations   = 2  ! GG skewness-correction sweeps (GREEN_GAUSS only)
    integer(ip)   :: n_nonortho_corr = 1  ! Jasak non-orthogonal corrector outer iterations
    logical     :: use_rotational_projection = .true.  ! Brown-Cortez-Minion rotational correction (-nudiv u*)
    logical     :: use_sgs = .true.   ! enable SGS model; set false for no-model/ILES runs
    character(16) :: sgs_model = 'VREMAN'  ! VREMAN (default), SMAGORINSKY, or DYNAMIC_SMAG
    !> Wall-model sampling height; <= 0 samples the wall-adjacent cell (Larsson: set h independent of the grid).
    real(rp)    :: wall_model_height = 0.0_rp
    !> Wall-model velocity filter time constant; <= 0 disables filtering. Try ~ delta/u_tau.
    real(rp)    :: wall_model_filter_time = 0.0_rp
    !> Gates verbose runtime diagnostics (projection idempotency, interface KE budget, worst-divergence-cell report).
    logical     :: debug_diagnostics = .false.
    character(len=:), allocatable :: petsc_options(:)
  end type solver_params_t

  !> Field count for solver_params_t, kept in lockstep with copy_solver_params below.
  !> Bump this whenever a field is added to/removed from the type, and add/remove the
  !> matching dst%x = src%x line: copy_solver_params asserts the two stay in sync, since
  !> storage_size() can't be used as a portable completeness check (compiler-dependent
  !> layout: gfortran/nvfortran report different bit sizes for the identical type).
  integer(ip), parameter :: N_SOLVER_PARAMS_FIELDS = 49

  type, public :: FVSolver
    type(Mesh),    pointer :: mesh => null()
    type(solver_params_t)  :: params

    ! -- primary fields
    real(rp), allocatable :: U(:,:)  ! velocity (3, nc+ng)
    real(rp), allocatable :: p(:)  ! kinematic pressure (nc+ng)
    real(rp), allocatable :: p_old(:)  ! (nc) pressure at step start for incremental projection
    real(rp), allocatable :: T(:)  ! temperature (nc+ng)
    real(rp), allocatable :: C(:,:)  ! passive scalars (n_scalars, nc+ng)

    ! -- SGS viscosity
    real(rp), allocatable :: nu_sgs(:)  ! (nc + ng) - ghost halo populated after compute


    ! -- boundary values
    real(rp), allocatable :: U_bnd(:,:)  ! (3, nbnd)
    real(rp), allocatable :: p_bnd(:)  ! (nbnd)
    real(rp), allocatable :: phi_prev(:)      ! (nc+ng) stage-start p^{k-1/2}, fixed across non-ortho iterations for RC
    real(rp), allocatable :: phi_prev_bnd(:)  ! (nbnd) BCs for phi_prev
    real(rp), allocatable :: T_bnd(:)  ! (nbnd)
    real(rp), allocatable :: C_bnd(:,:)  ! (n_scalars, nbnd)
    real(rp), allocatable :: tau_w_bnd(:,:)  ! (3, nbnd) wall shear stress
    real(rp), allocatable :: sem_fluc_bnd(:,:)  ! (3, nbnd) SEM inflow fluctuation; zero when SEM is off
    real(rp), allocatable :: abl_mean_bnd(:,:)  ! (3, nbnd) prescribed inlet mean; only when a profile file is loaded
    real(rp), allocatable :: cell_z(:)  ! (nc) cached cell-centre height, used by nudging
    real(rp), allocatable :: q_w_bnd(:)  ! (nbnd) wall heat flux for ABL
    real(rp), allocatable :: L_bnd(:)  ! (nbnd) per-face Obukhov length for ABL

    ! -- ABL module state
    type(ABLConfig)   :: abl_config
    type(SEMState)    :: sem_state
    type(InflowPlanes) :: inflow_planes   ! database-driven inlet fluctuation (inflow_source = 'planes')
    logical           :: inflow_planes_on = .false.
    type(NudgingData) :: nudging_data

    ! -- gradients (recomputed each step / stage; pressure gradient is built
    ! -- on the fly from the Poisson stencil in correct_velocity_psolver_dual)
    real(rp), allocatable :: grad_U(:,:,:)  ! (3,3,nc+ng) - ghost halo populated after compute
    real(rp), allocatable :: grad_T(:,:)  ! (3,nc+ng)
    real(rp), allocatable :: grad_C(:,:,:)  ! (3,nc+ng,n_scalars)
    real(rp), allocatable :: face_flux(:)  ! (n_faces) working face-normal flux (predictor during assembly,
                                           ! divergence-free after correct_face_fluxes)
    real(rp), allocatable :: U_face(:)  ! (n_faces) divergence-free face flux from last projection; used for advection
    logical     :: u_face_ready = .false.  ! true once U_face holds a real projection flux (avoids an any() scan in compute_dt)
    real(rp), allocatable :: inv_vol(:)  ! (nc) cached 1/volume

    ! -- persistent work buffers (avoid per-step allocations)
    real(rp), allocatable :: w_H_U(:,:), w_H_T(:), w_H_C(:,:), w_rhs_p(:), w_div_diag(:)
    real(rp), allocatable :: w_U0(:,:), w_T0(:)
    real(rp), allocatable :: w_C0(:,:)
    real(rp), allocatable :: w_conv_U(:,:), w_visc_U(:,:), w_conv_T(:), w_visc_T(:)
    real(rp), allocatable :: w_conv_C(:), w_visc_C(:)
    logical, allocatable :: w_flux_built(:)
    real(rp), allocatable :: w_gp(:,:)
    real(rp), allocatable :: w_gp_rc(:,:)   ! grad(p^{k-1/2}) as applied to U; Rhie-Chow add-back reference
    real(rp), allocatable :: w_gp_cell(:,:) ! (3,nc+ng) pressure gradient last applied by the velocity correction
    real(rp), allocatable :: w_p_tot(:), w_p_tot_bnd(:)  ! p^{k-1/2} + current increment, for the lagged cross term
    real(rp), allocatable :: w_gp_prev(:,:) ! scratch for GG pressure gradient sweeps
    !> (3,3,nc) inverse of sum_f |Sf| n_f (x) n_f; reconstructs a cell vector from face-normal components.
    real(rp), allocatable :: recon_Tinv(:,:,:)
    ! -- gradient halo pack buffer: (9+3+3*n_scalars, nc+ng), one alloc per run.
    real(rp), allocatable :: w_pack_grad(:,:)
    ! -- previous-gradient scratch for GG skewness sweeps (GREEN_GAUSS + gg_iterations>1).
    real(rp), allocatable :: w_g_prev_U(:,:,:)
    real(rp), allocatable :: w_g_prev_T(:,:)

    ! -- BC configs per patch
    type(BCConfig), allocatable :: bc_configs(:)
    integer(ip), allocatable :: periodic_patches(:)
    integer(ip), allocatable :: nonperiodic_patches(:)
    integer(ip), allocatable :: wallmodel_patches(:)
    integer(ip), allocatable :: periodic_bfaces(:)
    integer(ip), allocatable :: nonperiodic_bfaces(:)
    integer(ip) :: n_periodic_patches = 0
    integer(ip) :: n_nonperiodic_patches = 0
    integer(ip) :: n_wallmodel_patches = 0
    integer(ip) :: n_periodic_bfaces = 0
    integer(ip) :: n_nonperiodic_bfaces = 0
    logical,  allocatable :: w_outlet_cell(:)   ! (nc) precomputed outlet cell mask
    real(rp), allocatable :: w_outlet_pval(:)   ! (nc) reference pressure per outlet cell
    real(rp), allocatable :: patch_p_val(:)     ! (n_patches) bc_configs(:)%p_val, flat for device kernels
    !> Static projection operator arrays (psolver indices/coefficients, non-orthogonal caches,
    !> recon_Tinv, boundary-face lists) have been staged on the device (stage_static_on_device).
    logical :: static_on_device = .false.
    !> Every field and work array used by the time step lives on the device (to_device): the device
    !> copy is authoritative and the host copy is only refreshed by sync_host.
    logical :: on_device = .false.
    !> Contiguous single-scalar scratch: C(s,:)/C_bnd(s,:) are strided rows, which would force a
    !> host temporary when handed to a device kernel.
    real(rp), allocatable :: w_Cs(:), w_Cs_bnd(:)

    ! -- Pressure solver backend
    type(pressure_solver_handle_t) :: psolver
    real(rp)    :: V_min_global  = 0.0_rp  ! global minimum cell volume
    real(rp)    :: V_total_global = 0.0_rp  ! global domain volume, for bulk-velocity averaging
    real(rp)    :: p_atol_factor = 1.0e-10_rp  ! set to a machine-precision floor in solver_init; atol = factor/dt
    !> Recomputed each step from the current bulk velocity when constant_mass_flux_forcing is on.
    real(rp)    :: mass_flux_force(3) = [0.0_rp, 0.0_rp, 0.0_rp]

    !> .true. when no face anywhere passes the NONORTHO_SIN2_TOL test, so every Jasak correction is provably zero.
    logical     :: skip_nonortho = .false.
    !> .true. when skip_nonortho holds and no patch needs grad(p) for O(h^2) Neumann extrapolation.
    logical     :: skip_grad_p = .false.

    ! -- Non-orthogonal k-vector cache (Jasak min-correction): static mesh geometry, built once in
    ! -- build_nonortho_cache instead of recomputed (sqrt + divides) on every RK stage/corrector.
    real(rp), allocatable :: nonortho_k_int(:,:)       ! (3, psolver%n_faces) interior/periodic faces
    logical,  allocatable :: nonortho_active_int(:)    ! (psolver%n_faces)
    real(rp), allocatable :: nonortho_k_bnd(:,:)       ! (3, n_nonperiodic_bfaces) BC_OUTLET faces
    logical,  allocatable :: nonortho_active_bnd(:)    ! (n_nonperiodic_bfaces)

    ! -- time tracking
    real(rp)    :: time = 0.0_rp
    integer(ip) :: step = 0
    !> Active when temperature feeds back to momentum (Boussinesq only; neutral ABL is isothermal).
    logical     :: solve_temperature = .false.

    ! -- lightweight runtime performance counters (wall time, seconds)
    logical     :: perf_enabled = .true.
    integer(ip) :: perf_steps = 0_ip
    real(rp)    :: perf_t_step = 0.0_rp
    real(rp)    :: perf_t_halo_bc = 0.0_rp
    real(rp)    :: perf_t_rhs = 0.0_rp
    real(rp)    :: perf_t_rhs_sgs = 0.0_rp
    real(rp)    :: perf_t_rhs_grad = 0.0_rp
    real(rp)    :: perf_t_rhs_transport = 0.0_rp
    real(rp)    :: perf_t_rhs_sources = 0.0_rp
    !> Predictor->pre-correction gap (velocity update, BC/halo, pre-correction projection, grad_U refresh); previously untimed.
    real(rp)    :: perf_t_predictor = 0.0_rp
    !> Sub-split of perf_t_predictor: _bc = update+BC/halo, _precorrect = projection, _regrad = the second compute_all_gradients call.
    real(rp)    :: perf_t_predictor_bc = 0.0_rp
    real(rp)    :: perf_t_predictor_precorrect = 0.0_rp
    real(rp)    :: perf_t_predictor_regrad = 0.0_rp
    real(rp)    :: perf_t_pressure = 0.0_rp
    real(rp)    :: perf_t_projection = 0.0_rp
    integer(ip) :: perf_pressure_solves = 0_ip
    integer(ip) :: perf_pressure_iters = 0_ip
    real(rp)    :: perf_pressure_residual_max = 0.0_rp

    ! -- P0 interface diagnostics: KE budget and divergence split by cell class
    logical                  :: diag_interface = .false.
    logical, allocatable     :: cell_is_refine(:)  ! (n_cells) touches a refinement-interface face
    real(rp)                 :: diag_div_refine = 0.0_rp
    real(rp)                 :: diag_div_smooth = 0.0_rp
    integer(ip)              :: diag_nonortho_iters = 0_ip
    integer(ip)              :: diag_last_step = 0_ip
    real(rp)                 :: diag_nonortho_curr = 0.0_rp   ! last divU_norm_curr in the corrector loop
    real(rp)                 :: diag_nonortho_prev = 0.0_rp   ! divU_norm_prev at that exit
    real(rp)                 :: diag_nonortho_atol = 0.0_rp   ! atol_dynamic at that exit

  contains
    procedure :: init         => solver_init
    procedure :: refresh_inflow => solver_refresh_inflow
    procedure :: advance      => solver_advance
    procedure :: project_once => solver_project_once
    procedure :: projection_idempotency => solver_projection_idempotency
    procedure :: divergence_metrics => solver_divergence_metrics
    procedure :: destroy      => solver_destroy
    procedure :: compute_dt   => solver_compute_dt
    procedure :: compute_cfl_rate => solver_compute_cfl_rate
    procedure :: to_device    => solver_to_device
    procedure :: sync_host    => solver_sync_host
  end type FVSolver

  real(rp), allocatable, save :: dt_flux_sum_buf(:)  ! solver_compute_dt scratch, grown once
  real(rp), allocatable, save :: dt_diff_sum_buf(:)  ! solver_compute_dt scratch, grown once
  !$acc declare create(dt_flux_sum_buf, dt_diff_sum_buf)

contains

  !> Explicit component-wise copy of solver_params_t, replacing plain `dst = src`: nvfortran's
  !> compiler-generated intrinsic assignment for this type corrupts the heap (observed as a
  !> stray malloc() abort later in solver_init), likely from the deferred-length allocatable
  !> character array component mixed among the others.
  subroutine copy_solver_params(dst, src)
    type(solver_params_t), intent(inout) :: dst
    type(solver_params_t), intent(in)    :: src
    integer(ip) :: i
    integer(ip) :: n_copied

    n_copied = 0
    dst%nu                          = src%nu
    dst%alpha_T                     = src%alpha_T
    dst%beta                        = src%beta
    dst%T_ref                       = src%T_ref
    dst%g                           = src%g
    dst%f_const                     = src%f_const
    dst%constant_mass_flux_forcing  = src%constant_mass_flux_forcing
    dst%U_bulk_target                = src%U_bulk_target
    dst%use_rotation                = src%use_rotation
    dst%omega_rot                   = src%omega_rot
    dst%use_harmonic_forcing        = src%use_harmonic_forcing
    dst%forcing_amp                 = src%forcing_amp
    dst%forcing_omega               = src%forcing_omega
    dst%forcing_phase               = src%forcing_phase
    dst%forcing_dir                 = src%forcing_dir
    dst%Pr_t                        = src%Pr_t
    dst%dt                          = src%dt
    dst%t_end                       = src%t_end
    dst%nstep_max                   = src%nstep_max
    dst%cfl_max                     = src%cfl_max
    dst%scheme                      = src%scheme
    dst%time_integrator             = src%time_integrator
    dst%use_skew_symmetric_advection = src%use_skew_symmetric_advection
    dst%write_interval              = src%write_interval
    dst%write_time                  = src%write_time
    dst%p_tol                       = src%p_tol
    dst%p_max_iter                  = src%p_max_iter
    dst%p_target_divU               = src%p_target_divU
    dst%pressure_solver             = src%pressure_solver
    dst%boussinesq                  = src%boussinesq
    dst%adapt_dt                    = src%adapt_dt
    dst%n_scalars                   = src%n_scalars
    dst%scalar_scheme               = src%scalar_scheme
    dst%temperature_scheme          = src%temperature_scheme
    dst%debug_diagnostics           = src%debug_diagnostics
    n_copied = n_copied + 35  ! nu .. debug_diagnostics above

    if (allocated(dst%scalar_names)) deallocate(dst%scalar_names)
    if (allocated(src%scalar_names)) then
      allocate(dst%scalar_names(size(src%scalar_names)))
      dst%scalar_names = src%scalar_names
    end if
    n_copied = n_copied + 1

    if (allocated(dst%scalar_diffusivity)) deallocate(dst%scalar_diffusivity)
    if (allocated(src%scalar_diffusivity)) then
      allocate(dst%scalar_diffusivity(size(src%scalar_diffusivity)))
      dst%scalar_diffusivity = src%scalar_diffusivity
    end if
    n_copied = n_copied + 1

    if (allocated(dst%scalar_pr_t)) deallocate(dst%scalar_pr_t)
    if (allocated(src%scalar_pr_t)) then
      allocate(dst%scalar_pr_t(size(src%scalar_pr_t)))
      dst%scalar_pr_t = src%scalar_pr_t
    end if
    n_copied = n_copied + 1

    dst%n_scalar_sources = src%n_scalar_sources
    n_copied = n_copied + 1
    if (allocated(dst%scalar_sources)) deallocate(dst%scalar_sources)
    if (allocated(src%scalar_sources)) then
      allocate(dst%scalar_sources(size(src%scalar_sources)))
      do i = 1, size(src%scalar_sources)
        dst%scalar_sources(i) = src%scalar_sources(i)
      end do
    end if
    n_copied = n_copied + 1

    dst%gradient_scheme             = src%gradient_scheme
    dst%gg_iterations               = src%gg_iterations
    dst%n_nonortho_corr             = src%n_nonortho_corr
    dst%use_rotational_projection   = src%use_rotational_projection
    dst%use_sgs                     = src%use_sgs
    dst%sgs_model                   = src%sgs_model
    dst%wall_model_height           = src%wall_model_height
    dst%wall_model_filter_time      = src%wall_model_filter_time
    n_copied = n_copied + 8  ! gradient_scheme .. wall_model_filter_time above

    if (allocated(dst%petsc_options)) deallocate(dst%petsc_options)
    if (allocated(src%petsc_options)) then
      allocate(character(len=len(src%petsc_options)) :: dst%petsc_options(size(src%petsc_options)))
      dst%petsc_options = src%petsc_options
    end if
    n_copied = n_copied + 1

    ! Catches the common oversight this routine exists to guard against: a field added to
    ! solver_params_t without a matching copy line here (see N_SOLVER_PARAMS_FIELDS comment).
    if (n_copied /= N_SOLVER_PARAMS_FIELDS) &
      error stop 'copy_solver_params: field count mismatch vs N_SOLVER_PARAMS_FIELDS'
  end subroutine copy_solver_params

  ! Initialisation
  subroutine solver_init(self, msh, params, bc_configs, abl_cfg, inflow_cfg)
    class(FVSolver),    intent(inout) :: self
    type(Mesh), target, intent(inout) :: msh
    type(solver_params_t), intent(in) :: params
    type(BCConfig),    intent(in)     :: bc_configs(:)
    type(ABLConfig),   intent(in), optional :: abl_cfg
    type(InflowConfig), intent(in), optional :: inflow_cfg

    integer(ip) :: nc, ng, nbnd, s, c
    logical     :: want_planes
    integer(ip) :: q_out, f_loc_out, f_bnd_out, f_glob_out, ocell_out

    self%mesh    => msh
    call copy_solver_params(self%params, params)
    if (present(abl_cfg)) self%abl_config = abl_cfg
    nc   = msh%n_cells
    ng   = msh%n_ghost
    nbnd = msh%n_bnd_faces

    ! Allocate and initialize primary fields
    allocate(self%U(3, nc+ng),        source=0.0_rp)
    allocate(self%p(nc+ng),           source=0.0_rp)  ! pressure initialized to 0
    allocate(self%p_old(nc),           source=0.0_rp)
    allocate(self%T(nc+ng),           source=params%T_ref)
    if (params%n_scalars > 0) allocate(self%C(params%n_scalars, nc+ng), source=0.0_rp)

    ! Ensure boundary/ghost values are properly initialized
    self%U(:, nc+1:nc+ng) = 0.0_rp  ! ghost cells zero initially
    self%T(nc+1:nc+ng)    = params%T_ref  ! ghost temps at reference
    self%p(nc+1:nc+ng)    = 0.0_rp  ! ghost pressure zero
    allocate(self%nu_sgs(nc+ng),      source=0.0_rp)
    ! Boundary face arrays - explicitly zero
    allocate(self%U_bnd(3, nbnd),      source=0.0_rp)
    allocate(self%p_bnd(nbnd),         source=0.0_rp)
    allocate(self%phi_prev(nc+ng),      source=0.0_rp)
    allocate(self%phi_prev_bnd(nbnd),   source=0.0_rp)
    allocate(self%T_bnd(nbnd),         source=params%T_ref)
    if (params%n_scalars > 0) allocate(self%C_bnd(params%n_scalars, nbnd), source=0.0_rp)
    allocate(self%tau_w_bnd(3, nbnd),  source=0.0_rp)
    allocate(self%sem_fluc_bnd(3, nbnd), source=0.0_rp)
    allocate(self%cell_z(nc), source=0.0_rp)
    do c = 1, nc
      self%cell_z(c) = msh%cells(c)%centroid(3)
    end do
    allocate(self%grad_U(3,3,nc+ng),  source=0.0_rp)
    allocate(self%grad_T(3,nc+ng),    source=0.0_rp)
    if (params%n_scalars > 0) allocate(self%grad_C(3,nc+ng,params%n_scalars), source=0.0_rp)
    allocate(self%face_flux(msh%n_faces), source=0.0_rp)
    allocate(self%U_face(msh%n_faces),    source=0.0_rp)
    allocate(self%inv_vol(nc), source=0.0_rp)
    do c = 1, nc
      self%inv_vol(c) = 1.0_rp / max(self%mesh%cells(c)%volume, SMALL)
    end do
    ! Persistent device residency: inv_vol never changes after setup, U_face is refreshed via
    ! !$acc update device at each write site instead of re-copyin on every advection_rhs_* call.
    ! Both are rank-local with no ghost/halo concept (inv_vol sized nc, U_face sized n_faces,
    ! neither ever exchanged across ranks), so this is safe for any nranks, not just nranks<=1.
#ifdef _OPENACC
    !$acc enter data copyin(self%inv_vol, self%U_face)
#endif

    allocate(self%w_H_U(3, nc), self%w_H_T(nc), self%w_rhs_p(nc), self%w_div_diag(nc), source=0.0_rp)
    ! RK3W stage storage: w_U0/w_T0/w_C0 hold previous-stage H buffers.
    allocate(self%w_U0(3, nc), self%w_T0(nc), source=0.0_rp)
    ! w_C0/w_H_C keep max(1,...) so the advance_rk3w associate block can always reference them.
    allocate(self%w_C0(max(1_ip, params%n_scalars), nc), source=0.0_rp)
    allocate(self%w_conv_U(3, nc), self%w_visc_U(3, nc), self%w_conv_T(nc), self%w_visc_T(nc), source=0.0_rp)
    allocate(self%w_H_C(max(1_ip, params%n_scalars), nc), source=0.0_rp)
    if (params%n_scalars > 0) &
      allocate(self%w_conv_C(nc), self%w_visc_C(nc), self%w_Cs(nc+ng), self%w_Cs_bnd(nbnd), source=0.0_rp)
    allocate(self%w_flux_built(msh%n_faces))
    self%w_flux_built = .false.
    allocate(self%w_gp(3, nc+ng), source=0.0_rp)
    if (msh%nranks > 1) then
      ! F2: pre-allocate the gradient halo pack buffer once. 9 components
      ! for grad_U (3x3) + 3 for grad_T + 3*n_scalars for grad_C.
      allocate(self%w_pack_grad(9 + 3 + 3*params%n_scalars, nc+ng))
    end if
    if (trim(params%gradient_scheme) == 'GREEN_GAUSS' .and. params%gg_iterations > 1) then
      ! F3: persistent previous-gradient scratch for GG skewness sweeps.
      allocate(self%w_g_prev_U(3, 3, nc+ng))
      allocate(self%w_g_prev_T(3, nc+ng))
      allocate(self%w_gp_prev(3, nc+ng), source=0.0_rp)
    end if
    ! Always allocate w_gp_rc for frozen pressure gradient (used in Rhie-Chow for all cases)
    allocate(self%w_gp_rc(3, nc+ng), source=0.0_rp)
    allocate(self%w_gp_cell(3, nc+ng), source=0.0_rp)
    allocate(self%w_p_tot(nc+ng), source=0.0_rp)
    allocate(self%w_p_tot_bnd(nbnd), source=0.0_rp)
    ! Identity default so the velocity correction stays defined when no Poisson solve runs.
    allocate(self%recon_Tinv(3, 3, nc), source=0.0_rp)
    do q_out = 1, 3
      self%recon_Tinv(q_out, q_out, :) = 1.0_rp
    end do
    allocate(self%bc_configs(msh%n_patches))
    self%bc_configs = bc_configs
    allocate(self%patch_p_val(max(1_ip, msh%n_patches)), source=0.0_rp)
    do q_out = 1, msh%n_patches
      self%patch_p_val(q_out) = bc_configs(q_out)%p_val
    end do
    call build_patch_classification(self)

    ! Precompute outlet cell mask once so anchor_pressure_to_outlet avoids per-call allocation.
    allocate(self%w_outlet_cell(nc), source=.false.)
    allocate(self%w_outlet_pval(nc), source=0.0_rp)
    do q_out = 1, msh%n_patches
      if (bc_configs(q_out)%bc_id /= BC_OUTLET) cycle
      associate(patch_out => msh%patches(q_out))
        do f_loc_out = 1, patch_out%n_faces
          f_bnd_out  = patch_out%start_face + f_loc_out - 1
          f_glob_out = msh%n_int_faces + f_bnd_out
          ocell_out  = msh%faces(f_glob_out)%owner
          if (ocell_out >= 1 .and. ocell_out <= nc) then
            if (.not. self%w_outlet_cell(ocell_out)) then
              self%w_outlet_cell(ocell_out) = .true.
              self%w_outlet_pval(ocell_out) = bc_configs(q_out)%p_val
            end if
          end if
        end do
      end associate
    end do

    ! Initialize boundary face values after BC configurations are available.
    call apply_all_bc(self)

    if (self%params%n_scalar_sources > 0) then
      do s = 1, self%params%n_scalar_sources
        call scalar_source_prepare(self%params%scalar_sources(s))
      end do
    end if

    call initialize_velocity_from_bc(self)

    ! Only initialize pressure solver if needed (p_max_iter > 0)
    if (self%params%p_max_iter > 0) then
      call mesh_assign_compact_ids(msh)  ! Assign compact IDs for PETSc matrix assembly
      call pressure_solver_init(self%psolver, msh, bc_configs, &
                                solver_type=self%params%pressure_solver, &
                                petsc_options=self%params%petsc_options)
      ! atol = p_atol_factor/dt; use V_mean (not V_min) to avoid unachievable tolerances on graded meshes.
      block
        real(rp) :: vmin_local, vsum_local, vsum_global, vmean_global
        integer  :: ierr
        if (nc > 0) then
          vmin_local = minval(msh%cells(1:nc)%volume)
          vsum_local = sum(msh%cells(1:nc)%volume)
        else
          vmin_local = huge(1.0_rp)
          vsum_local = 0.0_rp
        end if
        call MPI_Allreduce(vmin_local, self%V_min_global, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MIN, msh%comm, ierr)
        call MPI_Allreduce(vsum_local, vsum_global, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, msh%comm, ierr)
        self%V_total_global = vsum_global
        if (msh%n_cells_global > 0) then
          vmean_global = vsum_global / real(msh%n_cells_global, rp)
        else
          vmean_global = self%V_min_global
        end if
        ! Pressure solver atol = max(eps_div_rel * V_min / dt, 1e-15) to target velocity-relative divergence
        self%p_atol_factor = 1.0e-15_rp  ! Machine precision floor
        if (msh%rank == 0) then
          write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
            ' [psolver] V_min_global=', self%V_min_global, &
            '  V_mean_global=', vmean_global, &
            '  divergence_tolerance (machine-precision, velocity-relative)=', self%p_atol_factor
        end if
      end block
      call detect_nonortho_skip(self)
      call build_nonortho_cache(self)
      call build_reconstruct_operator(self)
      call stage_static_on_device(self)
    end if

    ! ABL module initialization
    ! SEM inflow, MOST wall treatment and nudging are independent; enable any combination.
    if (self%abl_config%wall_model_id == WM_MOST) then
      allocate(self%q_w_bnd(nbnd), source=0.0_rp)
      ! Zero means "unseeded"; the wall model falls back to abl_config%L_obukhov on the first call.
      allocate(self%L_bnd(nbnd), source=0.0_rp)
    end if

    want_planes = .false.
    if (present(inflow_cfg)) then
      if (len_trim(inflow_cfg%export_grid) > 0) call export_inflow_grid(self, inflow_cfg%export_grid)
      want_planes = (trim(inflow_cfg%source) == 'planes')
    end if

    if (want_planes) then
      call setup_inflow_planes(self, inflow_cfg)
      call update_sem_inflow(self, self%time)
    else if (self%abl_config%use_sem) then
      ! sem_init derives the inlet frame from the ABL inlet patch and reduces it across ranks.
      call sem_init(self%sem_state, self%abl_config, self%mesh)
      ! A prescribed profile or power-law mean also supplies the inlet mean, not just the stresses.
      if (self%sem_state%override_mean) allocate(self%abl_mean_bnd(3, nbnd), source=0.0_rp)
      call update_sem_inflow(self, self%time)
    end if

    if (self%abl_config%use_nudging) then
      call nudging_load(self%nudging_data, self%abl_config%nudging_file)
      if (.not. self%nudging_data%loaded) &
        error stop 'solver_init: use_nudging is set but the nudging file could not be loaded'
    end if

    call report_wall_model_setup(self)

    ! Temperature transport only runs under Boussinesq; a neutral ABL needs no stratification.
    self%solve_temperature = self%params%boussinesq

    call build_interface_diagnostics(self)

  end subroutine solver_init

  !> Measure ||P(P(u)) - P(u)||: a consistent projection is idempotent, so this must be solver noise.
  subroutine solver_projection_idempotency(self, dt_proj)
    class(FVSolver), intent(inout) :: self
    real(rp), optional, intent(in) :: dt_proj

    integer(ip) :: nc, ng, c
    integer     :: ierr
    real(rp)    :: loc(6), glb(6), d2, rel_l2, dt_use
    real(rp), allocatable :: U_save(:,:), p_save(:), pb_save(:), Ub_save(:,:)
    real(rp), allocatable :: uf_save(:), ff_save(:), pold_save(:)
    real(rp), allocatable :: pp_save(:), ppb_save(:), gp_save(:,:), gprc_save(:,:)
    real(rp), allocatable :: U1(:,:)

    nc = self%mesh%n_cells
    ng = self%mesh%n_ghost

    ! Measure at the operating dt: the lagged cross-term gradient carries the previous
    ! dt's pressure scaling, so mixing dt between calls corrupts the comparison.
    dt_use = self%params%dt
    if (present(dt_proj)) then
      dt_use = dt_proj
    else if (self%params%adapt_dt) then
      dt_use = self%compute_dt()
    end if

    ! Snapshot every field project_once writes, so the run continues from an untouched state.
    U_save    = self%U
    p_save    = self%p
    pb_save   = self%p_bnd
    Ub_save   = self%U_bnd
    uf_save   = self%U_face
    ff_save   = self%face_flux
    pold_save = self%p_old
    pp_save   = self%phi_prev
    ppb_save  = self%phi_prev_bnd
    gp_save   = self%w_gp
    gprc_save = self%w_gp_rc

    ! Warm-up pass so p and w_gp are both consistent with dt_use before comparing.
    call solver_project_once(self, dt_use)
    call solver_project_once(self, dt_use)
    U1 = self%U
    call solver_project_once(self, dt_use)

    loc = 0.0_rp
    do c = 1, nc
      d2 = sum((self%U(:, c) - U1(:, c))**2)
      if (self%cell_is_refine(c)) then
        loc(1) = max(loc(1), sqrt(d2))
      else
        loc(2) = max(loc(2), sqrt(d2))
      end if
      loc(5) = max(loc(5), sqrt(sum(U1(:, c)**2)))
      loc(6) = max(loc(6), sqrt(sum(self%U(:, c)**2)))
      loc(3) = loc(3) + self%mesh%cells(c)%volume * d2
      loc(4) = loc(4) + self%mesh%cells(c)%volume * sum(U1(:, c)**2)
    end do

    call MPI_Allreduce(loc(1:2), glb(1:2), 2, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
    call MPI_Allreduce(loc(3:4), glb(3:4), 2, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)
    call MPI_Allreduce(loc(5:6), glb(5:6), 2, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)

    rel_l2 = 0.0_rp
    if (glb(4) > SMALL) rel_l2 = sqrt(glb(3) / glb(4))

    if (self%mesh%rank == 0) then
      write(*,'(A)') '--- Projection Idempotency ---'
      write(*,'(A,1X,ES12.4)') 'dt used                ', dt_use
      write(*,'(A,1X,ES12.4)') 'max |Pu|               ', glb(5)
      write(*,'(A,1X,ES12.4)') 'max |P(Pu)|            ', glb(6)
      write(*,'(A,1X,ES12.4)') 'max |P(Pu)-Pu|  refine ', glb(1)
      write(*,'(A,1X,ES12.4)') 'max |P(Pu)-Pu|  smooth ', glb(2)
      write(*,'(A,1X,ES12.4)') 'relative L2 change     ', rel_l2
      if (rel_l2 > 1.0e-6_rp) &
        write(*,'(A)') ' [projection] NOT idempotent: cell correction is inconsistent with the face flux'
      write(*,'(A)') '--- End Projection Idempotency ---'
    end if

    self%U            = U_save
    self%p            = p_save
    self%p_bnd        = pb_save
    self%U_bnd        = Ub_save
    self%U_face       = uf_save
    self%face_flux    = ff_save
    self%p_old        = pold_save
    self%phi_prev     = pp_save
    self%phi_prev_bnd = ppb_save
    self%w_gp         = gp_save
    self%w_gp_rc      = gprc_save
  end subroutine solver_projection_idempotency

  !> Record post-projection max |div u| by cell class.
  subroutine record_divergence_split(self, jno)
    class(FVSolver), intent(inout) :: self
    integer(ip),     intent(in)    :: jno

    integer(ip) :: c, nc
    real(rp)    :: loc(2), glb(2)
    integer     :: ierr

    nc  = self%mesh%n_cells
    ! Divergence of the corrected face fluxes: what the projection actually achieved. Uses its
    ! own buffer (not w_rhs_p) since the caller still needs w_rhs_p = div(u*) for the rotational
    ! projection correction after this diagnostic runs.
    call face_flux_divergence(self, self%w_div_diag(1:nc))

    loc = 0.0_rp
    do c = 1, nc
      if (self%cell_is_refine(c)) then
        loc(1) = max(loc(1), abs(self%w_div_diag(c)))
      else
        loc(2) = max(loc(2), abs(self%w_div_diag(c)))
      end if
    end do

    call MPI_Allreduce(loc, glb, 2, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
    self%diag_div_refine     = glb(1)
    self%diag_div_smooth     = glb(2)
    self%diag_nonortho_iters = jno
  end subroutine record_divergence_split

  !> Mark cells touching a refinement-interface face and enable the KE budget when any exist.
  subroutine build_interface_diagnostics(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: f, o, nb, nc, n_local, n_global, ierr

    nc = self%mesh%n_cells
    allocate(self%cell_is_refine(nc), source=.false.)

    do f = 1, self%mesh%n_int_faces
      if (.not. self%mesh%faces(f)%is_refine) cycle
      o  = self%mesh%faces(f)%owner
      nb = self%mesh%faces(f)%neighbour
      if (o  >= 1 .and. o  <= nc) self%cell_is_refine(o)  = .true.
      if (nb >= 1 .and. nb <= nc) self%cell_is_refine(nb) = .true.
    end do

    n_local = count(self%cell_is_refine)
    call MPI_Allreduce(n_local, n_global, 1, MPI_INTEGER, MPI_SUM, self%mesh%comm, ierr)

    self%diag_interface = (n_global > 0)
    ! The KE budget accumulates in the host array-of-structs advection loops only. Under OpenACC the
    ! fields live on the device, so enabling it would route KEP/CDS through a host path that reads
    ! stale data; the budget is a debug diagnostic, so it stays off in device builds.
#ifdef _OPENACC
    call advection_set_ke_budget(.false.)
#else
    call advection_set_ke_budget(self%diag_interface)
#endif

    if (self%mesh%rank == 0 .and. self%diag_interface) &
      write(*,'(A,I0,A)') ' [interface-diag] ', n_global, ' cells adjacent to a refinement interface'

    ! Green-Gauss loses linear-exactness where one cell side is split into several faces.
    if (self%diag_interface .and. trim(self%params%gradient_scheme) == 'GREEN_GAUSS') &
      call warn_log_message('[gradient]  ', 'GREEN_GAUSS is inconsistent at refinement interfaces; use least-squares')
  end subroutine build_interface_diagnostics

  !> One-time report and sanity check of the wall-model setup; rank 0 prints.
  subroutine report_wall_model_setup(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: p, f, f_bnd, f_global, o, n_wall_local, ierr
    real(rp)    :: y_w, n_hat(3)
    real(rp)    :: ymin_l, ymax_l, buf(4), gbuf(4)
    real(rp)    :: z0_min, ratio, z0_p, ratio_min_l
    character(len=16) :: wm_name

    n_wall_local = 0
    ymin_l =  huge(1.0_rp)
    ymax_l = -huge(1.0_rp)
    z0_min =  huge(1.0_rp)
    ratio_min_l = huge(1.0_rp)

    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_WALL_MODEL) cycle
      ! Compare y_w against this patch's own z0: a global min/min pairing can hide a bad patch.
      z0_p = max(self%bc_configs(p)%abl_z0_m, SMALL)
      if (z0_p <= 0.0_rp) z0_p = max(self%abl_config%z0_m, SMALL)
      z0_min = min(z0_min, z0_p)
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        o        = self%mesh%faces(f_global)%owner
        if (o < 1) cycle
        n_hat = self%mesh%faces(f_global)%normal
        y_w = abs(dot_product(self%mesh%faces(f_global)%centre - self%mesh%cells(o)%centroid, n_hat))
        n_wall_local = n_wall_local + 1
        ymin_l = min(ymin_l, y_w)
        ymax_l = max(ymax_l, y_w)
        ratio_min_l = min(ratio_min_l, y_w / z0_p)
      end do
    end do

    ! Reduce as maxima of negated minima so one call covers both ends.
    buf(1) = -huge(1.0_rp)
    buf(2) = merge(-ymin_l, -huge(1.0_rp), n_wall_local > 0)
    buf(3) = merge(ymax_l,  -huge(1.0_rp), n_wall_local > 0)
    buf(4) = merge(-ratio_min_l, -huge(1.0_rp), n_wall_local > 0)
    call MPI_Allreduce(buf, gbuf, 4, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
    buf(1) = real(n_wall_local, rp)
    call MPI_Allreduce(buf(1), gbuf(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)
    ratio = -gbuf(4)

    ! Every rank evaluates the fatal check so the stop is collective, not a deadlock.
    if (nint(gbuf(1)) > 0 .and. self%abl_config%wall_model_id == WM_MOST) then
      if (ratio <= 1.0_rp) then
        if (self%mesh%rank == 0) then
          write(*,'(A,F0.3)') ' [wallmodel] ERROR: min y_w/z0_m = ', ratio
          write(*,'(A)') ' [wallmodel] The first cell lies inside the roughness layer.'
          write(*,'(A)') ' [wallmodel] Lower z0_m/patch_z0_m, or set wall_model = loglaw for a smooth wall.'
        end if
        error stop 'wall model: first cell inside roughness layer'
      end if
    end if

    if (self%mesh%rank /= 0) return

    select case (self%abl_config%wall_model_id)
    case (WM_MOST);   wm_name = 'MOST'
    case (WM_NONE);   wm_name = 'none'
    case default;     wm_name = 'loglaw'
    end select

    if (nint(gbuf(1)) == 0) then
      if (self%abl_config%wall_model_id /= WM_LOGLAW) &
        write(*,'(A)') ' [wallmodel] wall_model set but no wallModel patches exist'
      return
    end if

    write(*,'(A,A,A,I0,A)') ' [wallmodel] ', trim(wm_name), ' on ', nint(gbuf(1)), ' faces'
    write(*,'(A,ES10.3,A,ES10.3)') ' [wallmodel] sampling height y_w: ', -gbuf(2), ' to ', gbuf(3)

    if (self%abl_config%wall_model_id == WM_NONE) then
      write(*,'(A)') ' [wallmodel] WARNING: walls carry no stress with wall_model = none'
      return
    end if

    if (self%abl_config%wall_model_id == WM_MOST) then
      write(*,'(A,ES10.3,A,F0.2)') ' [wallmodel] smallest z0_m = ', z0_min, ', min y_w/z0_m = ', ratio
      if (ratio < 3.0_rp) then
        write(*,'(A)') ' [wallmodel] WARNING: y_w/z0_m < 3, the log law is barely valid here'
      end if
    else
      ! Smooth-wall log law is meaningless inside the viscous sublayer.
      if (self%params%nu > SMALL) &
        write(*,'(A,A)') ' [wallmodel] check the first cell sits in the log layer', &
                         ' (y+ >~ 30); see the y_w range above'
    end if
  end subroutine report_wall_model_setup

  subroutine initialize_velocity_from_bc(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: p
    real(rp) :: u0(3)

    u0 = 0.0_rp
    do p = 1, size(self%bc_configs)
      if (self%bc_configs(p)%bc_id == BC_INLET) then
        u0 = self%bc_configs(p)%U_val
        exit
      end if
    end do

    self%U(1, :) = u0(1)
    self%U(2, :) = u0(2)
    self%U(3, :) = u0(3)

    if (self%mesh%rank == 0 .and. any(abs(u0) > 0.0_rp)) then
      write(*, '(A,3F8.4)') ' IC: initialised velocity field to inlet U = ', u0
    end if
  end subroutine initialize_velocity_from_bc

  ! One-time pressure projection for externally seeded initial conditions.
  subroutine solver_project_once(self, dt_proj)
    class(FVSolver), intent(inout) :: self
    real(rp), optional, intent(in) :: dt_proj

    integer(ip) :: nc, c
    real(rp)    :: dt_use

    nc = self%mesh%n_cells
    dt_use = self%params%dt
    if (present(dt_proj)) dt_use = dt_proj
    if (dt_use <= 0.0_rp) dt_use = self%params%dt

    self%w_rhs_p(1:nc) = 0.0_rp

    ! No pre-correction happens here, so Rhie-Chow has nothing to add back.
    self%phi_prev     = 0.0_rp
    self%phi_prev_bnd = 0.0_rp
    self%w_gp_rc      = 0.0_rp

    call refresh_halos_and_bc(self)

    call compute_face_fluxes_and_rhs(self, dt_use, self%w_rhs_p(1:nc))

    ! Only solve pressure if solver was initialized (p_max_iter > 0)
    if (self%params%p_max_iter > 0) then
      block
        real(rp) :: u_bulk_local(2), u_bulk_global(2), u_scale, atol_dynamic
        integer :: ierr
        ! Velocity-relative divergence tolerance: |div(U)| ~ eps_rel * |U|.
        ! Volume-weighted mean speed, not the domain max: see the matching block above.
        u_bulk_local(1) = sum(sqrt(self%U(1,1:nc)**2 + self%U(2,1:nc)**2 + self%U(3,1:nc)**2) * &
                               self%mesh%cells(1:nc)%volume)
        u_bulk_local(2) = sum(self%mesh%cells(1:nc)%volume)
        call MPI_Allreduce(u_bulk_local, u_bulk_global, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, self%mesh%comm, ierr)
        u_scale = u_bulk_global(1) / max(u_bulk_global(2), SMALL)
        atol_dynamic = max(1.0e-12_rp * max(u_scale, 1.0_rp) * self%V_min_global / dt_use, &
                           self%p_atol_factor)
        call pressure_solver_solve(self%psolver, self%w_rhs_p(1:nc), self%p(1:nc), &
                                   tol=self%params%p_tol, &
                                   max_iter=self%params%p_max_iter, &
                                   atol_in=atol_dynamic)
      end block
    end if

    call anchor_pressure_to_outlet(self%p(1:nc), self%mesh, self%w_outlet_cell, self%w_outlet_pval, self%on_device)
    call self%mesh%exchange_halos_scalar(self%p)
    call correct_face_fluxes(self, dt_use)  ! uses phi -> divergence-free
    self%U_face(:) = self%face_flux(:)  ! persist divergence-free flux for advection
    self%u_face_ready = .true.
#ifdef _OPENACC
    !$acc update device(self%U_face)
#endif
    call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh)
    call correct_velocity_psolver_dual(self, dt_use, with_cross=.true.)  ! uses phi
    call apply_velocity_bc(self%U, self%U_bnd, &
                           self%bc_configs, self%mesh, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)

    ! -- Rotational correction: p = phi - (nu+nu_sgs).div u* (Brown, Cortez & Minion 2001)
    if (self%params%use_rotational_projection) then
      !$omp parallel do schedule(static)
      do c = 1, nc
        self%p(c) = self%p(c) + (self%params%nu + self%nu_sgs(c)) * self%w_rhs_p(c) * dt_use * self%inv_vol(c)
      end do
      !$omp end parallel do
      call self%mesh%exchange_halos_scalar(self%p)
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh)
    end if

    self%p_old(1:nc) = self%p(1:nc)
    self%phi_prev(1:nc+self%mesh%n_ghost) = self%p(1:nc+self%mesh%n_ghost)
    self%phi_prev_bnd(:) = self%p_bnd(:)
    call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh)
    ! Seed w_gp for the first stage's pre-correction.
    call compute_pressure_gradient(self)

  end subroutine solver_project_once

  ! Dispatch time integration
  subroutine solver_advance(self, dt_in)
    class(FVSolver), intent(inout) :: self
    real(rp), optional, intent(in) :: dt_in

    real(rp) :: dt, dt_safe, t_step0
    dt = self%params%dt
    if (present(dt_in)) dt = dt_in
    t_step0 = wall_time_now()

    ! CFL guard: for fixed-dt runs only; adapt_dt handles this via global reduction.
    ! RK schemes are not TVD/SSP; unchecked CFL will drive explicit terms unstable.
    if (.not. self%params%adapt_dt) then
      dt_safe = self%compute_dt()
      if (dt > dt_safe) then
        call warn_log_bump(WARN_CFL_CLAMP, 1_ip, dt / max(dt_safe, tiny(1.0_rp)))
        dt = dt_safe
      end if
    end if

    if (self%params%constant_mass_flux_forcing) call update_mass_flux_forcing(self, dt)

    call prof_push('step')
    call advance_rk3w(self, dt)
    call prof_pop()

    ! U/p ghosts are already current from the RK3 stage-3 tail; only T/C need a post-RK refresh.
    call apply_scalar_bc(self)

    self%time = self%time + dt
    self%step = self%step + 1

    if (self%perf_enabled) then
      self%perf_t_step = self%perf_t_step + (wall_time_now() - t_step0)
      self%perf_steps = self%perf_steps + 1_ip
      ! emit periodic performance summary (self-gates on PERF_REPORT_INTERVAL).
      call maybe_report_perf(self, .false.)
    end if
  end subroutine solver_advance

  !> Constant mass flux forcing: recompute the uniform body force each step so the volume-averaged
  !> velocity tracks U_bulk_target. Uniform forces survive the pressure projection on a periodic
  !> domain (grad of a periodic potential has zero mean), so f = (U_target - <U>) / dt drives the
  !> bulk velocity to target over successive steps.
  subroutine update_mass_flux_forcing(self, dt)
    class(FVSolver), intent(inout), target :: self
    real(rp), intent(in) :: dt

    integer(ip) :: nc, c, ierr
    real(rp)    :: mom_local(3), mom_global(3), U_bulk(3), m1, m2, m3
    real(rp), pointer :: U(:,:), cvol(:)

    nc = self%mesh%n_cells
    U => self%U
    cvol => self%mesh%cell_volume_soa
    m1 = 0.0_rp;  m2 = 0.0_rp;  m3 = 0.0_rp
    !$acc parallel loop if(self%on_device) default(present) reduction(+:m1, m2, m3)
    do c = 1, nc
      m1 = m1 + U(1, c) * cvol(c)
      m2 = m2 + U(2, c) * cvol(c)
      m3 = m3 + U(3, c) * cvol(c)
    end do
    mom_local = [m1, m2, m3]
    call MPI_Allreduce(mom_local, mom_global, 3, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)

    U_bulk = mom_global / max(self%V_total_global, SMALL)
    self%mass_flux_force = (self%params%U_bulk_target - U_bulk) / dt
  end subroutine update_mass_flux_forcing

  subroutine refresh_halos_and_bc(self)
    class(FVSolver), intent(inout) :: self

    if (self%params%n_scalars > 0) then
      call self%mesh%begin_halos_state(self%U, self%p, self%T, self%C)
    else
      call self%mesh%begin_halos_state(self%U, self%p, self%T)
    end if
    call apply_all_bc_subset(self, .false.)
    if (self%params%n_scalars > 0) then
      call self%mesh%end_halos_state(self%U, self%p, self%T, self%C)
    else
      call self%mesh%end_halos_state(self%U, self%p, self%T)
    end if
    call apply_all_bc_subset(self, .true.)
  end subroutine refresh_halos_and_bc

  !> Compute and halo-exchange velocity, temperature, and scalar cell gradients for explicit RHS.
  subroutine compute_all_gradients(self)
    class(FVSolver), intent(inout), target :: self
    integer(ip) :: s, iter, c, i, j, nc, ntot
    logical     :: gg
    real(rp), pointer :: gU(:,:,:), gT(:,:), gC(:,:,:), gpU(:,:,:), gpT(:,:)

    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    gg = (trim(self%params%gradient_scheme) == 'GREEN_GAUSS')
    gU => self%grad_U
    gT => self%grad_T

    if (gg) then
      call green_gauss_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd)
      if (self%solve_temperature) &
        call green_gauss_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd)
    else
      call ls_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd)
      if (self%solve_temperature) &
        call ls_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd)
    end if
    do s = 1, self%params%n_scalars
      call load_scalar_row(self, s)
      if (gg) then
        call green_gauss_gradient_scalar(self%w_Cs, self%mesh, self%grad_C(:, :, s), self%w_Cs_bnd)
      else
        call ls_gradient_scalar(self%w_Cs, self%mesh, self%grad_C(:, :, s), self%w_Cs_bnd)
      end if
    end do

    if (gg .and. self%params%gg_iterations > 1) then
      gpU => self%w_g_prev_U
      gpT => self%w_g_prev_T
      do iter = 2, self%params%gg_iterations
        ! Exchange halos so ghost-cell gradients from the previous sweep are
        ! valid before building the skewness correction for this sweep.
        call exchange_grad_halos(self)
        ! F3: reuse persistent buffers (allocated in solver_init).
        !$acc parallel loop collapse(3) default(present)
        do c = 1, ntot
          do j = 1, 3
            do i = 1, 3
              gpU(i, j, c) = gU(i, j, c)
            end do
          end do
        end do
        call green_gauss_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd, self%w_g_prev_U)
        if (self%solve_temperature) then
          !$acc parallel loop collapse(2) default(present)
          do c = 1, ntot
            do i = 1, 3
              gpT(i, c) = gT(i, c)
            end do
          end do
          call green_gauss_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd, self%w_g_prev_T)
        end if
        do s = 1, self%params%n_scalars
          ! Reload prev-gradient scratch for each scalar field.
          gC => self%grad_C
          !$acc parallel loop collapse(2) default(present)
          do c = 1, ntot
            do i = 1, 3
              gpT(i, c) = gC(i, c, s)
            end do
          end do
          call load_scalar_row(self, s)
          call green_gauss_gradient_scalar(self%w_Cs, self%mesh, self%grad_C(:,:,s), &
                                           self%w_Cs_bnd, self%w_g_prev_T)
        end do
      end do
    end if

    ! Halo-exchange gradients so ghost cells carry valid values for viscous/advection stencils.
    call exchange_grad_halos(self)
  end subroutine compute_all_gradients

  !> Copy row s of C/C_bnd into the contiguous w_Cs/w_Cs_bnd scratch (on the device when the
  !> fields are resident, on the host otherwise).
  subroutine load_scalar_row(self, s)
    class(FVSolver), intent(inout), target :: self
    integer(ip),     intent(in)    :: s
    integer(ip) :: c, ntot, nbnd
    real(rp), pointer :: Cf(:,:), Cb(:,:), Cs(:), Csb(:)

    ntot = self%mesh%n_cells + self%mesh%n_ghost
    nbnd = self%mesh%n_bnd_faces
    Cf => self%C;  Cb => self%C_bnd;  Cs => self%w_Cs;  Csb => self%w_Cs_bnd
    !$acc parallel loop if(self%on_device) default(present)
    do c = 1, ntot
      Cs(c) = Cf(s, c)
    end do
    !$acc parallel loop if(self%on_device) default(present)
    do c = 1, nbnd
      Csb(c) = Cb(s, c)
    end do
  end subroutine load_scalar_row

  !> Pack, halo-exchange, and unpack grad_U/grad_T/grad_C across MPI ranks (on the device when
  !> the fields are resident).
  subroutine exchange_grad_halos(self)
    class(FVSolver), intent(inout), target :: self

    integer(ip) :: nc, ng, ntot, n_scal, n_fields, c, s, k, i, j
    real(rp), pointer :: pk(:,:), gU(:,:,:), gT(:,:), gC(:,:,:)
    logical :: dev

    if (self%mesh%nranks <= 1) return

    nc     = self%mesh%n_cells
    ng     = self%mesh%n_ghost
    ntot   = nc + ng
    n_scal = self%params%n_scalars
    ! grad_U: 9 components (3x3), grad_T: 3 components, grad_C: 3*n_scal
    n_fields = 9 + 3 + 3 * n_scal
    if (n_fields == 0 .or. ntot == 0) return
    if (.not. allocated(self%w_pack_grad)) then
      allocate(self%w_pack_grad(n_fields, ntot))
    end if
    pk => self%w_pack_grad
    gU => self%grad_U
    gT => self%grad_T
    dev = self%on_device

    ! Pack owned cells; ghost slots will be filled by the exchange.
    !$acc parallel loop if(dev) default(present)
    do c = 1, nc
      do j = 1, 3
        do i = 1, 3
          pk(3*(j-1) + i, c) = gU(i, j, c)
        end do
      end do
      pk(10, c) = gT(1, c);  pk(11, c) = gT(2, c);  pk(12, c) = gT(3, c)
    end do
    if (n_scal > 0) then
      gC => self%grad_C
      !$acc parallel loop collapse(2) if(dev) default(present) private(k)
      do c = 1, nc
        do s = 1, n_scal
          k = 12 + 3*(s-1)
          pk(k+1, c) = gC(1, c, s);  pk(k+2, c) = gC(2, c, s);  pk(k+3, c) = gC(3, c, s)
        end do
      end do
    end if

    call self%mesh%exchange_halos_scalars(self%w_pack_grad)

    ! Unpack ghost cells (no need to overwrite owned cells).
    !$acc parallel loop if(dev) default(present)
    do c = nc+1, ntot
      do j = 1, 3
        do i = 1, 3
          gU(i, j, c) = pk(3*(j-1) + i, c)
        end do
      end do
      gT(1, c) = pk(10, c);  gT(2, c) = pk(11, c);  gT(3, c) = pk(12, c)
    end do
    if (n_scal > 0) then
      !$acc parallel loop collapse(2) if(dev) default(present) private(k)
      do c = nc+1, ntot
        do s = 1, n_scal
          k = 12 + 3*(s-1)
          gC(1, c, s) = pk(k+1, c);  gC(2, c, s) = pk(k+2, c);  gC(3, c, s) = pk(k+3, c)
        end do
      end do
    end if
  end subroutine exchange_grad_halos

  !> Halo-exchange the SGS viscosity so that processor-face viscous flux
  !> uses the same nu_sgs on both sides.
  subroutine exchange_nu_sgs_halo(self)
    class(FVSolver), intent(inout) :: self
    if (self%mesh%nranks <= 1) return
    if (.not. allocated(self%nu_sgs)) return
    call self%mesh%exchange_halos_scalar(self%nu_sgs)
  end subroutine exchange_nu_sgs_halo

  ! Wray RK3 (Williamson 1980), incremental projection each stage; alpha={8/15,5/12,3/4}, beta={0,-17/60,-5/12}, gamma={8/15,2/15,1/3}.
  ! Stage k: u*=u^k+dt*(alpha_k*H^k+beta_k*H_prev)-gamma_k*dt*grad(p^{k-1/2}); phi=solve(div(u*)); u+=...; p+=phi. w_U0/w_T0/w_C0=H_prev.
  subroutine advance_rk3w(self, dt)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(in)    :: dt

    integer(ip) :: nc, ntot, nbnd, nfaces, ns, nsc, c, k, s, jno, p_iter
    real(rp)    :: gk_dt, time_n, a_k, b_k, nu_mol
    real(rp)    :: t_perf, t_sub, p_resid
    ! Plain-array aliases for the device loops below (the polymorphic self can't be mapped).
    real(rp), pointer :: U(:,:), T(:), Cf(:,:), p(:), p_old(:), nu_sgs(:), gp_rc(:,:), gp_cell(:,:)
    real(rp), pointer :: phi_prev(:), phi_prev_bnd(:), p_bnd(:), U_face(:), face_flux(:), inv_vol(:)
    real(rp), pointer :: cvol(:), tau(:,:), hU(:,:), hT(:), hC(:,:), hU0(:,:), hT0(:), hC0(:,:), rhsp(:)

    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nbnd = self%mesh%n_bnd_faces
    nfaces = self%mesh%n_faces
    ns = self%params%n_scalars
    nsc = max(1_ip, ns)
    nu_mol = self%params%nu
    U => self%U;  T => self%T;  p => self%p;  p_old => self%p_old;  nu_sgs => self%nu_sgs
    gp_rc => self%w_gp_rc;  gp_cell => self%w_gp_cell;  phi_prev => self%phi_prev
    phi_prev_bnd => self%phi_prev_bnd;  p_bnd => self%p_bnd;  U_face => self%U_face
    face_flux => self%face_flux;  inv_vol => self%inv_vol;  cvol => self%mesh%cell_volume_soa
    tau => self%tau_w_bnd
    hU => self%w_H_U;  hT => self%w_H_T;  hC => self%w_H_C;  rhsp => self%w_rhs_p
    hU0 => self%w_U0;  hT0 => self%w_T0;  hC0 => self%w_C0
    Cf => self%w_H_C  ! placeholder when no scalars (never touched)
    if (ns > 0) Cf => self%C

    associate(H_U => self%w_H_U(:, 1:nc), H_T => self%w_H_T(1:nc), &
          rhs_p => self%w_rhs_p(1:nc), &
          H_C => self%w_H_C(1:nsc, 1:nc))

    !$acc parallel loop default(present)
    do c = 1, nc
      hU(1, c) = 0.0_rp;  hU(2, c) = 0.0_rp;  hU(3, c) = 0.0_rp
      hT(c) = 0.0_rp
      rhsp(c) = 0.0_rp
      ! H_prev: beta(1)=0 so this is never read, but zero for safety
      hU0(1, c) = 0.0_rp;  hU0(2, c) = 0.0_rp;  hU0(3, c) = 0.0_rp
      hT0(c) = 0.0_rp
    end do
    !$acc parallel loop collapse(2) default(present)
    do c = 1, nc
      do s = 1, nsc
        hC(s, c) = 0.0_rp
        hC0(s, c) = 0.0_rp
      end do
    end do

    time_n = self%time

    ! Counters accumulate over all three stages and are harvested once per step.
    call advection_reset_tvd_fallback()
    call wall_model_reset_cap_count()

    ! -- Three-stage Wray RK3 loop
    do k = 1, 3
      gk_dt = RK3W_GAMMA(k) * dt
      self%time = time_n + RK3W_C(k) * dt

      ! -- BCs + halos
      t_perf = wall_time_now()
      call prof_push('halo_bc')
      call refresh_halos_and_bc(self)
      call prof_pop()
      self%perf_t_halo_bc = self%perf_t_halo_bc + (wall_time_now() - t_perf)

      ! -- SGS viscosity, gradients, wall model, explicit RHS.
      t_perf = wall_time_now()
      t_sub = wall_time_now()
      call prof_push('gradients')
      call compute_all_gradients(self)
      call prof_pop()
      self%perf_t_rhs_grad = self%perf_t_rhs_grad + (wall_time_now() - t_sub)

      t_sub = wall_time_now()
      call prof_push('sgs_wall_model')
      call prof_push('sgs_model')
      if (self%params%use_sgs) then
        select case (trim(self%params%sgs_model))
        case ('SMAGORINSKY')
          call smagorinsky_sgs(self%grad_U, self%mesh, self%nu_sgs, dt_stab=gk_dt)
        case ('DYNAMIC_SMAG')  ! Germano-Lilly dynamic procedure, see mod_dynamic_smagorinsky
          call dynamic_smagorinsky_sgs(self%U, self%grad_U, self%mesh, self%nu_sgs, &
                                       U_bnd=self%U_bnd, dt_stab=gk_dt)
        case default  ! 'VREMAN' or unknown
          call vreman_sgs_viscosity(self%U, self%mesh, self%nu_sgs, U_bnd=self%U_bnd, &
                                    grad_U_ext=self%grad_U, dt_stab=gk_dt)
        end select
      else
        !$acc parallel loop default(present)
        do c = 1, nc
          nu_sgs(c) = 0.0_rp
        end do
      end if

      call prof_pop()
      call prof_push('nu_sgs_halo')
      call exchange_nu_sgs_halo(self)
      call prof_pop()
      call prof_push('wall_model')
      select case (self%abl_config%wall_model_id)
      case (WM_MOST)
        ! MOST runs on the host: refresh its inputs there and push its outputs back.
        !$acc update host(self%U, self%T, self%tau_w_bnd, self%q_w_bnd, self%L_bnd) if(self%on_device)
        call abl_wall_model(self%U, self%T, self%mesh, self%abl_config, self%bc_configs, &
                            self%tau_w_bnd, self%q_w_bnd, self%L_bnd, update_L=(k == 1))
        !$acc update device(self%tau_w_bnd, self%q_w_bnd, self%L_bnd) if(self%on_device)
      case (WM_NONE)
        !$acc parallel loop default(present)
        do c = 1, nbnd
          tau(1, c) = 0.0_rp;  tau(2, c) = 0.0_rp;  tau(3, c) = 0.0_rp
        end do
      case default
        ! Filter advances once per step, so stages 2 and 3 pass dt = 0.
        call apply_wall_model_bc(self%U, self%params%nu, self%mesh, self%tau_w_bnd, &
                                 h_wm=self%params%wall_model_height, &
                                 t_filt=self%params%wall_model_filter_time, &
                                 dt=merge(dt, 0.0_rp, k == 1))
      end select
      call prof_pop()
      call prof_pop()
      self%perf_t_rhs_sgs = self%perf_t_rhs_sgs + (wall_time_now() - t_sub)

      ! -- Explicit RHS H(u^k)
      call prof_push('explicit_rhs')
      call compute_explicit_rhs(self, H_U, H_T, H_C, self%time)
      call prof_pop()
      self%perf_t_rhs = self%perf_t_rhs + (wall_time_now() - t_perf)

      ! -- Predictor: u* = u^k + dt*(alpha_k*H + beta_k*H_prev)
      t_perf = wall_time_now()
      call prof_push('predictor')
      t_sub = wall_time_now()
      a_k = RK3W_ALPHA(k)
      b_k = RK3W_BETA(k)
#ifdef _OPENACC
      !$acc parallel loop default(present)
#else
      !$omp parallel do schedule(static)
#endif
      do c = 1, nc
        U(1, c) = U(1, c) + dt * (a_k * hU(1, c) + b_k * hU0(1, c))
        U(2, c) = U(2, c) + dt * (a_k * hU(2, c) + b_k * hU0(2, c))
        U(3, c) = U(3, c) + dt * (a_k * hU(3, c) + b_k * hU0(3, c))
        T(c)    = T(c)    + dt * (a_k * hT(c)    + b_k * hT0(c))
      end do
#ifndef _OPENACC
      !$omp end parallel do
#endif
      if (ns > 0) then
#ifdef _OPENACC
        !$acc parallel loop collapse(2) default(present)
#else
        !$omp parallel do schedule(static) private(s)
#endif
        do c = 1, nc
          do s = 1, ns
            Cf(s, c) = Cf(s, c) + dt * (a_k * hC(s, c) + b_k * hC0(s, c))
          end do
        end do
#ifndef _OPENACC
        !$omp end parallel do
#endif
      end if

      call self%mesh%begin_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, &
                             periodic_only=.false., grad_U=self%grad_U, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)
      call self%mesh%end_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      self%perf_t_predictor_bc = self%perf_t_predictor_bc + (wall_time_now() - t_sub)

      ! -- Pre-correction: u_hat = u* - gk_dt*grad(p^{k-1/2}); net gives p^{k+1/2}=p^{k-1/2}+phi.
      t_sub = wall_time_now()
      !$acc parallel loop default(present)
      do c = 1, nc
        p_old(c) = p(c)
      end do
      call correct_velocity_psolver_dual(self, gk_dt, with_cross=.false.)
      ! Rhie-Chow must add back exactly the gradient just removed, so record it here.
      !$acc parallel loop default(present)
      do c = 1, nc
        gp_rc(1, c) = gp_cell(1, c);  gp_rc(2, c) = gp_cell(2, c);  gp_rc(3, c) = gp_cell(3, c)
      end do
      if (self%mesh%nranks > 1 .and. self%mesh%n_ghost > 0) &
        call self%mesh%exchange_halos_vec(self%w_gp_rc)
      self%perf_t_predictor_precorrect = self%perf_t_predictor_precorrect + (wall_time_now() - t_sub)

      t_sub = wall_time_now()
      call self%mesh%begin_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, &
                             periodic_only=.false., grad_U=self%grad_U, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)
      call self%mesh%end_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      self%perf_t_predictor_bc = self%perf_t_predictor_bc + (wall_time_now() - t_sub)

      ! Refresh grad_U on the predicted field (mass-flux skewness correction is O(1) at refinement faces).
      t_sub = wall_time_now()
      call compute_all_gradients(self)
      self%perf_t_predictor_regrad = self%perf_t_predictor_regrad + (wall_time_now() - t_sub)
      call prof_pop()
      self%perf_t_predictor = self%perf_t_predictor + (wall_time_now() - t_perf)

      ! -- Save p^{k-1/2} into phi_prev for the Rhie-Chow correction.
      ! phi_prev is frozen across jno iterations; RC always sees the same pressure.
      !$acc parallel loop default(present)
      do c = 1, ntot
        phi_prev(c) = p(c)
      end do
      !$acc parallel loop default(present)
      do c = 1, nbnd
        phi_prev_bnd(c) = p_bnd(c)
      end do

      ! -- Pressure Poisson with effective sub-step gk_dt
      t_perf = wall_time_now()
      call prof_push('pressure')
      block
        real(rp) :: u_bulk_local(2), u_bulk_global(2), u_scale, atol_dynamic
        integer :: ierr_tol
        integer(ip) :: n_nonortho_max
        real(rp) :: divU_norm_prev, divU_norm_curr, divU_local, stagnate_tol, ub1, ub2
        ! Volume-weighted mean speed, not the domain max: max(|U|) is a single-cell
        ! quantity and loosens atol precisely when one cell is diverging.
        ub1 = 0.0_rp
        ub2 = 0.0_rp
        !$acc parallel loop default(present) reduction(+:ub1, ub2)
        do c = 1, nc
          ub1 = ub1 + sqrt(U(1, c)**2 + U(2, c)**2 + U(3, c)**2) * cvol(c)
          ub2 = ub2 + cvol(c)
        end do
        u_bulk_local(1) = ub1
        u_bulk_local(2) = ub2
        call MPI_Allreduce(u_bulk_local, u_bulk_global, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, self%mesh%comm, ierr_tol)
        u_scale = u_bulk_global(1) / max(u_bulk_global(2), SMALL)
        atol_dynamic = max(1.0e-12_rp * max(u_scale, 1.0_rp) * self%V_min_global / gk_dt, &
                           self%p_atol_factor)
        ! A 2:1 face contracts at ~0.47/sweep, so the 3-sweep ceiling stops well short of div-free.
        if (self%diag_interface) then
          n_nonortho_max = min(self%params%n_nonortho_corr, 8_ip)
          stagnate_tol   = 0.99_rp
        else
          n_nonortho_max = min(self%params%n_nonortho_corr, 3_ip)
          stagnate_tol   = 0.9_rp
        end if
        divU_norm_prev = huge(1.0_rp)
        do jno = 1, n_nonortho_max
          call compute_face_fluxes_and_rhs(self, gk_dt, rhs_p, p_is_increment=(jno > 1))
          ! Global reduction: the exit test below must fire on every rank at the
          ! same jno or the ranks desynchronise on the collective solve.
          divU_local = 0.0_rp
          !$acc parallel loop default(present) reduction(max:divU_local)
          do c = 1, nc
            divU_local = max(divU_local, abs(rhsp(c)))
          end do
          call MPI_Allreduce(divU_local, divU_norm_curr, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, self%mesh%comm, ierr_tol)
          if (self%diag_interface) self%diag_nonortho_iters = jno
          ! self%p enters this loop holding p^{k-1/2} (needed by compute_face_fluxes_and_rhs
          ! above for the jno==1 gradient), but the KSP solves for the increment phi. Seeding
          ! its nonzero initial guess with O(1) p^{k-1/2} instead of the O(dt*div u*) increment
          ! wastes iterations and loosens the effective (relative) convergence tolerance.
          if (jno == 1) then
            !$acc parallel loop default(present)
            do c = 1, nc
              p(c) = 0.0_rp
            end do
          end if
          call prof_push('ksp_solve')
          ! With resident fields the RHS and solution go to PETSc as device pointers.
          call pressure_solver_solve(self%psolver, rhs_p, self%p(1:nc), tol=self%params%p_tol, &
                         max_iter=self%params%p_max_iter, n_iter=p_iter, residual=p_resid, &
                         atol_in=atol_dynamic, on_device=self%on_device)
          call prof_pop()
          call self%mesh%begin_halos_scalar(self%p)
          if (self%skip_grad_p) then
            call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, &
                                   periodic_only=.false.)
          else
            call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, &
                                   periodic_only=.false., grad_p=self%w_gp)
          end if
          call self%mesh%end_halos_scalar(self%p)
          call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
          self%perf_pressure_solves = self%perf_pressure_solves + 1_ip
          self%perf_pressure_iters  = self%perf_pressure_iters  + p_iter
          self%perf_pressure_residual_max = max(self%perf_pressure_residual_max, p_resid)
          if (self%diag_interface) then
            self%diag_nonortho_curr = divU_norm_curr
            self%diag_nonortho_prev = divU_norm_prev
            self%diag_nonortho_atol = atol_dynamic
          end if
          ! Exit early if converged or not improving (orthogonal mesh converges in 1 iter).
          if (jno > 1 .and. (divU_norm_curr < atol_dynamic .or. divU_norm_curr > divU_norm_prev * stagnate_tol)) exit
          ! phi_prev stays at p^{k-1/2} for every corrector: it is what the pre-correction removed.
          divU_norm_prev = divU_norm_curr
        end do
      end block
      call prof_pop()
      self%perf_t_pressure = self%perf_t_pressure + (wall_time_now() - t_perf)

      ! -- Project: correct face fluxes and velocity (uses final phi from outer loop)
      t_perf = wall_time_now()
      call prof_push('projection')

      call correct_face_fluxes(self, gk_dt)            ! phi -> divergence-free
      ! Persist the divergence-free flux (U_face is resident for the whole run).
      !$acc parallel loop default(present)
      do c = 1, nfaces
        U_face(c) = face_flux(c)
      end do
      self%u_face_ready = .true.
      if (self%diag_interface .and. self%params%debug_diagnostics) then
        call solver_sync_host(self)
        call record_divergence_split(self, jno)
      end if
      call correct_velocity_psolver_dual(self, gk_dt, with_cross=.true.)  ! u^{k+1} = u* - gk_dt*grad(phi)
      call self%mesh%begin_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, &
                             periodic_only=.false., grad_U=self%grad_U, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)
      call self%mesh%end_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      call prof_pop()
      self%perf_t_projection = self%perf_t_projection + (wall_time_now() - t_perf)

      ! -- Rotational correction to phi: phi = phi - (nu+nu_sgs)*div(u*) (Brown, Cortez & Minion 2001)
      if (self%params%use_rotational_projection) then
#ifdef _OPENACC
        !$acc parallel loop default(present)
#else
        !$omp parallel do schedule(static)
#endif
        do c = 1, nc
          p(c) = p(c) + (nu_mol + nu_sgs(c)) * rhsp(c) * gk_dt * inv_vol(c)
        end do
#ifndef _OPENACC
        !$omp end parallel do
#endif
        call self%mesh%begin_halos_scalar(self%p)
        call apply_pressure_bc_o2(self)
        call self%mesh%end_halos_scalar(self%p)
        call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      end if

      ! -- Accumulate: p^{k+1/2} = p^{k-1/2} + phi (CaNS Eq. 2d).
      !$acc parallel loop default(present)
      do c = 1, nc
        p(c) = p(c) + p_old(c)
      end do
      call anchor_pressure_to_outlet(self%p(1:nc), self%mesh, self%w_outlet_cell, self%w_outlet_pval, self%on_device)
      call self%mesh%begin_halos_scalar(self%p)
      call apply_pressure_bc_o2(self)
      call self%mesh%end_halos_scalar(self%p)
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.true.)

      ! Leave w_gp = grad(p^{k+1/2}) for the next stage's Neumann-p extrapolation and
      ! pre-correction; without this they would see the gradient of the increment.
      call compute_pressure_gradient(self)

      ! -- Save current H as H_prev for next stage
      !$acc parallel loop default(present)
      do c = 1, nc
        hU0(1, c) = hU(1, c);  hU0(2, c) = hU(2, c);  hU0(3, c) = hU(3, c)
        hT0(c) = hT(c)
      end do
      if (ns > 0) then
        !$acc parallel loop collapse(2) default(present)
        do c = 1, nc
          do s = 1, ns
            hC0(s, c) = hC(s, c)
          end do
        end do
      end if

    end do  ! k = 1, 3

    ! Harvest rank-local diagnostic counters into the warning log (no I/O, no MPI here).
    call warn_log_bump(WARN_TVD_FALLBACK, advection_get_tvd_fallback_count())
    call warn_log_bump(WARN_WALL_UTAU_CAP, wall_model_get_cap_count())

    self%time = time_n

    if (self%abl_config%use_sem .or. self%inflow_planes_on) then
      if (.not. self%inflow_planes_on) call sem_advance(self%sem_state, self%abl_config, dt)
      call update_sem_inflow(self, time_n + dt)
      !$acc update device(self%sem_fluc_bnd) if(self%on_device)
      if (allocated(self%abl_mean_bnd)) then
        !$acc update device(self%abl_mean_bnd) if(self%on_device)
      end if
    end if

    end associate

  end subroutine advance_rk3w

  ! Explicit RHS: H = -conv + visc + buoy
  subroutine compute_explicit_rhs(self, H_U, H_T, H_C, time_eval)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(inout) :: H_U(:,:)  ! (3, nc)
    real(rp),        intent(inout) :: H_T(:)  ! (nc)
    real(rp),        intent(inout) :: H_C(:,:)  ! (n_scalars, nc)
    real(rp),        intent(in)    :: time_eval

    integer(ip) :: nc, c, s, ntot, nbnd
    real(rp)    :: force_bg(3), dir_n, hscale, hdir(3), om(3), gvec(3), beta_b, T_ref
    real(rp)    :: t_sub
    logical     :: rot, bous
    real(rp), pointer :: cU(:,:), vU(:,:), cT(:), vT(:), cC(:), vC(:), U(:,:), T(:), Cf(:,:), Cb(:,:)
    real(rp), pointer :: Cs(:), Csb(:)

    t_sub = wall_time_now()
    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nbnd = self%mesh%n_bnd_faces
    cU => self%w_conv_U;  vU => self%w_visc_U;  cT => self%w_conv_T;  vT => self%w_visc_T
    U => self%U;  T => self%T

    ! The advection/viscous kernels zero their own output on the device.
    ! Advection uses U_face (divergence-free flux from last projection) for discrete mass conservation.
    call advection_rhs_vector(self%U, self%mesh, self%params%scheme, self%w_conv_U(:, 1:nc), &
            U_bnd=self%U_bnd, face_flux=self%U_face, inv_vol=self%inv_vol, grad_U=self%grad_U)
    call viscous_rhs_vector(self%U, self%params%nu, self%nu_sgs, &
           self%grad_U, self%mesh, self%w_visc_U(:, 1:nc), U_bnd=self%U_bnd, inv_vol=self%inv_vol)

    if (self%solve_temperature) then
      call advection_rhs_scalar(self%T, self%U, self%mesh, self%params%temperature_scheme, self%w_conv_T(1:nc), &
                    phi_bnd=self%T_bnd, U_bnd=self%U_bnd, &
                    face_flux=self%U_face, inv_vol=self%inv_vol, grad_phi=self%grad_T)
      call viscous_rhs_scalar(self%T, self%params%alpha_T, self%nu_sgs, &
             self%params%Pr_t, self%grad_T, self%mesh, self%w_visc_T(1:nc), &
                   phi_bnd=self%T_bnd, inv_vol=self%inv_vol)
    else
      !$acc parallel loop default(present)
      do c = 1, nc
        cT(c) = 0.0_rp
        vT(c) = 0.0_rp
      end do
    end if

    ! Precompute cell-independent body force: constant + time-varying harmonic.
    ! Coriolis (-2 Omega x U) is U-dependent and remains inside the loop.
    if (self%params%constant_mass_flux_forcing) then
      force_bg = self%mass_flux_force
    else
      force_bg = self%params%f_const
    end if
    if (self%params%use_harmonic_forcing) then
      dir_n = sqrt(sum(self%params%forcing_dir**2))
      if (dir_n > SMALL) then
        hdir = self%params%forcing_dir / dir_n
      else
        hdir = [1.0_rp, 0.0_rp, 0.0_rp]
      end if
      hscale   = self%params%forcing_amp * cos(self%params%forcing_omega * time_eval + self%params%forcing_phase)
      force_bg = force_bg + hscale * hdir
    end if
    rot = self%params%use_rotation
    bous = self%params%boussinesq
    om = self%params%omega_rot
    gvec = self%params%g
    beta_b = self%params%beta
    T_ref = self%params%T_ref

    ! H = -conv + visc + Coriolis + buoyancy + forcing (grad p applied in projection)
#ifdef _OPENACC
    !$acc parallel loop default(present) copyin(force_bg, om, gvec)
#else
    !$omp parallel do schedule(static)
#endif
    do c = 1, nc
      H_U(1, c) = cU(1, c) + vU(1, c) + force_bg(1)
      H_U(2, c) = cU(2, c) + vU(2, c) + force_bg(2)
      H_U(3, c) = cU(3, c) + vU(3, c) + force_bg(3)
      if (rot) then
        H_U(1, c) = H_U(1, c) - 2.0_rp * (om(2) * U(3, c) - om(3) * U(2, c))
        H_U(2, c) = H_U(2, c) - 2.0_rp * (om(3) * U(1, c) - om(1) * U(3, c))
        H_U(3, c) = H_U(3, c) - 2.0_rp * (om(1) * U(2, c) - om(2) * U(1, c))
      end if
      if (bous) then
        H_U(1, c) = H_U(1, c) - beta_b * (T(c) - T_ref) * gvec(1)
        H_U(2, c) = H_U(2, c) - beta_b * (T(c) - T_ref) * gvec(2)
        H_U(3, c) = H_U(3, c) - beta_b * (T(c) - T_ref) * gvec(3)
      end if
      H_T(c) = cT(c) + vT(c)
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif

    ! Wall-model traction/heat flux enter H so RK3 weights them like every other term.
    call add_wall_stress_to_rhs(self, H_U)
    if (self%abl_config%wall_model_id == WM_MOST) call add_wall_heat_flux_to_rhs(self, H_T)

    if (self%params%n_scalars > 0) then
      cC => self%w_conv_C;  vC => self%w_visc_C;  Cf => self%C;  Cb => self%C_bnd
      Cs => self%w_Cs;  Csb => self%w_Cs_bnd
      do s = 1, self%params%n_scalars
        ! Rows of C/C_bnd are strided; hand the kernels contiguous copies.
        !$acc parallel loop default(present)
        do c = 1, ntot
          Cs(c) = Cf(s, c)
        end do
        !$acc parallel loop default(present)
        do c = 1, nbnd
          Csb(c) = Cb(s, c)
        end do
        call advection_rhs_scalar(self%w_Cs, self%U, self%mesh, self%params%scalar_scheme, self%w_conv_C(1:nc), &
                      phi_bnd=self%w_Cs_bnd, U_bnd=self%U_bnd, face_flux=self%U_face, inv_vol=self%inv_vol, &
                      grad_phi=self%grad_C(:, :, s))
        call viscous_rhs_scalar(self%w_Cs, self%params%scalar_diffusivity(s), self%nu_sgs, &
                     self%params%scalar_pr_t(s), self%grad_C(:, :, s), self%mesh, self%w_visc_C(1:nc), &
                     phi_bnd=self%w_Cs_bnd, inv_vol=self%inv_vol)
        !$acc parallel loop default(present)
        do c = 1, nc
          H_C(s, c) = cC(c) + vC(c)
        end do
      end do
    end if
    self%perf_t_rhs_transport = self%perf_t_rhs_transport + (wall_time_now() - t_sub)

    ! Scalar-source and large-scale nudging work: not transport. Both run on the host.
    t_sub = wall_time_now()
    if (self%params%n_scalars > 0 .and. self%params%n_scalar_sources > 0) then
      !$acc update host(H_C) if(self%on_device)
      call add_scalar_source_terms(self, H_C, time_eval)
      !$acc update device(H_C) if(self%on_device)
    end if

    ! Apply Davies nudging relaxation to large-scale profiles
    if (self%abl_config%use_nudging) then
      !$acc update host(self%U, self%T, H_U, H_T) if(self%on_device)
      call nudging_apply_rhs(self%nudging_data, self%U, self%T, self%cell_z, &
                             nc, time_eval, self%abl_config%nudging_tau, &
                             self%abl_config%nudging_z_start, H_U, H_T, &
                             z_full=self%abl_config%nudging_z_full, &
                             horizontal_only=self%abl_config%nudging_horizontal_only, &
                             dt=self%params%dt)
      !$acc update device(H_U, H_T) if(self%on_device)
    end if
    self%perf_t_rhs_sources = self%perf_t_rhs_sources + (wall_time_now() - t_sub)
  end subroutine compute_explicit_rhs

  subroutine add_scalar_source_terms(self, H_C, time_eval)
    class(FVSolver), intent(in) :: self
    real(rp), intent(inout) :: H_C(:,:)
    real(rp), intent(in) :: time_eval

    integer(ip) :: src_id, c, sid
    real(rp) :: src_val

    if (self%params%n_scalar_sources <= 0) return
    if (self%params%n_scalars <= 0) return

    do src_id = 1, self%params%n_scalar_sources
      if (.not. self%params%scalar_sources(src_id)%enabled) cycle
      sid = self%params%scalar_sources(src_id)%scalar_id
      if (sid < 1 .or. sid > self%params%n_scalars) cycle
      do c = 1, self%mesh%n_cells
        src_val = scalar_source_value_at(self%params%scalar_sources(src_id), self%mesh%cells(c)%centroid, time_eval)
        if (abs(src_val) > SMALL) H_C(sid, c) = H_C(sid, c) + src_val
      end do
    end do
  end subroutine add_scalar_source_terms

  pure function cross3(a, b) result(c)
    real(rp), intent(in) :: a(3), b(3)
    real(rp) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross3

  ! Shared helpers

  subroutine build_patch_classification(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: p, f, i_per, i_non, i_wall, i_per_bf, i_non_bf

    self%n_periodic_patches = 0
    self%n_nonperiodic_patches = 0
    self%n_wallmodel_patches = 0
    self%n_periodic_bfaces = 0
    self%n_nonperiodic_bfaces = 0
    do p = 1, self%mesh%n_patches
      if (self%bc_configs(p)%bc_id == BC_PERIODIC) then
        self%n_periodic_patches = self%n_periodic_patches + 1
        self%n_periodic_bfaces = self%n_periodic_bfaces + self%mesh%patches(p)%n_faces
      else
        self%n_nonperiodic_patches = self%n_nonperiodic_patches + 1
        self%n_nonperiodic_bfaces = self%n_nonperiodic_bfaces + self%mesh%patches(p)%n_faces
      end if
      if (self%bc_configs(p)%bc_id == BC_WALL_MODEL) self%n_wallmodel_patches = self%n_wallmodel_patches + 1
    end do

    if (allocated(self%periodic_patches)) deallocate(self%periodic_patches)
    if (allocated(self%nonperiodic_patches)) deallocate(self%nonperiodic_patches)
    if (allocated(self%wallmodel_patches)) deallocate(self%wallmodel_patches)
    if (allocated(self%periodic_bfaces)) deallocate(self%periodic_bfaces)
    if (allocated(self%nonperiodic_bfaces)) deallocate(self%nonperiodic_bfaces)

    allocate(self%periodic_patches(max(1_ip, self%n_periodic_patches)))
    allocate(self%nonperiodic_patches(max(1_ip, self%n_nonperiodic_patches)))
    allocate(self%wallmodel_patches(max(1_ip, self%n_wallmodel_patches)))
    allocate(self%periodic_bfaces(max(1_ip, self%n_periodic_bfaces)))
    allocate(self%nonperiodic_bfaces(max(1_ip, self%n_nonperiodic_bfaces)))

    i_per = 0
    i_non = 0
    i_wall = 0
    i_per_bf = 0
    i_non_bf = 0
    do p = 1, self%mesh%n_patches
      if (self%bc_configs(p)%bc_id == BC_PERIODIC) then
        i_per = i_per + 1
        self%periodic_patches(i_per) = p
        do f = 1, self%mesh%patches(p)%n_faces
          i_per_bf = i_per_bf + 1
          self%periodic_bfaces(i_per_bf) = self%mesh%patches(p)%start_face + f - 1
        end do
      else
        i_non = i_non + 1
        self%nonperiodic_patches(i_non) = p
        do f = 1, self%mesh%patches(p)%n_faces
          i_non_bf = i_non_bf + 1
          self%nonperiodic_bfaces(i_non_bf) = self%mesh%patches(p)%start_face + f - 1
        end do
      end if
      if (self%bc_configs(p)%bc_id == BC_WALL_MODEL) then
        i_wall = i_wall + 1
        self%wallmodel_patches(i_wall) = p
      end if
    end do

    !> Validate wallmodel_patches array for consistency.
    !> If BC_WALL_MODEL patches are marked but not populated correctly,
    if (self%n_wallmodel_patches > 0) then
      do i_wall = 1, self%n_wallmodel_patches
        p = self%wallmodel_patches(i_wall)
        if (p < 1 .or. p > self%mesh%n_patches) then
          error stop 'build_patch_classification: wallmodel_patches index out of bounds'
        end if
        if (self%bc_configs(p)%bc_id /= BC_WALL_MODEL) then
          error stop 'build_patch_classification: wallmodel_patches contains non-BC_WALL_MODEL patch'
        end if
      end do
    end if
  end subroutine build_patch_classification

  !> Set skip_nonortho/skip_grad_p by evaluating the exact runtime predicate on every face once at init.
  subroutine detect_nonortho_skip(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: idx, o, nb, f_mesh, ib, f_bnd, f_glob, q, nc, ierr
    integer(ip) :: n_active_local, n_active_global, n_gradbc_local, n_gradbc_global
    real(rp)    :: d_on(3), k_f(3), s_dot_dhat, d_on_len

    nc = self%mesh%n_cells
    n_active_local = 0
    n_gradbc_local = 0

    ! Mirrors the interior/periodic correction loop in compute_face_fluxes_and_rhs.
    do idx = 1, self%psolver%n_faces
      if (self%psolver%skip_face(idx) /= 0) cycle
      o      = int(self%psolver%face_owner(idx), ip)
      nb     = self%psolver%face_neighbor_ext(idx)
      f_mesh = self%psolver%face_mesh_idx(idx)
      if (o < 1 .or. o > nc) cycle
      if (nb < 1 .or. nb > nc + self%mesh%n_ghost) cycle
      if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle
      d_on = self%mesh%cells(nb)%centroid - self%mesh%cells(o)%centroid
      call nonortho_kfactor_overrelaxed(self%mesh%faces(f_mesh)%sf, d_on, k_f, s_dot_dhat, d_on_len)
      if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(self%mesh%faces(f_mesh)%sf**2)) &
        n_active_local = n_active_local + 1
    end do

    ! Mirrors the Dirichlet-pressure (outlet) correction loop.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd  = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      if (self%bc_configs(q)%bc_id /= BC_OUTLET) cycle
      d_on = self%mesh%faces(f_glob)%centre - self%mesh%cells(o)%centroid
      call nonortho_kfactor_overrelaxed(self%mesh%faces(f_glob)%sf, d_on, k_f, s_dot_dhat, d_on_len)
      if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(self%mesh%faces(f_glob)%sf**2)) &
        n_active_local = n_active_local + 1
    end do

    ! apply_pressure_bc also consumes grad(p) for O(h^2) Neumann extrapolation on any other patch type.
    do q = 1, self%mesh%n_patches
      select case (self%bc_configs(q)%bc_id)
      case (BC_OUTLET, BC_PERIODIC, BC_WALL_NOSLIP, BC_WALL_VELOCITY, BC_WALL_MODEL, BC_SLIP, BC_SYMMETRY)
        continue
      case default
        n_gradbc_local = n_gradbc_local + 1
      end select
    end do

    call MPI_Allreduce(n_active_local, n_active_global, 1, MPI_INTEGER, MPI_SUM, self%mesh%comm, ierr)
    call MPI_Allreduce(n_gradbc_local, n_gradbc_global, 1, MPI_INTEGER, MPI_SUM, self%mesh%comm, ierr)

    self%skip_nonortho = (n_active_global == 0)
    self%skip_grad_p   = self%skip_nonortho .and. (n_gradbc_global == 0)

    ! n_nonortho_corr <= 0 means auto: one corrector on orthogonal meshes, two otherwise.
    if (self%params%n_nonortho_corr <= 0_ip) then
      if (self%skip_nonortho) then
        self%params%n_nonortho_corr = 1_ip
      else
        self%params%n_nonortho_corr = 2_ip
      end if
    end if

    if (self%mesh%rank == 0) then
      if (self%skip_grad_p) then
        write(*,'(A)') ' [nonortho] mesh orthogonal within tolerance: skipping deferred correction and grad(p)'
      else if (self%skip_nonortho) then
        write(*,'(A)') ' [nonortho] mesh orthogonal within tolerance: skipping deferred correction'
      end if
      write(*,'(A,I0)') ' [nonortho] correctors per stage = ', self%params%n_nonortho_corr
    end if
  end subroutine detect_nonortho_skip

  !> Precompute the Jasak min-correction k-vector per face; static mesh geometry, so this replaces
  !> up to a dozen sqrt+divide recomputations per timestep in the hot RK/corrector loops below.
  subroutine build_nonortho_cache(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: idx, o, nb, f_mesh, ib, f_bnd, f_glob, q
    real(rp)    :: d_on(3), k_f(3), s_dot_dhat, d_on_len

    if (allocated(self%nonortho_k_int))      deallocate(self%nonortho_k_int)
    if (allocated(self%nonortho_active_int)) deallocate(self%nonortho_active_int)
    if (allocated(self%nonortho_k_bnd))      deallocate(self%nonortho_k_bnd)
    if (allocated(self%nonortho_active_bnd)) deallocate(self%nonortho_active_bnd)

    allocate(self%nonortho_k_int(3, self%psolver%n_faces), source=0.0_rp)
    allocate(self%nonortho_active_int(self%psolver%n_faces), source=.false.)
    allocate(self%nonortho_k_bnd(3, self%n_nonperiodic_bfaces), source=0.0_rp)
    allocate(self%nonortho_active_bnd(self%n_nonperiodic_bfaces), source=.false.)

    if (self%skip_nonortho) return  ! every face fails the tolerance test; leave the caches zeroed

    ! Mirrors the interior/periodic d_on used by compute_face_fluxes_and_rhs, correct_face_fluxes,
    ! and correct_velocity_psolver_dual (identical for all three: same psolver face list).
    do idx = 1, self%psolver%n_faces
      if (self%psolver%skip_face(idx) /= 0) cycle
      o      = int(self%psolver%face_owner(idx), ip)
      nb     = self%psolver%face_neighbor_ext(idx)
      f_mesh = self%psolver%face_mesh_idx(idx)
      if (o < 1 .or. o > self%mesh%n_cells) cycle
      if (nb < 1 .or. nb > self%mesh%n_cells + self%mesh%n_ghost) cycle
      if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle
      d_on = self%mesh%cells(nb)%centroid - self%mesh%cells(o)%centroid
      call nonortho_kfactor_overrelaxed(self%mesh%faces(f_mesh)%sf, d_on, k_f, s_dot_dhat, d_on_len)
      if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(self%mesh%faces(f_mesh)%sf**2)) then
        self%nonortho_k_int(:, idx) = k_f
        self%nonortho_active_int(idx) = .true.
      end if
    end do

    ! Mirrors the Dirichlet-pressure (BC_OUTLET) boundary-face d_on in the same three routines.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd  = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > self%mesh%n_cells) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      if (self%bc_configs(q)%bc_id /= BC_OUTLET) cycle
      d_on = self%mesh%faces(f_glob)%centre - self%mesh%cells(o)%centroid
      call nonortho_kfactor_overrelaxed(self%mesh%faces(f_glob)%sf, d_on, k_f, s_dot_dhat, d_on_len)
      if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(self%mesh%faces(f_glob)%sf**2)) then
        self%nonortho_k_bnd(:, ib) = k_f
        self%nonortho_active_bnd(ib) = .true.
      end if
    end do
  end subroutine build_nonortho_cache

  !> Apply all boundary conditions for current time step.
  subroutine apply_all_bc(self)
    class(FVSolver), intent(inout) :: self
    call apply_all_bc_subset(self)
  end subroutine apply_all_bc

  !> Refresh T/C ghost BCs only; U/p are assumed already current (RK3 stage tail).
  subroutine apply_scalar_bc(self)
    class(FVSolver), intent(inout) :: self
    integer(ip) :: s
    logical     :: has_diff

    call apply_temperature_bc(self%T, self%T_bnd, self%bc_configs, self%mesh, grad_T=self%grad_T)
    do s = 1, self%params%n_scalars
      has_diff = .false.
      if (allocated(self%params%scalar_diffusivity)) has_diff = (s <= size(self%params%scalar_diffusivity))
      if (has_diff) then
        call apply_passive_scalar_bc(self%C, self%C_bnd, self%bc_configs, self%mesh, s, &
                                     self%params%scalar_diffusivity(s), self%time)
      else
        call apply_passive_scalar_bc(self%C, self%C_bnd, self%bc_configs, self%mesh, s, 0.0_rp, self%time)
      end if
    end do
  end subroutine apply_scalar_bc

  subroutine apply_all_bc_subset(self, periodic_only)
    class(FVSolver), intent(inout) :: self
    logical, optional, intent(in) :: periodic_only
    integer(ip) :: s
    logical     :: has_diff

    ! Gradients held from previous sweep; zeros on first call = zeroth-order (no cold-start error).
    call apply_velocity_bc(self%U, self%U_bnd, &
                            self%bc_configs, self%mesh, periodic_only=periodic_only, &
                            grad_U=self%grad_U, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)
    ! w_gp is not a gradient when skip_grad_p is set; omit it so the BC falls back to zero-gradient.
    if (self%skip_grad_p) then
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, &
                              periodic_only=periodic_only)
    else
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, &
                              periodic_only=periodic_only, grad_p=self%w_gp)
    end if
    call apply_temperature_bc(self%T, self%T_bnd, self%bc_configs, self%mesh, &
                               periodic_only=periodic_only, grad_T=self%grad_T)
    do s = 1, self%params%n_scalars
      ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
      ! must never be evaluated in the same expression as the allocated() check.
      has_diff = .false.
      if (allocated(self%params%scalar_diffusivity)) has_diff = (s <= size(self%params%scalar_diffusivity))
      if (has_diff) then
        call apply_passive_scalar_bc(self%C, self%C_bnd, self%bc_configs, self%mesh, s, &
                                     self%params%scalar_diffusivity(s), self%time, periodic_only=periodic_only)
      else
        call apply_passive_scalar_bc(self%C, self%C_bnd, self%bc_configs, self%mesh, s, 0.0_rp, self%time, &
                                     periodic_only=periodic_only)
      end if
    end do
  end subroutine apply_all_bc_subset

  !> Refresh the inflow fluctuation (and, when prescribed, the mean) on every ABL inlet face for
  !> physical time t. The source is the online SEM or the inflow-plane database.
  !> The patch-mean normal component is removed so the inlet mass flux is unchanged.
  subroutine update_sem_inflow(self, t)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: t

    integer(ip) :: p, f, f_bnd, f_global, ierr, kf
    real(rp)    :: u_mean(3), u_fluc(3), buf(2), gbuf(2), delta, R(6), L(3)
    real(rp), allocatable :: fl_planes(:,:)

    if (.not. allocated(self%sem_fluc_bnd)) return
    self%sem_fluc_bnd = 0.0_rp
    if (.not. self%inflow_planes_on .and. .not. self%sem_state%initialized) return

    if (self%inflow_planes_on) then
      allocate(fl_planes(3, self%inflow_planes%nf))
      if (self%inflow_planes%nf > 0) call inflow_planes_eval(self%inflow_planes, t, fl_planes)
    end if

    buf = 0.0_rp
    kf = 0
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        if (self%inflow_planes_on) then
          kf = kf + 1
          if (self%inflow_planes%f_bnd(kf) /= f_bnd) error stop 'update_inflow: inlet face order changed since binding'
          u_fluc = fl_planes(:, kf)
          call sem_target_at(self%sem_state, self%abl_config, self%mesh%faces(f_global)%centre(3), u_mean, R, L)
        else
          call sem_velocity_at_face(self%sem_state, self%abl_config, &
                                    self%mesh%faces(f_global)%centre, u_mean, u_fluc)
        end if
        self%sem_fluc_bnd(:, f_bnd) = u_fluc
        if (allocated(self%abl_mean_bnd)) self%abl_mean_bnd(:, f_bnd) = u_mean
        buf(1) = buf(1) + dot_product(u_fluc, self%mesh%faces(f_global)%sf)
        buf(2) = buf(2) + self%mesh%faces(f_global)%area
      end do
    end do

    call MPI_Allreduce(buf, gbuf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)
    if (gbuf(2) <= SMALL) return
    delta = gbuf(1) / gbuf(2)

    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        self%sem_fluc_bnd(:, f_bnd) = self%sem_fluc_bnd(:, f_bnd) &
                                    - delta * self%mesh%faces(f_global)%normal
      end do
    end do
  end subroutine update_sem_inflow

  !> Open the inflow-plane database and map it onto this rank's ABL inlet faces. Stops unless the
  !> domain is periodic in the spanwise direction, which the replay (spanwise shifts) relies on.
  subroutine setup_inflow_planes(self, icfg)
    class(FVSolver),    intent(inout) :: self
    type(InflowConfig), intent(in)    :: icfg

    logical     :: found, ok
    integer(ip) :: p, f, nf, k, nbnd, f_bnd, f_global
    real(rp)    :: W
    real(rp), allocatable :: xf(:,:)
    integer(ip), allocatable :: fb(:)

    nbnd = self%mesh%n_bnd_faces
    call sem_init_targets(self%sem_state, self%abl_config, self%mesh, found)
    if (.not. found) error stop "inflow_source = 'planes' needs an ablInlet patch"
    call inflow_check_lateral_periodic(self%mesh, self%sem_state%t1, W, ok)
    if (.not. ok) then
      if (self%mesh%rank == 0) then
        write(*,'(A)') 'ERROR: inflow_source = "planes" needs the domain to be periodic in the spanwise direction'
        write(*,'(A)') '       (paired periodic patches with normals +-t1 across the ablInlet plane); none were found.'
      end if
      error stop 'inflow planes: spanwise periodicity required'
    end if

    nf = 0
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id == BC_ABL_INLET) nf = nf + self%mesh%patches(p)%n_faces
    end do
    allocate(xf(3, nf), fb(nf))
    k = 0
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, self%mesh%patches(p)%n_faces
        k = k + 1
        f_bnd = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        xf(:, k) = self%mesh%faces(f_global)%centre
        fb(k) = f_bnd
      end do
    end do

    call inflow_planes_open(self%inflow_planes, icfg, self%abl_config, self%mesh%rank)
    call inflow_planes_bind(self%inflow_planes, xf, fb, self%sem_state%n_in, self%sem_state%t1, W, self%mesh%rank)
    self%inflow_planes_on = .true.
    if (self%sem_state%override_mean) allocate(self%abl_mean_bnd(3, nbnd), source=0.0_rp)
  end subroutine setup_inflow_planes

  !> Write the inlet plane grid (uniform spanwise nodes, distinct heights) in the format read by
  !> dopamine-inflow-gen, so a database can be generated on exactly these face positions, then stop.
  subroutine export_inflow_grid(self, filename)
    class(FVSolver),  intent(inout) :: self
    character(len=*), intent(in)    :: filename

    type(SEMState) :: st
    logical     :: found, ok
    integer(ip) :: p, f, fg, nl, i, j, ierr, ny
    real(rp)    :: W, y0, dy, zscale
    real(rp), allocatable :: al(:), zl(:), ua(:), uz(:)

    call sem_init_targets(st, self%abl_config, self%mesh, found)
    if (.not. found) error stop 'inflow_export_grid needs an ablInlet patch'
    call inflow_check_lateral_periodic(self%mesh, st%t1, W, ok)
    if (.not. ok) error stop 'inflow_export_grid: the domain is not periodic in the spanwise direction'

    nl = 0
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id == BC_ABL_INLET) nl = nl + self%mesh%patches(p)%n_faces
    end do
    allocate(al(nl), zl(nl))
    i = 0
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, self%mesh%patches(p)%n_faces
        i = i + 1
        fg = self%mesh%n_int_faces + self%mesh%patches(p)%start_face + f - 1
        al(i) = dot_product(self%mesh%faces(fg)%centre, st%t1)
        zl(i) = self%mesh%faces(fg)%centre(3)
      end do
    end do
    zscale = 1.0_rp
    if (nl > 0) zscale = max(maxval(abs(zl)), 1.0_rp)
    call MPI_Allreduce(MPI_IN_PLACE, zscale, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
    call inflow_unique_sorted_global(self%mesh%comm, al, 1.0e-9_rp * W, ua)
    call inflow_unique_sorted_global(self%mesh%comm, zl, 1.0e-9_rp * zscale, uz)
    ny = size(ua)

    if (self%mesh%rank == 0) then
      y0 = ua(1); dy = W / real(ny, rp)
      do j = 1, ny
        if (abs(ua(j) - (y0 + real(j - 1, rp) * dy)) > 1.0e-6_rp * dy) then
          write(*,'(A,I0,A,ES11.4,A)') 'ERROR: inlet spanwise nodes are not uniform (node ', j, ' is off by ', &
            ua(j) - (y0 + real(j - 1, rp) * dy), ' m).'
          error stop 'inflow_export_grid: the inlet must be uniformly spaced across the periodic span'
        end if
      end do
      call inflow_write_grid_file(trim(filename), W, y0, int(ny, ip), st%n_in, uz)
      write(*,'(A,I0,A,I0,A,A)') ' Inflow grid exported: ', ny, ' spanwise nodes x ', size(uz), ' heights -> ', trim(filename)
      write(*,'(A)') ' Use it as gen_grid_file in dopamine-inflow-gen; stopping.'
    end if
    call MPI_Barrier(self%mesh%comm, ierr)
    call MPI_Finalize(ierr)
    stop
  end subroutine export_inflow_grid

  !> Re-evaluate the inflow at the solver's current time; used after a restart restores the clock.
  subroutine solver_refresh_inflow(self)
    class(FVSolver), intent(inout) :: self
    if (.not. self%inflow_planes_on) return
    call update_sem_inflow(self, self%time)
    !$acc update device(self%sem_fluc_bnd) if(self%on_device)
    if (allocated(self%abl_mean_bnd)) then
      !$acc update device(self%abl_mean_bnd) if(self%on_device)
    end if
  end subroutine solver_refresh_inflow

  !> Add the wall-model traction to H_U as an acceleration (thin-wall approximation).
  !> In H rather than applied post-predictor so it carries the RK3 alpha/beta weighting.
  subroutine add_wall_stress_to_rhs(self, H_U)
    class(FVSolver), intent(in), target :: self
    real(rp),        intent(inout) :: H_U(:,:)  ! (3, nc)

    integer(ip) :: f_bnd, f_global, o, q, nc, n_int, n_bnd, n_patches
    real(rp)    :: w
    integer(ip), pointer :: fown(:), fpid(:), pbc(:)
    real(rp),    pointer :: farea(:), cvol(:), tau(:,:)

    if (self%n_wallmodel_patches <= 0) return
    nc = self%mesh%n_cells
    n_int = self%mesh%n_int_faces
    n_bnd = self%mesh%n_bnd_faces
    n_patches = self%mesh%n_patches
    fown => self%mesh%face_owner_soa
    fpid => self%mesh%face_patch_id_soa
    pbc => self%mesh%patch_bc_id_soa
    farea => self%mesh%face_area_soa
    cvol => self%mesh%cell_volume_soa
    tau => self%tau_w_bnd

    ! Wall faces can share an owner cell (edges/corners), hence the atomics.
    !$acc parallel loop default(present) copyin(tau) copy(H_U) private(f_global, o, q, w)
    do f_bnd = 1, n_bnd
      f_global = n_int + f_bnd
      q = fpid(f_global)
      if (q < 1 .or. q > n_patches) cycle
      if (pbc(q) /= BC_WALL_MODEL) cycle
      o = fown(f_global)
      if (o < 1 .or. o > nc) cycle
      w = farea(f_global) / cvol(o)
      !$acc atomic update
      H_U(1, o) = H_U(1, o) + tau(1, f_bnd) * w
      !$acc atomic update
      H_U(2, o) = H_U(2, o) + tau(2, f_bnd) * w
      !$acc atomic update
      H_U(3, o) = H_U(3, o) + tau(3, f_bnd) * w
    end do
  end subroutine add_wall_stress_to_rhs

  !> Add the wall heat flux to H_T as a rate (mirrors add_wall_stress_to_rhs).
  subroutine add_wall_heat_flux_to_rhs(self, H_T)
    class(FVSolver), intent(in), target :: self
    real(rp),        intent(inout) :: H_T(:)  ! (nc)

    integer(ip) :: f_bnd, f_global, o, q, nc, n_int, n_bnd, n_patches
    integer(ip), pointer :: fown(:), fpid(:), pbc(:)
    real(rp),    pointer :: farea(:), cvol(:), qw(:)

    if (.not. allocated(self%q_w_bnd)) return
    if (self%n_wallmodel_patches <= 0) return
    nc = self%mesh%n_cells
    n_int = self%mesh%n_int_faces
    n_bnd = self%mesh%n_bnd_faces
    n_patches = self%mesh%n_patches
    fown => self%mesh%face_owner_soa
    fpid => self%mesh%face_patch_id_soa
    pbc => self%mesh%patch_bc_id_soa
    farea => self%mesh%face_area_soa
    cvol => self%mesh%cell_volume_soa
    qw => self%q_w_bnd

    ! q_w_bnd is kinematic heat flux [K m/s]; area/volume converts to a rate. Atomics: wall
    ! faces can share an owner cell.
    !$acc parallel loop default(present) copyin(qw) copy(H_T) private(f_global, o, q)
    do f_bnd = 1, n_bnd
      f_global = n_int + f_bnd
      q = fpid(f_global)
      if (q < 1 .or. q > n_patches) cycle
      if (pbc(q) /= BC_WALL_MODEL) cycle
      o = fown(f_global)
      if (o < 1 .or. o > nc) cycle
      !$acc atomic update
      H_T(o) = H_T(o) + qw(f_bnd) * farea(f_global) / cvol(o)
    end do
  end subroutine add_wall_heat_flux_to_rhs

  !> Fill w_gp with grad(p) using the configured gradient scheme and halo-exchange it.
  subroutine compute_pressure_gradient(self)
    class(FVSolver), intent(inout) :: self

    if (self%skip_grad_p) return
    call gradient_of_scalar(self, self%p, self%p_bnd, self%w_gp)
  end subroutine compute_pressure_gradient

  !> Non-periodic pressure BCs with the O(h^2) Neumann extrapolation whenever w_gp is a real gradient.
  subroutine apply_pressure_bc_o2(self)
    class(FVSolver), intent(inout) :: self

    if (self%skip_grad_p) then
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.false.)
    else
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, &
                             periodic_only=.false., grad_p=self%w_gp)
    end if
  end subroutine apply_pressure_bc_o2

  !> Gradient of an arbitrary cell scalar with its boundary values, halo-exchanged.
  subroutine gradient_of_scalar(self, fld, fld_bnd, grad)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(in)    :: fld(:)
    real(rp),        intent(in)    :: fld_bnd(:)
    real(rp),        intent(inout) :: grad(:,:)

    integer(ip) :: iter, c
    real(rp), pointer :: gprev(:,:)

    select case (trim(self%params%gradient_scheme))
    case ('GREEN_GAUSS')
      call green_gauss_gradient_scalar(fld, self%mesh, grad, fld_bnd)
      do iter = 2, self%params%gg_iterations
        if (self%mesh%nranks > 1 .and. self%mesh%n_ghost > 0) &
          call self%mesh%exchange_halos_vec(grad)
        gprev => self%w_gp_prev
        !$acc parallel loop if(self%on_device) default(present)
        do c = 1, size(grad, 2)
          gprev(1, c) = grad(1, c);  gprev(2, c) = grad(2, c);  gprev(3, c) = grad(3, c)
        end do
        call green_gauss_gradient_scalar(fld, self%mesh, grad, fld_bnd, self%w_gp_prev)
      end do
    case default
      call ls_gradient_scalar(fld, self%mesh, grad, fld_bnd)
    end select

    if (self%mesh%nranks > 1 .and. self%mesh%n_ghost > 0) &
      call self%mesh%exchange_halos_vec(grad)
  end subroutine gradient_of_scalar

  !> w_p_tot = phi_prev + p (and on the boundary): p^{k-1/2} plus the current increment.
  subroutine build_total_pressure(self)
    class(FVSolver), intent(inout), target :: self
    integer(ip) :: c
    real(rp), pointer :: pt(:), ptb(:), pp(:), ppb(:), p(:), pb(:)

    pt => self%w_p_tot;  ptb => self%w_p_tot_bnd;  pp => self%phi_prev;  ppb => self%phi_prev_bnd
    p => self%p;  pb => self%p_bnd
    !$acc parallel loop if(self%on_device) default(present)
    do c = 1, size(pt)
      pt(c) = pp(c) + p(c)
    end do
    !$acc parallel loop if(self%on_device) default(present)
    do c = 1, size(ptb)
      ptb(c) = ppb(c) + pb(c)
    end do
  end subroutine build_total_pressure

  !> Assemble predictor face fluxes and Poisson RHS: rhs(c) = -(1/dt)*sum_f q*_f.
  subroutine compute_face_fluxes_and_rhs(self, dt, rhs_p, p_is_increment)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(in)    :: dt
    real(rp),        intent(out)   :: rhs_p(:)
    logical, optional, intent(in)  :: p_is_increment  ! .true. once self%p holds phi, not p^{k-1/2}

    integer(ip) :: nc, ntot, nf, n_int, n_faces, n_bnd, n_patches, n_nper, n_per, c, f, q, ib, f_bnd, f_glob
    integer(ip) :: idx, o, nb, f_mesh, pair_bnd, pair_face, om, nm, j
    real(rp)    :: flux, coeff_rc, gc_f, d_o, d_nb, d_tot, w_blend, s_dot_dhat, d_on_len, cross, rdt
    real(rp)    :: sf(3), d_on(3), k_f(3), k_rc(3), Uf(3), Uf_plain(3), Uf_eps(3), gm(3)
    logical     :: incr, skip_no
    integer(ip) :: col, idx2, n_loops, lo, hi
    logical     :: use_color
    integer(c_int), pointer :: skip(:), fown_p(:)
    integer(ip),    pointer :: fnb_p(:), fmi(:), perm(:), nper(:), per(:), fown(:), fnbm(:), fpid(:), pbc(:)
    integer(ip),    pointer :: ppair(:), pnb(:)
    real(c_double), pointer :: fcoef(:)
    logical,        pointer :: act_int(:), act_bnd(:), built(:)
    real(rp),       pointer :: k_int(:,:), k_bnd(:,:), gcs(:), eps(:,:), skew(:), fsf(:,:), fdcf(:), pdcf(:)
    real(rp),       pointer :: fctr(:,:), cctr(:,:), U(:,:), gU(:,:,:), Ub(:,:), gprc(:,:), phi(:), gp(:,:), ff(:)

    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nf = self%psolver%n_faces
    n_int = self%mesh%n_int_faces
    n_faces = self%mesh%n_faces
    n_bnd = self%mesh%n_bnd_faces
    n_patches = self%mesh%n_patches
    n_nper = self%n_nonperiodic_bfaces
    n_per = self%n_periodic_bfaces
    skip_no = self%skip_nonortho

    ! The cross term must lag grad(p^{k+1/2}); after the first corrector self%p is only
    ! the increment, so the total pressure has to be rebuilt before differentiating.
    incr = .false.
    if (present(p_is_increment)) incr = p_is_increment
    if (.not. self%skip_grad_p) then
      if (incr) then
        call build_total_pressure(self)
        call gradient_of_scalar(self, self%w_p_tot, self%w_p_tot_bnd, self%w_gp)
      else
        call gradient_of_scalar(self, self%p, self%p_bnd, self%w_gp)
      end if
    end if

    skip => self%psolver%skip_face
    fown_p => self%psolver%face_owner
    fnb_p => self%psolver%face_neighbor_ext
    fmi => self%psolver%face_mesh_idx
    fcoef => self%psolver%face_coeff
    act_int => self%nonortho_active_int
    k_int => self%nonortho_k_int
    act_bnd => self%nonortho_active_bnd
    k_bnd => self%nonortho_k_bnd
    nper => self%nonperiodic_bfaces
    per => self%periodic_bfaces
    gcs => self%mesh%face_gc_soa
    eps => self%mesh%face_eps_soa
    skew => self%mesh%face_skew_soa
    fsf => self%mesh%face_sf_soa
    fdcf => self%mesh%face_dcf_soa
    pdcf => self%mesh%periodic_pair_dcf_soa
    ppair => self%mesh%periodic_pair_soa
    pnb => self%mesh%periodic_neighbor_soa
    fown => self%mesh%face_owner_soa
    fnbm => self%mesh%face_neighbour_soa
    fpid => self%mesh%face_patch_id_soa
    pbc => self%mesh%patch_bc_id_soa
    fctr => self%mesh%face_centre_soa
    cctr => self%mesh%cell_centroid_soa
    U => self%U
    gU => self%grad_U
    Ub => self%U_bnd
    gprc => self%w_gp_rc
    phi => self%phi_prev
    gp => self%w_gp
    ff => self%face_flux
    built => self%w_flux_built

    use_color = (self%psolver%n_colors > 0 .and. allocated(self%psolver%idx_color_perm) &
                 .and. allocated(self%psolver%color_offsets))
    n_loops = merge(self%psolver%n_colors, 1_ip, use_color)
    if (use_color) then
      perm => self%psolver%idx_color_perm
    else
      perm => fmi  ! placeholder, unused when not colored (idx is taken directly)
    end if

    !$acc data copyin(U, gU, Ub, gprc, phi, gp) copyout(ff, built, rhs_p)
    !$acc parallel loop default(present)
    do c = 1, nc
      rhs_p(c) = 0.0_rp
    end do
    !$acc parallel loop default(present)
    do f = 1, n_faces
      ff(f) = 0.0_rp
      built(f) = .false.
    end do

    do col = 1, n_loops
      if (use_color) then
        lo = self%psolver%color_offsets(col)
        hi = self%psolver%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) default(present) &
      !$acc      private(idx, o, nb, f_mesh, om, nm, j, flux, coeff_rc, gc_f, d_o, d_nb, d_tot, w_blend, &
      !$acc              s_dot_dhat, d_on_len, sf, d_on, k_f, k_rc, Uf, Uf_plain, Uf_eps, gm, f_bnd, pair_bnd, pair_face)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(idx, idx2, o, nb, f_mesh, om, nm, j, flux, coeff_rc, gc_f, d_o, d_nb, d_tot, w_blend, &
      !$omp&         s_dot_dhat, d_on_len, sf, d_on, k_f, k_rc, Uf, Uf_plain, Uf_eps, gm, f_bnd, pair_bnd, pair_face)
#endif
      do idx2 = lo, hi
        if (use_color) then
          idx = perm(idx2)
        else
          idx = idx2
        end if
        if (skip(idx) /= 0) cycle

        o      = int(fown_p(idx), ip)
        nb     = fnb_p(idx)
        f_mesh = fmi(idx)

        if (o < 1 .or. o > nc) cycle
        if (nb < 1 .or. nb > ntot) cycle
        if (f_mesh < 1 .or. f_mesh > n_faces) cycle

        ! Face velocity interpolation using geometric gc (FP&S Sect. 7.3, Kim & Choi Eq. 25)
        ! gc = distance ratio (no volume weighting): gc = d_cf_nb / (d_cf_o + d_cf_nb)
        if (f_mesh <= n_int) then
          ! Internal faces: use precomputed geometric gc
          gc_f = gcs(f_mesh)
        else
          ! Periodic boundary faces: compute gc from d_cf distances, preferring the precomputed
          ! partner distance, else the partner face's own d_cf.
          gc_f = 0.5_rp
          f_bnd = f_mesh - n_int
          d_o  = fdcf(f_mesh)
          d_nb = d_o
          if (f_bnd >= 1 .and. f_bnd <= n_bnd) then
            if (pdcf(f_bnd) > SMALL) d_nb = pdcf(f_bnd)
            if (d_nb <= SMALL) then
              pair_bnd = ppair(f_bnd)
              if (pair_bnd >= 1) then
                pair_face = n_int + pair_bnd
                if (pair_face >= 1 .and. pair_face <= n_faces) d_nb = fdcf(pair_face)
              end if
            end if
          end if
          d_tot = d_o + d_nb
          if (d_tot > SMALL) gc_f = d_nb / d_tot
        end if
        d_on = cctr(:, nb) - cctr(:, o)
        sf = fsf(:, f_mesh)

        ! Face velocity interpolation: eps.grad correction smoothly blended off by face skew, plain
        ! gc for periodic (gradient discontinuity across the seam breaks KE conservation).
        Uf_plain = gc_f * U(:, o) + (1.0_rp - gc_f) * U(:, nb)
        if (f_mesh <= n_int) then
          ! Blend weight -> 1 on conforming faces, -> 0 near/above the refinement-skew tolerance.
          w_blend = refine_blend_weight(skew(f_mesh))
          ! face_interpolate_vector with grad_U on an interior face (mesh owner/neighbour, cached eps).
          om = fown(f_mesh)
          nm = fnbm(f_mesh)
          Uf_eps = gcs(f_mesh) * U(:, om) + (1.0_rp - gcs(f_mesh)) * U(:, nm)
          do j = 1, 3
            gm = gcs(f_mesh) * gU(:, j, om) + (1.0_rp - gcs(f_mesh)) * gU(:, j, nm)
            Uf_eps(j) = Uf_eps(j) + gm(1)*eps(1, f_mesh) + gm(2)*eps(2, f_mesh) + gm(3)*eps(3, f_mesh)
          end do
          Uf = Uf_plain + w_blend * (Uf_eps - Uf_plain)
        else
          Uf = Uf_plain
        end if
        flux = Uf(1)*sf(1) + Uf(2)*sf(2) + Uf(3)*sf(3)

        ! Rhie-Chow: the compact jump must carry the cross part too, or it fails to vanish for a linear
        ! pressure field on non-orthogonal faces. gc_f unclamped for the pressure-gradient weight.
        coeff_rc = real(fcoef(idx), rp)
        k_rc = 0.0_rp
        if (.not. skip_no) then
          call nonortho_kfactor_overrelaxed(sf, d_on, k_f, s_dot_dhat, d_on_len)
          if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(sf**2)) k_rc = k_f
        end if
        flux = flux + dt * nonortho_rhie_chow(gprc(:, o), gprc(:, nb), gc_f, sf, k_rc, coeff_rc, &
                                              phi(nb) - phi(o))

        ff(f_mesh) = flux
        built(f_mesh) = .true.

        ! rhs = -(1/dt)*sum_f q_f : outward flux from owner -> subtract from owner,
        ! outward from owner = inward to nb -> add to nb.
        rhs_p(o) = rhs_p(o) - flux
        if (nb <= nc) rhs_p(nb) = rhs_p(nb) + flux
      end do
#ifndef _OPENACC
      !$omp end parallel do
#endif
    end do

    ! Non-periodic boundary faces contribute only to their owner cell (atomics: faces share owners).
    !$acc parallel loop default(present) private(f_bnd, f_glob, o, q, flux)
    do ib = 1, n_nper
      f_bnd = nper(ib)
      f_glob = n_int + f_bnd
      o = fown(f_glob)
      if (o < 1 .or. o > nc) cycle
      q = fpid(f_glob)
      if (q < 1 .or. q > n_patches) cycle
      select case (pbc(q))
      case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
        flux = 0.0_rp
      case default
        ! BC_OUTLET: extrapolated boundary velocity (O(h^2) on skewed meshes); BC_INLET/BC_ABL_INLET:
        ! prescribed Dirichlet velocity, so include its flux in the RHS.
        flux = Ub(1, f_bnd)*fsf(1, f_glob) + Ub(2, f_bnd)*fsf(2, f_glob) + Ub(3, f_bnd)*fsf(3, f_glob)
      end select
      ff(f_glob) = flux
      built(f_glob) = .true.
      !$acc atomic update
      rhs_p(o) = rhs_p(o) - flux
    end do

    ! Periodic boundary faces that are not explicitly represented in the pressure-face list must
    ! still receive a valid convective mass flux. Prefer paired-face anti-symmetry when the companion
    ! was built, else the midpoint flux from the resolved partner cell; for an unbuilt pair both
    ! routes give the same value (the partner's midpoint flux is exactly the negation), so the
    ! faces can be done in any order.
    !$acc parallel loop default(present) private(f_bnd, f_glob, o, nb, pair_bnd, pair_face)
    do ib = 1, n_per
      f_bnd = per(ib)
      f_glob = n_int + f_bnd
      if (f_glob < 1 .or. f_glob > n_faces) cycle
      if (built(f_glob)) cycle
      o = fown(f_glob)
      if (o < 1 .or. o > nc) cycle

      pair_bnd = ppair(f_bnd)
      if (pair_bnd >= 1 .and. pair_bnd <= n_bnd) then
        pair_face = n_int + pair_bnd
        if (built(pair_face)) then
          ff(f_glob) = -ff(pair_face)
          built(f_glob) = .true.
          cycle
        end if
      end if

      ! Fallback: midpoint flux from resolved periodic partner cell.
      nb = pnb(f_bnd)
      if (nb <= 0 .and. pair_bnd >= 1 .and. pair_bnd <= n_bnd) nb = fown(n_int + pair_bnd)
      if (nb >= 1 .and. nb <= ntot) then
        ff(f_glob) = 0.5_rp * ((U(1, o) + U(1, nb))*fsf(1, f_glob) + (U(2, o) + U(2, nb))*fsf(2, f_glob) &
                             + (U(3, o) + U(3, nb))*fsf(3, f_glob))
      else
        ff(f_glob) = U(1, o)*fsf(1, f_glob) + U(2, o)*fsf(2, f_glob) + U(3, o)*fsf(3, f_glob)
      end if
      built(f_glob) = .true.
    end do

    ! Divide by dt to complete rhs = -(1/dt)*sum_f q_f.
    rdt = 1.0_rp / max(dt, SMALL)
    !$acc parallel loop default(present)
    do c = 1, nc
      rhs_p(c) = rhs_p(c) * rdt
    end do

    ! correct_face_fluxes must reuse this same lagged w_gp for the divergence to cancel.
    ! skip_nonortho means every face below fails its tolerance test, so both loops are provably no-ops.
    if (.not. skip_no) then
      ! use_color/n_loops carried over from the main flux loop above (same psolver idx-space).
      do col = 1, n_loops
        if (use_color) then
          lo = self%psolver%color_offsets(col)
          hi = self%psolver%color_offsets(col+1) - 1_ip
        else
          lo = 1_ip
          hi = nf
        end if
#ifdef _OPENACC
        !$acc parallel loop if(use_color) default(present) private(idx, o, nb, f_mesh, cross)
#else
        !$omp parallel do if(use_color) default(shared) private(idx, idx2, o, nb, f_mesh, cross)
#endif
        do idx2 = lo, hi
          if (use_color) then
            idx = perm(idx2)
          else
            idx = idx2
          end if
          if (skip(idx) /= 0) cycle
          o      = int(fown_p(idx), ip)
          nb     = fnb_p(idx)
          f_mesh = fmi(idx)
          if (o  < 1 .or. o  > nc) cycle
          ! Allow ghost neighbours: the matrix includes processor faces, so the
          ! deferred-correction RHS must include them too.
          if (nb < 1 .or. nb > ntot) cycle
          if (f_mesh < 1 .or. f_mesh > n_faces) cycle
          if (.not. act_int(idx)) cycle

          ! Use mesh-precomputed gc (consistent with gradient and viscous operators).
          cross = nonortho_cross_flux(gp(:, o), gp(:, nb), gcs(f_mesh), k_int(:, idx))
          rhs_p(o)  = rhs_p(o)  + cross
          if (nb <= nc) rhs_p(nb) = rhs_p(nb) - cross
        end do
#ifndef _OPENACC
        !$omp end parallel do
#endif
      end do

      ! Jasak deferred non-orthogonal correction on Dirichlet pressure boundary faces (BC_OUTLET),
      ! which are not in psolver%n_faces. Owner-cell gradient is the best available estimate of
      ! grad_p at the face (no off-cell information across the BC).
      !$acc parallel loop default(present) private(f_bnd, f_glob, o, q, cross)
      do ib = 1, n_nper
        f_bnd = nper(ib)
        f_glob = n_int + f_bnd
        if (f_glob < 1 .or. f_glob > n_faces) cycle
        o = fown(f_glob)
        if (o < 1 .or. o > nc) cycle
        q = fpid(f_glob)
        if (q < 1 .or. q > n_patches) cycle
        if (pbc(q) /= BC_OUTLET) cycle
        if (.not. act_bnd(ib)) cycle
        cross = gp(1, o)*k_bnd(1, ib) + gp(2, o)*k_bnd(2, ib) + gp(3, o)*k_bnd(3, ib)
        !$acc atomic update
        rhs_p(o) = rhs_p(o) + cross
      end do
    end if
    !$acc end data
  end subroutine compute_face_fluxes_and_rhs

  !> Correct previously assembled face fluxes by pressure jump using the
  !> same coefficients used in the pressure operator.
  subroutine correct_face_fluxes(self, dt)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(in)    :: dt

    integer(ip) :: idx, o, nb, f_mesh, nc, ntot, nf, n_int, n_faces, n_patches, n_nper, n_per
    integer(ip) :: ib, f_bnd, f_glob, pair_bnd, pair_face, o_self, o_pair, q
    real(rp)    :: coeff, d_orth, sn, gc, cross, d1, d2, d3
    logical     :: have_pair
    integer(c_int), pointer :: skip(:), fown_p(:)
    integer(ip),    pointer :: fnb_p(:), fmi(:), nper(:), per(:), fown(:), fpid(:), pbc(:), ppair(:)
    real(c_double), pointer :: fcoef(:)
    logical,        pointer :: act_int(:), act_bnd(:)
    integer(lp),    pointer :: gid(:)
    real(rp),       pointer :: k_int(:,:), k_bnd(:,:), gcs(:), ff(:), p(:), gp(:,:), pval(:)
    real(rp),       pointer :: fctr(:,:), fnrm(:,:), farea(:), cctr(:,:)

    if (.not. allocated(self%face_flux)) return

    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nf = self%psolver%n_faces
    n_int = self%mesh%n_int_faces
    n_faces = self%mesh%n_faces
    n_patches = self%mesh%n_patches
    n_nper = self%n_nonperiodic_bfaces
    n_per = self%n_periodic_bfaces
    skip => self%psolver%skip_face
    fown_p => self%psolver%face_owner
    fnb_p => self%psolver%face_neighbor_ext
    fmi => self%psolver%face_mesh_idx
    fcoef => self%psolver%face_coeff
    act_int => self%nonortho_active_int
    k_int => self%nonortho_k_int
    act_bnd => self%nonortho_active_bnd
    k_bnd => self%nonortho_k_bnd
    nper => self%nonperiodic_bfaces
    per => self%periodic_bfaces
    pval => self%patch_p_val
    gcs => self%mesh%face_gc_soa
    fown => self%mesh%face_owner_soa
    fpid => self%mesh%face_patch_id_soa
    pbc => self%mesh%patch_bc_id_soa
    fctr => self%mesh%face_centre_soa
    fnrm => self%mesh%face_normal_soa
    farea => self%mesh%face_area_soa
    cctr => self%mesh%cell_centroid_soa
    ppair => self%mesh%periodic_pair_soa
    gid => self%mesh%global_id
    ff => self%face_flux
    p => self%p
    gp => self%w_gp
    have_pair = allocated(self%mesh%periodic_pair)

    !$acc data copyin(p, gp) copy(ff)
    ! idx <-> f_mesh is a bijection over the assembled faces: each idx writes a distinct
    ! face_flux(f_mesh), so unlike the rhs_p scatter-adds this loop is race-free as-is.
#ifdef _OPENACC
    !$acc parallel loop default(present) private(o, nb, f_mesh, coeff, gc, cross)
#else
    !$omp parallel do private(idx, o, nb, f_mesh, coeff, gc, cross) schedule(static)
#endif
    do idx = 1, nf
      if (skip(idx) /= 0) cycle
      o = int(fown_p(idx), ip)
      nb = fnb_p(idx)
      if (o < 1 .or. o > nc) cycle
      if (nb < 1 .or. nb > ntot) cycle
      f_mesh = fmi(idx)
      if (f_mesh < 1 .or. f_mesh > n_faces) cycle
      coeff = real(fcoef(idx), rp)
      ff(f_mesh) = ff(f_mesh) - dt * coeff * (p(nb) - p(o))
      ! Jasak (1996) Sect. 6.4.3: cross-diffusion using grad(p^old) from compute_face_fluxes_and_rhs.
      ! Same lagged gradient as in rhs_p -> algebraic cancellation (div = 0).
      if (act_int(idx)) then
        ! Use precomputed projection-based gc (FPS Sect. 7.2) - consistent with
        ! the gradient, advection, and viscous operators.
        gc = gcs(f_mesh)
        cross = nonortho_cross_flux(gp(:, o), gp(:, nb), gc, k_int(:, idx))
        ff(f_mesh) = ff(f_mesh) - dt * cross
      end if
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif

    ! Project non-periodic boundary face fluxes with the Dirichlet/Neumann p BC.
    !$acc parallel loop default(present) private(f_bnd, f_glob, o, q, sn, d_orth, coeff, d1, d2, d3)
    do ib = 1, n_nper
      f_bnd  = nper(ib)
      f_glob = n_int + f_bnd
      if (f_glob < 1 .or. f_glob > n_faces) cycle
      o = fown(f_glob)
      if (o < 1 .or. o > nc) cycle
      q = fpid(f_glob)
      if (q < 1 .or. q > n_patches) cycle
      if (pbc(q) /= BC_OUTLET) cycle
      d1 = fctr(1, f_glob) - cctr(1, o)
      d2 = fctr(2, f_glob) - cctr(2, o)
      d3 = fctr(3, f_glob) - cctr(3, o)
      sn = abs(d1*fnrm(1, f_glob) + d2*fnrm(2, f_glob) + d3*fnrm(3, f_glob))
      d_orth = max(sn, SMALL)
      coeff = farea(f_glob) / d_orth
      ff(f_glob) = ff(f_glob) - dt * coeff * (pval(q) - p(o))
      ! Cross-diffusion correction for Dirichlet p BC (outlet); owner-side gradient only (no off-cell info).
      if (act_bnd(ib)) then
        ff(f_glob) = ff(f_glob) - dt * (gp(1, o)*k_bnd(1, ib) + gp(2, o)*k_bnd(2, ib) + gp(3, o)*k_bnd(3, ib))
      end if
    end do

    ! Re-sync non-canonical same-rank periodic faces by anti-symmetry. Only non-canonical faces
    ! are written and only canonical ones read, so the faces are independent.
    if (n_per > 0 .and. have_pair) then
      !$acc parallel loop default(present) private(f_bnd, pair_bnd, f_glob, pair_face, o_self, o_pair)
      do ib = 1, n_per
        f_bnd = per(ib)
        if (f_bnd < 1 .or. f_bnd > size(ppair)) cycle
        pair_bnd = ppair(f_bnd)
        if (pair_bnd <= 0) cycle  ! cross-rank pair: handled by remote rank
        f_glob    = n_int + f_bnd
        pair_face = n_int + pair_bnd
        if (f_glob    < 1 .or. f_glob    > n_faces) cycle
        if (pair_face < 1 .or. pair_face > n_faces) cycle
        o_self = fown(f_glob)
        o_pair = fown(pair_face)
        if (o_self < 1 .or. o_self > nc) cycle
        if (o_pair < 1 .or. o_pair > nc) cycle
        if (o_self == o_pair) cycle
        ! Canonical face is the one whose owner has the smaller global_id;
        ! mirror the non-canonical face flux from it.
        if (gid(o_self) > gid(o_pair)) ff(f_glob) = -ff(pair_face)
      end do
    end if
    !$acc end data
  end subroutine correct_face_fluxes

  !> Correct cell velocities: u = u - dt * T^-1 . sum_f n_f g_f, g_f the jump the flux correction removes.
  ! with_cross = .false. for the pre-correction: Rhie-Chow replaces only the orthogonal part at the face.
  subroutine correct_velocity_psolver_dual(self, dt, with_cross)
    class(FVSolver), intent(inout), target :: self
    real(rp),        intent(in)    :: dt
    logical,         intent(in)    :: with_cross

    integer(ip) :: nc, ntot, c, idx, o, nb, f_mesh, nf, n_int, n_faces, n_patches, n_nper
    integer(ip) :: ib, f_bnd, f_glob, q
    real(rp)    :: coeff, g_f, gc, sn, d_orth, n1, n2, n3, g1, g2, g3, d1, d2, d3
    logical     :: do_cross
    integer(ip) :: col, idx2, n_loops, lo, hi
    logical     :: use_color
    integer(c_int), pointer :: skip(:), fown_p(:)
    integer(ip),    pointer :: fnb_p(:), fmi(:), perm(:), nper(:), fown(:), fpid(:), pbc(:)
    real(c_double), pointer :: fcoef(:)
    logical,        pointer :: act_int(:), act_bnd(:)
    real(rp),       pointer :: k_int(:,:), k_bnd(:,:), gcs(:), p(:), gp(:,:), gpc(:,:), U(:,:), Tinv(:,:,:)
    real(rp),       pointer :: pval(:), fctr(:,:), fnrm(:,:), farea(:), cctr(:,:)

    do_cross = with_cross .and. (.not. self%skip_nonortho)
    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nf = self%psolver%n_faces
    n_int = self%mesh%n_int_faces
    n_faces = self%mesh%n_faces
    n_patches = self%mesh%n_patches
    n_nper = self%n_nonperiodic_bfaces
    skip => self%psolver%skip_face
    fown_p => self%psolver%face_owner
    fnb_p => self%psolver%face_neighbor_ext
    fmi => self%psolver%face_mesh_idx
    fcoef => self%psolver%face_coeff
    act_int => self%nonortho_active_int
    k_int => self%nonortho_k_int
    act_bnd => self%nonortho_active_bnd
    k_bnd => self%nonortho_k_bnd
    nper => self%nonperiodic_bfaces
    pval => self%patch_p_val
    Tinv => self%recon_Tinv
    gcs => self%mesh%face_gc_soa
    fown => self%mesh%face_owner_soa
    fpid => self%mesh%face_patch_id_soa
    pbc => self%mesh%patch_bc_id_soa
    fctr => self%mesh%face_centre_soa
    fnrm => self%mesh%face_normal_soa
    farea => self%mesh%face_area_soa
    cctr => self%mesh%cell_centroid_soa
    p => self%p
    gp => self%w_gp
    gpc => self%w_gp_cell
    U => self%U

    use_color = (self%psolver%n_colors > 0 .and. allocated(self%psolver%idx_color_perm) &
                 .and. allocated(self%psolver%color_offsets))
    n_loops = merge(self%psolver%n_colors, 1_ip, use_color)
    if (use_color) then
      perm => self%psolver%idx_color_perm
    else
      perm => fmi  ! placeholder, unused when not colored (idx is taken directly)
    end if

    !$acc data copyin(p, gp) copy(gpc, U)
    !$acc parallel loop default(present)
    do c = 1, nc
      gpc(1, c) = 0.0_rp
      gpc(2, c) = 0.0_rp
      gpc(3, c) = 0.0_rp
    end do

    ! g_f is identical to the pressure jump correct_face_fluxes removes from the
    ! mass flux, so the cell and face projections share one discrete operator.
    do col = 1, n_loops
      if (use_color) then
        lo = self%psolver%color_offsets(col)
        hi = self%psolver%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) default(present) private(idx, o, nb, f_mesh, coeff, g_f, gc, n1, n2, n3)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(idx, idx2, o, nb, f_mesh, coeff, g_f, gc, n1, n2, n3)
#endif
      do idx2 = lo, hi
        if (use_color) then
          idx = perm(idx2)
        else
          idx = idx2
        end if
        if (skip(idx) /= 0) cycle
        o      = int(fown_p(idx), ip)
        nb     = fnb_p(idx)
        f_mesh = fmi(idx)
        if (o      < 1 .or. o      > nc) cycle
        if (nb     < 1 .or. nb     > ntot) cycle
        if (f_mesh < 1 .or. f_mesh > n_faces) cycle

        coeff = real(fcoef(idx), rp)
        g_f   = coeff * (p(nb) - p(o))
        if (do_cross) then
          if (act_int(idx)) then
            gc  = gcs(f_mesh)
            g_f = g_f + nonortho_cross_flux(gp(:, o), gp(:, nb), gc, k_int(:, idx))
          end if
        end if

        ! Sf points out of the owner, so the neighbour sees (-n_f)*(-g_f): same sign.
        n1 = fnrm(1, f_mesh)
        n2 = fnrm(2, f_mesh)
        n3 = fnrm(3, f_mesh)
        gpc(1, o) = gpc(1, o) + g_f * n1
        gpc(2, o) = gpc(2, o) + g_f * n2
        gpc(3, o) = gpc(3, o) + g_f * n3
        if (nb <= nc) then
          gpc(1, nb) = gpc(1, nb) + g_f * n1
          gpc(2, nb) = gpc(2, nb) + g_f * n2
          gpc(3, nb) = gpc(3, nb) + g_f * n3
        end if
      end do
#ifndef _OPENACC
      !$omp end parallel do
#endif
    end do

    ! Fixed-flux patches (walls, symmetry, inlets) get no flux correction, so their
    ! normal pressure gradient is zero; only Dirichlet-p outlets contribute. Outlet faces
    ! can share an owner cell, hence the atomics.
    !$acc parallel loop default(present) private(f_bnd, f_glob, o, q, sn, d_orth, coeff, g_f, d1, d2, d3, g1, g2, g3)
    do ib = 1, n_nper
      f_bnd  = nper(ib)
      f_glob = n_int + f_bnd
      if (f_glob < 1 .or. f_glob > n_faces) cycle
      o = fown(f_glob)
      if (o < 1 .or. o > nc) cycle
      q = fpid(f_glob)
      if (q < 1 .or. q > n_patches) cycle
      if (pbc(q) /= BC_OUTLET) cycle

      d1 = fctr(1, f_glob) - cctr(1, o)
      d2 = fctr(2, f_glob) - cctr(2, o)
      d3 = fctr(3, f_glob) - cctr(3, o)
      sn     = abs(d1*fnrm(1, f_glob) + d2*fnrm(2, f_glob) + d3*fnrm(3, f_glob))
      d_orth = max(sn, SMALL)
      coeff  = farea(f_glob) / d_orth
      g_f    = coeff * (pval(q) - p(o))
      if (do_cross) then
        if (act_bnd(ib)) g_f = g_f + gp(1, o)*k_bnd(1, ib) + gp(2, o)*k_bnd(2, ib) + gp(3, o)*k_bnd(3, ib)
      end if
      g1 = g_f * fnrm(1, f_glob)
      g2 = g_f * fnrm(2, f_glob)
      g3 = g_f * fnrm(3, f_glob)
      !$acc atomic update
      gpc(1, o) = gpc(1, o) + g1
      !$acc atomic update
      gpc(2, o) = gpc(2, o) + g2
      !$acc atomic update
      gpc(3, o) = gpc(3, o) + g3
    end do

#ifdef _OPENACC
    !$acc parallel loop default(present) private(g1, g2, g3)
#else
    !$omp parallel do schedule(static) private(g1, g2, g3)
#endif
    do c = 1, nc
      g1 = Tinv(1, 1, c)*gpc(1, c) + Tinv(1, 2, c)*gpc(2, c) + Tinv(1, 3, c)*gpc(3, c)
      g2 = Tinv(2, 1, c)*gpc(1, c) + Tinv(2, 2, c)*gpc(2, c) + Tinv(2, 3, c)*gpc(3, c)
      g3 = Tinv(3, 1, c)*gpc(1, c) + Tinv(3, 2, c)*gpc(2, c) + Tinv(3, 3, c)*gpc(3, c)
      gpc(1, c) = g1
      gpc(2, c) = g2
      gpc(3, c) = g3
      U(1, c) = U(1, c) - dt * g1
      U(2, c) = U(2, c) - dt * g2
      U(3, c) = U(3, c) - dt * g3
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif
    !$acc end data
  end subroutine correct_velocity_psolver_dual

  !> Move every field and work array of the time step to the device (OpenACC builds; a no-op
  !> otherwise). Call once, after initial/restart conditions are set and before the time loop:
  !> from here on the device copy is authoritative and the host copy is refreshed only by
  !> sync_host (output, statistics, sampling, checkpoints).
  subroutine solver_to_device(self, gpu_aware_mpi)
    class(FVSolver), intent(inout) :: self
    logical,         intent(in)    :: gpu_aware_mpi
#ifdef _OPENACC
    if (self%on_device) return
    if (self%params%p_max_iter <= 0) &
      error stop 'solver_to_device: OpenACC builds need the pressure solver (p_max_iter > 0)'
    if (self%mesh%n_colors <= 0 .or. self%psolver%n_colors <= 0) &
      error stop 'solver_to_device: face colouring is missing (a rank with no faces?); the device kernels need it'
    call stage_static_on_device(self)
    !$acc enter data copyin(self%U, self%p, self%p_old, self%T, self%nu_sgs, self%U_bnd, self%p_bnd, &
    !$acc      self%phi_prev, self%phi_prev_bnd, self%T_bnd, self%tau_w_bnd, self%sem_fluc_bnd, &
    !$acc      self%grad_U, self%grad_T, self%face_flux, self%w_H_U, self%w_H_T, self%w_H_C, &
    !$acc      self%w_rhs_p, self%w_U0, self%w_T0, self%w_C0, self%w_conv_U, self%w_visc_U, &
    !$acc      self%w_conv_T, self%w_visc_T, self%w_flux_built, self%w_gp, self%w_gp_rc, &
    !$acc      self%w_gp_cell, self%w_p_tot, self%w_p_tot_bnd, self%w_outlet_cell, self%w_outlet_pval)
    !$acc update device(self%U_face, self%inv_vol)
    if (self%params%n_scalars > 0) then
      !$acc enter data copyin(self%C, self%C_bnd, self%grad_C, self%w_conv_C, self%w_visc_C, &
      !$acc      self%w_Cs, self%w_Cs_bnd)
    end if
    if (allocated(self%w_pack_grad)) then
      !$acc enter data copyin(self%w_pack_grad)
    end if
    if (allocated(self%w_g_prev_U)) then
      !$acc enter data copyin(self%w_g_prev_U, self%w_g_prev_T, self%w_gp_prev)
    end if
    if (allocated(self%abl_mean_bnd)) then
      !$acc enter data copyin(self%abl_mean_bnd)
    end if
    if (allocated(self%q_w_bnd)) then
      !$acc enter data copyin(self%q_w_bnd, self%L_bnd)
    end if
    call self%mesh%set_halo_device(.true., gpu_aware_mpi)
    self%on_device = .true.
#endif
  end subroutine solver_to_device

  !> Refresh the host copy of everything output, statistics, sampling and checkpoints read.
  subroutine solver_sync_host(self)
    class(FVSolver), intent(inout) :: self
#ifdef _OPENACC
    if (.not. self%on_device) return
    !$acc update host(self%U, self%p, self%T, self%nu_sgs, self%U_bnd, self%p_bnd, self%T_bnd, &
    !$acc      self%tau_w_bnd, self%grad_U, self%grad_T, self%face_flux, self%U_face, self%w_gp)
    if (self%params%n_scalars > 0) then
      !$acc update host(self%C, self%C_bnd, self%grad_C)
    end if
    if (allocated(self%q_w_bnd)) then
      !$acc update host(self%q_w_bnd, self%L_bnd)
    end if
#endif
  end subroutine solver_sync_host

  !> Release what solver_to_device staged.
  subroutine solver_leave_device(self)
    class(FVSolver), intent(inout) :: self
#ifdef _OPENACC
    if (.not. self%on_device) return
    !$acc exit data delete(self%U, self%p, self%p_old, self%T, self%nu_sgs, self%U_bnd, self%p_bnd, &
    !$acc      self%phi_prev, self%phi_prev_bnd, self%T_bnd, self%tau_w_bnd, self%sem_fluc_bnd, &
    !$acc      self%grad_U, self%grad_T, self%face_flux, self%w_H_U, self%w_H_T, self%w_H_C, &
    !$acc      self%w_rhs_p, self%w_U0, self%w_T0, self%w_C0, self%w_conv_U, self%w_visc_U, &
    !$acc      self%w_conv_T, self%w_visc_T, self%w_flux_built, self%w_gp, self%w_gp_rc, &
    !$acc      self%w_gp_cell, self%w_p_tot, self%w_p_tot_bnd, self%w_outlet_cell, self%w_outlet_pval)
    if (self%params%n_scalars > 0) then
      !$acc exit data delete(self%C, self%C_bnd, self%grad_C, self%w_conv_C, self%w_visc_C, &
      !$acc      self%w_Cs, self%w_Cs_bnd)
    end if
    if (allocated(self%w_pack_grad)) then
      !$acc exit data delete(self%w_pack_grad)
    end if
    if (allocated(self%w_g_prev_U)) then
      !$acc exit data delete(self%w_g_prev_U, self%w_g_prev_T, self%w_gp_prev)
    end if
    if (allocated(self%abl_mean_bnd)) then
      !$acc exit data delete(self%abl_mean_bnd)
    end if
    if (allocated(self%q_w_bnd)) then
      !$acc exit data delete(self%q_w_bnd, self%L_bnd)
    end if
    self%on_device = .false.
#endif
  end subroutine solver_leave_device

  !> Stage the static projection-operator arrays on the device once (built in solver_init, never
  !> modified afterwards), so the projection kernels find them present.
  subroutine stage_static_on_device(self)
    class(FVSolver), intent(inout) :: self
#ifdef _OPENACC
    if (self%static_on_device) return
    !$acc enter data copyin(self%psolver%skip_face, self%psolver%face_owner, self%psolver%face_neighbor_ext, &
    !$acc      self%psolver%face_mesh_idx, self%psolver%face_coeff, &
    !$acc      self%nonortho_k_int, self%nonortho_active_int, self%nonortho_k_bnd, self%nonortho_active_bnd, &
    !$acc      self%recon_Tinv, self%nonperiodic_bfaces, self%periodic_bfaces, self%patch_p_val, &
    !$acc      self%mesh%global_id)
    if (allocated(self%psolver%idx_color_perm)) then
      !$acc enter data copyin(self%psolver%idx_color_perm, self%psolver%color_offsets)
    end if
    self%static_on_device = .true.
#endif
  end subroutine stage_static_on_device

  subroutine release_static_on_device(self)
    class(FVSolver), intent(inout) :: self
#ifdef _OPENACC
    if (.not. self%static_on_device) return
    !$acc exit data delete(self%psolver%skip_face, self%psolver%face_owner, self%psolver%face_neighbor_ext, &
    !$acc      self%psolver%face_mesh_idx, self%psolver%face_coeff, &
    !$acc      self%nonortho_k_int, self%nonortho_active_int, self%nonortho_k_bnd, self%nonortho_active_bnd, &
    !$acc      self%recon_Tinv, self%nonperiodic_bfaces, self%periodic_bfaces, self%patch_p_val, &
    !$acc      self%mesh%global_id)
    if (allocated(self%psolver%idx_color_perm)) then
      !$acc exit data delete(self%psolver%idx_color_perm, self%psolver%color_offsets)
    end if
    self%static_on_device = .false.
#endif
  end subroutine release_static_on_device

  !> Build T^-1, T = sum_f |Sf| n_f (x) n_f: recovers a cell vector from its face-normal components.
  subroutine build_reconstruct_operator(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: nc, c, i, j, idx, o, nb, f_mesh, ib, f_bnd, f_glob
    real(rp)    :: area, n_hat(3)
    real(rp), allocatable :: T(:,:,:)

    nc = self%mesh%n_cells
    allocate(T(3, 3, nc), source=0.0_rp)

    do idx = 1, self%psolver%n_faces
      if (self%psolver%skip_face(idx) /= 0) cycle
      o      = int(self%psolver%face_owner(idx), ip)
      nb     = self%psolver%face_neighbor_ext(idx)
      f_mesh = self%psolver%face_mesh_idx(idx)
      if (o      < 1 .or. o      > nc) cycle
      if (nb     < 1 .or. nb     > nc + self%mesh%n_ghost) cycle
      if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle

      area  = self%mesh%faces(f_mesh)%area
      n_hat = self%mesh%faces(f_mesh)%normal
      ! n (x) n is sign-independent, and a periodic partner face has the same area.
      do j = 1, 3
        do i = 1, 3
          T(i, j, o) = T(i, j, o) + area * n_hat(i) * n_hat(j)
          if (nb <= nc) T(i, j, nb) = T(i, j, nb) + area * n_hat(i) * n_hat(j)
        end do
      end do
    end do

    ! Boundary faces close the cell: they carry a zero normal gradient on fixed-flux
    ! patches and the Dirichlet jump on outlets, so both belong in the fit.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd  = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      area  = self%mesh%faces(f_glob)%area
      n_hat = self%mesh%faces(f_glob)%normal
      do j = 1, 3
        do i = 1, 3
          T(i, j, o) = T(i, j, o) + area * n_hat(i) * n_hat(j)
        end do
      end do
    end do

    do c = 1, nc
      call invert_sym3(T(:, :, c), self%recon_Tinv(:, :, c))
    end do

    deallocate(T)
  end subroutine build_reconstruct_operator

  !> Invert a symmetric positive-definite 3x3; falls back to the isotropic part if degenerate.
  pure subroutine invert_sym3(A, Ainv)
    real(rp), intent(in)  :: A(3,3)
    real(rp), intent(out) :: Ainv(3,3)

    real(rp) :: c11, c12, c13, c22, c23, c33, det, tr
    integer(ip) :: i

    c11 = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    c12 = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    c13 = A(1,2)*A(2,3) - A(1,3)*A(2,2)
    c22 = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    c23 = A(1,3)*A(2,1) - A(1,1)*A(2,3)
    c33 = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    det = A(1,1)*c11 + A(1,2)*(A(2,3)*A(3,1) - A(2,1)*A(3,3)) + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

    tr = A(1,1) + A(2,2) + A(3,3)
    if (abs(det) <= 1.0e-12_rp * max(tr, SMALL)**3) then
      Ainv = 0.0_rp
      do i = 1, 3
        Ainv(i,i) = 3.0_rp / max(tr, SMALL)
      end do
      return
    end if

    Ainv(1,1) = c11 / det;  Ainv(1,2) = c12 / det;  Ainv(1,3) = c13 / det
    Ainv(2,1) = c12 / det;  Ainv(2,2) = c22 / det;  Ainv(2,3) = c23 / det
    Ainv(3,1) = c13 / det;  Ainv(3,2) = c23 / det;  Ainv(3,3) = c33 / det
  end subroutine invert_sym3

  !> Closed-form eigenvalues of a symmetric 3x3 matrix (Smith 1961), ascending order.
  pure function eig3_sym(A) result(eigs)
    real(rp), intent(in) :: A(3,3)
    real(rp) :: eigs(3)

    real(rp) :: q, p1, p2, p, r, phi, B(3,3)
    integer(ip) :: i, j

    p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
    if (p1 <= SMALL) then
      eigs = [A(1,1), A(2,2), A(3,3)]
    else
      q  = (A(1,1) + A(2,2) + A(3,3)) / 3.0_rp
      p2 = (A(1,1)-q)**2 + (A(2,2)-q)**2 + (A(3,3)-q)**2 + 2.0_rp*p1
      p  = sqrt(max(p2, SMALL) / 6.0_rp)
      do j = 1, 3
        do i = 1, 3
          B(i,j) = A(i,j) / p
        end do
        B(j,j) = B(j,j) - q/p
      end do
      r = 0.5_rp * (B(1,1)*(B(2,2)*B(3,3) - B(2,3)*B(3,2)) &
                  - B(1,2)*(B(2,1)*B(3,3) - B(2,3)*B(3,1)) &
                  + B(1,3)*(B(2,1)*B(3,2) - B(2,2)*B(3,1)))
      r   = max(-1.0_rp, min(1.0_rp, r))
      phi = acos(r) / 3.0_rp
      eigs(3) = q + 2.0_rp*p*cos(phi)
      eigs(1) = q + 2.0_rp*p*cos(phi + 2.0_rp*3.14159265358979_rp/3.0_rp)
      eigs(2) = 3.0_rp*q - eigs(1) - eigs(3)
    end if
    ! Ascending sort (3 elements: fixed compare-swap network).
    if (eigs(1) > eigs(2)) eigs([1,2]) = eigs([2,1])
    if (eigs(2) > eigs(3)) eigs([2,3]) = eigs([3,2])
    if (eigs(1) > eigs(2)) eigs([1,2]) = eigs([2,1])
  end function eig3_sym

  !> Compute divergence from the stored face fluxes with conservative
  !> owner/neighbour accumulation, plus non-periodic boundary fluxes.
  subroutine face_flux_divergence(self, div_u)
    class(FVSolver), intent(in)  :: self
    real(rp),        intent(out) :: div_u(:)

    integer(ip) :: nc, q, ib, f_bnd, f_glob
    integer(ip) :: idx, o, nb, f_mesh
    real(rp)    :: flux

    nc = self%mesh%n_cells
    div_u = 0.0_rp

    if (allocated(self%face_flux)) then
      do idx = 1, self%psolver%n_faces
        if (self%psolver%skip_face(idx) /= 0) cycle
        o = int(self%psolver%face_owner(idx), ip)
        nb = self%psolver%face_neighbor_ext(idx)
        if (o < 1 .or. o > nc) cycle
        f_mesh = self%psolver%face_mesh_idx(idx)
        if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle
        flux = self%face_flux(f_mesh)

        div_u(o) = div_u(o) + flux * self%inv_vol(o)
        if (nb >= 1 .and. nb <= nc) div_u(nb) = div_u(nb) - flux * self%inv_vol(nb)
      end do
    end if

    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (allocated(self%face_flux)) then
        flux = self%face_flux(f_glob)
      else
        if (q < 1 .or. q > self%mesh%n_patches) cycle
        select case (self%bc_configs(q)%bc_id)
        case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
          flux = 0.0_rp
        case (BC_OUTLET)
          flux = 0.0_rp
        case default
          flux = dot_product(self%U_bnd(:, f_bnd), self%mesh%faces(f_glob)%sf)
        end select
      end if
      div_u(o) = div_u(o) + flux * self%inv_vol(o)
    end do
  end subroutine face_flux_divergence

  !> Green-Gauss divergence of the cell-centred velocity field (not the face fluxes).
  subroutine cell_velocity_divergence(self, div_u)
    class(FVSolver), intent(in)  :: self
    real(rp),        intent(out) :: div_u(:)

    integer(ip) :: nc, f, o, nb, f_bnd
    real(rp)    :: Uf(3), gc, flux

    nc = self%mesh%n_cells
    div_u = 0.0_rp

    do f = 1, self%mesh%n_faces
      o  = self%mesh%faces(f)%owner
      nb = self%mesh%faces(f)%neighbour
      if (nb >= 1) then
        gc = self%mesh%faces(f)%gc
        Uf = gc * self%U(:, o) + (1.0_rp - gc) * self%U(:, nb)
      else
        f_bnd = f - self%mesh%n_int_faces
        if (f_bnd < 1 .or. f_bnd > size(self%U_bnd, 2)) cycle
        Uf = self%U_bnd(:, f_bnd)
      end if
      flux = dot_product(Uf, self%mesh%faces(f)%sf)
      if (o  >= 1 .and. o  <= nc) div_u(o)  = div_u(o)  + flux * self%inv_vol(o)
      if (nb >= 1 .and. nb <= nc) div_u(nb) = div_u(nb) - flux * self%inv_vol(nb)
    end do
  end subroutine cell_velocity_divergence

  !> Shift pressure so the mean at outlet cells equals p_ref (MPI-consistent).
  subroutine anchor_pressure_to_outlet(p, msh, outlet_cell, outlet_pval, dev)
    real(rp),    intent(inout) :: p(:)  ! pressure array (n_cells)
    type(Mesh),  intent(in), target :: msh
    logical,     intent(in)    :: outlet_cell(:)  ! (n_cells) precomputed outlet mask
    real(rp),    intent(in)    :: outlet_pval(:)  ! (n_cells) reference pressure per outlet cell
    logical,     intent(in)    :: dev  ! p/outlet arrays are device-resident

    integer(ip) :: ocell, nc
    real(rp)    :: p_sum, p_ref_sum, p_ref, p_shift
    integer(ip) :: n_out, n_ref
    integer(ip) :: ierr
    real(rp)    :: buf(2)
    real(rp), pointer :: cvol(:)

    nc = msh%n_cells
    cvol => msh%cell_volume_soa
    p_sum = 0.0_rp;  n_out = 0
    p_ref_sum = 0.0_rp
    p_ref = 0.0_rp

    !$acc parallel loop if(dev) default(present) reduction(+:p_sum, p_ref_sum, n_out)
    do ocell = 1, nc
      if (outlet_cell(ocell)) then
        p_sum     = p_sum     + p(ocell)
        p_ref_sum = p_ref_sum + outlet_pval(ocell)
        n_out     = n_out     + 1
      end if
    end do

    ! Global sum across MPI ranks; n_out reduced as integer to avoid FP precision issues.
    buf(1) = p_sum
    buf(2) = p_ref_sum
    call MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, n_out, 1, MPI_INTEGER, &
                       MPI_SUM, msh%comm, ierr)
    p_sum = buf(1)
    p_ref_sum = buf(2)
    n_ref = n_out

    if (n_out > 0 .and. n_ref > 0) then
      p_ref = p_ref_sum / real(n_ref, rp)
      p_shift = p_sum / real(n_out, rp) - p_ref
    else
      ! pure-Neumann case (no outlet anywhere) - pressure is determined
      ! only up to a constant. Anchor by subtracting the global volume-weighted mean.
      block
        real(rp) :: vol_local, vol_global, pmean_local, pmean_global, gbuf(2)
        integer(ip) :: c
        vol_local   = 0.0_rp
        pmean_local = 0.0_rp
        !$acc parallel loop if(dev) default(present) reduction(+:vol_local, pmean_local)
        do c = 1, nc
          vol_local   = vol_local   + cvol(c)
          pmean_local = pmean_local + p(c) * cvol(c)
        end do
        gbuf(1) = vol_local
        gbuf(2) = pmean_local
        call MPI_Allreduce(MPI_IN_PLACE, gbuf, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, msh%comm, ierr)
        vol_global   = gbuf(1)
        pmean_global = gbuf(2)
        p_shift = 0.0_rp
        if (vol_global > 0.0_rp) p_shift = pmean_global / vol_global
      end block
    end if
    !$acc parallel loop if(dev) default(present)
    do ocell = 1, size(p)
      p(ocell) = p(ocell) - p_shift
    end do
  end subroutine anchor_pressure_to_outlet


  !> Face-flux divergence, plus optionally the cell-centred one (a different operator, not zeroed by the projection).
  subroutine solver_divergence_metrics(self, max_abs_div, mean_abs_div, &
                                       max_abs_div_cell, mean_abs_div_cell)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(out)   :: max_abs_div
    real(rp),        intent(out)   :: mean_abs_div
    real(rp), optional, intent(out) :: max_abs_div_cell
    real(rp), optional, intent(out) :: mean_abs_div_cell

    integer(ip) :: nc, c, ierr
    real(rp)    :: vol_sum_local, vol_sum_global
    real(rp)    :: abs_sum_local, abs_sum_global
    real(rp)    :: max_local, max_cell_global

    nc = self%mesh%n_cells
    self%w_rhs_p(1:nc) = 0.0_rp

    call self%mesh%exchange_halos_vec(self%U)

    call apply_velocity_bc(self%U, self%U_bnd, &
                 self%bc_configs, self%mesh, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)

    ! Use the corrected face fluxes stored by correct_face_fluxes during the
    ! last projection step. By F&P construction sum_f q_f = 0 exactly, so
    call face_flux_divergence(self, self%w_rhs_p(1:nc))

    max_local = maxval(abs(self%w_rhs_p(1:nc)))
    abs_sum_local = 0.0_rp
    vol_sum_local = 0.0_rp
    do c = 1, nc
      abs_sum_local = abs_sum_local + abs(self%w_rhs_p(c)) * self%mesh%cells(c)%volume
      vol_sum_local = vol_sum_local + self%mesh%cells(c)%volume
    end do

    call MPI_Allreduce(max_local, max_abs_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, self%mesh%comm, ierr)
    block
      real(rp) :: sum_local(2), sum_global(2)
      sum_local = [abs_sum_local, vol_sum_local]
      call MPI_Allreduce(sum_local, sum_global, 2, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, self%mesh%comm, ierr)
      abs_sum_global = sum_global(1)
      vol_sum_global = sum_global(2)
    end block

    if (vol_sum_global > SMALL) then
      mean_abs_div = abs_sum_global / vol_sum_global
    else
      mean_abs_div = 0.0_rp
    end if

    if (present(max_abs_div_cell) .or. present(mean_abs_div_cell)) then
      call cell_velocity_divergence(self, self%w_rhs_p(1:nc))
      max_local     = maxval(abs(self%w_rhs_p(1:nc)))
      abs_sum_local = 0.0_rp
      do c = 1, nc
        abs_sum_local = abs_sum_local + abs(self%w_rhs_p(c)) * self%mesh%cells(c)%volume
      end do
      call MPI_Allreduce(max_local, max_cell_global, 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, self%mesh%comm, ierr)
      if (present(max_abs_div_cell)) max_abs_div_cell = max_cell_global
      call MPI_Allreduce(abs_sum_local, abs_sum_global, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, self%mesh%comm, ierr)
      if (present(mean_abs_div_cell)) then
        if (vol_sum_global > SMALL) then
          mean_abs_div_cell = abs_sum_global / vol_sum_global
        else
          mean_abs_div_cell = 0.0_rp
        end if
      end if
    end if

  end subroutine solver_divergence_metrics

  !> CFL-based adaptive time step: per-cell sums of all face fluxes/diffusion, not per-face max.
  real(rp) function solver_compute_dt(self) result(dt)
    class(FVSolver), intent(in), target :: self
    integer(ip) :: f, o, nb, ierr, nc, ntot, nf, n_patches, patch_id, bc_id
    real(rp)    :: cfl_max_local, cfl_max_global
    real(rp)    :: diff_max_local, diff_max_global, nu_eff, area, w, gc
    real(rp)    :: dt_min, safety_factor, dt_conv, dt_diff
    real(rp)    :: u_face_use, sf1, sf2, sf3, nu_mol
    logical     :: have_uface, dev
    real(rp), parameter :: diff_cfl_limit = 1.0_rp
    integer(ip), pointer :: fown(:), fnb(:), fpid(:), pbc(:)
    real(rp),    pointer :: fsf(:,:), gcs(:), dcf(:), cvol(:), U(:,:), Uf(:), nus(:)

    nc = self%mesh%n_cells
    ntot = nc + self%mesh%n_ghost
    nf = self%mesh%n_faces
    n_patches = self%mesh%n_patches
    nu_mol = self%params%nu
    dev = self%on_device
    if (.not. allocated(dt_flux_sum_buf)) then
      allocate(dt_flux_sum_buf(nc))
      allocate(dt_diff_sum_buf(nc))
    else if (size(dt_flux_sum_buf) < nc) then
      deallocate(dt_flux_sum_buf, dt_diff_sum_buf)
      allocate(dt_flux_sum_buf(nc))
      allocate(dt_diff_sum_buf(nc))
    end if
    fown => self%mesh%face_owner_soa
    fnb  => self%mesh%face_neighbour_soa
    fpid => self%mesh%face_patch_id_soa
    pbc  => self%mesh%patch_bc_id_soa
    fsf  => self%mesh%face_sf_soa
    gcs  => self%mesh%face_gc_soa
    dcf  => self%mesh%face_dcf_soa
    cvol => self%mesh%cell_volume_soa
    U => self%U
    Uf => self%U_face
    nus => self%nu_sgs

    !$acc parallel loop if(dev) default(present)
    do o = 1, nc
      dt_flux_sum_buf(o) = 0.0_rp
      dt_diff_sum_buf(o) = 0.0_rp
    end do

    ! detect step-1 case where U_face has not been built yet. Fall back
    ! to a centred owner-side dot(U,sf) so the convective CFL reflects the
    have_uface = allocated(self%U_face) .and. self%u_face_ready

    ! Explicit stability depends on the SUM of every face's contribution into a cell
    ! (mod_advection/mod_viscous accumulate over the full face loop), not any single
    ! face alone: a per-face max undercounts polyhedral buffer cells at deep
    ! snappyHexMesh refinement interfaces, which can carry far more than 6 faces.
#ifdef _OPENACC
    !$acc parallel loop if(dev) default(present) &
    !$acc      private(o, nb, sf1, sf2, sf3, area, u_face_use, nu_eff, patch_id, bc_id, w, gc)
#else
    !$omp parallel do default(shared) &
    !$omp& private(f, o, nb, sf1, sf2, sf3, area, u_face_use, nu_eff, patch_id, bc_id, w, gc)
#endif
    do f = 1, nf
      o  = fown(f)
      nb = fnb(f)
      sf1 = fsf(1, f);  sf2 = fsf(2, f);  sf3 = fsf(3, f)
      area = sqrt(sf1**2 + sf2**2 + sf3**2)

      if (have_uface) then
        u_face_use = Uf(f)
      else
        u_face_use = U(1, o)*sf1 + U(2, o)*sf2 + U(3, o)*sf3
      end if
      if (o  >= 1 .and. o  <= nc) then
#ifdef _OPENACC
        !$acc atomic update
#else
        !$omp atomic
#endif
        dt_flux_sum_buf(o)  = dt_flux_sum_buf(o)  + abs(u_face_use)
      end if
      if (nb >= 1 .and. nb <= nc) then
#ifdef _OPENACC
        !$acc atomic update
#else
        !$omp atomic
#endif
        dt_flux_sum_buf(nb) = dt_flux_sum_buf(nb) + abs(u_face_use)
      end if

      ! Diffusive limit: gc-weighted nu_eff (mirrors mod_viscous face assembly).
      w = area / max(dcf(f), SMALL)
      if (nb >= 1 .and. nb <= ntot) then
        if (nb <= nc) then
          gc = gcs(f)
          nu_eff = nu_mol + gc * max(0.0_rp, nus(o)) + (1.0_rp - gc) * max(0.0_rp, nus(nb))
        else
          nu_eff = nu_mol + max(0.0_rp, nus(o))
        end if
#ifdef _OPENACC
        !$acc atomic update
#else
        !$omp atomic
#endif
        dt_diff_sum_buf(o) = dt_diff_sum_buf(o) + nu_eff * w
        if (nb >= 1 .and. nb <= nc) then
#ifdef _OPENACC
          !$acc atomic update
#else
          !$omp atomic
#endif
          dt_diff_sum_buf(nb) = dt_diff_sum_buf(nb) + nu_eff * w
        end if
      else if (nb < 1 .and. o >= 1 .and. o <= nc) then
        ! No-slip/Dirichlet boundary faces apply a real explicit viscous stress at the
        ! owner-face distance (mod_viscous.f90); excluding them let near-wall cells with a
        ! thin wall-normal spacing (thinner than volume^(1/3)) run at an unstable dt.
        ! Wall-model faces are excluded: their stress is applied via a separate limiter.
        patch_id = fpid(f)
        bc_id = -1
        if (patch_id >= 1 .and. patch_id <= n_patches) bc_id = pbc(patch_id)
        if (bc_id /= BC_WALL_MODEL) then
          nu_eff = nu_mol + max(0.0_rp, nus(o))
#ifdef _OPENACC
          !$acc atomic update
#else
          !$omp atomic
#endif
          dt_diff_sum_buf(o) = dt_diff_sum_buf(o) + nu_eff * w
        end if
      end if
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif

    cfl_max_local  = 0.0_rp
    diff_max_local = 0.0_rp
    !$acc parallel loop if(dev) default(present) reduction(max:cfl_max_local, diff_max_local)
    do o = 1, nc
      cfl_max_local  = max(cfl_max_local,  dt_flux_sum_buf(o) / (cvol(o) + SMALL))
      diff_max_local = max(diff_max_local, dt_diff_sum_buf(o) / (cvol(o) + SMALL))
    end do

    block
      real(rp) :: red_local(2), red_global(2)
      red_local = [cfl_max_local, diff_max_local]
      call MPI_Allreduce(red_local, red_global, 2, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, self%mesh%comm, ierr)
      cfl_max_global  = red_global(1)
      diff_max_global = red_global(2)
    end block

    ! Compute convective and diffusive dt limits, then apply safety margin.
    safety_factor = 0.8_rp
    dt_conv = huge(1.0_rp)
    dt_diff = huge(1.0_rp)
    if (cfl_max_global > SMALL) dt_conv = self%params%cfl_max / cfl_max_global
    if (diff_max_global > SMALL) dt_diff = diff_cfl_limit / diff_max_global

    dt = safety_factor * min(dt_conv, dt_diff)
    if (.not. (dt > 0.0_rp)) dt = self%params%dt

    ! Keep adaptive dt in a numerically safe range (both upper and lower bounds)
    ! Upper bound: never exceed user-specified maximum dt
    dt = min(dt, self%params%dt)

    ! Lower bound: abort if dt becomes dangerously small (instability indicator)
    dt_min = 1.0e-6_rp * self%params%dt
    if (.not. (dt > 0.0_rp)) then
      dt = self%params%dt
    else if (dt < dt_min) then
      ! dt this far below the requested value means the explicit limits are violated (blow-up);
      ! clamping up would run an unstable step, so stop. Every rank sees the same reduced dt.
      if (self%mesh%rank == 0) write(*,'(A,ES12.4,A,ES12.4)') &
        ' ERROR: stable dt =', dt, ' fell below 1e-6*dt =', dt_min
      error stop 'solver_compute_dt: stable time step collapsed (solution diverging)'
    else
      dt = max(dt, dt_min)
    end if
  end function solver_compute_dt

  !> Clamp dt down (never up) so solver_time+dt lands exactly on the next write_time
  !> boundary; only ever shrinks dt, so it cannot make the step less CFL-stable.
  pure subroutine snap_dt_to_write_time(solver_time, write_time, dt, hit_boundary)
    real(rp), intent(in)    :: solver_time, write_time
    real(rp), intent(inout) :: dt
    logical,  intent(out)   :: hit_boundary

    real(rp) :: t_next
    real(rp), parameter :: REL_TOL = 1.0e-9_rp

    hit_boundary = .false.
    if (write_time <= 0.0_rp .or. dt <= 0.0_rp) return

    ! The +REL_TOL guards against a prior exact snap leaving solver_time a hair below
    ! n*write_time (floating-point residue), which would otherwise re-target the same
    ! boundary instead of advancing to the next one.
    ! floor(...) uses int64 (lp): int32 would overflow for a long run with a small write_time.
    t_next = write_time * real(floor(solver_time / write_time + REL_TOL, lp) + 1_lp, rp)

    if (solver_time + dt >= t_next - REL_TOL * write_time) then
      dt = max(t_next - solver_time, 0.0_rp)
      hit_boundary = .true.
    end if
  end subroutine snap_dt_to_write_time

  !> Global max convective CFL rate (1/s): CFL number = compute_cfl_rate() * dt.
  !> MPI-collective; must be called by all ranks.
  real(rp) function solver_compute_cfl_rate(self) result(cfl_max)
    class(FVSolver), intent(in) :: self
    integer(ip) :: f, o, nb, ierr
    real(rp)    :: cfl_loc, cfl_max_local, u_face_use, sf(3), vol_min
    logical     :: have_uface

    cfl_max_local = 0.0_rp
    have_uface = allocated(self%U_face) .and. self%u_face_ready

    do f = 1, self%mesh%n_faces
      o  = self%mesh%faces(f)%owner
      nb = self%mesh%faces(f)%neighbour
      if (nb >= 1 .and. nb <= self%mesh%n_cells) then
        vol_min = min(self%mesh%cells(o)%volume, self%mesh%cells(nb)%volume)
      else
        vol_min = self%mesh%cells(o)%volume
      end if
      if (have_uface) then
        u_face_use = self%U_face(f)
      else
        sf = self%mesh%faces(f)%sf
        u_face_use = dot_product(self%U(:, o), sf)
      end if
      cfl_loc = abs(u_face_use) / (vol_min + SMALL)
      cfl_max_local = max(cfl_max_local, cfl_loc)
    end do

    call MPI_Allreduce(cfl_max_local, cfl_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, self%mesh%comm, ierr)
  end function solver_compute_cfl_rate

  real(rp) function wall_time_now() result(t)
    integer :: count, rate

    call system_clock(count, rate)
    if (rate > 0) then
      t = real(count, rp) / real(rate, rp)
    else
      t = 0.0_rp
    end if
  end function wall_time_now

  subroutine maybe_report_perf(self, force)
    class(FVSolver), intent(inout) :: self
    logical,         intent(in)    :: force

    integer(ip) :: ierr
    integer(ip) :: steps_local, steps_global
    integer(ip) :: p_solves_local, p_solves_sum
    real(rp) :: t_step_max, t_halo_max, t_rhs_max, t_predictor_max, t_pressure_max, t_proj_max
    real(rp) :: t_step_avg_sum, t_halo_avg_sum, t_rhs_avg_sum, t_predictor_avg_sum, t_pressure_avg_sum, t_proj_avg_sum
    real(rp) :: denom, p_iters_per_solve_local, p_iters_per_solve_max, p_iters_per_solve_avg_sum
    real(rp) :: p_resid_max_global
    real(rp) :: t_rhs_sgs_max, t_rhs_grad_max, t_rhs_transport_max, t_rhs_sources_max
    real(rp) :: t_rhs_sgs_avg_sum, t_rhs_grad_avg_sum, t_rhs_transport_avg_sum, t_rhs_sources_avg_sum
    ! Sub-split of predictor (see perf_t_predictor_bc/precorrect/regrad).
    real(rp) :: t_predictor_bc_max, t_predictor_precorrect_max, t_predictor_regrad_max
    real(rp) :: t_predictor_bc_avg_sum, t_predictor_precorrect_avg_sum, t_predictor_regrad_avg_sum

    if (.not. self%perf_enabled) return
    if (.not. force) then
      if (self%step <= 0) return
      if (mod(self%step, PERF_REPORT_INTERVAL) /= 0) return
    end if
    if (.not. associated(self%mesh)) return

    steps_local = self%perf_steps
    call MPI_Allreduce(steps_local, steps_global, 1, MPI_INTEGER, MPI_MAX, self%mesh%comm, ierr)
    if (steps_global <= 0_ip) return

    p_solves_local = self%perf_pressure_solves
    if (p_solves_local > 0) then
      p_iters_per_solve_local = real(self%perf_pressure_iters, rp) / real(p_solves_local, rp)
    else
      p_iters_per_solve_local = 0.0_rp
    end if

    ! F7: pack 7 MAX-reduce values + 6 SUM-reduce values into 2 Allreduces
    ! (was 13 separate Allreduce calls).
    block
      real(rp) :: max_in(8), max_out(8), sum_in(7), sum_out(7)
      integer(ip) :: int_in(1), int_out(1)

      max_in = [ self%perf_t_step, self%perf_t_halo_bc, self%perf_t_rhs, self%perf_t_predictor, &
                 self%perf_t_pressure, self%perf_t_projection,           &
                 p_iters_per_solve_local, self%perf_pressure_residual_max ]
      sum_in = [ self%perf_t_step, self%perf_t_halo_bc, self%perf_t_rhs, self%perf_t_predictor, &
                 self%perf_t_pressure, self%perf_t_projection,           &
                 p_iters_per_solve_local ]

      call MPI_Allreduce(max_in, max_out, 8, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
      call MPI_Allreduce(sum_in, sum_out, 7, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)

      int_in(1) = p_solves_local
      call MPI_Allreduce(int_in, int_out, 1, MPI_INTEGER, MPI_SUM, self%mesh%comm, ierr)

      t_step_max      = max_out(1)
      t_halo_max      = max_out(2)
      t_rhs_max       = max_out(3)
      t_predictor_max = max_out(4)
      t_pressure_max  = max_out(5)
      t_proj_max      = max_out(6)
      p_iters_per_solve_max = max_out(7)
      p_resid_max_global    = max_out(8)

      t_step_avg_sum      = sum_out(1)
      t_halo_avg_sum      = sum_out(2)
      t_rhs_avg_sum       = sum_out(3)
      t_predictor_avg_sum = sum_out(4)
      t_pressure_avg_sum  = sum_out(5)
      t_proj_avg_sum      = sum_out(6)
      p_iters_per_solve_avg_sum = sum_out(7)

      p_solves_sum = int_out(1)
    end block

    block
      real(rp) :: max_in(7), max_out(7), sum_in(7), sum_out(7)

      max_in = [ self%perf_t_rhs_sgs, self%perf_t_rhs_grad, &
                 self%perf_t_rhs_transport, self%perf_t_rhs_sources, &
                 self%perf_t_predictor_bc, self%perf_t_predictor_precorrect, &
                 self%perf_t_predictor_regrad ]
      sum_in = max_in

      call MPI_Allreduce(max_in, max_out, 7, MPI_DOUBLE_PRECISION, MPI_MAX, self%mesh%comm, ierr)
      call MPI_Allreduce(sum_in, sum_out, 7, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)

      t_rhs_sgs_max       = max_out(1)
      t_rhs_grad_max       = max_out(2)
      t_rhs_transport_max = max_out(3)
      t_rhs_sources_max   = max_out(4)
      t_predictor_bc_max         = max_out(5)
      t_predictor_precorrect_max = max_out(6)
      t_predictor_regrad_max     = max_out(7)

      t_rhs_sgs_avg_sum       = sum_out(1)
      t_rhs_grad_avg_sum       = sum_out(2)
      t_rhs_transport_avg_sum = sum_out(3)
      t_rhs_sources_avg_sum   = sum_out(4)
      t_predictor_bc_avg_sum         = sum_out(5)
      t_predictor_precorrect_avg_sum = sum_out(6)
      t_predictor_regrad_avg_sum     = sum_out(7)
    end block

    if (self%mesh%rank == 0) then
      denom = real(max(1_ip, steps_global), rp)
      write(*,'(A)') '--- Performance Snapshot (sec/step; max-rank and rank-average) ---'
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'step_total   max', t_step_max/denom, 'avg', (t_step_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'halo+bc      max', t_halo_max/denom, 'avg', (t_halo_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'rhs          max', t_rhs_max/denom, 'avg', (t_rhs_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  rhs_sgs       max', t_rhs_sgs_max/denom, 'avg', (t_rhs_sgs_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  rhs_grad      max', t_rhs_grad_max/denom, 'avg', (t_rhs_grad_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  rhs_transport max', t_rhs_transport_max/denom, 'avg', (t_rhs_transport_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  rhs_sources   max', t_rhs_sources_max/denom, 'avg', (t_rhs_sources_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'predictor    max', t_predictor_max/denom, 'avg', (t_predictor_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  pred_bc         max', t_predictor_bc_max/denom, 'avg', (t_predictor_bc_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  pred_precorrect max', t_predictor_precorrect_max/denom, 'avg', (t_predictor_precorrect_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') '  pred_regrad     max', t_predictor_regrad_max/denom, 'avg', (t_predictor_regrad_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'pressure     max', t_pressure_max/denom, 'avg', (t_pressure_avg_sum/real(self%mesh%nranks, rp))/denom
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'projection   max', t_proj_max/denom, 'avg', (t_proj_avg_sum/real(self%mesh%nranks, rp))/denom
      if (p_solves_sum > 0) then
        write(*,'(A,1X,F10.3,1X,A,1X,F10.3)') 'p_iters/solve max', p_iters_per_solve_max, 'avg', &
          p_iters_per_solve_avg_sum/real(self%mesh%nranks, rp)
        write(*,'(A,1X,ES12.4)') 'p_resid_max   max', p_resid_max_global
      end if
      write(*,'(A)') '--- End Performance Snapshot ---'
    end if

    call report_interface_diagnostics(self)

  end subroutine maybe_report_perf

  !> Report the convective KE budget and divergence split; smooth-face production is the KEP control.
  subroutine report_interface_diagnostics(self)
    class(FVSolver), intent(inout) :: self

    real(rp) :: loc(5), glb(5), denom
    integer  :: ierr

    if (.not. self%diag_interface) return

    call advection_get_ke_budget(loc(1), loc(2), loc(3), loc(4), loc(5))
    call MPI_Allreduce(loc, glb, 5, MPI_DOUBLE_PRECISION, MPI_SUM, self%mesh%comm, ierr)
    call advection_reset_ke_budget()

    ! The budget resets every report, so it must be divided by steps since the last one.
    denom = real(max(1_ip, self%step - self%diag_last_step), rp)
    self%diag_last_step = self%step

    if (.not. self%params%debug_diagnostics) return

    ! The worst-cell report below runs on the host.
    call solver_sync_host(self)

    if (self%mesh%rank == 0) then
#ifdef _OPENACC
      write(*,'(A)') ' [interface-diag] KE budget is not accumulated in device builds (values below are zero)'
#endif
      write(*,'(A)') '--- Interface Diagnostics ---'
      write(*,'(A,1X,ES12.4)') 'KE production  smooth  ', glb(1)/denom
      write(*,'(A,1X,ES12.4)') 'KE production  refine  ', glb(2)/denom
      write(*,'(A,1X,ES12.4)') 'KE split-form  cancel  ', glb(3)/denom
      write(*,'(A,1X,ES12.4)') 'KE boundary    flux    ', glb(4)/denom
      write(*,'(A,1X,ES12.4)') 'net dE/dt      convect ', (glb(1)+glb(2)+glb(3)+glb(4))/denom
      ! Interior conservation error: everything except the physical boundary throughflow.
      write(*,'(A,1X,ES12.4)') 'CONSERVATION   error   ', &
        (glb(1)+glb(2)+glb(3)+glb(4)-glb(5))/denom
      write(*,'(A,1X,ES12.4)') 'post-proj |div| refine ', self%diag_div_refine
      write(*,'(A,1X,ES12.4)') 'post-proj |div| smooth ', self%diag_div_smooth
      write(*,'(A,1X,I0)')     'nonortho correctors used', self%diag_nonortho_iters
      write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4,1X,A,1X,ES12.4)') &
        'nonortho exit: curr', self%diag_nonortho_curr, 'prev', self%diag_nonortho_prev, &
        'atol', self%diag_nonortho_atol
      write(*,'(A)') '--- End Interface Diagnostics ---'
    end if

    call report_worst_divergence_cell(self)
  end subroutine report_interface_diagnostics

  !> Localize the cell with the largest Green-Gauss cell-velocity divergence (never zeroed by projection).
  subroutine report_worst_divergence_cell(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: nc, c, ic_max, f, o, nb
    real(rp)    :: div_u(self%mesh%n_cells)
    real(rp)    :: val_local, in2(2), out2(2)
    logical, allocatable :: cell_touches_refine(:)
    integer(ip) :: ic_max_rf
    real(rp)    :: val_local_rf, in2_rf(2), out2_rf(2)
    integer     :: ierr, owner_rank, owner_rank_rf

    nc = self%mesh%n_cells
    call cell_velocity_divergence(self, div_u(1:nc))

    ! Flag every cell that has at least one refinement-interface face, to separate
    ! this session's original suspect (refine faces) from an unrelated bad cell.
    allocate(cell_touches_refine(nc), source=.false.)
    do f = 1, self%mesh%n_faces
      if (.not. self%mesh%faces(f)%is_refine) cycle
      o  = self%mesh%faces(f)%owner
      nb = self%mesh%faces(f)%neighbour
      if (o  >= 1 .and. o  <= nc) cell_touches_refine(o)  = .true.
      if (nb >= 1 .and. nb <= nc) cell_touches_refine(nb) = .true.
    end do

    val_local    = 0.0_rp
    ic_max       = 0
    val_local_rf = 0.0_rp
    ic_max_rf    = 0
    do c = 1, nc
      if (abs(div_u(c)) > val_local) then
        val_local = abs(div_u(c))
        ic_max = c
      end if
      if (cell_touches_refine(c) .and. abs(div_u(c)) > val_local_rf) then
        val_local_rf = abs(div_u(c))
        ic_max_rf = c
      end if
    end do
    deallocate(cell_touches_refine)

    in2(1) = val_local
    in2(2) = real(self%mesh%rank, rp)
    call MPI_Allreduce(in2, out2, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, self%mesh%comm, ierr)
    owner_rank = nint(out2(2))

    in2_rf(1) = val_local_rf
    in2_rf(2) = real(self%mesh%rank, rp)
    call MPI_Allreduce(in2_rf, out2_rf, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, self%mesh%comm, ierr)
    owner_rank_rf = nint(out2_rf(2))

    call gather_and_print(self, ic_max, owner_rank, out2(1), 'Worst Divergence Cell (global)')
    call gather_and_print(self, ic_max_rf, owner_rank_rf, out2_rf(1), &
                          'Worst Divergence Cell (refine-adjacent)')

  contains

    !> Gather this cell's geometry on its owning rank and print it from rank 0.
    subroutine gather_and_print(slv, ic, rank_owner, worst_val, label)
      class(FVSolver), intent(in) :: slv
      integer(ip),      intent(in) :: ic
      integer,          intent(in) :: rank_owner
      real(rp),         intent(in) :: worst_val
      character(len=*), intent(in) :: label

      integer(ip) :: fi2, f2, o2, nb2, other2, nfc, f_worst
      real(rp)    :: minv, mskew, buf(14)
      real(rp)    :: d_on2(3), k_f2(3), s_dot_dhat2, d_on_len2, sf_dot_don
      real(rp)    :: t_recon(3,3), t_eigs(3)
      logical     :: rf
      integer     :: ierr2

      buf = 0.0_rp
      if (worst_val <= 0.0_rp) then
        if (slv%mesh%rank == 0) then
          write(*,'(A,1X,A,1X,A)') '---', trim(label), '--- (no cells matched) ---'
        end if
        return
      end if

      if (slv%mesh%rank == rank_owner .and. ic >= 1) then
        rf      = .false.
        minv    = huge(1.0_rp)
        mskew   = 0.0_rp
        f_worst = 0
        nfc     = slv%mesh%cells(ic)%n_faces
        do fi2 = slv%mesh%cells(ic)%face_offset, slv%mesh%cells(ic)%face_offset + nfc - 1
          f2 = slv%mesh%cell_faces(fi2)
          if (f2 < 1 .or. f2 > size(slv%mesh%faces)) cycle
          if (slv%mesh%faces(f2)%is_refine) rf = .true.
          if (slv%mesh%faces(f2)%skew > mskew) then
            mskew   = slv%mesh%faces(f2)%skew
            f_worst = f2
          end if
          o2  = slv%mesh%faces(f2)%owner
          nb2 = slv%mesh%faces(f2)%neighbour
          if (nb2 < 1) cycle
          other2 = merge(nb2, o2, o2 == ic)
          if (other2 >= 1 .and. other2 <= size(slv%mesh%cells) .and. other2 /= ic) &
            minv = min(minv, slv%mesh%cells(other2)%volume)
        end do
        ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
        ! must never be evaluated in the same expression as the allocated() check.
        buf(1) = real(ic, rp)
        if (allocated(slv%mesh%global_id)) then
          if (ic <= size(slv%mesh%global_id)) buf(1) = real(slv%mesh%global_id(ic), rp)
        end if
        buf(2)   = slv%mesh%cells(ic)%volume
        buf(3:5) = slv%mesh%cells(ic)%centroid
        buf(6)   = merge(1.0_rp, 0.0_rp, rf)
        buf(7)   = real(nfc, rp)
        buf(8)   = merge(minv, -1.0_rp, minv < huge(1.0_rp))
        buf(9)   = mskew

        ! Non-orthogonality diagnostic on the worst face: sin2theta = |k|^2/|Sf|^2 must be
        ! < 1 and Sf.d_on > 0 for the over-relaxed correction's contraction guarantee to hold.
        if (f_worst >= 1) then
          o2  = slv%mesh%faces(f_worst)%owner
          nb2 = slv%mesh%faces(f_worst)%neighbour
          if (nb2 >= 1) then
            d_on2 = slv%mesh%cells(nb2)%centroid - slv%mesh%cells(o2)%centroid
            call nonortho_kfactor_overrelaxed(slv%mesh%faces(f_worst)%sf, d_on2, k_f2, &
                                              s_dot_dhat2, d_on_len2)
            sf_dot_don = dot_product(slv%mesh%faces(f_worst)%sf, d_on2)
            buf(10) = sum(k_f2**2) / max(sum(slv%mesh%faces(f_worst)%sf**2), SMALL)
            buf(11) = sf_dot_don
            buf(12) = refine_blend_weight(mskew)
          end if
        end if

        ! T = sum_f |Sf| n_f(x)n_f recovered from its stored inverse; its condition number
        ! (lambda_max/lambda_min) measures how ill-posed the cell-velocity reconstruction is.
        ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
        ! must never be evaluated in the same expression as the allocated() check.
        if (allocated(slv%nu_sgs)) then
          if (ic <= size(slv%nu_sgs)) buf(13) = slv%nu_sgs(ic)
        end if
        if (allocated(slv%recon_Tinv)) then
          if (ic <= size(slv%recon_Tinv, 3)) then
            call invert_sym3(slv%recon_Tinv(:, :, ic), t_recon)
            t_eigs  = eig3_sym(t_recon)
            buf(14) = t_eigs(3) / max(t_eigs(1), SMALL)
          end if
        end if
      end if

      call MPI_Bcast(buf, 14, MPI_DOUBLE_PRECISION, rank_owner, slv%mesh%comm, ierr2)

      if (slv%mesh%rank == 0) then
        write(*,'(A,1X,A,1X,A)') '---', trim(label), '---'
        write(*,'(A,1X,ES12.4,1X,A,1X,I0)') 'max |div(U_cell)|', worst_val, 'on rank', rank_owner
        write(*,'(A,1X,I0,1X,A,1X,ES12.4,1X,A,1X,I0)') 'global_id', nint(buf(1),ip), &
          'volume', buf(2), 'n_faces', nint(buf(7),ip)
        write(*,'(A,1X,ES12.4,1X,ES12.4,1X,ES12.4)') 'centroid', buf(3), buf(4), buf(5)
        write(*,'(A,1X,L1)') 'touches refine-interface face', (buf(6) > 0.5_rp)
        write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'min neighbour volume', buf(8), &
          'volume ratio (this/min_nb)', buf(2) / max(buf(8), SMALL)
        write(*,'(A,1X,ES12.4)') 'max face skew on this cell', buf(9)
        write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4,1X,A,1X,ES12.4)') 'worst face: sin2theta', buf(10), &
          'Sf.d_on', buf(11), 'w_blend', buf(12)
        write(*,'(A,1X,ES12.4,1X,A,1X,ES12.4)') 'nu_sgs', buf(13), &
          'T cond number (lambda_max/lambda_min)', buf(14)
        write(*,'(A)') '--- End ' // trim(label) // ' ---'
      end if
    end subroutine gather_and_print

  end subroutine report_worst_divergence_cell

  subroutine solver_destroy(self)
    class(FVSolver), intent(inout) :: self
    call maybe_report_perf(self, .true.)
    call solver_leave_device(self)
    call release_static_on_device(self)
#ifdef _OPENACC
    if (allocated(self%inv_vol)) then
      !$acc exit data delete(self%inv_vol, self%U_face)
    end if
#endif
    call pressure_solver_free(self%psolver)
    if (allocated(self%U))          deallocate(self%U)
    if (allocated(self%p))          deallocate(self%p)
    if (allocated(self%p_old))      deallocate(self%p_old)
    if (allocated(self%T))          deallocate(self%T)
    if (allocated(self%C))          deallocate(self%C)
    if (allocated(self%nu_sgs))     deallocate(self%nu_sgs)
    if (allocated(self%U_bnd))      deallocate(self%U_bnd)
    if (allocated(self%p_bnd))          deallocate(self%p_bnd)
    if (allocated(self%phi_prev))        deallocate(self%phi_prev)
    if (allocated(self%phi_prev_bnd))    deallocate(self%phi_prev_bnd)
    if (allocated(self%T_bnd))           deallocate(self%T_bnd)
    if (allocated(self%C_bnd))      deallocate(self%C_bnd)
    if (allocated(self%tau_w_bnd))  deallocate(self%tau_w_bnd)
    if (allocated(self%sem_fluc_bnd)) deallocate(self%sem_fluc_bnd)
    if (allocated(self%cell_z))       deallocate(self%cell_z)
    if (allocated(self%grad_U))     deallocate(self%grad_U)
    if (allocated(self%grad_T))     deallocate(self%grad_T)
    if (allocated(self%grad_C))     deallocate(self%grad_C)
    if (allocated(self%face_flux))  deallocate(self%face_flux)
    if (allocated(self%U_face))     deallocate(self%U_face)
    if (allocated(self%inv_vol))    deallocate(self%inv_vol)
    if (allocated(self%w_H_U))      deallocate(self%w_H_U)
    if (allocated(self%w_H_T))      deallocate(self%w_H_T)
    if (allocated(self%w_H_C))      deallocate(self%w_H_C)
    if (allocated(self%w_rhs_p))    deallocate(self%w_rhs_p)
    if (allocated(self%w_div_diag)) deallocate(self%w_div_diag)
    if (allocated(self%w_U0))       deallocate(self%w_U0)
    if (allocated(self%w_T0))       deallocate(self%w_T0)
    if (allocated(self%w_C0))       deallocate(self%w_C0)
    if (allocated(self%w_conv_U))   deallocate(self%w_conv_U)
    if (allocated(self%w_visc_U))   deallocate(self%w_visc_U)
    if (allocated(self%w_conv_T))   deallocate(self%w_conv_T)
    if (allocated(self%w_visc_T))   deallocate(self%w_visc_T)
    if (allocated(self%w_conv_C))   deallocate(self%w_conv_C)
    if (allocated(self%w_visc_C))   deallocate(self%w_visc_C)
    if (allocated(self%w_flux_built)) deallocate(self%w_flux_built)
    if (allocated(self%w_gp))       deallocate(self%w_gp)
    if (allocated(self%w_gp_rc))    deallocate(self%w_gp_rc)
    if (allocated(self%w_gp_cell))  deallocate(self%w_gp_cell)
    if (allocated(self%w_p_tot))     deallocate(self%w_p_tot)
    if (allocated(self%w_p_tot_bnd)) deallocate(self%w_p_tot_bnd)
    if (allocated(self%recon_Tinv)) deallocate(self%recon_Tinv)
    if (allocated(self%cell_is_refine)) deallocate(self%cell_is_refine)
    if (allocated(self%w_gp_prev))  deallocate(self%w_gp_prev)
    if (allocated(self%w_outlet_cell)) deallocate(self%w_outlet_cell)
    if (allocated(self%w_outlet_pval)) deallocate(self%w_outlet_pval)
    if (allocated(self%patch_p_val))   deallocate(self%patch_p_val)
    if (allocated(self%w_Cs))          deallocate(self%w_Cs)
    if (allocated(self%w_Cs_bnd))      deallocate(self%w_Cs_bnd)
    if (allocated(self%w_pack_grad)) deallocate(self%w_pack_grad)
    if (allocated(self%w_g_prev_U))  deallocate(self%w_g_prev_U)
    if (allocated(self%w_g_prev_T))  deallocate(self%w_g_prev_T)
    if (allocated(self%bc_configs)) deallocate(self%bc_configs)
    if (allocated(self%periodic_patches)) deallocate(self%periodic_patches)
    if (allocated(self%nonperiodic_patches)) deallocate(self%nonperiodic_patches)
    if (allocated(self%wallmodel_patches)) deallocate(self%wallmodel_patches)
    if (allocated(self%periodic_bfaces)) deallocate(self%periodic_bfaces)
    if (allocated(self%nonperiodic_bfaces)) deallocate(self%nonperiodic_bfaces)
    if (allocated(self%nonortho_k_int))      deallocate(self%nonortho_k_int)
    if (allocated(self%nonortho_active_int)) deallocate(self%nonortho_active_int)
    if (allocated(self%nonortho_k_bnd))      deallocate(self%nonortho_k_bnd)
    if (allocated(self%nonortho_active_bnd)) deallocate(self%nonortho_active_bnd)
    if (allocated(self%params%scalar_names)) deallocate(self%params%scalar_names)
    if (allocated(self%params%scalar_diffusivity)) deallocate(self%params%scalar_diffusivity)
    if (allocated(self%params%scalar_pr_t)) deallocate(self%params%scalar_pr_t)
    if (allocated(self%params%scalar_sources)) deallocate(self%params%scalar_sources)
    if (allocated(self%params%petsc_options))  deallocate(self%params%petsc_options)
    if (allocated(self%q_w_bnd))               deallocate(self%q_w_bnd)
    if (allocated(self%abl_mean_bnd))          deallocate(self%abl_mean_bnd)
    if (allocated(self%L_bnd))                 deallocate(self%L_bnd)
    if (allocated(self%sem_state%pos))         deallocate(self%sem_state%pos)
    if (allocated(self%sem_state%amp))         deallocate(self%sem_state%amp)
    if (allocated(self%sem_state%sig))         deallocate(self%sem_state%sig)
    if (allocated(self%sem_state%uconv))       deallocate(self%sem_state%uconv)
    if (self%inflow_planes_on) call inflow_planes_close(self%inflow_planes)
    call vreman_free()
    call dynamic_smagorinsky_free()
    call wall_model_free()
    self%mesh => null()
  end subroutine solver_destroy

end module mod_fv_solver
