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
  use mod_abl_weather
  use mod_abl_wall
  use mod_pressure_solver
  use mod_bc
  use mod_warn_log
  use mod_profiling, only: prof_push, prof_pop
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
    procedure :: advance      => solver_advance
    procedure :: project_once => solver_project_once
    procedure :: projection_idempotency => solver_projection_idempotency
    procedure :: divergence_metrics => solver_divergence_metrics
    procedure :: destroy      => solver_destroy
    procedure :: compute_dt   => solver_compute_dt
    procedure :: compute_cfl_rate => solver_compute_cfl_rate
  end type FVSolver

  real(rp), allocatable, save :: dt_flux_sum_buf(:)  ! solver_compute_dt scratch, grown once
  real(rp), allocatable, save :: dt_diff_sum_buf(:)  ! solver_compute_dt scratch, grown once

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
  subroutine solver_init(self, msh, params, bc_configs, abl_cfg)
    class(FVSolver),    intent(inout) :: self
    type(Mesh), target, intent(inout) :: msh
    type(solver_params_t), intent(in) :: params
    type(BCConfig),    intent(in)     :: bc_configs(:)
    type(ABLConfig),   intent(in), optional :: abl_cfg

    integer(ip) :: nc, ng, nbnd, s, c
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
      allocate(self%w_conv_C(nc), self%w_visc_C(nc), source=0.0_rp)
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
    end if

    ! ABL module initialization
    ! SEM inflow, MOST wall treatment and nudging are independent; enable any combination.
    if (self%abl_config%wall_model_id == WM_MOST) then
      allocate(self%q_w_bnd(nbnd), source=0.0_rp)
      ! Zero means "unseeded"; the wall model falls back to abl_config%L_obukhov on the first call.
      allocate(self%L_bnd(nbnd), source=0.0_rp)
    end if

    if (self%abl_config%use_sem) then
      ! sem_init derives the inlet frame from the ABL inlet patch and reduces it across ranks.
      call sem_init(self%sem_state, self%abl_config, self%mesh)
      ! A prescribed profile also supplies the inlet mean, not just the stresses.
      if (self%sem_state%has_profile) allocate(self%abl_mean_bnd(3, nbnd), source=0.0_rp)
      call update_sem_inflow(self)
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
    call advection_set_ke_budget(self%diag_interface)

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

    call anchor_pressure_to_outlet(self%p(1:nc), self%mesh, self%w_outlet_cell, self%w_outlet_pval)
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
    class(FVSolver), intent(inout) :: self
    real(rp), intent(in) :: dt

    integer(ip) :: nc, c, ierr
    real(rp) :: mom_local(3), mom_global(3), U_bulk(3)

    nc = self%mesh%n_cells
    mom_local = 0.0_rp
    do c = 1, nc
      mom_local = mom_local + self%U(:, c) * self%mesh%cells(c)%volume
    end do
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
    class(FVSolver), intent(inout) :: self
    integer(ip) :: s

    select case (trim(self%params%gradient_scheme))
    case ('GREEN_GAUSS')
      call green_gauss_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd)
      if (self%solve_temperature) &
        call green_gauss_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd)
      do s = 1, self%params%n_scalars
        call green_gauss_gradient_scalar(self%C(s, :), self%mesh, self%grad_C(:, :, s), self%C_bnd(s, :))
      end do
      if (self%params%gg_iterations > 1) then
        block
          integer(ip) :: iter, si
          do iter = 2, self%params%gg_iterations
            ! Exchange halos so ghost-cell gradients from the previous sweep are
            ! valid before building the skewness correction for this sweep.
            call exchange_grad_halos(self)
            ! F3: reuse persistent buffers (allocated in solver_init).
            self%w_g_prev_U = self%grad_U
            call green_gauss_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd, self%w_g_prev_U)
            if (self%solve_temperature) then
              self%w_g_prev_T = self%grad_T
              call green_gauss_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd, self%w_g_prev_T)
            end if
            do si = 1, self%params%n_scalars
              ! Reload prev-gradient scratch for each scalar field.
              self%w_g_prev_T = self%grad_C(:, :, si)
              call green_gauss_gradient_scalar(self%C(si,:), self%mesh, self%grad_C(:,:,si), &
                                               self%C_bnd(si,:), self%w_g_prev_T)
            end do
          end do
        end block
      end if
    case default
      call ls_gradient_vector(self%U, self%mesh, self%grad_U, self%U_bnd)
      if (self%solve_temperature) &
        call ls_gradient_scalar(self%T, self%mesh, self%grad_T, self%T_bnd)
      do s = 1, self%params%n_scalars
        call ls_gradient_scalar(self%C(s, :), self%mesh, self%grad_C(:, :, s), self%C_bnd(s, :))
      end do
    end select

    ! Halo-exchange gradients so ghost cells carry valid values for viscous/advection stencils.
    call exchange_grad_halos(self)
  end subroutine compute_all_gradients

  !> Pack, halo-exchange, and unpack grad_U/grad_T/grad_C across MPI ranks.
  subroutine exchange_grad_halos(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: nc, ng, ntot, n_scal, n_fields, c, s, k

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

    associate(pack => self%w_pack_grad)
    ! Pack owned cells; ghost slots will be filled by the exchange.
    if (n_scal == 0) then
      ! F8: lighter U+T-only path (no scalar transport configured).
      do c = 1, nc
        pack(1:9, c)   = reshape(self%grad_U(:, :, c), [9])
        pack(10:12, c) = self%grad_T(:, c)
      end do
    else
      do c = 1, nc
        pack(1:9, c)   = reshape(self%grad_U(:, :, c), [9])
        pack(10:12, c) = self%grad_T(:, c)
        do s = 1, n_scal
          k = 12 + 3*(s-1)
          pack(k+1:k+3, c) = self%grad_C(:, c, s)
        end do
      end do
    end if

    call self%mesh%exchange_halos_scalars(pack)

    ! Unpack ghost cells (no need to overwrite owned cells).
    if (n_scal == 0) then
      do c = nc+1, ntot
        self%grad_U(:, :, c) = reshape(pack(1:9, c), [3, 3])
        self%grad_T(:, c)    = pack(10:12, c)
      end do
    else
      do c = nc+1, ntot
        self%grad_U(:, :, c) = reshape(pack(1:9, c), [3, 3])
        self%grad_T(:, c)    = pack(10:12, c)
        do s = 1, n_scal
          k = 12 + 3*(s-1)
          self%grad_C(:, c, s) = pack(k+1:k+3, c)
        end do
      end do
    end if
    end associate

    ! pack's ghost range is already device-synced by exchange_halos_scalars, but the unpack above
    ! into grad_U/grad_T/grad_C is a separate host step, so their device ghost cells need their own
    ! push. A no-op if these aren't currently device-resident.
    !$acc update device(self%grad_U(:,:,nc+1:ntot), self%grad_T(:,nc+1:ntot)) if(ng > 0) if_present
    if (n_scal > 0) then
      !$acc update device(self%grad_C(:,nc+1:ntot,:)) if(ng > 0) if_present
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
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: dt

    integer(ip) :: nc, c, k, jno, p_iter
    real(rp)    :: gk_dt, time_n
    real(rp)    :: t_perf, t_sub, p_resid

    nc = self%mesh%n_cells
    associate(H_U => self%w_H_U(:, 1:nc), H_T => self%w_H_T(1:nc), &
          rhs_p => self%w_rhs_p(1:nc), &
          H_C => self%w_H_C(1:max(1_ip,self%params%n_scalars), 1:nc), &
          H_U_prev => self%w_U0(:, 1:nc), H_T_prev => self%w_T0(1:nc), &
          H_C_prev => self%w_C0(1:max(1_ip,self%params%n_scalars), 1:nc))

    H_U = 0.0_rp
    H_T = 0.0_rp
    rhs_p = 0.0_rp
    H_C = 0.0_rp
    H_U_prev = 0.0_rp  ! beta(1)=0 so this is never read, but zero for safety
    H_T_prev = 0.0_rp
    H_C_prev = 0.0_rp

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
        self%nu_sgs(1:nc) = 0.0_rp
      end if

      call exchange_nu_sgs_halo(self)
      select case (self%abl_config%wall_model_id)
      case (WM_MOST)
        call abl_wall_model(self%U, self%T, self%mesh, self%abl_config, self%bc_configs, &
                            self%tau_w_bnd, self%q_w_bnd, self%L_bnd, update_L=(k == 1))
      case (WM_NONE)
        self%tau_w_bnd = 0.0_rp
      case default
        ! Filter advances once per step, so stages 2 and 3 pass dt = 0.
        call apply_wall_model_bc(self%U, self%params%nu, self%mesh, self%tau_w_bnd, &
                                 h_wm=self%params%wall_model_height, &
                                 t_filt=self%params%wall_model_filter_time, &
                                 dt=merge(dt, 0.0_rp, k == 1))
      end select
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
      !$omp parallel do schedule(static)
      do c = 1, nc
        self%U(:, c) = self%U(:, c) + dt * (RK3W_ALPHA(k) * H_U(:, c) + RK3W_BETA(k) * H_U_prev(:, c))
        self%T(c)    = self%T(c)    + dt * (RK3W_ALPHA(k) * H_T(c)     + RK3W_BETA(k) * H_T_prev(c))
      end do
      !$omp end parallel do
      if (self%params%n_scalars > 0) then
        !$omp parallel do schedule(static)
        do c = 1, nc
          self%C(:, c) = self%C(:, c) + dt * (RK3W_ALPHA(k) * H_C(:, c) + RK3W_BETA(k) * H_C_prev(:, c))
        end do
        !$omp end parallel do
      end if

      call self%mesh%begin_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, &
                             periodic_only=.false., grad_U=self%grad_U, u_fluc_bnd=self%sem_fluc_bnd, u_mean_bnd=self%abl_mean_bnd)
      call self%mesh%end_halos_vec(self%U)
      call apply_velocity_bc(self%U, self%U_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      self%perf_t_predictor_bc = self%perf_t_predictor_bc + (wall_time_now() - t_sub)

      ! -- Pre-correction: u_hat = u* - gk_dt*grad(p^{k-1/2}); net gives p^{k+1/2}=p^{k-1/2}+phi.
      t_sub = wall_time_now()
      self%p_old(1:nc) = self%p(1:nc)
      call correct_velocity_psolver_dual(self, gk_dt, with_cross=.false.)
      ! Rhie-Chow must add back exactly the gradient just removed, so record it here.
      self%w_gp_rc(:, 1:nc) = self%w_gp_cell(:, 1:nc)
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
      self%phi_prev(1:nc+self%mesh%n_ghost) = self%p(1:nc+self%mesh%n_ghost)
      self%phi_prev_bnd(:) = self%p_bnd(:)

      ! -- Pressure Poisson with effective sub-step gk_dt
      t_perf = wall_time_now()
      call prof_push('pressure')
      block
        real(rp) :: u_bulk_local(2), u_bulk_global(2), u_scale, atol_dynamic
        integer :: ierr_tol
        integer(ip) :: n_nonortho_max
        real(rp) :: divU_norm_prev, divU_norm_curr, divU_local, stagnate_tol
        ! Volume-weighted mean speed, not the domain max: max(|U|) is a single-cell
        ! quantity and loosens atol precisely when one cell is diverging.
        u_bulk_local(1) = sum(sqrt(self%U(1,1:nc)**2 + self%U(2,1:nc)**2 + self%U(3,1:nc)**2) * &
                               self%mesh%cells(1:nc)%volume)
        u_bulk_local(2) = sum(self%mesh%cells(1:nc)%volume)
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
          divU_local = maxval(abs(rhs_p(1:nc)))
          call MPI_Allreduce(divU_local, divU_norm_curr, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, self%mesh%comm, ierr_tol)
          if (self%diag_interface) self%diag_nonortho_iters = jno
          ! self%p enters this loop holding p^{k-1/2} (needed by compute_face_fluxes_and_rhs
          ! above for the jno==1 gradient), but the KSP solves for the increment phi. Seeding
          ! its nonzero initial guess with O(1) p^{k-1/2} instead of the O(dt*div u*) increment
          ! wastes iterations and loosens the effective (relative) convergence tolerance.
          if (jno == 1) self%p(1:nc) = 0.0_rp
          call prof_push('ksp_solve')
          call pressure_solver_solve(self%psolver, rhs_p, self%p(1:nc), tol=self%params%p_tol, &
                         max_iter=self%params%p_max_iter, n_iter=p_iter, residual=p_resid, &
                         atol_in=atol_dynamic)
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
      self%U_face(:) = self%face_flux(:)               ! persist divergence-free flux
      self%u_face_ready = .true.
#ifdef _OPENACC
      !$acc update device(self%U_face)
#endif
      if (self%diag_interface) call record_divergence_split(self, jno)
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
        !$omp parallel do schedule(static)
        do c = 1, nc
          self%p(c) = self%p(c) + (self%params%nu + self%nu_sgs(c)) * rhs_p(c) * gk_dt * self%inv_vol(c)
        end do
        !$omp end parallel do
        call self%mesh%begin_halos_scalar(self%p)
        call apply_pressure_bc_o2(self)
        call self%mesh%end_halos_scalar(self%p)
        call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.true.)
      end if

      ! -- Accumulate: p^{k+1/2} = p^{k-1/2} + phi (CaNS Eq. 2d).
      self%p(1:nc) = self%p(1:nc) + self%p_old(1:nc)
      call anchor_pressure_to_outlet(self%p(1:nc), self%mesh, self%w_outlet_cell, self%w_outlet_pval)
      call self%mesh%begin_halos_scalar(self%p)
      call apply_pressure_bc_o2(self)
      call self%mesh%end_halos_scalar(self%p)
      call apply_pressure_bc(self%p, self%p_bnd, self%bc_configs, self%mesh, periodic_only=.true.)

      ! Leave w_gp = grad(p^{k+1/2}) for the next stage's Neumann-p extrapolation and
      ! pre-correction; without this they would see the gradient of the increment.
      call compute_pressure_gradient(self)

      ! -- Save current H as H_prev for next stage
      H_U_prev = H_U
      H_T_prev = H_T
      if (self%params%n_scalars > 0) H_C_prev = H_C(1:self%params%n_scalars, :)

    end do  ! k = 1, 3

    ! Harvest rank-local diagnostic counters into the warning log (no I/O, no MPI here).
    call warn_log_bump(WARN_TVD_FALLBACK, advection_get_tvd_fallback_count())
    call warn_log_bump(WARN_WALL_UTAU_CAP, wall_model_get_cap_count())

    self%time = time_n

    if (self%abl_config%use_sem) then
      call sem_advance(self%sem_state, dt)
      call update_sem_inflow(self)
    end if

    end associate

  end subroutine advance_rk3w

  ! Explicit RHS: H = -conv + visc + buoy
  subroutine compute_explicit_rhs(self, H_U, H_T, H_C, time_eval)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(inout) :: H_U(:,:)  ! (3, nc)
    real(rp),        intent(inout) :: H_T(:)  ! (nc)
    real(rp),        intent(inout) :: H_C(:,:)  ! (n_scalars, nc)
    real(rp),        intent(in)    :: time_eval

    integer(ip) :: nc, c, s
    real(rp)    :: force_bg(3), dir_n, hscale, hdir(3)
    real(rp)    :: t_sub

    t_sub = wall_time_now()
    nc = self%mesh%n_cells
    self%w_conv_U(:, 1:nc) = 0.0_rp
    self%w_visc_U(:, 1:nc) = 0.0_rp
    self%w_conv_T(1:nc) = 0.0_rp
    self%w_visc_T(1:nc) = 0.0_rp

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

    !$omp parallel do schedule(static)
    do c = 1, nc
      ! H = -conv + visc + Coriolis + buoyancy + forcing (grad p applied in projection)
      H_U(:, c) = self%w_conv_U(:, c) + self%w_visc_U(:, c) + force_bg
      if (self%params%use_rotation) &
        H_U(:, c) = H_U(:, c) - 2.0_rp * cross3(self%params%omega_rot, self%U(:, c))
      if (self%params%boussinesq) &
        H_U(:, c) = H_U(:, c) - self%params%beta * (self%T(c) - self%params%T_ref) * self%params%g
      H_T(c) = self%w_conv_T(c) + self%w_visc_T(c)
    end do
    !$omp end parallel do

    ! Wall-model traction/heat flux enter H so RK3 weights them like every other term.
    call add_wall_stress_to_rhs(self, H_U)
    if (self%abl_config%wall_model_id == WM_MOST) call add_wall_heat_flux_to_rhs(self, H_T)

    if (self%params%n_scalars > 0) then
      self%w_conv_C(1:nc) = 0.0_rp
      self%w_visc_C(1:nc) = 0.0_rp
      do s = 1, self%params%n_scalars
        call advection_rhs_scalar(self%C(s, :), self%U, self%mesh, self%params%scalar_scheme, self%w_conv_C(1:nc), &
                      phi_bnd=self%C_bnd(s, :), U_bnd=self%U_bnd, face_flux=self%U_face, inv_vol=self%inv_vol, &
                      grad_phi=self%grad_C(:, :, s))
        call viscous_rhs_scalar(self%C(s, :), self%params%scalar_diffusivity(s), self%nu_sgs, &
                     self%params%scalar_pr_t(s), self%grad_C(:, :, s), self%mesh, self%w_visc_C(1:nc), &
                     phi_bnd=self%C_bnd(s, :), inv_vol=self%inv_vol)
        H_C(s, :) = self%w_conv_C(1:nc) + self%w_visc_C(1:nc)
      end do
    end if
    self%perf_t_rhs_transport = self%perf_t_rhs_transport + (wall_time_now() - t_sub)

    ! Scalar-source and large-scale nudging work: not transport
    t_sub = wall_time_now()
    if (self%params%n_scalars > 0) call add_scalar_source_terms(self, H_C, time_eval)

    ! Apply Davies nudging relaxation to large-scale profiles
    if (self%abl_config%use_nudging) then
      call nudging_apply_rhs(self%nudging_data, self%U, self%T, self%cell_z, &
                             nc, time_eval, self%abl_config%nudging_tau, &
                             self%abl_config%nudging_z_start, H_U, H_T, &
                             z_full=self%abl_config%nudging_z_full, &
                             horizontal_only=self%abl_config%nudging_horizontal_only, &
                             dt=self%params%dt)
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

  !> Refresh the SEM inflow fluctuation on every ABL inlet face.
  !> The patch-mean normal component is removed so the inlet mass flux is unchanged.
  subroutine update_sem_inflow(self)
    class(FVSolver), intent(inout) :: self

    integer(ip) :: p, f, f_bnd, f_global, ierr
    real(rp)    :: u_mean(3), u_fluc(3), buf(2), gbuf(2), delta

    if (.not. allocated(self%sem_fluc_bnd)) return
    self%sem_fluc_bnd = 0.0_rp
    if (.not. self%sem_state%initialized) return

    buf = 0.0_rp
    do p = 1, self%mesh%n_patches
      if (self%mesh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        call sem_velocity_at_face(self%sem_state, self%abl_config, &
                                  self%mesh%faces(f_global)%centre, u_mean, u_fluc)
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

  !> Add the wall-model traction to H_U as an acceleration (thin-wall approximation).
  !> In H rather than applied post-predictor so it carries the RK3 alpha/beta weighting.
  subroutine add_wall_stress_to_rhs(self, H_U)
    class(FVSolver), intent(in)    :: self
    real(rp),        intent(inout) :: H_U(:,:)  ! (3, nc)

    integer(ip) :: p, qidx, f, f_bnd, f_global, o

    do qidx = 1, self%n_wallmodel_patches
      p = self%wallmodel_patches(qidx)
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        o        = self%mesh%faces(f_global)%owner
        if (o < 1 .or. o > self%mesh%n_cells) cycle
        H_U(:, o) = H_U(:, o) + self%tau_w_bnd(:, f_bnd) * &
                    self%mesh%faces(f_global)%area / self%mesh%cells(o)%volume
      end do
    end do
  end subroutine add_wall_stress_to_rhs

  !> Add the wall heat flux to H_T as a rate (mirrors add_wall_stress_to_rhs).
  subroutine add_wall_heat_flux_to_rhs(self, H_T)
    class(FVSolver), intent(in)    :: self
    real(rp),        intent(inout) :: H_T(:)  ! (nc)

    integer(ip) :: p, qidx, f, f_bnd, f_global, o

    if (.not. allocated(self%q_w_bnd)) return

    do qidx = 1, self%n_wallmodel_patches
      p = self%wallmodel_patches(qidx)
      do f = 1, self%mesh%patches(p)%n_faces
        f_bnd    = self%mesh%patches(p)%start_face + f - 1
        f_global = self%mesh%n_int_faces + f_bnd
        o        = self%mesh%faces(f_global)%owner
        if (o < 1 .or. o > self%mesh%n_cells) cycle
        ! q_w_bnd is kinematic heat flux [K m/s]; area/volume converts to a rate.
        H_T(o) = H_T(o) + self%q_w_bnd(f_bnd) * &
                 self%mesh%faces(f_global)%area / self%mesh%cells(o)%volume
      end do
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
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: fld(:)
    real(rp),        intent(in)    :: fld_bnd(:)
    real(rp),        intent(inout) :: grad(:,:)

    integer(ip) :: iter

    select case (trim(self%params%gradient_scheme))
    case ('GREEN_GAUSS')
      call green_gauss_gradient_scalar(fld, self%mesh, grad, fld_bnd)
      do iter = 2, self%params%gg_iterations
        if (self%mesh%nranks > 1 .and. self%mesh%n_ghost > 0) &
          call self%mesh%exchange_halos_vec(grad)
        self%w_gp_prev = grad
        call green_gauss_gradient_scalar(fld, self%mesh, grad, fld_bnd, self%w_gp_prev)
      end do
    case default
      call ls_gradient_scalar(fld, self%mesh, grad, fld_bnd)
    end select

    if (self%mesh%nranks > 1 .and. self%mesh%n_ghost > 0) &
      call self%mesh%exchange_halos_vec(grad)
  end subroutine gradient_of_scalar

  !> Assemble predictor face fluxes and Poisson RHS: rhs(c) = -(1/dt)*sum_f q*_f.
  subroutine compute_face_fluxes_and_rhs(self, dt, rhs_p, p_is_increment)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: dt
    real(rp),        intent(out)   :: rhs_p(:)
    logical, optional, intent(in)  :: p_is_increment  ! .true. once self%p holds phi, not p^{k-1/2}

    integer(ip) :: nc, nf, q, ib, f_bnd, f_glob
    integer(ip) :: idx, o, nb, f_mesh, pair_bnd, pair_face
    real(rp)    :: Uf(3), Uf_plain(3), flux, coeff_rc
    integer(ip) :: f_bnd_gc, pair_bnd_gc, pair_face_gc
    real(rp)    :: gc_f, d_cf_o_gc, d_cf_nb_gc, d_tot_gc, w_blend
    real(rp)    :: d_on(3), d_on_len, k_f(3), s_dot_dhat, gc_cross, cross, k_rc(3)
    real(rp)    :: eps_f(3), eps_proj
    logical     :: incr
    integer(ip) :: col, idx2, n_loops, lo, hi
    logical     :: use_color

    nc = self%mesh%n_cells
    nf = self%psolver%n_faces
    rhs_p = 0.0_rp

    self%face_flux = 0.0_rp
    self%w_flux_built = .false.

    ! The cross term must lag grad(p^{k+1/2}); after the first corrector self%p is only
    ! the increment, so the total pressure has to be rebuilt before differentiating.
    incr = .false.
    if (present(p_is_increment)) incr = p_is_increment
    if (.not. self%skip_grad_p) then
      if (incr) then
        self%w_p_tot     = self%phi_prev     + self%p
        self%w_p_tot_bnd = self%phi_prev_bnd + self%p_bnd
        call gradient_of_scalar(self, self%w_p_tot, self%w_p_tot_bnd, self%w_gp)
      else
        call gradient_of_scalar(self, self%p, self%p_bnd, self%w_gp)
      end if
    end if

    use_color = (self%psolver%n_colors > 0 .and. allocated(self%psolver%idx_color_perm) &
                 .and. allocated(self%psolver%color_offsets))
    n_loops = merge(self%psolver%n_colors, 1_ip, use_color)

    do col = 1, n_loops
      if (use_color) then
        lo = self%psolver%color_offsets(col)
        hi = self%psolver%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(idx, idx2, o, nb, f_mesh, Uf, Uf_plain, flux, coeff_rc, &
      !$omp&         f_bnd_gc, pair_bnd_gc, pair_face_gc, gc_f, d_cf_o_gc, d_cf_nb_gc, d_tot_gc, w_blend, &
      !$omp&         d_on, d_on_len, k_f, s_dot_dhat, k_rc, eps_f, eps_proj)
      do idx2 = lo, hi
        if (use_color) then
          idx = self%psolver%idx_color_perm(idx2)
        else
          idx = idx2
        end if
        if (self%psolver%skip_face(idx) /= 0) cycle

        o      = int(self%psolver%face_owner(idx), ip)
        nb     = self%psolver%face_neighbor_ext(idx)
        f_mesh = self%psolver%face_mesh_idx(idx)

        if (o < 1 .or. o > nc) cycle
        if (nb < 1 .or. nb > nc + self%mesh%n_ghost) cycle
        if (f_mesh < 1 .or. f_mesh > size(self%mesh%faces)) cycle

      ! Face velocity interpolation using geometric gc (FP&S Sect. 7.3, Kim & Choi Eq. 25)
      ! gc = distance ratio (no volume weighting): gc = d_cf_nb / (d_cf_o + d_cf_nb)
      if (f_mesh <= self%mesh%n_int_faces) then
        ! Internal faces: use precomputed geometric gc
        gc_f = self%mesh%faces(f_mesh)%gc
      else
        ! Periodic boundary faces: compute gc from d_cf distances
        gc_f = 0.5_rp
        f_bnd_gc = f_mesh - self%mesh%n_int_faces
        d_cf_o_gc = self%mesh%faces(f_mesh)%d_cf
        d_cf_nb_gc = d_cf_o_gc  ! default fallback

        ! Try periodic_pair_dcf first (precomputed distance from partner face)
        if (allocated(self%mesh%periodic_pair_dcf)) then
          if (f_bnd_gc >= 1 .and. f_bnd_gc <= size(self%mesh%periodic_pair_dcf)) then
            if (self%mesh%periodic_pair_dcf(f_bnd_gc) > SMALL) &
              d_cf_nb_gc = self%mesh%periodic_pair_dcf(f_bnd_gc)
          end if
        end if

        ! Fallback: find partner face and extract d_cf from it
        if (d_cf_nb_gc <= SMALL .and. allocated(self%mesh%periodic_pair)) then
          if (f_bnd_gc >= 1 .and. f_bnd_gc <= size(self%mesh%periodic_pair)) then
            if (self%mesh%periodic_pair(f_bnd_gc) >= 1) then
              pair_bnd_gc = self%mesh%periodic_pair(f_bnd_gc)
              pair_face_gc = self%mesh%n_int_faces + pair_bnd_gc
              if (pair_face_gc >= 1 .and. pair_face_gc <= size(self%mesh%faces)) then
                d_cf_nb_gc = self%mesh%faces(pair_face_gc)%d_cf
              end if
            end if
          end if
        end if

        ! Compute geometric gc from distance ratio
        d_tot_gc = d_cf_o_gc + d_cf_nb_gc
        if (d_tot_gc > SMALL) gc_f = d_cf_nb_gc / d_tot_gc
      end if
      ! O(h^2) skewness correction on interior faces only; eps_f=0 on periodic faces (breaks energy conservation).
      d_on = self%mesh%cells(nb)%centroid - self%mesh%cells(o)%centroid
      if (f_mesh <= self%mesh%n_int_faces) then
        ! gc_f is the owner weight, so the gc-interpolation point is x_o + (1-gc_f)*d_on.
        eps_f    = self%mesh%faces(f_mesh)%centre - (self%mesh%cells(o)%centroid + (1.0_rp - gc_f) * d_on)
        eps_proj = dot_product(eps_f, d_on) / (sum(d_on**2) + SMALL)
      else
        eps_f    = 0.0_rp
        eps_proj = 0.0_rp
      end if

      ! Face velocity interpolation: eps.grad correction smoothly blended off by face skew, plain gc for periodic.
      Uf_plain = gc_f * self%U(:, o) + (1.0_rp - gc_f) * self%U(:, nb)
      if (f_mesh <= self%mesh%n_int_faces) then
        ! Blend weight -> 1 on conforming faces, -> 0 near/above the refinement-skew tolerance.
        w_blend = refine_blend_weight(self%mesh%faces(f_mesh)%skew)
        Uf = Uf_plain + w_blend * (face_interpolate_vector(self%U, f_mesh, self%mesh, grad_U=self%grad_U) - Uf_plain)
      else
        ! Periodic: plain gc interpolation (no eps correction). Reason: gradient discontinuity across
        ! the seam means averaging gradients violates kinetic energy conservation.
        Uf = Uf_plain
      end if
      flux = dot_product(Uf, self%mesh%faces(f_mesh)%sf)

      ! Rhie-Chow: the compact jump must carry the cross part too, or it fails to vanish for a linear pressure field on non-orthogonal faces.
      ! Use gc_f directly (unclamped) for pressure gradient weight; clamping breaks O(h^2).
      coeff_rc = real(self%psolver%face_coeff(idx), rp)
      k_rc = 0.0_rp
      if (.not. self%skip_nonortho) then
        ! d_on already holds cells(nb)%centroid - cells(o)%centroid from the eps-skewness block above.
        call nonortho_kfactor_overrelaxed(self%mesh%faces(f_mesh)%sf, d_on, k_f, s_dot_dhat, d_on_len)
        if (sum(k_f**2) >= NONORTHO_SIN2_TOL * sum(self%mesh%faces(f_mesh)%sf**2)) k_rc = k_f
      end if
      flux = flux + dt * nonortho_rhie_chow(self%w_gp_rc(:, o), self%w_gp_rc(:, nb), gc_f, &
                                            self%mesh%faces(f_mesh)%sf, k_rc, coeff_rc, &
                                            self%phi_prev(nb) - self%phi_prev(o))

      self%face_flux(f_mesh) = flux
      self%w_flux_built(f_mesh) = .true.

      ! rhs = -(1/dt)*sum_f q_f : outward flux from owner -> subtract from owner,
      ! outward from owner = inward to nb -> add to nb.
      rhs_p(o) = rhs_p(o) - flux
      if (nb >= 1 .and. nb <= nc) rhs_p(nb) = rhs_p(nb) + flux
      end do
      !$omp end parallel do
    end do

    ! Non-periodic boundary faces contribute only to their owner cell.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      select case (self%bc_configs(q)%bc_id)
      case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
        flux = 0.0_rp
      case (BC_OUTLET)
        ! Dirichlet p BC: use extrapolated boundary velocity (O(h^2) on skewed meshes).
        flux = dot_product(self%U_bnd(:, f_bnd), self%mesh%faces(f_glob)%sf)
      case default
        ! BC_INLET/BC_ABL_INLET: prescribed Dirichlet velocity, so include its flux in the RHS.
        flux = dot_product(self%U_bnd(:, f_bnd), self%mesh%faces(f_glob)%sf)
      end select
      self%face_flux(f_glob) = flux
      self%w_flux_built(f_glob) = .true.
      rhs_p(o) = rhs_p(o) - flux
    end do

    ! Periodic boundary faces that are not explicitly represented in the
    ! pressure-face list must still receive a valid convective mass flux.
    do ib = 1, self%n_periodic_bfaces
      f_bnd = self%periodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      if (self%w_flux_built(f_glob)) cycle

      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle

        ! Prefer paired-face anti-symmetry when the companion face was built.
        if (allocated(self%mesh%periodic_pair)) then
          if (f_bnd >= 1 .and. f_bnd <= size(self%mesh%periodic_pair)) then
            pair_bnd = self%mesh%periodic_pair(f_bnd)
            if (pair_bnd >= 1 .and. pair_bnd <= self%mesh%n_bnd_faces) then
              pair_face = self%mesh%n_int_faces + pair_bnd
              if (pair_face >= 1 .and. pair_face <= self%mesh%n_faces) then
                if (self%w_flux_built(pair_face)) then
                  self%face_flux(f_glob) = -self%face_flux(pair_face)
                  self%w_flux_built(f_glob) = .true.
                  cycle
                end if
              end if
            end if
          end if
        end if

        ! Fallback: midpoint flux from resolved periodic partner cell.
        nb = 0
        if (allocated(self%mesh%periodic_neighbor)) then
          if (f_bnd >= 1 .and. f_bnd <= size(self%mesh%periodic_neighbor)) nb = self%mesh%periodic_neighbor(f_bnd)
        end if
        if (nb <= 0 .and. allocated(self%mesh%periodic_pair)) then
          if (f_bnd >= 1 .and. f_bnd <= size(self%mesh%periodic_pair)) then
            pair_bnd = self%mesh%periodic_pair(f_bnd)
            if (pair_bnd >= 1 .and. pair_bnd <= self%mesh%n_bnd_faces) then
              pair_face = self%mesh%n_int_faces + pair_bnd
              if (pair_face >= 1 .and. pair_face <= self%mesh%n_faces) nb = self%mesh%faces(pair_face)%owner
            end if
          end if
        end if

        if (nb >= 1 .and. nb <= nc + self%mesh%n_ghost) then
          Uf = 0.5_rp * (self%U(:, o) + self%U(:, nb))
          self%face_flux(f_glob) = dot_product(Uf, self%mesh%faces(f_glob)%sf)
        else
          self%face_flux(f_glob) = dot_product(self%U(:, o), self%mesh%faces(f_glob)%sf)
        end if
      self%w_flux_built(f_glob) = .true.
    end do

    ! Divide by dt to complete rhs = -(1/dt)*sum_f q_f.
    rhs_p = rhs_p / max(dt, SMALL)

    ! correct_face_fluxes must reuse this same lagged w_gp for the divergence to cancel.
    ! skip_nonortho means every face below fails its tolerance test, so both loops are provably no-ops.
    if (self%skip_nonortho) return

    ! use_color/n_loops carried over from the main flux loop above (same psolver idx-space).
    do col = 1, n_loops
      if (use_color) then
        lo = self%psolver%color_offsets(col)
        hi = self%psolver%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(idx, idx2, o, nb, f_mesh, k_f, gc_cross, cross)
      do idx2 = lo, hi
        if (use_color) then
          idx = self%psolver%idx_color_perm(idx2)
        else
          idx = idx2
        end if
        if (self%psolver%skip_face(idx) /= 0) cycle
        o      = int(self%psolver%face_owner(idx), ip)
        nb     = self%psolver%face_neighbor_ext(idx)
        f_mesh = self%psolver%face_mesh_idx(idx)
        if (o  < 1 .or. o  > nc) cycle
        ! Allow ghost neighbours: the matrix includes processor faces, so the
        ! deferred-correction RHS must include them too. Owner-side rhs always
        if (nb < 1 .or. nb > nc + self%mesh%n_ghost) cycle
        if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle

        if (.not. self%nonortho_active_int(idx)) cycle
        k_f = self%nonortho_k_int(:, idx)

        ! Use mesh-precomputed gc (consistent with gradient and viscous operators).
        gc_cross = self%mesh%faces(f_mesh)%gc
        cross    = nonortho_cross_flux(self%w_gp(:, o), self%w_gp(:, nb), gc_cross, k_f)

        rhs_p(o)  = rhs_p(o)  + cross
        if (nb <= nc) rhs_p(nb) = rhs_p(nb) - cross
      end do
      !$omp end parallel do
    end do

    ! Jasak deferred non-orthogonal correction on Dirichlet pressure
    ! boundary faces (BC_OUTLET). These faces are not in psolver%n_faces
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      if (self%bc_configs(q)%bc_id /= BC_OUTLET) cycle

      if (.not. self%nonortho_active_bnd(ib)) cycle
      k_f = self%nonortho_k_bnd(:, ib)

      ! Owner-cell gradient is the best available estimate of grad_p at the
      ! Dirichlet face (no off-cell information available across the BC).
      cross = dot_product(self%w_gp(:, o), k_f)
      rhs_p(o) = rhs_p(o) + cross
    end do
  end subroutine compute_face_fluxes_and_rhs

  !> Correct previously assembled face fluxes by pressure jump using the
  !> same coefficients used in the pressure operator.
  subroutine correct_face_fluxes(self, dt)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: dt

    integer(ip) :: idx, o, nb, f_mesh
    integer(ip) :: ib, f_bnd, f_glob, pair_bnd, pair_face, o_self, o_pair, q
    real(rp)    :: coeff, d_orth, sn
    real(rp)    :: k_f(3), gc, cross

    if (.not. allocated(self%face_flux)) return

    ! idx <-> f_mesh is a bijection over the assembled faces: each idx writes a distinct
    ! face_flux(f_mesh), so unlike the rhs_p scatter-adds this loop is race-free as-is.
    !$omp parallel do private(idx, o, nb, f_mesh, coeff, k_f, gc, cross) schedule(static)
    do idx = 1, self%psolver%n_faces
      if (self%psolver%skip_face(idx) /= 0) cycle
      o = int(self%psolver%face_owner(idx), ip)
      nb = self%psolver%face_neighbor_ext(idx)
      if (o < 1 .or. o > self%mesh%n_cells) cycle
      if (nb < 1 .or. nb > self%mesh%n_cells + self%mesh%n_ghost) cycle
      f_mesh = self%psolver%face_mesh_idx(idx)
      if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle
      coeff = real(self%psolver%face_coeff(idx), rp)
      self%face_flux(f_mesh) = self%face_flux(f_mesh) - dt * coeff * (self%p(nb) - self%p(o))
      ! Jasak (1996) Sect. 6.4.3: cross-diffusion using grad(p^old) from compute_face_fluxes_and_rhs.
      ! Same lagged gradient as in rhs_p -> algebraic cancellation (div = 0).
      if (self%nonortho_active_int(idx)) then
        k_f = self%nonortho_k_int(:, idx)
        ! Use precomputed projection-based gc (FPS Sect. 7.2) - consistent with
        ! the gradient, advection, and viscous operators.
        gc = self%mesh%faces(f_mesh)%gc
        cross = nonortho_cross_flux(self%w_gp(:, o), self%w_gp(:, nb), gc, k_f)
        self%face_flux(f_mesh) = self%face_flux(f_mesh) - dt * cross
      end if
    end do
    !$omp end parallel do

    ! Project non-periodic boundary face fluxes with the Dirichlet/Neumann p BC.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd  = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > self%mesh%n_cells) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      select case (self%bc_configs(q)%bc_id)
      case (BC_OUTLET)
        sn = abs(dot_product(self%mesh%faces(f_glob)%centre - self%mesh%cells(o)%centroid, &
                             self%mesh%faces(f_glob)%normal))
        d_orth = max(sn, SMALL)
        coeff = self%mesh%faces(f_glob)%area / d_orth
        self%face_flux(f_glob) = self%face_flux(f_glob) &
                                 - dt * coeff * (self%bc_configs(q)%p_val - self%p(o))
        ! Cross-diffusion correction for Dirichlet p BC (outlet); owner-side gradient only (no off-cell info).
        if (self%nonortho_active_bnd(ib)) then
          k_f = self%nonortho_k_bnd(:, ib)
          cross = dot_product(self%w_gp(:, o), k_f)
          self%face_flux(f_glob) = self%face_flux(f_glob) - dt * cross
        end if
      end select
    end do

    ! Re-sync non-canonical same-rank periodic faces by anti-symmetry.
    if (self%n_periodic_bfaces > 0 .and. allocated(self%mesh%periodic_pair)) then
      do ib = 1, self%n_periodic_bfaces
        f_bnd = self%periodic_bfaces(ib)
        if (f_bnd < 1 .or. f_bnd > size(self%mesh%periodic_pair)) cycle
        pair_bnd = self%mesh%periodic_pair(f_bnd)
        if (pair_bnd <= 0) cycle  ! cross-rank pair: handled by remote rank
        f_glob    = self%mesh%n_int_faces + f_bnd
        pair_face = self%mesh%n_int_faces + pair_bnd
        if (f_glob    < 1 .or. f_glob    > self%mesh%n_faces) cycle
        if (pair_face < 1 .or. pair_face > self%mesh%n_faces) cycle
        o_self = self%mesh%faces(f_glob   )%owner
        o_pair = self%mesh%faces(pair_face)%owner
        if (o_self < 1 .or. o_self > self%mesh%n_cells) cycle
        if (o_pair < 1 .or. o_pair > self%mesh%n_cells) cycle
        if (o_self == o_pair) cycle
        ! Canonical face is the one whose owner has the smaller global_id;
        ! mirror the non-canonical face flux from it.
        if (self%mesh%global_id(o_self) > self%mesh%global_id(o_pair)) then
          self%face_flux(f_glob) = -self%face_flux(pair_face)
        end if
      end do
    end if
  end subroutine correct_face_fluxes

  !> Correct cell velocities: u = u - dt * T^-1 . sum_f n_f g_f, g_f the jump the flux correction removes.
  ! with_cross = .false. for the pre-correction: Rhie-Chow replaces only the orthogonal part at the face.
  subroutine correct_velocity_psolver_dual(self, dt, with_cross)
    class(FVSolver), intent(inout) :: self
    real(rp),        intent(in)    :: dt
    logical,         intent(in)    :: with_cross

    integer(ip) :: nc, c, idx, o, nb, f_mesh
    integer(ip) :: ib, f_bnd, f_glob, q
    real(rp)    :: coeff, g_f, gc, sn, d_orth
    real(rp)    :: k_f(3)
    logical     :: do_cross
    integer(ip) :: col, idx2, n_loops, lo, hi
    logical     :: use_color

    do_cross = with_cross .and. (.not. self%skip_nonortho)
    nc = self%mesh%n_cells
    self%w_gp_cell(:, 1:nc) = 0.0_rp

    use_color = (self%psolver%n_colors > 0 .and. allocated(self%psolver%idx_color_perm) &
                 .and. allocated(self%psolver%color_offsets))
    n_loops = merge(self%psolver%n_colors, 1_ip, use_color)

    ! g_f is identical to the pressure jump correct_face_fluxes removes from the
    ! mass flux, so the cell and face projections share one discrete operator.
    do col = 1, n_loops
      if (use_color) then
        lo = self%psolver%color_offsets(col)
        hi = self%psolver%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = self%psolver%n_faces
      end if
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(idx, idx2, o, nb, f_mesh, coeff, g_f, gc, k_f)
      do idx2 = lo, hi
        if (use_color) then
          idx = self%psolver%idx_color_perm(idx2)
        else
          idx = idx2
        end if
        if (self%psolver%skip_face(idx) /= 0) cycle
        o      = int(self%psolver%face_owner(idx), ip)
        nb     = self%psolver%face_neighbor_ext(idx)
        f_mesh = self%psolver%face_mesh_idx(idx)
        if (o      < 1 .or. o      > nc) cycle
        if (nb     < 1 .or. nb     > nc + self%mesh%n_ghost) cycle
        if (f_mesh < 1 .or. f_mesh > self%mesh%n_faces) cycle

        coeff = real(self%psolver%face_coeff(idx), rp)
        g_f   = coeff * (self%p(nb) - self%p(o))
        if (do_cross .and. self%nonortho_active_int(idx)) then
          k_f = self%nonortho_k_int(:, idx)
          gc  = self%mesh%faces(f_mesh)%gc
          g_f = g_f + nonortho_cross_flux(self%w_gp(:, o), self%w_gp(:, nb), gc, k_f)
        end if

      ! Sf points out of the owner, so the neighbour sees (-n_f)*(-g_f): same sign.
        self%w_gp_cell(:, o) = self%w_gp_cell(:, o) + g_f * self%mesh%faces(f_mesh)%normal
        if (nb <= nc) &
          self%w_gp_cell(:, nb) = self%w_gp_cell(:, nb) + g_f * self%mesh%faces(f_mesh)%normal
      end do
      !$omp end parallel do
    end do

    ! Fixed-flux patches (walls, symmetry, inlets) get no flux correction, so their
    ! normal pressure gradient is zero; only Dirichlet-p outlets contribute.
    do ib = 1, self%n_nonperiodic_bfaces
      f_bnd  = self%nonperiodic_bfaces(ib)
      f_glob = self%mesh%n_int_faces + f_bnd
      if (f_glob < 1 .or. f_glob > self%mesh%n_faces) cycle
      o = self%mesh%faces(f_glob)%owner
      if (o < 1 .or. o > nc) cycle
      q = self%mesh%faces(f_glob)%patch_id
      if (q < 1 .or. q > self%mesh%n_patches) cycle
      if (self%bc_configs(q)%bc_id /= BC_OUTLET) cycle

      sn     = abs(dot_product(self%mesh%faces(f_glob)%centre - self%mesh%cells(o)%centroid, &
                               self%mesh%faces(f_glob)%normal))
      d_orth = max(sn, SMALL)
      coeff  = self%mesh%faces(f_glob)%area / d_orth
      g_f    = coeff * (self%bc_configs(q)%p_val - self%p(o))
      if (do_cross .and. self%nonortho_active_bnd(ib)) then
        k_f = self%nonortho_k_bnd(:, ib)
        g_f = g_f + dot_product(self%w_gp(:, o), k_f)
      end if
      self%w_gp_cell(:, o) = self%w_gp_cell(:, o) + g_f * self%mesh%faces(f_glob)%normal
    end do

    !$omp parallel do schedule(static)
    do c = 1, nc
      self%w_gp_cell(:, c) = matmul(self%recon_Tinv(:, :, c), self%w_gp_cell(:, c))
      self%U(:, c) = self%U(:, c) - dt * self%w_gp_cell(:, c)
    end do
    !$omp end parallel do
  end subroutine correct_velocity_psolver_dual

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
  subroutine anchor_pressure_to_outlet(p, msh, outlet_cell, outlet_pval)
    real(rp),    intent(inout) :: p(:)  ! pressure array (n_cells)
    type(Mesh),  intent(in)    :: msh
    logical,     intent(in)    :: outlet_cell(:)  ! (n_cells) precomputed outlet mask
    real(rp),    intent(in)    :: outlet_pval(:)  ! (n_cells) reference pressure per outlet cell

    integer(ip) :: ocell
    real(rp)    :: p_sum, p_ref_sum, p_ref, p_shift
    integer(ip) :: n_out, n_ref
    integer(ip) :: ierr
    real(rp)    :: buf(2)

    p_sum = 0.0_rp;  n_out = 0
    p_ref_sum = 0.0_rp; n_ref = 0
    p_ref = 0.0_rp

    do ocell = 1, msh%n_cells
      if (outlet_cell(ocell)) then
        p_sum     = p_sum     + p(ocell)
        p_ref_sum = p_ref_sum + outlet_pval(ocell)
        n_out     = n_out     + 1
        n_ref     = n_ref     + 1
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
      p(:) = p(:) - p_shift
    else
      ! pure-Neumann case (no outlet anywhere) - pressure is determined
      ! only up to a constant. Anchor by subtracting the global volume-
      block
        real(rp) :: vol_local, vol_global, pmean_local, pmean_global, gbuf(2)
        integer(ip) :: c
        vol_local   = 0.0_rp
        pmean_local = 0.0_rp
        do c = 1, msh%n_cells
          vol_local   = vol_local   + msh%cells(c)%volume
          pmean_local = pmean_local + p(c) * msh%cells(c)%volume
        end do
        gbuf(1) = vol_local
        gbuf(2) = pmean_local
        call MPI_Allreduce(MPI_IN_PLACE, gbuf, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, msh%comm, ierr)
        vol_global   = gbuf(1)
        pmean_global = gbuf(2)
        if (vol_global > 0.0_rp) then
          p_shift = pmean_global / vol_global
          p(:) = p(:) - p_shift
        end if
      end block
    end if
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
    class(FVSolver), intent(in) :: self
    integer(ip) :: f, o, nb, ierr, nc, patch_id, bc_id
    real(rp)    :: cfl_max_local, cfl_max_global
    real(rp)    :: diff_max_local, diff_max_global, nu_eff, area
    real(rp)    :: dt_min, safety_factor, dt_conv, dt_diff
    real(rp)    :: u_face_use, sf(3)
    logical     :: have_uface
    real(rp), parameter :: diff_cfl_limit = 1.0_rp

    nc = self%mesh%n_cells
    if (.not. allocated(dt_flux_sum_buf)) then
      allocate(dt_flux_sum_buf(nc))
      allocate(dt_diff_sum_buf(nc))
    else if (size(dt_flux_sum_buf) < nc) then
      deallocate(dt_flux_sum_buf, dt_diff_sum_buf)
      allocate(dt_flux_sum_buf(nc))
      allocate(dt_diff_sum_buf(nc))
    end if
    dt_flux_sum_buf(1:nc) = 0.0_rp
    dt_diff_sum_buf(1:nc) = 0.0_rp

    ! detect step-1 case where U_face has not been built yet. Fall back
    ! to a centred owner-side dot(U,sf) so the convective CFL reflects the
    have_uface = allocated(self%U_face) .and. self%u_face_ready

    ! Explicit stability depends on the SUM of every face's contribution into a cell
    ! (mod_advection/mod_viscous accumulate over the full face loop), not any single
    ! face alone: a per-face max undercounts polyhedral buffer cells at deep
    ! snappyHexMesh refinement interfaces, which can carry far more than 6 faces.
    !$omp parallel do default(shared) private(f, o, nb, sf, area, u_face_use, nu_eff, patch_id, bc_id)
    do f = 1, self%mesh%n_faces
      o  = self%mesh%faces(f)%owner
      nb = self%mesh%faces(f)%neighbour
      sf = self%mesh%faces(f)%sf
      area = sqrt(sf(1)**2 + sf(2)**2 + sf(3)**2)

      if (have_uface) then
        u_face_use = self%U_face(f)
      else
        u_face_use = dot_product(self%U(:, o), sf)
      end if
      if (o  >= 1 .and. o  <= nc) then
        !$omp atomic
        dt_flux_sum_buf(o)  = dt_flux_sum_buf(o)  + abs(u_face_use)
      end if
      if (nb >= 1 .and. nb <= nc) then
        !$omp atomic
        dt_flux_sum_buf(nb) = dt_flux_sum_buf(nb) + abs(u_face_use)
      end if

      ! Diffusive limit: gc-weighted nu_eff (mirrors mod_viscous face assembly).
      if (nb >= 1 .and. nb <= self%mesh%n_cells + self%mesh%n_ghost) then
        if (nb <= nc) then
          if (allocated(self%nu_sgs)) then
            nu_eff = self%params%nu &
                   + self%mesh%faces(f)%gc       * max(0.0_rp, self%nu_sgs(o))  &
                   + (1.0_rp - self%mesh%faces(f)%gc) * max(0.0_rp, self%nu_sgs(nb))
          else
            nu_eff = self%params%nu
          end if
        else
          nu_eff = self%params%nu
          if (allocated(self%nu_sgs)) nu_eff = nu_eff + max(0.0_rp, self%nu_sgs(o))
        end if
        !$omp atomic
        dt_diff_sum_buf(o) = dt_diff_sum_buf(o) + nu_eff * area / max(self%mesh%faces(f)%d_cf, SMALL)
        if (nb >= 1 .and. nb <= nc) then
          !$omp atomic
          dt_diff_sum_buf(nb) = dt_diff_sum_buf(nb) + nu_eff * area / max(self%mesh%faces(f)%d_cf, SMALL)
        end if
      else if (nb < 1 .and. o >= 1 .and. o <= nc) then
        ! No-slip/Dirichlet boundary faces apply a real explicit viscous stress at the
        ! owner-face distance (mod_viscous.f90); excluding them let near-wall cells with a
        ! thin wall-normal spacing (thinner than volume^(1/3)) run at an unstable dt.
        ! Wall-model faces are excluded: their stress is applied via a separate limiter.
        patch_id = self%mesh%faces(f)%patch_id
        bc_id = -1
        if (patch_id >= 1 .and. patch_id <= self%mesh%n_patches) bc_id = self%mesh%patches(patch_id)%bc_id
        if (bc_id /= BC_WALL_MODEL) then
          nu_eff = self%params%nu
          if (allocated(self%nu_sgs)) nu_eff = nu_eff + max(0.0_rp, self%nu_sgs(o))
          !$omp atomic
          dt_diff_sum_buf(o) = dt_diff_sum_buf(o) + nu_eff * area / max(self%mesh%faces(f)%d_cf, SMALL)
        end if
      end if
    end do
    !$omp end parallel do

    cfl_max_local  = 0.0_rp
    diff_max_local = 0.0_rp
    do o = 1, nc
      cfl_max_local  = max(cfl_max_local,  dt_flux_sum_buf(o) / (self%mesh%cells(o)%volume + SMALL))
      diff_max_local = max(diff_max_local, dt_diff_sum_buf(o) / (self%mesh%cells(o)%volume + SMALL))
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

    if (self%mesh%rank == 0) then
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
    call vreman_free()
    call dynamic_smagorinsky_free()
    call wall_model_free()
    self%mesh => null()
  end subroutine solver_destroy

end module mod_fv_solver
