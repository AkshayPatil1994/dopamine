! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Namelist-based input reader; all parameters and defaults documented in reference.input.
module mod_input_parser
  use mod_kinds
  use mod_parameters
  use mod_scalar
  use mod_bc
  use mod_fv_solver, only: solver_params_t
  use mod_sampling_types
  use mod_abl_types, only: ABLConfig, InflowGenParams, InflowConfig, InflowRecordConfig, WM_NONE, WM_LOGLAW, WM_MOST
  implicit none
  private

  public :: InputConfig
  public :: parse_input_file
  public :: resolve_wall_model

  integer(ip), parameter :: MAX_BC_PATCHES = 32
  integer(ip), parameter :: MAX_NAME_LEN   = 128
  integer(ip), parameter :: max_opts = 32

  type, public :: InputConfig
    character(256)  :: mesh_file  = ''
    character(32)   :: mesh_type  = 'gmsh'
    type(solver_params_t)  :: params
    integer(ip)            :: n_bcs = 0
    character(MAX_NAME_LEN):: bc_patch_names(MAX_BC_PATCHES) = ''
    type(BCConfig)         :: bc_configs(MAX_BC_PATCHES)
    character(32)          :: ic_type = 'uniform'
    real(rp)               :: ic_uref = 0.0_rp
    real(rp)               :: ic_lref = 0.0_rp
    real(rp)               :: ic_noise_level = 0.0_rp
    integer(ip)            :: ic_noise_seed = 0  ! 0 = seed from the clock (not reproducible)
    integer(ip)            :: wall_normal_direction = 3
    character(32)          :: output_prefix = 'dopamine'
    character(16)          :: gradient_scheme = 'LEAST_SQUARES'
    logical                :: hdf5_use_file_locking = .false.
    character(256)         :: warn_log_file = 'solver_warnings.txt'
    integer(ip)            :: warn_flush_interval = 100
    integer(ip)            :: stats_global_frequency = 0
    integer(ip)            :: stats_slice_frequency = 0
    integer(ip)            :: stats_probe_frequency = 0
    integer(ip)            :: stats_line_frequency = 0
    integer(ip)            :: stats_curve_frequency = 0
    real(rp)               :: stats_intersection_decimate_spacing = 0.0_rp
    logical                :: stats_sort_intersection_points = .true.
    integer(ip)            :: n_planes = 0
    integer(ip)            :: n_lines = 0
    integer(ip)            :: n_curves = 0
    integer(ip)            :: n_surfaces = 0
    type(PlaneSampleDef)   :: planes(MAX_SAMPLE_OBJECTS)
    type(LineSampleDef)    :: lines(MAX_SAMPLE_OBJECTS)
    type(CurveSampleDef)   :: curves(MAX_SAMPLE_OBJECTS)
    type(SurfaceSampleDef) :: surfaces(MAX_SAMPLE_OBJECTS)
    integer(ip)            :: n_scalars = 0
    integer(ip)            :: n_scalar_sources = 0
    character(64)          :: scalar_names(MAX_PASSIVE_SCALARS) = ''
    real(rp)               :: scalar_diffusivity(MAX_PASSIVE_SCALARS) = 0.0_rp
    real(rp)               :: scalar_pr_t(MAX_PASSIVE_SCALARS) = 1.0_rp
    type(ScalarSourceDef)  :: scalar_sources(MAX_SCALAR_SOURCES)
    character(len=:), allocatable :: petsc_options(:)
    ! Checkpoint / Restart
    logical        :: restart_enabled      = .false.
    character(256) :: restart_dir          = ''
    integer(ip)    :: restart_step         = -1
    integer(ip)    :: checkpoint_interval  = 0
    ! ABL parameters
    logical    :: use_abl           = .false.
    real(rp)   :: abl_u_ref         = 8.0_rp
    real(rp)   :: abl_z_ref         = 100.0_rp
    real(rp)   :: abl_z0_m          = 0.1_rp
    real(rp)   :: abl_z0_h          = 0.01_rp
    real(rp)   :: abl_T_surface     = 300.0_rp
    real(rp)   :: abl_T_ref         = 300.0_rp
    real(rp)   :: abl_L_obukhov     = 0.0_rp
    logical    :: abl_compute_L     = .true.
    real(rp)   :: abl_wind_dir_deg  = 270.0_rp
    logical    :: abl_use_sem       = .true.
    integer(ip) :: abl_n_eddies     = 500
    real(rp)   :: abl_sem_sigma     = 0.0_rp
    character(256) :: abl_sem_profile_file = ''
    logical    :: abl_use_nudging   = .false.
    character(256) :: abl_nudging_file = ''
    real(rp)   :: abl_nudging_tau   = 1800.0_rp
    real(rp)   :: abl_nudging_z_start = 0.0_rp
    real(rp)   :: abl_nudging_z_full  = 0.0_rp
    logical    :: abl_nudging_horizontal_only = .true.
    type(ABLConfig) :: abl_config  ! packed form passed to the solver
    type(InflowConfig)       :: inflow        ! successor inflow source
    type(InflowGenParams)    :: inflow_gen    ! dopamine-inflow-gen
    type(InflowRecordConfig) :: inflow_rec    ! precursor plane recorder

  end type InputConfig


contains
  !> Read all namelists from the input file.
  !> Missing namelists are silently skipped (defaults remain in effect).
  subroutine parse_input_file(filename, config)
    character(len=*), intent(in)    :: filename
    type(InputConfig), intent(inout) :: config

    integer(ip) :: iunit, ios

    open(newunit=iunit, file=trim(filename), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) error stop 'parse_input_file: cannot open '//trim(filename)

    call read_domain(iunit, config)
    rewind(iunit)
    call read_numerics(iunit, config)
    rewind(iunit)
    call read_stats(iunit, config)
    rewind(iunit)
    call read_physics(iunit, config)
    rewind(iunit)
    call read_abl(iunit, config)
    rewind(iunit)
    call read_inflow(iunit, config)
    rewind(iunit)
    call read_scalars(iunit, config)
    rewind(iunit)
    call read_bc(iunit, config)
    rewind(iunit)
    call read_scalar_sources(iunit, config)
    rewind(iunit)
    call read_initial_condition(iunit, config)
    rewind(iunit)
    call read_sampling_planes(iunit, config)
    rewind(iunit)
    call read_sampling_lines(iunit, config)
    rewind(iunit)
    call read_sampling_curves(iunit, config)
    rewind(iunit)
    call read_sampling_surfaces(iunit, config)
    rewind(iunit)
    call read_output(iunit, config)
    rewind(iunit)
    call read_checkpoint(iunit, config)
    rewind(iunit)
    call read_PetScParams(iunit, config)

  end subroutine parse_input_file

  subroutine read_domain(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    character(256) :: mesh_file = ''
    character(32)  :: mesh_type = 'gmsh'
    integer(ip)    :: ios

    namelist /dopamine_domain/ mesh_file, mesh_type

    read(iunit, nml=dopamine_domain, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_domain, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    if (len_trim(mesh_file) > 0) config%mesh_file = mesh_file
    config%mesh_type = mesh_type

    ! Validate mesh_type
    select case (trim(config%mesh_type))
    case ('gmsh', 'openfoam', 'foam', 'preprocessed')
      ! Valid types
    case default
      write(*, '(A,A)') 'WARNING: unrecognized mesh_type: ', trim(config%mesh_type)
    end select
  end subroutine read_domain

  subroutine read_numerics(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    ! Namelist variables with defaults
    character(16)  :: advection_scheme      = 'CENTRAL'
    character(16)  :: temperature_scheme    = 'VANLEER_TVD'
    character(16)  :: gradient_scheme  = 'LEAST_SQUARES'
    real(rp)       :: dt               = 1.0e-3_rp
    real(rp)       :: t_end            = 100.0_rp
    integer(ip)    :: nstep_max        = 0_ip
    real(rp)       :: cfl_max          = 0.4_rp
    logical        :: adapt_dt         = .false.
    integer(ip)    :: solver_stats_frequency = 0
    integer(ip)    :: calc_stats_frequency = 0
    real(rp)       :: p_tol            = 1.0e-10_rp
    integer(ip)    :: p_max_iter       = 1000
    real(rp)       :: p_target_divU    = 1.0e-08_rp
    character(32)  :: pressure_solver  = 'boomeramg'  ! boomeramg | gamg
    logical        :: use_rotational_projection = .true.
    logical        :: use_sgs = .true.
    integer(ip)    :: gg_iterations = 2
    integer(ip)    :: n_nonortho_corr = 0  ! 0 = auto (1 if mesh orthogonal, else 2)
    real(rp)       :: pressure_min_cos_theta = 0.0_rp  ! > 0 floors cos(theta) in the Poisson operator
    logical        :: rhie_chow_dual = .false.  ! consistent Rhie-Chow on non-orthogonal faces
    character(16)  :: time_scheme = 'RK3W'  ! legacy; use time_integrator instead
    character(16)  :: time_integrator = 'RK3W'  ! RK3W is the only integrator
    logical        :: use_skew_symmetric_advection = .false.
    character(16)  :: sgs_model = 'VREMAN'  ! VREMAN | SMAGORINSKY | DYNAMIC_SMAG | AMD
    real(rp)       :: c_amd = AMD_C_DEFAULT         ! AMD constant
    logical        :: amd_wall_gradient = .true.  ! AMD: log-law wall-normal gradient in wall-model cells
    real(rp)       :: wall_model_filter_tf_fac = 0.0_rp  ! > 0: adaptive filter T_f = f*y/(kappa*u_tau)
    real(rp)       :: wall_model_height = 0.0_rp       ! <= 0 uses the wall-adjacent cell
    real(rp)       :: wall_model_filter_time = 0.0_rp  ! <= 0 disables time filtering
    logical        :: debug_diagnostics = .false.  ! verbose runtime diagnostics (off by default)
    integer(ip)    :: ios

        namelist /dopamine_numerics/ advection_scheme, temperature_scheme, dt, t_end, &
           cfl_max, adapt_dt, nstep_max, &
             solver_stats_frequency, calc_stats_frequency, &
                               p_tol, p_max_iter, p_target_divU, gradient_scheme, pressure_solver, &
                               use_rotational_projection, gg_iterations, use_sgs, time_scheme, n_nonortho_corr, pressure_min_cos_theta, rhie_chow_dual, &
                               time_integrator, use_skew_symmetric_advection, sgs_model, &
                               wall_model_height, wall_model_filter_time, debug_diagnostics, &
                               c_amd, amd_wall_gradient, wall_model_filter_tf_fac

    read(iunit, nml=dopamine_numerics, iostat=ios)
    if (ios > 0) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_numerics, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    ! Names are matched case-insensitively and an unknown one stops the run: each used to fall back to a
    ! default without a word (or with a warning that scrolls past), so a typo silently changed the physics.
    advection_scheme   = upper_case(adjustl(advection_scheme))
    temperature_scheme = upper_case(adjustl(temperature_scheme))
    gradient_scheme    = upper_case(adjustl(gradient_scheme))
    sgs_model          = upper_case(adjustl(sgs_model))
    pressure_solver    = lower_case(adjustl(pressure_solver))
    select case (trim(gradient_scheme))
    case ('LEAST_SQUARES', 'GREEN_GAUSS')
    case default
      call input_fatal('gradient_scheme', gradient_scheme, 'LEAST_SQUARES, GREEN_GAUSS')
    end select
    select case (trim(sgs_model))
    case ('VREMAN', 'SMAGORINSKY', 'DYNAMIC_SMAG', 'AMD')
    case default
      call input_fatal('sgs_model', sgs_model, 'VREMAN, SMAGORINSKY, DYNAMIC_SMAG, AMD')
    end select
    select case (trim(pressure_solver))
    case ('boomeramg', 'gamg')
    case default
      call input_fatal('pressure_solver', pressure_solver, 'boomeramg, gamg')
    end select

    associate(p => config%params)
      select case (trim(advection_scheme))
      case ('KOREN_TVD', 'KOREN')
        p%scheme = SCHEME_KOREN_TVD
      case ('VANLEER_TVD', 'VANLEER')
        p%scheme = SCHEME_VANLEER_TVD
      case ('QUICK')
        p%scheme = SCHEME_QUICK
      case ('KEP', 'SKEW');    p%scheme = SCHEME_KEP
      case ('CDS', 'CENTRAL'); p%scheme = SCHEME_CDS
      case default
        call input_fatal('advection_scheme', advection_scheme, 'KOREN_TVD, VANLEER_TVD, QUICK, KEP, CDS')
      end select

      select case (trim(temperature_scheme))
      case ('KOREN_TVD', 'KOREN')
        p%temperature_scheme = SCHEME_KOREN_TVD
      case ('VANLEER_TVD', 'VANLEER')
        p%temperature_scheme = SCHEME_VANLEER_TVD
      case ('QUICK')
        p%temperature_scheme = SCHEME_QUICK
      case ('CDS', 'CENTRAL')
        p%temperature_scheme = SCHEME_CDS
      case default
        call input_fatal('temperature_scheme', temperature_scheme, 'KOREN_TVD, VANLEER_TVD, QUICK, CDS')
      end select

      p%dt             = dt
      p%gradient_scheme = gradient_scheme
      config%gradient_scheme = gradient_scheme
      p%t_end          = t_end
      p%nstep_max      = nstep_max
      p%cfl_max        = cfl_max
      p%adapt_dt       = adapt_dt
      if (config%stats_global_frequency <= 0) then
        if (solver_stats_frequency > 0) then
          config%stats_global_frequency = solver_stats_frequency
        else
          config%stats_global_frequency = max(0_ip, calc_stats_frequency)
        end if
      end if
      p%p_tol          = p_tol
      p%p_max_iter     = p_max_iter
      p%p_target_divU  = p_target_divU
      p%pressure_solver = pressure_solver
      p%use_rotational_projection = use_rotational_projection
      p%use_sgs = use_sgs
      p%gg_iterations = gg_iterations
      p%n_nonortho_corr = max(0_ip, n_nonortho_corr)  ! 0 kept as the auto sentinel
      p%pressure_min_cos_theta = min(max(0.0_rp, pressure_min_cos_theta), 1.0_rp)
      p%rhie_chow_dual = rhie_chow_dual
      p%time_integrator = time_integrator  ! RK3W 
      ! Parsed for input compatibility only; KEP already supplies the energy-conserving form.
      p%use_skew_symmetric_advection = use_skew_symmetric_advection
      if (use_skew_symmetric_advection) &
        write(*,'(A)') ' [input] use_skew_symmetric_advection is deprecated and ignored; use advection_scheme = KEP'
      p%sgs_model = sgs_model  ! VREMAN (default), SMAGORINSKY, DYNAMIC_SMAG, or AMD
      p%c_amd = c_amd
      p%amd_wall_gradient = amd_wall_gradient
      p%wall_model_height = wall_model_height
      p%wall_model_filter_time = wall_model_filter_time
      p%wall_model_filter_tf_fac = wall_model_filter_tf_fac
      p%debug_diagnostics = debug_diagnostics
    end associate
    ! Propagate PETSc options if present
    if (allocated(config%petsc_options)) then
      if (allocated(config%params%petsc_options)) deallocate(config%params%petsc_options)
      allocate(character(len=256) :: config%params%petsc_options(size(config%petsc_options)))
      config%params%petsc_options = config%petsc_options
    end if
  end subroutine read_numerics

  subroutine read_stats(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: solver_stats_frequency
    integer(ip) :: calc_stats_frequency
    integer(ip) :: slice_write_frequency
    integer(ip) :: probe_write_frequency
    integer(ip) :: line_write_frequency
    integer(ip) :: curve_write_frequency
    real(rp)    :: intersection_decimate_spacing
    logical     :: sort_intersection_points
    integer(ip) :: ios

    namelist /dopamine_stats/ solver_stats_frequency, calc_stats_frequency, &
                  slice_write_frequency, probe_write_frequency, &
                  line_write_frequency, curve_write_frequency, &
                  intersection_decimate_spacing, sort_intersection_points

    solver_stats_frequency = config%stats_global_frequency
    calc_stats_frequency = config%stats_global_frequency
    slice_write_frequency = config%stats_slice_frequency
    probe_write_frequency = config%stats_probe_frequency
    line_write_frequency = config%stats_line_frequency
    curve_write_frequency = config%stats_curve_frequency
    intersection_decimate_spacing = config%stats_intersection_decimate_spacing
    sort_intersection_points = config%stats_sort_intersection_points

    read(iunit, nml=dopamine_stats, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_stats, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    if (solver_stats_frequency > 0) then
      config%stats_global_frequency = solver_stats_frequency
    else
      config%stats_global_frequency = max(0_ip, calc_stats_frequency)
    end if

    config%stats_slice_frequency = max(0_ip, slice_write_frequency)
    config%stats_probe_frequency = max(0_ip, probe_write_frequency)
    config%stats_line_frequency = max(0_ip, line_write_frequency)
    config%stats_curve_frequency = max(0_ip, curve_write_frequency)
    config%stats_intersection_decimate_spacing = max(0.0_rp, intersection_decimate_spacing)
    config%stats_sort_intersection_points = sort_intersection_points

  end subroutine read_stats

  subroutine read_physics(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    real(rp)    :: nu         = 1.0e-5_rp
    real(rp)    :: alpha_T    = 1.4e-5_rp
    real(rp)    :: beta       = 3.4e-3_rp
    real(rp)    :: T_ref      = 300.0_rp
    real(rp)    :: Pr_t       = 0.85_rp
    real(rp)    :: gx         = 0.0_rp
    real(rp)    :: gy         = 0.0_rp
    real(rp)    :: gz         = 0.0_rp
    real(rp)    :: fx         = 0.0_rp
    real(rp)    :: fy         = 0.0_rp
    real(rp)    :: fz         = 0.0_rp
    logical     :: constant_mass_flux_forcing = .false.
    real(rp)    :: Ux_constant = 0.0_rp
    real(rp)    :: Uy_constant = 0.0_rp
    real(rp)    :: Uz_constant = 0.0_rp
    logical     :: use_rotation = .false.
    real(rp)    :: omega_x    = 0.0_rp
    real(rp)    :: omega_y    = 0.0_rp
    real(rp)    :: omega_z    = 0.0_rp
    logical     :: use_harmonic_forcing = .false.
    real(rp)    :: forcing_amp   = 0.0_rp
    real(rp)    :: forcing_omega = 0.0_rp
    real(rp)    :: forcing_phase = 0.0_rp
    real(rp)    :: forcing_dir_x = 1.0_rp
    real(rp)    :: forcing_dir_y = 0.0_rp
    real(rp)    :: forcing_dir_z = 0.0_rp
    logical     :: boussinesq = .false.
    integer(ip) :: ios

    namelist /dopamine_physics/ nu, alpha_T, beta, T_ref, Pr_t, &
                  gx, gy, gz, fx, fy, fz, &
                  constant_mass_flux_forcing, Ux_constant, Uy_constant, Uz_constant, &
                  use_rotation, omega_x, omega_y, omega_z, &
                  use_harmonic_forcing, forcing_amp, forcing_omega, forcing_phase, &
                  forcing_dir_x, forcing_dir_y, forcing_dir_z, boussinesq

    read(iunit, nml=dopamine_physics, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_physics, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    ! Constant pressure-gradient (fx/fy/fz) and constant mass-flux forcing are mutually exclusive.
    if (constant_mass_flux_forcing .and. &
        (abs(fx) > SMALL .or. abs(fy) > SMALL .or. abs(fz) > SMALL)) &
      error stop 'physics: fx/fy/fz must be 0 with constant_mass_flux_forcing'

    associate(p => config%params)
      p%nu         = nu
      p%alpha_T    = alpha_T
      p%beta       = beta
      p%T_ref      = T_ref
      p%Pr_t       = Pr_t
      p%g          = [gx, gy, gz]
      p%f_const    = [fx, fy, fz]
      p%constant_mass_flux_forcing = constant_mass_flux_forcing
      p%U_bulk_target = [Ux_constant, Uy_constant, Uz_constant]
      p%use_rotation = use_rotation
      p%omega_rot  = [omega_x, omega_y, omega_z]
      p%use_harmonic_forcing = use_harmonic_forcing
      p%forcing_amp   = forcing_amp
      p%forcing_omega = forcing_omega
      p%forcing_phase = forcing_phase
      p%forcing_dir = [forcing_dir_x, forcing_dir_y, forcing_dir_z]
      p%boussinesq = boussinesq
    end associate
  end subroutine read_physics

  subroutine read_bc(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip)    :: n_patches = 0
    integer(ip)    :: ios, p

    character(MAX_NAME_LEN) :: patch_name(MAX_BC_PATCHES)
    character(32)           :: patch_type(MAX_BC_PATCHES)
    real(rp)                :: patch_U(MAX_BC_PATCHES, 3)
    real(rp)                :: patch_T(MAX_BC_PATCHES)
    real(rp)                :: patch_p(MAX_BC_PATCHES)
    real(rp)                :: patch_z0_m(MAX_BC_PATCHES)
    real(rp)                :: patch_z0_h(MAX_BC_PATCHES)
    real(rp)                :: patch_scalar_value(MAX_BC_PATCHES, MAX_PASSIVE_SCALARS)
    real(rp)                :: patch_scalar_flux(MAX_BC_PATCHES, MAX_PASSIVE_SCALARS)
    real(rp)                :: patch_scalar_flux_omega(MAX_BC_PATCHES, MAX_PASSIVE_SCALARS)
    real(rp)                :: patch_scalar_flux_phase(MAX_BC_PATCHES, MAX_PASSIVE_SCALARS)
    character(16)           :: patch_scalar_mode(MAX_BC_PATCHES, MAX_PASSIVE_SCALARS)

    namelist /dopamine_bc/ n_patches, patch_name, patch_type, &
               patch_U, patch_T, patch_p, patch_z0_m, patch_z0_h, &
               patch_scalar_value, patch_scalar_flux, &
               patch_scalar_flux_omega, patch_scalar_flux_phase, patch_scalar_mode

    patch_name = ''
    patch_type = 'zeroGradient'
    patch_U    = 0.0_rp
    patch_T    = config%params%T_ref
    patch_p    = 0.0_rp
    ! Zero means "inherit the global &dopamine_abl roughness".
    patch_z0_m = 0.0_rp
    patch_z0_h = 0.0_rp
    patch_scalar_value = 0.0_rp
    patch_scalar_flux = 0.0_rp
    patch_scalar_flux_omega = 0.0_rp
    patch_scalar_flux_phase = 0.0_rp
    patch_scalar_mode = 'zeroGradient'

    read(iunit, nml=dopamine_bc, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_bc, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    if (n_patches > MAX_BC_PATCHES) then
      write(*, '(A,I0,A,I0)') 'ERROR: n_patches = ', n_patches, ' exceeds the supported maximum ', MAX_BC_PATCHES
      error stop 'input: too many boundary patches'
    end if
    config%n_bcs = max(n_patches, 0_ip)
    do p = 1, config%n_bcs
      config%bc_patch_names(p)     = trim(patch_name(p))
      config%bc_configs(p)%bc_id   = bc_name_to_id(trim(patch_type(p)))
      config%bc_configs(p)%U_val   = patch_U(p, :)
      config%bc_configs(p)%T_val   = patch_T(p)
      config%bc_configs(p)%p_val   = patch_p(p)
      ! &dopamine_abl is parsed before &dopamine_bc, so the ABL inlet sees the user's values.
      config%bc_configs(p)%abl_u_ref     = config%abl_u_ref
      config%bc_configs(p)%abl_z_ref     = config%abl_z_ref
      ! Per-patch roughness wins; 0 inherits the global value, <0 defers to kB^-1 in the wall model.
      config%bc_configs(p)%abl_z0_m = config%abl_z0_m
      if (patch_z0_m(p) > 0.0_rp) config%bc_configs(p)%abl_z0_m = patch_z0_m(p)
      config%bc_configs(p)%abl_z0_h = config%abl_z0_h
      if (patch_z0_h(p) < 0.0_rp .or. patch_z0_h(p) > 0.0_rp) &
        config%bc_configs(p)%abl_z0_h = patch_z0_h(p)
      config%bc_configs(p)%abl_T_surface = config%abl_T_surface
      config%bc_configs(p)%abl_T_ref     = config%abl_T_ref
      config%bc_configs(p)%abl_wind_dir  = config%abl_wind_dir_deg
      if (config%n_scalars > 0) then
        allocate(config%bc_configs(p)%scalar_bc_mode(config%n_scalars))
        allocate(config%bc_configs(p)%scalar_value(config%n_scalars))
        allocate(config%bc_configs(p)%scalar_flux(config%n_scalars))
        allocate(config%bc_configs(p)%scalar_flux_omega(config%n_scalars))
        allocate(config%bc_configs(p)%scalar_flux_phase(config%n_scalars))
        config%bc_configs(p)%scalar_value = patch_scalar_value(p, 1:config%n_scalars)
        config%bc_configs(p)%scalar_flux = patch_scalar_flux(p, 1:config%n_scalars)
        config%bc_configs(p)%scalar_flux_omega = patch_scalar_flux_omega(p, 1:config%n_scalars)
        config%bc_configs(p)%scalar_flux_phase = patch_scalar_flux_phase(p, 1:config%n_scalars)
        config%bc_configs(p)%scalar_bc_mode = SCALAR_BC_ZERO_GRADIENT
        config%bc_configs(p)%scalar_bc_mode(:) = scalar_bc_mode_from_name(patch_scalar_mode(p, 1))
        call map_scalar_bc_modes(config%bc_configs(p)%scalar_bc_mode, patch_scalar_mode(p, :), config%n_scalars)
      end if
    end do
  end subroutine read_bc

  subroutine read_scalars(iunit, config)
    integer(ip), intent(in) :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_scalars, ios, s
    character(16) :: scalar_tvd_scheme
    character(64) :: scalar_name(MAX_PASSIVE_SCALARS)
    real(rp) :: scalar_diffusivity(MAX_PASSIVE_SCALARS)
    real(rp) :: scalar_pr_t(MAX_PASSIVE_SCALARS)

    namelist /dopamine_scalars/ n_scalars, scalar_name, scalar_diffusivity, scalar_pr_t, scalar_tvd_scheme

    n_scalars = 0
    scalar_name = ''
    scalar_diffusivity = 0.0_rp
    scalar_pr_t = 1.0_rp
    scalar_tvd_scheme = 'VANLEER_TVD'

    read(iunit, nml=dopamine_scalars, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_scalars, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_scalars = min(max(0_ip, n_scalars), MAX_PASSIVE_SCALARS)
    if (config%n_scalars <= 0) then
      config%params%n_scalars = 0
      return
    end if

    config%params%n_scalars = config%n_scalars
    select case (trim(scalar_tvd_scheme))
    case ('KOREN_TVD', 'KOREN')
      config%params%scalar_scheme = SCHEME_KOREN_TVD
    case ('VANLEER_TVD', 'VANLEER')
      config%params%scalar_scheme = SCHEME_VANLEER_TVD
    case ('QUICK')
      config%params%scalar_scheme = SCHEME_QUICK
    case ('CDS', 'CENTRAL')
      config%params%scalar_scheme = SCHEME_CDS
    case default
      config%params%scalar_scheme = SCHEME_VANLEER_TVD
      write(*,'(A,A,A)') 'WARNING: unsupported scalar_tvd_scheme=', trim(scalar_tvd_scheme), '; using VANLEER_TVD'
    end select
    allocate(config%params%scalar_names(config%n_scalars))
    allocate(config%params%scalar_diffusivity(config%n_scalars))
    allocate(config%params%scalar_pr_t(config%n_scalars))

    do s = 1, config%n_scalars
      config%scalar_names(s) = trim(scalar_name(s))
      if (len_trim(config%scalar_names(s)) == 0) write(config%scalar_names(s), '(A,I0)') 'scalar_', s
      config%scalar_diffusivity(s) = max(0.0_rp, scalar_diffusivity(s))
      config%scalar_pr_t(s) = max(1.0e-6_rp, scalar_pr_t(s))
      config%params%scalar_names(s) = config%scalar_names(s)
      config%params%scalar_diffusivity(s) = config%scalar_diffusivity(s)
      config%params%scalar_pr_t(s) = config%scalar_pr_t(s)
    end do
  end subroutine read_scalars

  subroutine read_abl(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    logical    :: use_abl         = .false.
    character(16) :: wall_model   = ''
    real(rp)   :: u_ref           = 8.0_rp
    real(rp)   :: z_ref           = 100.0_rp
    real(rp)   :: z0_m            = 0.1_rp
    real(rp)   :: z0_h            = 0.01_rp
    real(rp)   :: T_surface       = 300.0_rp
    real(rp)   :: T_ref           = 300.0_rp
    real(rp)   :: L_obukhov       = 0.0_rp
    logical    :: compute_L       = .true.
    real(rp)   :: wind_dir_deg    = 270.0_rp
    logical    :: use_sem         = .true.
    integer(ip) :: n_eddies       = 500
    real(rp)   :: sem_sigma       = 0.0_rp
    character(16) :: sem_type     = 'dfsem'
    character(256) :: sem_profile_file = ''
    character(8)   :: sem_profile_format = 'stress'
    character(8)   :: sem_mean_law    = 'log'
    real(rp)   :: sem_alpha       = 0.2_rp
    real(rp)   :: sem_Iu_ref      = 0.0_rp
    real(rp)   :: sem_Iv_ratio    = 0.75_rp
    real(rp)   :: sem_Iw_ratio    = 0.5_rp
    real(rp)   :: sem_rho_uw      = -0.3_rp
    integer(ip) :: sem_seed       = 0
    logical    :: use_nudging     = .false.
    character(256) :: nudging_file = ''
    real(rp)   :: nudging_tau     = 1800.0_rp
    real(rp)   :: nudging_z_start = 0.0_rp
    real(rp)   :: nudging_z_full  = 0.0_rp
    logical    :: nudging_horizontal_only = .true.
    integer(ip) :: ios

    namelist /dopamine_abl/ use_abl, wall_model, u_ref, z_ref, z0_m, z0_h, T_surface, T_ref, &
                            L_obukhov, compute_L, wind_dir_deg, use_sem, n_eddies, &
                            sem_sigma, sem_type, sem_profile_file, sem_profile_format, sem_mean_law, sem_alpha, &
                            sem_Iu_ref, sem_Iv_ratio, sem_Iw_ratio, sem_rho_uw, sem_seed, &
                            use_nudging, nudging_file, nudging_tau, nudging_z_start, &
                            nudging_z_full, nudging_horizontal_only

    read(iunit, nml=dopamine_abl, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_abl, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%use_abl           = use_abl
    config%abl_u_ref         = u_ref
    config%abl_z_ref         = z_ref
    config%abl_z0_m          = z0_m
    config%abl_z0_h          = z0_h
    config%abl_T_surface     = T_surface
    config%abl_T_ref         = T_ref
    config%abl_L_obukhov     = L_obukhov
    config%abl_compute_L     = compute_L
    config%abl_wind_dir_deg  = wind_dir_deg
    config%abl_use_sem       = use_sem
    config%abl_n_eddies      = n_eddies
    config%abl_sem_sigma     = sem_sigma
    config%abl_sem_profile_file = sem_profile_file
    config%abl_use_nudging   = use_nudging
    config%abl_nudging_file  = nudging_file
    config%abl_nudging_tau   = nudging_tau
    config%abl_nudging_z_start = nudging_z_start
    config%abl_nudging_z_full  = nudging_z_full
    config%abl_nudging_horizontal_only = nudging_horizontal_only

    ! Pack into the ABLConfig the solver actually uses.
    config%abl_config%use_abl          = use_abl
    config%abl_config%wall_model       = wall_model
    config%abl_config%wall_model_id    = resolve_wall_model(wall_model, use_abl)
    config%abl_config%u_ref            = u_ref
    config%abl_config%z_ref            = z_ref
    config%abl_config%z0_m             = z0_m
    config%abl_config%z0_h             = z0_h
    config%abl_config%T_surface        = T_surface
    config%abl_config%T_ref            = T_ref
    config%abl_config%L_obukhov        = L_obukhov
    config%abl_config%compute_L        = compute_L
    config%abl_config%wind_dir_deg     = wind_dir_deg
    config%abl_config%use_sem          = use_sem
    config%abl_config%n_eddies         = n_eddies
    config%abl_config%sem_sigma        = sem_sigma
    config%abl_config%sem_type         = lower_case(adjustl(sem_type))
    if (trim(config%abl_config%sem_type) /= 'dfsem' .and. trim(config%abl_config%sem_type) /= 'jarrin') then
      write(*,'(A)') 'ERROR: sem_type must be ''dfsem'' or ''jarrin'', got '''//trim(sem_type)//''''
      error stop 'read_abl: unknown sem_type'
    end if
    config%abl_config%sem_profile_file = sem_profile_file
    config%abl_config%sem_profile_format = sem_profile_format
    config%abl_config%sem_mean_law     = sem_mean_law
    config%abl_config%sem_alpha        = sem_alpha
    config%abl_config%sem_Iu_ref       = sem_Iu_ref
    config%abl_config%sem_Iv_ratio     = sem_Iv_ratio
    config%abl_config%sem_Iw_ratio     = sem_Iw_ratio
    config%abl_config%sem_rho_uw       = sem_rho_uw
    config%abl_config%sem_seed         = sem_seed
    call validate_sem_inputs(config%abl_config)
    config%abl_config%use_nudging      = use_nudging
    config%abl_config%nudging_file     = nudging_file
    config%abl_config%nudging_tau      = nudging_tau
    config%abl_config%nudging_z_start  = nudging_z_start
    config%abl_config%nudging_z_full   = nudging_z_full
    config%abl_config%nudging_horizontal_only = nudging_horizontal_only

  end subroutine read_abl

  !> Reject SEM inputs that would silently produce nonsense.
  subroutine validate_sem_inputs(c)
    type(ABLConfig), intent(in) :: c

    if (trim(c%sem_profile_format) /= 'stress' .and. trim(c%sem_profile_format) /= 'ti') &
      error stop "sem_profile_format must be 'stress' or 'ti'"
    if (trim(c%sem_mean_law) /= 'log' .and. trim(c%sem_mean_law) /= 'power') &
      error stop "sem_mean_law must be 'log' or 'power'"
    if (c%sem_Iu_ref < 0.0_rp .or. c%sem_Iv_ratio < 0.0_rp .or. c%sem_Iw_ratio < 0.0_rp) &
      error stop 'sem intensities must be non-negative'
    if (abs(c%sem_rho_uw) >= 1.0_rp) error stop 'sem_rho_uw must lie in (-1, 1)'
    if (c%sem_sigma < 0.0_rp) error stop 'sem_sigma must be >= 0'
  end subroutine validate_sem_inputs

  subroutine read_inflow(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    ! &dopamine_inflow: successor source
    character(16)  :: inflow_source = 'sem'
    character(256) :: inflow_file = ''
    real(rp)       :: inflow_t_offset = 0.0_rp
    real(rp)       :: inflow_crossfade_time = 0.0_rp
    integer(ip)    :: inflow_seed = 1
    character(256) :: inflow_export_grid = ''
    ! &dopamine_inflow_gen
    character(256) :: gen_file = 'inflow_planes.nc', gen_grid_file = ''
    real(rp)       :: gen_W = 0.0_rp, gen_y0 = 0.0_rp, gen_z_min = 0.0_rp, gen_z_max = 0.0_rp
    real(rp)       :: gen_dt = 0.0_rp, gen_n_turnovers = 2.0_rp, gen_T_turnover = 0.0_rp, gen_n_sigma = 3.0_rp
    integer(ip)    :: gen_ny = 0, gen_nz = 0, gen_n_segments = 4, gen_seed = 1
    logical        :: gen_verify = .true., gen_decorrelate = .true.
    ! &dopamine_inflow_record
    logical        :: record_enabled = .false.
    character(256) :: record_file = 'inflow_record.nc'
    real(rp)       :: record_point(3) = 0.0_rp, record_normal(3) = [1.0_rp, 0.0_rp, 0.0_rp]
    real(rp)       :: record_dt = 0.0_rp, record_t_start = 0.0_rp, record_segment_time = 0.0_rp
    integer(ip)    :: record_n_segments = 1
    integer(ip)    :: ios

    namelist /dopamine_inflow/ inflow_source, inflow_file, inflow_t_offset, inflow_crossfade_time, &
                               inflow_seed, inflow_export_grid
    namelist /dopamine_inflow_gen/ gen_file, gen_grid_file, gen_W, gen_y0, gen_ny, gen_nz, gen_z_min, gen_z_max, &
                               gen_dt, gen_n_turnovers, gen_T_turnover, gen_n_segments, gen_n_sigma, gen_seed, &
                                   gen_verify, gen_decorrelate
    namelist /dopamine_inflow_record/ record_enabled, record_file, record_point, record_normal, record_dt, &
                               record_t_start, record_segment_time, record_n_segments

    read(iunit, nml=dopamine_inflow, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_inflow, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if
    rewind(iunit)
    read(iunit, nml=dopamine_inflow_gen, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_inflow_gen, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if
    rewind(iunit)
    read(iunit, nml=dopamine_inflow_record, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_inflow_record, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    if (trim(inflow_source) /= 'sem' .and. trim(inflow_source) /= 'planes') &
      error stop "inflow_source must be 'sem' or 'planes'"
    if (trim(inflow_source) == 'planes' .and. len_trim(inflow_file) == 0) &
      error stop "inflow_source = 'planes' needs inflow_file"
    if (inflow_crossfade_time < 0.0_rp) error stop 'inflow_crossfade_time must be >= 0'
    if (record_enabled) then
      if (record_dt <= 0.0_rp .or. record_segment_time <= 0.0_rp .or. record_n_segments < 1) &
        error stop 'inflow_record: need record_dt, segment_time > 0, n_segments >= 1'
      if (sqrt(sum(record_normal**2)) <= 0.0_rp) error stop 'record_normal must be non-zero'
    end if

    config%inflow%source = inflow_source
    config%inflow%file = inflow_file
    config%inflow%t_offset = inflow_t_offset
    config%inflow%crossfade_time = inflow_crossfade_time
    config%inflow%seed = inflow_seed
    config%inflow%export_grid = inflow_export_grid
    config%inflow_gen%file = gen_file
    config%inflow_gen%grid_file = gen_grid_file
    config%inflow_gen%W = gen_W
    config%inflow_gen%y0 = gen_y0
    config%inflow_gen%ny = gen_ny
    config%inflow_gen%nz = gen_nz
    config%inflow_gen%z_min = gen_z_min
    config%inflow_gen%z_max = gen_z_max
    config%inflow_gen%dt = gen_dt
    config%inflow_gen%n_turnovers = gen_n_turnovers
    config%inflow_gen%T_turnover = gen_T_turnover
    config%inflow_gen%n_segments = gen_n_segments
    config%inflow_gen%n_sigma = gen_n_sigma
    config%inflow_gen%seed = gen_seed
    config%inflow_gen%verify = gen_verify
    config%inflow_gen%decorrelate = gen_decorrelate
    config%inflow_rec%enabled = record_enabled
    config%inflow_rec%file = record_file
    config%inflow_rec%point = record_point
    config%inflow_rec%normal = record_normal / sqrt(sum(record_normal**2))
    config%inflow_rec%dt = record_dt
    config%inflow_rec%t_start = record_t_start
    config%inflow_rec%segment_time = record_segment_time
    config%inflow_rec%n_segments = record_n_segments
  end subroutine read_inflow

  subroutine read_scalar_sources(iunit, config)
    integer(ip), intent(in) :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_scalar_sources, source_scalar_id(MAX_SCALAR_SOURCES), ios, i
    character(64) :: source_name(MAX_SCALAR_SOURCES), source_kind(MAX_SCALAR_SOURCES)
    character(32) :: source_shape(MAX_SCALAR_SOURCES), source_time_mode(MAX_SCALAR_SOURCES)
    character(256) :: source_file(MAX_SCALAR_SOURCES)
    real(rp) :: source_magnitude(MAX_SCALAR_SOURCES), source_omega(MAX_SCALAR_SOURCES), source_phase(MAX_SCALAR_SOURCES)
    real(rp) :: source_center(MAX_SCALAR_SOURCES,3), source_p0(MAX_SCALAR_SOURCES,3), source_p1(MAX_SCALAR_SOURCES,3)
    real(rp) :: source_normal(MAX_SCALAR_SOURCES,3), source_half_thickness(MAX_SCALAR_SOURCES)
    real(rp) :: source_radius(MAX_SCALAR_SOURCES), source_box_min(MAX_SCALAR_SOURCES,3), source_box_max(MAX_SCALAR_SOURCES,3)
    real(rp) :: source_tag_radius(MAX_SCALAR_SOURCES)

    namelist /dopamine_scalar_sources/ n_scalar_sources, source_scalar_id, source_name, source_kind, source_shape, &
      source_time_mode, source_file, source_magnitude, source_omega, source_phase, source_center, source_p0, source_p1, &
      source_normal, source_half_thickness, source_radius, source_box_min, source_box_max, source_tag_radius

    n_scalar_sources = 0
    source_scalar_id = 1
    source_name = ''
    source_kind = 'volume'
    source_shape = 'sphere'
    source_time_mode = 'constant'
    source_file = ''
    source_magnitude = 0.0_rp
    source_omega = 0.0_rp
    source_phase = 0.0_rp
    source_center = 0.0_rp
    source_p0 = 0.0_rp
    source_p1 = 0.0_rp
    source_normal = 0.0_rp
    source_half_thickness = 0.0_rp
    source_radius = 0.0_rp
    source_box_min = 0.0_rp
    source_box_max = 0.0_rp
    source_tag_radius = 0.0_rp

    read(iunit, nml=dopamine_scalar_sources, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_scalar_sources, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_scalar_sources = min(max(0_ip, n_scalar_sources), MAX_SCALAR_SOURCES)
    if (config%n_scalar_sources <= 0) then
      config%params%n_scalar_sources = 0
      return
    end if

    config%params%n_scalar_sources = config%n_scalar_sources
    allocate(config%params%scalar_sources(config%n_scalar_sources))
    do i = 1, config%n_scalar_sources
      config%scalar_sources(i)%enabled = .true.
      config%scalar_sources(i)%scalar_id = max(1_ip, min(config%n_scalars, source_scalar_id(i)))
      config%scalar_sources(i)%name = trim(source_name(i))
      config%scalar_sources(i)%src_kind = scalar_source_kind_from_name(source_kind(i))
      config%scalar_sources(i)%shape = scalar_source_shape_from_name(source_shape(i))
      config%scalar_sources(i)%time_mode = scalar_time_mode_from_name(source_time_mode(i))
      config%scalar_sources(i)%file = trim(source_file(i))
      config%scalar_sources(i)%magnitude = source_magnitude(i)
      config%scalar_sources(i)%omega = source_omega(i)
      config%scalar_sources(i)%phase = source_phase(i)
      config%scalar_sources(i)%center = source_center(i,:)
      config%scalar_sources(i)%p0 = source_p0(i,:)
      config%scalar_sources(i)%p1 = source_p1(i,:)
      config%scalar_sources(i)%normal = source_normal(i,:)
      config%scalar_sources(i)%half_thickness = max(0.0_rp, source_half_thickness(i))
      config%scalar_sources(i)%radius = max(0.0_rp, source_radius(i))
      config%scalar_sources(i)%box_min = source_box_min(i,:)
      config%scalar_sources(i)%box_max = source_box_max(i,:)
      config%scalar_sources(i)%tag_radius = max(0.0_rp, source_tag_radius(i))
      ! A sphere source with radius 0 selects no cell at all; say so instead of silently doing nothing.
      if (config%scalar_sources(i)%shape == SCALAR_SHAPE_SPHERE .and. config%scalar_sources(i)%radius <= 0.0_rp &
          .and. (config%scalar_sources(i)%src_kind == SCALAR_SRC_POINT .or. &
                 config%scalar_sources(i)%src_kind == SCALAR_SRC_LINE .or. &
                 config%scalar_sources(i)%src_kind == SCALAR_SRC_VOLUME)) &
        write(*, '(A,I0,A)') ' WARNING: scalar source ', i, ' has source_radius = 0 and will match no cells'
      config%params%scalar_sources(i) = config%scalar_sources(i)
    end do
  end subroutine read_scalar_sources

  subroutine map_scalar_bc_modes(mode_out, mode_names, n_scalars)
    integer(ip), intent(out) :: mode_out(:)
    character(16), intent(in) :: mode_names(:)
    integer(ip), intent(in) :: n_scalars
    integer(ip) :: s

    do s = 1, n_scalars
      mode_out(s) = scalar_bc_mode_from_name(mode_names(s))
    end do
  end subroutine map_scalar_bc_modes

  subroutine read_initial_condition(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    character(32) :: ic_type = 'uniform'
    real(rp)      :: ic_uref = 1.0_rp
    real(rp)      :: ic_lref = 0.0_rp
    real(rp)      :: noise_level = 0.0_rp
    integer(ip)   :: noise_seed = 0
    integer(ip)   :: wall_normal_direction = 3
    integer(ip)   :: ios

    namelist /dopamine_initial_condition/ ic_type, ic_uref, ic_lref, noise_level, noise_seed, wall_normal_direction

    read(iunit, nml=dopamine_initial_condition, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_initial_condition, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%ic_type = trim(ic_type)
    config%ic_uref = ic_uref
    config%ic_lref = ic_lref
    config%ic_noise_level = noise_level
    config%ic_noise_seed = noise_seed
    config%wall_normal_direction = max(1, min(3, wall_normal_direction))
  end subroutine read_initial_condition

  subroutine read_output(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    character(256) :: prefix         = 'dopamine'
    integer(ip)    :: write_interval = 100
    real(rp)       :: write_time     = 0.0_rp
    logical        :: hdf5_use_file_locking = .false.
    character(256) :: warn_log_file = 'solver_warnings.txt'
    integer(ip)    :: warn_flush_interval = 100
    integer(ip)    :: ios

    namelist /dopamine_output/ prefix, write_interval, write_time, hdf5_use_file_locking, &
                               warn_log_file, warn_flush_interval

    read(iunit, nml=dopamine_output, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_output, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%output_prefix             = trim(prefix)
    config%params%write_interval     = write_interval
    config%params%write_time         = write_time
    config%hdf5_use_file_locking     = hdf5_use_file_locking
    config%warn_log_file             = trim(warn_log_file)
    config%warn_flush_interval       = max(1_ip, warn_flush_interval)
  end subroutine read_output

  subroutine read_checkpoint(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    logical        :: restart_enabled      = .false.
    character(256) :: restart_dir          = ''
    integer(ip)    :: restart_step         = -1
    integer(ip)    :: checkpoint_interval  = 0
    integer(ip)    :: ios

    namelist /dopamine_checkpoint/ restart_enabled, restart_dir, restart_step, checkpoint_interval

    read(iunit, nml=dopamine_checkpoint, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_checkpoint, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%restart_enabled      = restart_enabled
    config%restart_dir          = trim(restart_dir)
    config%restart_step         = restart_step
    config%checkpoint_interval  = checkpoint_interval
  end subroutine read_checkpoint

  subroutine read_sampling_planes(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_planes, plane_nx(MAX_SAMPLE_OBJECTS), plane_ny(MAX_SAMPLE_OBJECTS)
    integer(ip) :: plane_write_frequency(MAX_SAMPLE_OBJECTS), p, ios
    character(64) :: plane_name(MAX_SAMPLE_OBJECTS)
    real(rp) :: plane_origin(MAX_SAMPLE_OBJECTS,3), plane_normal(MAX_SAMPLE_OBJECTS,3)
    real(rp) :: plane_tangent_hint(MAX_SAMPLE_OBJECTS,3)
    real(rp) :: plane_half_width(MAX_SAMPLE_OBJECTS), plane_half_height(MAX_SAMPLE_OBJECTS)

    namelist /dopamine_sampling_planes/ n_planes, plane_name, plane_origin, plane_normal, &
      plane_tangent_hint, plane_half_width, plane_half_height, plane_nx, plane_ny, plane_write_frequency

    n_planes = 0
    plane_name = ''
    plane_origin = 0.0_rp
    plane_normal = 0.0_rp
    plane_tangent_hint = 0.0_rp
    plane_half_width = 0.5_rp
    plane_half_height = 0.5_rp
    plane_nx = 64
    plane_ny = 64
    plane_write_frequency = config%stats_slice_frequency

    read(iunit, nml=dopamine_sampling_planes, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_sampling_planes, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_planes = min(max(0_ip, n_planes), MAX_SAMPLE_OBJECTS)
    do p = 1, config%n_planes
      config%planes(p)%enabled = .true.
      config%planes(p)%name = trim(plane_name(p))
      config%planes(p)%origin = plane_origin(p, :)
      config%planes(p)%normal = plane_normal(p, :)
      config%planes(p)%tangent_hint = plane_tangent_hint(p, :)
      config%planes(p)%half_width = max(1.0e-8_rp, plane_half_width(p))
      config%planes(p)%half_height = max(1.0e-8_rp, plane_half_height(p))
      config%planes(p)%nx = max(2_ip, plane_nx(p))
      config%planes(p)%ny = max(2_ip, plane_ny(p))
      config%planes(p)%write_frequency = max(0_ip, plane_write_frequency(p))
      if (config%planes(p)%write_frequency == 0) config%planes(p)%write_frequency = config%stats_slice_frequency
      if (len_trim(config%planes(p)%name) == 0) write(config%planes(p)%name, '(A,I0)') 'plane_', p
    end do
  end subroutine read_sampling_planes

  subroutine read_sampling_lines(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_lines, line_npts(MAX_SAMPLE_OBJECTS), line_write_frequency(MAX_SAMPLE_OBJECTS)
    integer(ip) :: p, ios
    character(64) :: line_name(MAX_SAMPLE_OBJECTS)
    real(rp) :: line_p0(MAX_SAMPLE_OBJECTS,3), line_p1(MAX_SAMPLE_OBJECTS,3)

    namelist /dopamine_sampling_lines/ n_lines, line_name, line_p0, line_p1, line_npts, line_write_frequency

    n_lines = 0
    line_name = ''
    line_p0 = 0.0_rp
    line_p1 = 0.0_rp
    line_npts = 128
    line_write_frequency = config%stats_line_frequency

    read(iunit, nml=dopamine_sampling_lines, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_sampling_lines, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_lines = min(max(0_ip, n_lines), MAX_SAMPLE_OBJECTS)
    do p = 1, config%n_lines
      config%lines(p)%enabled = .true.
      config%lines(p)%name = trim(line_name(p))
      config%lines(p)%p0 = line_p0(p, :)
      config%lines(p)%p1 = line_p1(p, :)
      config%lines(p)%npts = max(2_ip, line_npts(p))
      config%lines(p)%write_frequency = max(0_ip, line_write_frequency(p))
      if (config%lines(p)%write_frequency == 0) then
        config%lines(p)%write_frequency = max(config%stats_line_frequency, config%stats_probe_frequency)
      end if
      if (len_trim(config%lines(p)%name) == 0) write(config%lines(p)%name, '(A,I0)') 'line_', p
    end do
  end subroutine read_sampling_lines

  subroutine read_sampling_curves(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_curves, curve_write_frequency(MAX_SAMPLE_OBJECTS), p, ios
    character(64) :: curve_name(MAX_SAMPLE_OBJECTS)
    character(256) :: curve_file(MAX_SAMPLE_OBJECTS)

    namelist /dopamine_sampling_curves/ n_curves, curve_name, curve_file, curve_write_frequency

    n_curves = 0
    curve_name = ''
    curve_file = ''
    curve_write_frequency = config%stats_curve_frequency

    read(iunit, nml=dopamine_sampling_curves, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_sampling_curves, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_curves = min(max(0_ip, n_curves), MAX_SAMPLE_OBJECTS)
    do p = 1, config%n_curves
      config%curves(p)%enabled = .true.
      config%curves(p)%name = trim(curve_name(p))
      config%curves(p)%file = trim(curve_file(p))
      config%curves(p)%write_frequency = max(0_ip, curve_write_frequency(p))
      if (config%curves(p)%write_frequency == 0) then
        config%curves(p)%write_frequency = max(config%stats_curve_frequency, config%stats_probe_frequency)
      end if
      if (len_trim(config%curves(p)%name) == 0) write(config%curves(p)%name, '(A,I0)') 'curve_', p
    end do
  end subroutine read_sampling_curves

  subroutine read_sampling_surfaces(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: n_surfaces, surface_write_frequency(MAX_SAMPLE_OBJECTS), p, ios
    character(64) :: surface_name(MAX_SAMPLE_OBJECTS)
    character(256) :: surface_file(MAX_SAMPLE_OBJECTS)
    logical :: surface_intersect_with_planes(MAX_SAMPLE_OBJECTS)
    logical :: surface_sort_intersection_points(MAX_SAMPLE_OBJECTS)
    real(rp) :: surface_intersection_decimate_spacing(MAX_SAMPLE_OBJECTS)

    namelist /dopamine_sampling_surfaces/ n_surfaces, surface_name, surface_file, surface_write_frequency, &
      surface_intersect_with_planes, surface_intersection_decimate_spacing, surface_sort_intersection_points

    n_surfaces = 0
    surface_name = ''
    surface_file = ''
    surface_write_frequency = config%stats_probe_frequency
    surface_intersect_with_planes = .true.
    surface_sort_intersection_points = config%stats_sort_intersection_points
    surface_intersection_decimate_spacing = config%stats_intersection_decimate_spacing

    read(iunit, nml=dopamine_sampling_surfaces, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_sampling_surfaces, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if

    config%n_surfaces = min(max(0_ip, n_surfaces), MAX_SAMPLE_OBJECTS)
    do p = 1, config%n_surfaces
      config%surfaces(p)%enabled = .true.
      config%surfaces(p)%name = trim(surface_name(p))
      config%surfaces(p)%file = trim(surface_file(p))
      config%surfaces(p)%write_frequency = max(0_ip, surface_write_frequency(p))
      if (config%surfaces(p)%write_frequency == 0) config%surfaces(p)%write_frequency = config%stats_probe_frequency
      config%surfaces(p)%intersect_with_planes = surface_intersect_with_planes(p)
      config%surfaces(p)%intersection_decimate_spacing = max(0.0_rp, surface_intersection_decimate_spacing(p))
      config%surfaces(p)%sort_intersection_points = surface_sort_intersection_points(p)
      if (len_trim(config%surfaces(p)%name) == 0) write(config%surfaces(p)%name, '(A,I0)') 'surface_', p
    end do
  end subroutine read_sampling_surfaces

  subroutine read_PetScParams(iunit, config)
    integer(ip),       intent(in)    :: iunit
    type(InputConfig), intent(inout) :: config

    integer(ip) :: ios, i
    character(256)        :: petsc_options_tmp(max_opts)
    integer(ip)           :: n_petsc_options

    namelist /dopamine_PetscParams/ petsc_options_tmp, n_petsc_options

    petsc_options_tmp = ''
    n_petsc_options = 0

    read(iunit, nml=dopamine_PetscParams, iostat=ios)
    if (ios /= 0 .and. ios /= -1) then
      write(*,'(A,I0)') 'ERROR: cannot parse &dopamine_PetscParams, ios=', ios
      error stop 'input namelist parse error (check variable names and syntax)'
    end if
    if (n_petsc_options > max_opts) then
      write(*, '(A,I0,A,I0)') 'ERROR: n_petsc_options = ', n_petsc_options, ' exceeds the supported maximum ', max_opts
      error stop 'input: too many PETSc options'
    end if
    if (n_petsc_options > 0) then
      if (allocated(config%petsc_options)) deallocate(config%petsc_options)
      allocate(character(len=256) :: config%petsc_options(n_petsc_options))
      do i = 1, n_petsc_options
        config%petsc_options(i) = trim(petsc_options_tmp(i))
      end do
      ! Propagate PETSc options to params for solver access
      if (allocated(config%params%petsc_options)) deallocate(config%params%petsc_options)
      allocate(character(len=256) :: config%params%petsc_options(n_petsc_options))
      config%params%petsc_options = config%petsc_options
    end if

    close(iunit)

  end subroutine read_PetScParams

  integer(ip) function bc_name_to_id(name)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(lower_case(name))))
    case ('velocity', 'inlet', 'fixedvalue')
      bc_name_to_id = BC_INLET
    case ('outlet', 'fixedpressure', 'outflow')
      bc_name_to_id = BC_OUTLET
    case ('zerogradient')
      bc_name_to_id = BC_ZEROGRADIENT
    case ('noslip')
      bc_name_to_id = BC_WALL_NOSLIP
    case ('movingwall', 'wallvelocity')
      bc_name_to_id = BC_WALL_VELOCITY
    case ('wallmodel', 'wall')
      bc_name_to_id = BC_WALL_MODEL
    case ('slip')
      bc_name_to_id = BC_SLIP
    case ('symmetry')
      bc_name_to_id = BC_SYMMETRY
    case ('periodic', 'cyclic')
      bc_name_to_id = BC_PERIODIC
    case ('fixedtemperature')
      bc_name_to_id = BC_FIXEDTEMP
    case ('ablinlet')
      bc_name_to_id = BC_ABL_INLET
    case default
      ! A typo (e.g. a wrong-case 'noslip') used to become an open boundary without a word.
      write(*, '(A,A,A)') 'ERROR: unrecognised patch_type ''', trim(name), ''''
      error stop 'bc_name_to_id: unknown patch_type'
    end select
  end function bc_name_to_id


  !> Resolve the wall-model selector; use_abl is the legacy way of asking for MOST.
  integer(ip) function resolve_wall_model(name, use_abl) result(id)
    character(len=*), intent(in) :: name
    logical,          intent(in) :: use_abl

    select case (trim(adjustl(lower_case(name))))
    case ('most', 'abl')
      id = WM_MOST
    case ('loglaw', 'log', 'smooth')
      id = WM_LOGLAW
    case ('none', 'off')
      id = WM_NONE
    case ('')
      id = merge(WM_MOST, WM_LOGLAW, use_abl)
    case default
      write(*,'(A)') 'WARNING: unknown wall_model "'//trim(name)//'", using loglaw'
      id = WM_LOGLAW
    end select
  end function resolve_wall_model

  subroutine input_fatal(name, value, allowed)
    character(len=*), intent(in) :: name, value, allowed
    write(*, '(A)') 'ERROR: unknown '//name//' = '''//trim(value)//''' (allowed: '//allowed//')'
    error stop 'input: unknown option value'
  end subroutine input_fatal

  pure function upper_case(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
      c = iachar(t(i:i))
      if (c >= 97 .and. c <= 122) t(i:i) = achar(c - 32)
    end do
  end function upper_case

  pure function lower_case(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
      c = iachar(t(i:i))
      if (c >= 65 .and. c <= 90) t(i:i) = achar(c + 32)
    end do
  end function lower_case

end module mod_input_parser
