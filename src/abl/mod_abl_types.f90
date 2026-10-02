! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Atmospheric Boundary Layer (ABL) types and configuration for urban flow simulations
module mod_abl_types
  use mod_kinds
  implicit none
  private

  ! Wall-model selector; resolved once at init from the input strings.
  integer(ip), parameter, public :: WM_NONE   = 0_ip
  integer(ip), parameter, public :: WM_LOGLAW = 1_ip  ! smooth-wall log law
  integer(ip), parameter, public :: WM_MOST   = 2_ip  ! Monin-Obukhov, rough/stratified

  !> ABL configuration - input parameters and runtime state
  type, public :: ABLConfig
    ! The three ABL features are independent; use_abl is a legacy alias for wall_model='most'.
    logical     :: use_abl        = .false.
    character(len=16) :: wall_model = ''  ! '', 'none', 'loglaw', 'most'
    integer(ip) :: wall_model_id  = WM_LOGLAW  ! resolved from wall_model / use_abl
    real(rp)    :: u_ref          = 8.0_rp  ! reference speed at z_ref [m/s]
    real(rp)    :: z_ref          = 100.0_rp  ! reference height [m]
    real(rp)    :: z0_m           = 0.1_rp  ! momentum roughness length [m]
    real(rp)    :: z0_h           = 0.01_rp  ! thermal roughness length [m]
    real(rp)    :: T_surface      = 300.0_rp  ! surface temperature [K]
    real(rp)    :: T_ref          = 300.0_rp  ! free-stream reference temperature [K]
    real(rp)    :: L_obukhov      = 0.0_rp  ! Obukhov length [m]: 0=neutral, >0=stable, <0=unstable
    logical     :: compute_L      = .true.  ! update L dynamically from surface fluxes
    real(rp)    :: wind_dir_deg   = 270.0_rp  ! wind direction (meteorological, degrees)
    logical     :: use_sem        = .true.  ! enable Synthetic Eddy Method
    integer(ip) :: n_eddies       = 500  ! number of eddies
    real(rp)    :: sem_sigma      = 0.0_rp  ! streamwise integral length scale Lx at z_ref [m]; 0 = 0.05*inlet span
    character(len=256) :: sem_profile_file = ''  ! CSV of U(z) plus stresses or intensities; '' = analytic
    character(len=8)   :: sem_profile_format = 'stress'  ! 'stress' (z,U,R.. ) or 'ti' (z,U,Iu,Iv,Iw..)
    character(len=8)   :: sem_mean_law   = 'log'  ! analytic mean profile: 'log' or 'power'
    real(rp)    :: sem_alpha      = 0.2_rp  ! power-law exponent when sem_mean_law = 'power'
    real(rp)    :: sem_Iu_ref     = 0.0_rp  ! streamwise turbulence intensity at z_ref; 0 = Panofsky-Dutton stresses
    real(rp)    :: sem_Iv_ratio   = 0.75_rp  ! Iv / Iu
    real(rp)    :: sem_Iw_ratio   = 0.5_rp   ! Iw / Iu
    real(rp)    :: sem_rho_uw     = -0.3_rp  ! u'w' correlation coefficient where the input has none
    real(rp)    :: sem_L_exponent = 0.0_rp  ! Lx(z) = sem_sigma*(z/z_ref)**exponent; 0 = height-independent
    real(rp)    :: sem_Ly_ratio   = 0.5_rp   ! lateral length scale / Lx
    real(rp)    :: sem_Lz_ratio   = 0.35_rp  ! vertical length scale / Lx
    real(rp)    :: sem_L_min      = 0.0_rp  ! floor on every length scale [m] (e.g. a few grid cells)
    integer(ip) :: sem_seed       = 0       ! 0 = built-in seed
    logical     :: use_nudging    = .false.  ! enable Davies nudging
    character(len=256) :: nudging_file = ''  ! netCDF/CSV file with reference profiles
    real(rp)    :: nudging_tau    = 1800.0_rp  ! nudging timescale [s]
    real(rp)    :: nudging_z_start = 0.0_rp  ! height above which nudging is active [m]
    real(rp)    :: nudging_z_full  = 0.0_rp  ! height of full nudging strength; <= z_start = no ramp [m]
    logical     :: nudging_horizontal_only = .true.  ! leave w unforced so resolved motion survives
  end type ABLConfig

  !> dopamine-inflow-gen: offline generation of an inflow-plane database from the SEM.
  type, public :: InflowGenParams
    character(len=256) :: file = 'inflow_planes.nc'
    character(len=256) :: grid_file = ''     ! plane grid exported by the solver (inflow_export_grid)
    real(rp)    :: W = 0.0_rp                ! spanwise period [m]; required unless grid_file is given
    real(rp)    :: y0 = 0.0_rp               ! t1 coordinate of the first node [m]
    integer(ip) :: ny = 0, nz = 0            ! plane grid; z uniform in [z_min, z_max]
    real(rp)    :: z_min = 0.0_rp, z_max = 0.0_rp
    real(rp)    :: dt = 0.0_rp               ! frame spacing [s]; 0 = from the smallest eddy time
    real(rp)    :: n_turnovers = 2.0_rp      ! segment length in eddy turnover times
    real(rp)    :: T_turnover = 0.0_rp       ! eddy turnover time [s]; 0 = Lx/sigma_u at z_ref
    integer(ip) :: n_segments = 4            ! independent realisations
    real(rp)    :: n_sigma = 3.0_rp          ! eddies per kernel volume (sets the eddy count)
    integer(ip) :: seed = 1
    logical     :: decorrelate = .true.      ! remove the segment-mean spanwise-mean waveform (no line at 1/P)
    logical     :: verify = .true.
  end type InflowGenParams

  !> Successor-run inflow source.
  type, public :: InflowConfig
    character(len=16)  :: source = 'sem'     ! 'sem' (online eddies) or 'planes' (database)
    character(len=256) :: file = ''          ! database for source = 'planes'
    real(rp)    :: t_offset = 0.0_rp         ! start offset into the replay sequence [s]
    real(rp)    :: crossfade_time = 0.0_rp   ! overlap between replayed segments [s]; 0 = 10% of a segment
    integer(ip) :: seed = 1                  ! seeds the segment order and spanwise shifts
    character(len=256) :: export_grid = ''   ! write the inlet plane grid for dopamine-inflow-gen and stop
  end type InflowConfig

  !> Precursor run: record fluctuation planes into a database.
  type, public :: InflowRecordConfig
    logical     :: enabled = .false.
    character(len=256) :: file = 'inflow_record.nc'
    real(rp)    :: point(3)  = 0.0_rp        ! a point on the recording plane
    real(rp)    :: normal(3) = [1.0_rp, 0.0_rp, 0.0_rp]  ! flow direction through the plane
    real(rp)    :: dt = 0.0_rp               ! frame spacing [s]
    real(rp)    :: t_start = 0.0_rp          ! discard spin-up before this time [s]
    real(rp)    :: segment_time = 0.0_rp     ! length of one segment [s]
    integer(ip) :: n_segments = 1
  end type InflowRecordConfig

  !> Synthetic Eddy Method state - persists between timesteps.
  !> Each eddy carries a vector-potential strength and its own kernel semi-axes, so its
  !> velocity (the curl of the potential) is exactly solenoidal for any scales.
  type, public :: SEMState
    logical     :: initialized    = .false.
    integer(ip) :: n_eddies       = 0
    real(rp), allocatable :: pos(:,:)    ! (3, n_eddies) eddy centres [m]
    real(rp), allocatable :: amp(:,:)    ! (3, n_eddies) potential strength in the profile frame
    real(rp), allocatable :: sig(:,:)    ! (3, n_eddies) kernel semi-axes in the profile frame [m]
    real(rp), allocatable :: uconv(:)    ! (n_eddies) convection speed along n_in [m/s]
    real(rp)    :: n_sigma         = 0.0_rp  ! eddies per kernel volume sx*sy*sz (sets the stress normalisation)
    real(rp)    :: L_ref           = 0.0_rp  ! Lx at z_ref actually used [m]
    real(rp)    :: z_floor         = 0.0_rp  ! heights below this use the profile value at z_floor [m]
    real(rp)    :: u_conv_min      = 0.0_rp  ! floor on the eddy convection speed [m/s]
    real(rp)    :: realiz_err      = 0.0_rp  ! worst relative stress error from the realisability projection
    ! Inlet-plane frame: origin, inward normal (convection direction), in-plane tangents.
    real(rp)    :: x_inlet(3)      = 0.0_rp
    real(rp)    :: n_in(3)         = [1.0_rp, 0.0_rp, 0.0_rp]
    real(rp)    :: t1(3)           = [0.0_rp, 1.0_rp, 0.0_rp]
    real(rp)    :: t2(3)           = [0.0_rp, 0.0_rp, 1.0_rp]
    real(rp)    :: t1_min          = 0.0_rp  ! eddy box bounds along t1 [m]
    real(rp)    :: t1_max          = 0.0_rp
    real(rp)    :: t2_min          = 0.0_rp  ! eddy box bounds along t2 [m]
    real(rp)    :: t2_max          = 0.0_rp
    real(rp)    :: s_half          = 0.0_rp  ! half depth of the recycle box along n_in [m]
    ! Spanwise (t1) periodicity: eddies wrap with minimum-image distances over t1_period.
    logical     :: periodic_t1     = .false.
    real(rp)    :: t1_period       = 0.0_rp
    ! Largest kernel half-extent along t1, t2, n_in over all heights (sets the bin size) [m].
    real(rp)    :: ext(3)          = 0.0_rp
    ! Spatial bins over (t1, t2, n_in); eddies are re-binned after every advance.
    logical     :: bins_ready      = .false.
    integer(ip) :: nbin(3)         = 1
    real(rp)    :: bin_lo(3)       = 0.0_rp
    real(rp)    :: bin_h(3)        = 1.0_rp
    integer(ip), allocatable :: bin_start(:)  ! (nbin_total+1) CSR offsets into bin_list
    integer(ip), allocatable :: bin_list(:)   ! eddy ids grouped by bin
    ! Profile frame: columns are the streamwise, lateral and vertical axes of the stresses and scales.
    real(rp)    :: frame(3,3)      = reshape([1.0_rp,0.0_rp,0.0_rp, 0.0_rp,1.0_rp,0.0_rp, 0.0_rp,0.0_rp,1.0_rp], [3,3])
    ! Eddy heights are drawn from the 1/(sx*sy*sz) density tabulated along t2.
    integer(ip) :: n_cdf           = 0
    real(rp)    :: cdf_b0          = 0.0_rp
    real(rp)    :: cdf_db          = 0.0_rp
    real(rp), allocatable :: cdf_seed(:)  ! seeding: uniform in s, density ~ 1/kernel volume
    real(rp), allocatable :: cdf_rec(:)   ! recycling: the same, weighted by convection speed
    ! Prescribed inflow profile; empty means the analytic law from the namelist.
    logical     :: has_profile     = .false.
    logical     :: override_mean   = .false.  ! the solver must take the inlet mean from the SEM
    logical     :: prof_global_axes = .false.  ! tables give U and R in x,y,z (frame = identity)
    integer(ip) :: n_levels        = 0
    real(rp), allocatable :: prof_z(:)    ! (n_levels) ascending heights [m]
    real(rp), allocatable :: prof_U(:,:)  ! (3, n_levels) mean velocity [m/s]
    real(rp), allocatable :: prof_R(:,:)  ! (6, n_levels) uu,vv,ww,uv,uw,vw [m2/s2]
    real(rp), allocatable :: prof_L(:,:)  ! (3, n_levels) integral scales Lx,Ly,Lz [m]; 0 = from the law
  end type SEMState

  !> Weather data for nudging - loaded from file at initialization
  type, public :: NudgingData
    logical     :: loaded = .false.
    integer(ip) :: n_times = 0
    integer(ip) :: n_levels = 0
    real(rp), allocatable :: times(:)  ! (n_times) [s]
    real(rp), allocatable :: z_levels(:)  ! (n_levels) [m]
    real(rp), allocatable :: U_profile(:,:,:)  ! (3, n_levels, n_times) [m/s]
    real(rp), allocatable :: T_profile(:,:)  ! (n_levels, n_times) [K]
  end type NudgingData

end module mod_abl_types
