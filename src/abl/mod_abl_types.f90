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
    real(rp)    :: sem_sigma      = 0.0_rp  ! eddy length scale [m]; 0 = compute from BL height
    character(len=256) :: sem_profile_file = ''  ! CSV of U(z) and Reynolds stresses; '' = analytic
    logical     :: use_nudging    = .false.  ! enable Davies nudging
    character(len=256) :: nudging_file = ''  ! netCDF/CSV file with reference profiles
    real(rp)    :: nudging_tau    = 1800.0_rp  ! nudging timescale [s]
    real(rp)    :: nudging_z_start = 0.0_rp  ! height above which nudging is active [m]
    real(rp)    :: nudging_z_full  = 0.0_rp  ! height of full nudging strength; <= z_start = no ramp [m]
    logical     :: nudging_horizontal_only = .true.  ! leave w unforced so resolved motion survives
  end type ABLConfig

  !> Synthetic Eddy Method state - persists between timesteps
  type, public :: SEMState
    logical     :: initialized    = .false.
    integer(ip) :: n_eddies       = 0
    real(rp), allocatable :: pos(:,:)  ! (3, n_eddies) - eddy positions [m]
    real(rp), allocatable :: amp(:,:)  ! (3, n_eddies) - eddy intensity signs (+1/-1)
    real(rp)    :: sigma           = 0.0_rp  ! eddy length scale [m]
    real(rp)    :: u_convect       = 0.0_rp  ! bulk advection speed [m/s]
    ! Inlet-plane frame: origin, inward normal (convection direction), in-plane tangents.
    real(rp)    :: x_inlet(3)      = 0.0_rp
    real(rp)    :: n_in(3)         = [1.0_rp, 0.0_rp, 0.0_rp]
    real(rp)    :: t1(3)           = [0.0_rp, 1.0_rp, 0.0_rp]
    real(rp)    :: t2(3)           = [0.0_rp, 0.0_rp, 1.0_rp]
    real(rp)    :: t1_min          = 0.0_rp  ! eddy box bounds along t1 [m]
    real(rp)    :: t1_max          = 0.0_rp
    real(rp)    :: t2_min          = 0.0_rp  ! eddy box bounds along t2 [m]
    real(rp)    :: t2_max          = 0.0_rp
    real(rp)    :: box_volume      = 0.0_rp  ! recycle box volume [m^3]
    ! Periodic in-plane directions (periodic patches parallel to the inlet normal): the eddy box spans
    ! exactly one period there and eddies are wrapped, so the fluctuation is periodic across the seam.
    logical     :: per1            = .false.
    logical     :: per2            = .false.
    real(rp)    :: period1         = 0.0_rp  ! period along t1 [m]
    real(rp)    :: period2         = 0.0_rp  ! period along t2 [m]
    real(rp)    :: per1_lo         = 0.0_rp  ! t1 of the lower periodic plane, relative to x_inlet [m]
    real(rp)    :: per2_lo         = 0.0_rp
    ! Prescribed inflow profiles; empty means the analytic MOST/Panofsky defaults.
    logical     :: has_profile     = .false.
    logical     :: prof_global_axes = .false.  ! 10-column files give U and R in x,y,z
    integer(ip) :: n_levels        = 0
    real(rp), allocatable :: prof_z(:)  ! (n_levels) ascending heights [m]
    real(rp), allocatable :: prof_U(:,:)  ! (3, n_levels) mean velocity [m/s]
    real(rp), allocatable :: prof_R(:,:)  ! (6, n_levels) uu,vv,ww,uv,uw,vw [m2/s2]
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
