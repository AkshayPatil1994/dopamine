! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Physical and numerical constants
module mod_parameters
  use mod_kinds
  implicit none
  private

  ! Mathematical constants
  real(rp), parameter, public :: PI    = 3.14159265358979323846_rp
  real(rp), parameter, public :: TWOPI = 6.28318530717958647692_rp

  ! Von Karman constant
  real(rp), parameter, public :: KAPPA = 0.41_rp

  ! Log-law additive constant
  real(rp), parameter, public :: BLOG  = 5.2_rp

  ! Vreman constant c = 2.5*Cs^2; 0.025 is the wall-bounded value (Cs = 0.1).
  real(rp), parameter, public :: C_VREMAN = 0.025_rp

  ! AMD (Rozema et al. 2015) constant: the Poincare constant of second-order central differences.
  real(rp), parameter, public :: AMD_C_DEFAULT = 1.0_rp / 3.0_rp

  ! Standard gravity
  real(rp), parameter, public :: G_STD = 9.81_rp

  ! Numerical tolerances
  real(rp), parameter, public :: SMALL  = 1.0e-30_rp
  real(rp), parameter, public :: VSMALL = 1.0e-200_rp
  real(rp), parameter, public :: GREAT  = 1.0e15_rp

  ! Builds mesh face coloring so mod_advection/mod_viscous OpenMP loops run race-free in parallel.
  logical, parameter, public :: THREADED_FACE_LOOPS = .true.

  ! Face type flags
  integer(ip), parameter, public :: FACE_INTERNAL  = 0
  integer(ip), parameter, public :: FACE_BOUNDARY  = 1

  ! Boundary condition type flags
  integer(ip), parameter, public :: BC_INLET        = 1
  integer(ip), parameter, public :: BC_OUTLET       = 2
  integer(ip), parameter, public :: BC_WALL_NOSLIP  = 3
  integer(ip), parameter, public :: BC_WALL_MODEL   = 4
  integer(ip), parameter, public :: BC_SLIP         = 5
  integer(ip), parameter, public :: BC_PERIODIC     = 6
  integer(ip), parameter, public :: BC_SYMMETRY     = 7
  integer(ip), parameter, public :: BC_FIXEDTEMP    = 8
  integer(ip), parameter, public :: BC_ZEROGRADIENT = 9
  integer(ip), parameter, public :: BC_WALL_VELOCITY = 10
  integer(ip), parameter, public :: BC_ABL_INLET    = 11

  ! Advection scheme flags
  integer(ip), parameter, public :: SCHEME_VANLEER_TVD = 1
  integer(ip), parameter, public :: SCHEME_KOREN_TVD   = 2
  integer(ip), parameter, public :: SCHEME_QUICK = 5
  integer(ip), parameter, public :: SCHEME_CDS   = 3
  integer(ip), parameter, public :: SCHEME_CENTRAL = SCHEME_CDS
  integer(ip), parameter, public :: SCHEME_KEP   = 4
  integer(ip), parameter, public :: SCHEME_SKEW  = SCHEME_KEP

  ! Wray low-storage RK3 coefficients (Williamson 1980, Table 1); abscissas c=(0, 8/15, 2/3).
  ! Stage k: u*=u^k+dt*(alpha_k*H^k+beta_k*H_prev)-gamma_k*dt*grad(p^{k-1/2}); p+=phi (incremental).
  real(rp), parameter, public :: RK3W_ALPHA(3) = &
      [8.0_rp/15.0_rp,  5.0_rp/12.0_rp,  3.0_rp/4.0_rp]
  real(rp), parameter, public :: RK3W_BETA(3)  = &
      [0.0_rp,         -17.0_rp/60.0_rp, -5.0_rp/12.0_rp]
  real(rp), parameter, public :: RK3W_GAMMA(3) = &
      [8.0_rp/15.0_rp,  2.0_rp/15.0_rp,  1.0_rp/3.0_rp]
  real(rp), parameter, public :: RK3W_C(3) = &
      [0.0_rp,           8.0_rp/15.0_rp,  2.0_rp/3.0_rp]

end module mod_parameters
