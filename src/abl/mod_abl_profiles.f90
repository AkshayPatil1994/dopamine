! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Monin-Obukhov Similarity Theory (MOST) profiles and Businger-Dyer stability functions
module mod_abl_profiles
  use mod_kinds
  use mod_parameters, only: SMALL, KAPPA, PI
  implicit none
  private

  public :: abl_u_profile, abl_T_profile
  public :: abl_psi_m, abl_psi_h
  public :: abl_ustar_from_ref, abl_tstar_from_ref, abl_wind_vector

contains

  !> Neutral friction velocity matching u_ref at z_ref: u* = kappa*u_ref/ln(z_ref/z0_m).
  pure real(rp) function abl_ustar_from_ref(u_ref, z_ref, z0_m) result(u_star)
    !$acc routine seq
    real(rp), intent(in) :: u_ref, z_ref, z0_m
    if (z_ref > z0_m .and. z0_m > SMALL .and. u_ref > SMALL) then
      u_star = KAPPA * u_ref / log(z_ref / z0_m)
    else
      u_star = 0.05_rp * u_ref  ! fallback for a degenerate roughness setup
    end if
  end function abl_ustar_from_ref

  !> Neutral temperature scale matching T_ref at z_ref over a surface at T_surface.
  pure real(rp) function abl_tstar_from_ref(T_ref, T_surface, z_ref, z0_h) result(T_star)
    !$acc routine seq
    real(rp), intent(in) :: T_ref, T_surface, z_ref, z0_h
    if (z_ref > z0_h .and. z0_h > SMALL) then
      T_star = KAPPA * (T_ref - T_surface) / log(z_ref / z0_h)
    else
      T_star = 0.0_rp
    end if
  end function abl_tstar_from_ref

  !> Horizontal wind vector from speed and meteorological direction (degrees the wind blows from).
  pure function abl_wind_vector(u_mag, wind_dir_deg) result(uh)
    !$acc routine seq
    real(rp), intent(in) :: u_mag, wind_dir_deg
    real(rp) :: uh(3), rad
    rad   = wind_dir_deg * PI / 180.0_rp
    uh(1) = -u_mag * sin(rad)
    uh(2) = -u_mag * cos(rad)
    uh(3) = 0.0_rp
  end function abl_wind_vector

  !> Integrated stability correction for momentum: psi_m(zeta)
  !> Used in the log-law profile: u(z) = (u*/kappa) * [ln(z/z0) - psi_m(z/L) + psi_m(z0/L)]
  real(rp) function abl_psi_m(zeta) result(psi)
    !$acc routine seq
    real(rp), intent(in) :: zeta
    real(rp) :: x

    if (abs(zeta) < 1.0e-6_rp) then
      psi = 0.0_rp
    else if (zeta > 0.0_rp) then
      ! Stable: psi_m = -5*zeta
      psi = -5.0_rp * zeta
    else
      ! Unstable: complex formula from Paulson (1970)
      x = (1.0_rp - 16.0_rp * zeta) ** 0.25_rp
      psi = 2.0_rp * log((1.0_rp + x) / 2.0_rp) &
          + log((1.0_rp + x*x) / 2.0_rp) &
          - 2.0_rp * atan(x) + PI / 2.0_rp
    end if
  end function abl_psi_m

  !> Integrated stability correction for heat: psi_h(zeta)
  real(rp) function abl_psi_h(zeta) result(psi)
    !$acc routine seq
    real(rp), intent(in) :: zeta
    real(rp) :: y

    if (abs(zeta) < 1.0e-6_rp) then
      psi = 0.0_rp
    else if (zeta > 0.0_rp) then
      ! Stable: psi_h = -5*zeta
      psi = -5.0_rp * zeta
    else
      ! Unstable: simplified formula
      y = (1.0_rp - 16.0_rp * zeta) ** 0.5_rp
      psi = 2.0_rp * log((1.0_rp + y) / 2.0_rp)
    end if
  end function abl_psi_h

  !> MOST velocity profile: u(z) given friction velocity and roughness
  !> u(z) = (u*/kappa) * [ln(z/z0_m) - psi_m(z/L) + psi_m(z0_m/L)]
  real(rp) function abl_u_profile(z, u_star, z0_m, L) result(u_z)
    !$acc routine seq
    real(rp), intent(in) :: z, u_star, z0_m, L
    real(rp) :: zeta, zeta0

    u_z = 0.0_rp

    if (z <= z0_m .or. u_star <= SMALL .or. z0_m <= SMALL) return

    if (abs(L) < SMALL) then
      ! Neutral: no stability correction
      u_z = (u_star / KAPPA) * log(z / z0_m)
    else
      zeta = z / L
      zeta0 = z0_m / L
      u_z = (u_star / KAPPA) * (log(z / z0_m) - abl_psi_m(zeta) + abl_psi_m(zeta0))
    end if
  end function abl_u_profile

  !> MOST temperature profile: T(z) given surface temperature and temperature scale
  !> T(z) = T_surface + (T*/kappa) * [ln(z/z0_h) - psi_h(z/L) + psi_h(z0_h/L)]
  real(rp) function abl_T_profile(z, T_surface, T_star, z0_h, L) result(T_z)
    !$acc routine seq
    real(rp), intent(in) :: z, T_surface, T_star, z0_h, L
    real(rp) :: zeta, zeta0

    T_z = T_surface

    if (z <= z0_h .or. abs(T_star) <= SMALL .or. z0_h <= SMALL) return

    if (abs(L) < SMALL) then
      ! Neutral: no stability correction
      T_z = T_surface + (T_star / KAPPA) * log(z / z0_h)
    else
      zeta = z / L
      zeta0 = z0_h / L
      T_z = T_surface + (T_star / KAPPA) * (log(z / z0_h) - abl_psi_h(zeta) + abl_psi_h(zeta0))
    end if
  end function abl_T_profile

end module mod_abl_profiles
