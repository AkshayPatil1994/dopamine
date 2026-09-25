! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Non-orthogonal face correction (Jasak 1996, Sect. 3.3.2 minimum-correction split).
module mod_nonortho
  use mod_kinds
  implicit none
  private

  !> sin^2(theta) threshold below which a face is treated as orthogonal.
  !> 1e-10 -> theta ~= 6x10^-4 degrees; set larger (e.g. 1e-6) for looser skipping.
  real(rp), parameter, public :: NONORTHO_SIN2_TOL = 1.0e-10_rp

  public :: nonortho_kfactor
  public :: nonortho_kfactor_overrelaxed
  public :: nonortho_cross_flux
  public :: nonortho_rhie_chow

contains

  !> Compute the minimum-correction k-vector and auxiliary scalars.
  !> @param sf face area vector (outward from owner) [m^2]
  pure subroutine nonortho_kfactor(sf, d_on, k, s_dot_dhat, d_on_len)
    !$acc routine seq
    real(rp), intent(in)  :: sf(3)
    real(rp), intent(in)  :: d_on(3)
    real(rp), intent(out) :: k(3)
    real(rp), intent(out) :: s_dot_dhat
    real(rp), intent(out) :: d_on_len

    real(rp) :: inv_len, k_mag2, e_mag

    d_on_len = sqrt(sum(d_on**2))
    if (d_on_len < 1.0e-30_rp) then
      k          = 0.0_rp
      s_dot_dhat = 0.0_rp
      d_on_len   = 1.0e-30_rp
      return
    end if

    inv_len    = 1.0_rp / d_on_len
    s_dot_dhat = dot_product(sf, d_on) * inv_len
    k          = sf - (s_dot_dhat * inv_len) * d_on

    k_mag2 = sum(k**2)
    e_mag  = abs(s_dot_dhat)
    if (k_mag2 > e_mag * e_mag .and. k_mag2 > 0.0_rp) then
      k = k * (e_mag / sqrt(k_mag2))
    end if
  end subroutine nonortho_kfactor

  !> Over-relaxed k-vector (Jasak 1996 Sect. 3.3.3): k = Sf - (|Sf|^2/(Sf.d_on))*d_on; no cap (spectral radius = sin theta < 1).
  pure subroutine nonortho_kfactor_overrelaxed(sf, d_on, k, s_dot_dhat, d_on_len)
    !$acc routine seq
    real(rp), intent(in)  :: sf(3)
    real(rp), intent(in)  :: d_on(3)
    real(rp), intent(out) :: k(3)
    real(rp), intent(out) :: s_dot_dhat
    real(rp), intent(out) :: d_on_len

    real(rp) :: d_len, sf_sq, inv_len

    d_on_len = sqrt(sum(d_on**2))
    if (d_on_len < 1.0e-30_rp) then
      k          = 0.0_rp
      s_dot_dhat = 0.0_rp
      d_on_len   = 1.0e-30_rp
      return
    end if

    inv_len    = 1.0_rp / d_on_len
    s_dot_dhat = dot_product(sf, d_on) * inv_len
    d_len      = s_dot_dhat * d_on_len  ! = dot_product(sf, d_on)
    sf_sq      = sum(sf**2)

    if (d_len > 1.0e-30_rp) then
      ! Over-relaxed k; no cap - spectral radius sin theta < 1 ensures convergence.
      k = sf - (sf_sq / d_len) * d_on
    else
      ! Degenerate / near-perpendicular face (Sf.d_on ~= 0): minimum-correction.
      ! |k_min| = |Sf| sin theta <= |Sf| always; no cap needed.
      k = sf - (s_dot_dhat * inv_len) * d_on
    end if
  end subroutine nonortho_kfactor_overrelaxed

  !> Rhie-Chow face term; subtracting the cross part makes it vanish for a linear pressure field.
  pure function nonortho_rhie_chow(grad_o, grad_nb, gc, sf, k, coeff, dp) result(rc)
    !$acc routine seq
    real(rp), intent(in) :: grad_o(3)
    real(rp), intent(in) :: grad_nb(3)
    real(rp), intent(in) :: gc
    real(rp), intent(in) :: sf(3)
    real(rp), intent(in) :: k(3)  ! over-relaxed non-orthogonal vector; zero on an orthogonal face
    real(rp), intent(in) :: coeff  ! face Laplacian coefficient |Sf|^2/(Sf.d_on)
    real(rp), intent(in) :: dp  ! p_nb - p_o
    real(rp)             :: rc

    real(rp) :: grad_mean(3)

    ! sf - k is the orthogonal projection the compact jump represents, so a uniform gradient cancels.
    grad_mean = gc * grad_o + (1.0_rp - gc) * grad_nb
    rc = dot_product(grad_mean, sf - k) - coeff * dp
  end function nonortho_rhie_chow

  !> Scalar cross-diffusion flux: (gc.grad phi_o + (1-gc).grad phi_nb) . k  !> @param grad_o gradient at owner cell [phi/m]
  pure function nonortho_cross_flux(grad_o, grad_nb, gc, k) result(flux)
    !$acc routine seq
    real(rp), intent(in) :: grad_o(3)
    real(rp), intent(in) :: grad_nb(3)
    real(rp), intent(in) :: gc
    real(rp), intent(in) :: k(3)
    real(rp)             :: flux
    flux = dot_product(gc * grad_o + (1.0_rp - gc) * grad_nb, k)
  end function nonortho_cross_flux

end module mod_nonortho
