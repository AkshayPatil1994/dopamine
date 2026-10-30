! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Checks the cos(theta) floor of the over-relaxed k-vector (serial, no mesh).
program test_nonortho_limiter
  use mod_kinds
  use mod_nonortho
  implicit none

  real(rp) :: sf(3), d(3), k0(3), k1(3), g(3), sd, dl, a, a0
  integer  :: nfail

  nfail = 0

  ! Orthogonal face: floor inactive, k = 0.
  sf = [2.0_rp, 0.0_rp, 0.0_rp]; d = [0.5_rp, 0.0_rp, 0.0_rp]
  call nonortho_kfactor_limited(sf, d, 0.3_rp, k1, sd, dl)
  call check('orthogonal k', maxval(abs(k1)), 0.0_rp)

  ! Floor 0 (default) is bitwise the over-relaxed routine, even at 89.4 degrees.
  sf = [1.0_rp, 0.0_rp, 0.0_rp]; d = [0.01_rp, 1.0_rp, 0.0_rp]
  call nonortho_kfactor_overrelaxed(sf, d, k0, sd, dl)
  call nonortho_kfactor_limited(sf, d, 0.0_rp, k1, sd, dl)
  call check('floor off identical', maxval(abs(k1 - k0)), 0.0_rp)

  ! Floor below the face cos(theta) changes nothing.
  call nonortho_kfactor_limited(sf, d, 0.005_rp, k1, sd, dl)
  call check('inactive floor identical', maxval(abs(k1 - k0)), 0.0_rp)

  ! Active floor: implicit coefficient a = (Sf - k).d / (d.d) is bounded by |Sf|/(|d| cmin).
  a0 = (sf(1) - k0(1)) / d(1)
  call nonortho_kfactor_limited(sf, d, 0.2_rp, k1, sd, dl)
  a = (sf(1) - k1(1)) / d(1)
  call check('unlimited coefficient large', merge(1.0_rp, 0.0_rp, a0 > 50.0_rp), 1.0_rp)
  call check('limited coefficient bound', a, 1.0_rp / (sqrt(sum(d**2)) * 0.2_rp))

  ! Consistency: a (d.g) + k.g = Sf.g for any gradient, so linear fields stay exact.
  g = [0.7_rp, -1.3_rp, 2.1_rp]
  call check('linear exactness', a * dot_product(d, g) + dot_product(k1, g), dot_product(sf, g))

  if (nfail > 0) error stop 'test_nonortho_limiter FAILED'
  print *, 'test_nonortho_limiter PASSED'

contains

  subroutine check(name, got, want)
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: got, want
    if (abs(got - want) > 1.0e-12_rp * max(1.0_rp, abs(want))) then
      nfail = nfail + 1
      print '(A,A,2ES14.6)', 'FAIL ', name, got, want
    end if
  end subroutine check

end program test_nonortho_limiter
