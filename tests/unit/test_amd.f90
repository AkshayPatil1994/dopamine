! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Point checks of the AMD closure against hand-evaluated values (no mesh, serial).
program test_amd
  use mod_kinds
  use mod_amd, only: amd_point_viscosity
  implicit none

  real(rp), parameter :: C = 1.0_rp / 3.0_rp
  real(rp) :: a(3,3), g(3,3), nu
  integer  :: nfail

  nfail = 0
  g = 0.0_rp
  g(1,1) = 1.0_rp; g(2,2) = 1.0_rp; g(3,3) = 1.0_rp

  ! Pure shear du/dy: the model switches off.
  a = 0.0_rp; a(2,1) = 3.0_rp
  call amd_point_viscosity(a, g, C, nu)
  call check('pure shear', nu, 0.0_rp)

  ! Solid-body rotation: alpha antisymmetric, no dissipation.
  a = 0.0_rp; a(2,1) = -2.0_rp; a(1,2) = 2.0_rp
  call amd_point_viscosity(a, g, C, nu)
  call check('solid-body rotation', nu, 0.0_rp)

  ! Axisymmetric strain diag(1,1,-2), G = I: beta = diag(1,1,4), num = 6, alpha_sq = 6, nu = C.
  a = 0.0_rp; a(1,1) = 1.0_rp; a(2,2) = 1.0_rp; a(3,3) = -2.0_rp
  call amd_point_viscosity(a, g, C, nu)
  call check('isotropic filter strain', nu, C)

  ! Same strain, G = diag(4,1,1): beta = diag(4,1,4), num = 3, nu = C/2.
  g(1,1) = 4.0_rp
  call amd_point_viscosity(a, g, C, nu)
  call check('anisotropic filter strain', nu, 0.5_rp * C)

  ! Reversed strain (backscatter): the numerator is negative and nu_t is clipped to zero.
  g(1,1) = 1.0_rp
  call amd_point_viscosity(-a, g, C, nu)
  call check('backscatter clipped', nu, 0.0_rp)

  ! Quiescent cell.
  a = 0.0_rp
  call amd_point_viscosity(a, g, C, nu)
  call check('zero gradient', nu, 0.0_rp)

  ! The viscosity scales linearly with both the constant and the gradient magnitude (5x strain: num/alpha_sq = 5).
  a = 0.0_rp; a(1,1) = 5.0_rp; a(2,2) = 5.0_rp; a(3,3) = -10.0_rp
  call amd_point_viscosity(a, g, 2.0_rp * C, nu)
  call check('constant scaling', nu, 2.0_rp * C * 5.0_rp)

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'test_amd: ', nfail, ' FAILED'
    error stop 1
  end if
  write(*,'(A)') 'test_amd: all checks passed'

contains

  subroutine check(name, got, want)
    character(*), intent(in) :: name
    real(rp),     intent(in) :: got, want
    if (abs(got - want) > 1.0e-12_rp * max(1.0_rp, abs(want))) then
      write(*,'(A,A,2ES16.8)') 'FAIL ', name, got, want
      nfail = nfail + 1
    else
      write(*,'(A,A)') 'ok   ', name
    end if
  end subroutine check

end program test_amd
