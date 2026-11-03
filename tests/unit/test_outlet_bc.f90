! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Order of the convective-outlet face derivative in space and of its RK3W integration in time (serial, no mesh).
program test_outlet_bc
  use mod_kinds
  use mod_parameters
  use mod_bc, only: outlet_dn_u
  implicit none

  integer  :: nfail
  nfail = 0

  call spatial_order(skew=.false.)
  call spatial_order(skew=.true.)
  call temporal_order()

  if (nfail > 0) error stop 'test_outlet_bc FAILED'
  write(*,'(A)') 'test_outlet_bc passed'

contains

  real(rp) function fld(x)
    real(rp), intent(in) :: x(3)
    fld = sin(1.3_rp*x(1) + 0.4_rp*x(2)) + 0.5_rp*cos(0.7_rp*x(3) + x(1))
  end function fld

  function grd(x) result(g)
    real(rp), intent(in) :: x(3)
    real(rp) :: g(3)
    g(1) = 1.3_rp*cos(1.3_rp*x(1) + 0.4_rp*x(2)) - 0.5_rp*sin(0.7_rp*x(3) + x(1))
    g(2) = 0.4_rp*cos(1.3_rp*x(1) + 0.4_rp*x(2))
    g(3) = -0.35_rp*sin(0.7_rp*x(3) + x(1))
  end function grd

  !> Error of n.grad u at the face with the exact owner gradient, new formula against the first-order (u_b-u_o)/|delta|.
  subroutine spatial_order(skew)
    logical, intent(in) :: skew
    real(rp) :: xo(3), xb(3), dl(3), nrm(3), tn(3), h, len_d, nd, e_new(4), e_old(4), exact
    real(rp) :: ord_new, ord_old
    integer  :: i

    xo = [0.3_rp, 0.2_rp, 0.1_rp]
    nrm = [1.0_rp, 0.0_rp, 0.0_rp]
    do i = 1, 4
      h = 0.4_rp / 2.0_rp**(i - 1)
      xb = xo + h*[0.5_rp, 0.0_rp, 0.0_rp]
      if (skew) xb = xb + h*[0.0_rp, 0.15_rp, -0.1_rp]
      dl = xb - xo
      len_d = sqrt(sum(dl**2))
      dl = dl / len_d
      nd = dot_product(nrm, dl)
      tn = nrm - nd*dl
      exact = dot_product(nrm, grd(xb))
      e_new(i) = abs(outlet_dn_u(fld(xb), fld(xo), dot_product(grd(xo), dl), dot_product(grd(xo), tn), len_d, nd) - exact)
      e_old(i) = abs(nd*(fld(xb) - fld(xo))/len_d + dot_product(grd(xo), tn) - exact)
    end do
    ord_new = log(e_new(3)/e_new(4)) / log(2.0_rp)
    ord_old = log(e_old(3)/e_old(4)) / log(2.0_rp)
    write(*,'(A,L1,A,F6.2,A,F6.2)') 'spatial order skew=', skew, ': new ', ord_new, '  first-order formula ', ord_old
    if (skew) then
      call require('skewed offset: new formula at least first order and more accurate', &
                   ord_new > 0.9_rp .and. e_new(4) < e_old(4))
    else
      call require('orthogonal offset: new formula second order', ord_new > 1.9_rp)
      call require('orthogonal offset: first-order formula is first order', ord_old < 1.1_rp)
    end if
  end subroutine spatial_order

  !> dU_b/dt = -Uc n.grad u with the exact travelling-wave owner data, integrated with the Wray weights of the solver.
  subroutine temporal_order()
    real(rp) :: e(4), uref, ord
    integer  :: i

    uref = integrate(4096)
    do i = 1, 4
      e(i) = abs(integrate(16 * 2**(i - 1)) - uref)
    end do
    ord = log(e(3)/e(4)) / log(2.0_rp)
    write(*,'(A,F6.2)') 'temporal order of the boundary state (RK3W weights): ', ord
    call require('boundary state integrated at third order', ord > 2.8_rp .and. ord < 3.3_rp)
  end subroutine temporal_order

  real(rp) function integrate(nstep) result(ub)
    integer, intent(in) :: nstep
    real(rp), parameter :: uc = 1.3_rp, wn = 3.0_rp, tend = 1.0_rp, h = 0.2_rp
    real(rp) :: dt, t, g, gprev
    integer  :: istep, k

    dt = tend / nstep
    ub = sin(wn*(0.5_rp*h))
    gprev = 0.0_rp
    do istep = 1, nstep
      do k = 1, 3
        t = (istep - 1)*dt + RK3W_C(k)*dt
        g = -uc * outlet_dn_u(ub, sin(wn*(-uc*t)), wn*cos(wn*(-uc*t)), 0.0_rp, 0.5_rp*h, 1.0_rp)
        ub = ub + dt*(RK3W_ALPHA(k)*g + RK3W_BETA(k)*gprev)
        gprev = g
      end do
    end do
  end function integrate

  subroutine require(label, ok)
    character(*), intent(in) :: label
    logical,      intent(in) :: ok
    if (.not. ok) then
      write(*,'(2A)') 'FAIL: ', label
      nfail = nfail + 1
    end if
  end subroutine require

end program test_outlet_bc
