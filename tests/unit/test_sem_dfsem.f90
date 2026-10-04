! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Unit test of the divergence-free synthetic eddy kernel (no mesh, no MPI-dependent paths).
!>  1. one eddy is solenoidal to finite-difference accuracy and has zero volume-integrated velocity;
!>  2. the eddy ensemble on a doubly periodic inlet plane reproduces the target Reynolds tensor, with and
!>     without shear stress and for a rotated mean-wind direction (prescribed full-tensor tables);
!>  3. the Jacobi eigen-solver reconstructs its input;
!>  4. the binned evaluation equals the all-eddy sum (both kernels; periodic, partly periodic and open inlets,
!>     two-bin periodic widths, query points outside the box, time look-ahead, after the eddies have moved).
program test_sem_dfsem
  use mod_kinds
  use mod_abl_types
  use mod_sem_inflow
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none

  integer :: nfail = 0

  call check_eigen()
  call check_single_eddy()
  call check_statistics('analytic stresses, wind from 270', 270.0_rp, .false.)
  call check_statistics('analytic stresses, wind from 235', 235.0_rp, .false.)
  call check_statistics('prescribed tensor with shear uw', 270.0_rp, .true.)
  call check_bins('dfsem, doubly periodic', 'dfsem', .true., .true., 1.0_rp, 0.1_rp)
  call check_bins('dfsem, periodic in t1 only', 'dfsem', .true., .false., 1.0_rp, 0.1_rp)
  call check_bins('dfsem, open inlet', 'dfsem', .false., .false., 1.0_rp, 0.1_rp)
  call check_bins('dfsem, two bins across the period', 'dfsem', .true., .true., 0.45_rp, 0.2_rp)
  call check_bins('jarrin, doubly periodic', 'jarrin', .true., .true., 1.0_rp, 0.1_rp)
  call check_bins('jarrin, open inlet', 'jarrin', .false., .false., 1.0_rp, 0.1_rp)

  if (nfail > 0) then
    write(*,'(A,I0,A)') 'FAILED: ', nfail, ' check(s)'
    error stop 1
  end if
  write(*,'(A)') 'test_sem_dfsem: all checks passed'

contains

  subroutine report(name, ok, val, tol)
    character(len=*), intent(in) :: name
    logical,          intent(in) :: ok
    real(rp),         intent(in) :: val, tol
    write(*,'(A,A,A,ES10.3,A,ES9.2,A)') merge('PASS  ', 'FAIL  ', ok), name, '  value=', val, '  tol=', tol, ''
    if (.not. ok) nfail = nfail + 1
  end subroutine report

  real(rp) function urand(s) result(u)
    integer(int64), intent(inout) :: s
    s = ieor(s, ishft(s, 13)); s = ieor(s, ishft(s, -7)); s = ieor(s, ishft(s, 17))
    u = real(iand(s, 4503599627370495_int64), rp) / 4503599627370496.0_rp
  end function urand

  subroutine check_eigen()
    real(rp) :: A(3,3), lam(3), V(3,3), B(3,3), err
    integer  :: i
    A = reshape([5.0_rp, -1.0_rp, 0.7_rp, -1.0_rp, 3.0_rp, 0.2_rp, 0.7_rp, 0.2_rp, 1.5_rp], [3,3])
    call sem_eigen3(A, lam, V)
    B = 0.0_rp
    do i = 1, 3
      B = B + lam(i) * spread(V(:, i), 2, 3) * spread(V(:, i), 1, 3)
    end do
    err = maxval(abs(B - A)) + maxval(abs(matmul(transpose(V), V) - reshape([1,0,0,0,1,0,0,0,1], [3,3])))
    call report('Jacobi eigen-solver reconstructs A', err < 1.0e-13_rp, err, 1.0e-13_rp)
  end subroutine check_eigen

  subroutine check_single_eddy()
    real(rp) :: frame(9), sig(3), ampl(3), d(3), h, div, gradmax, flux_sum, dv
    real(rp) :: up(3), um(3), g(3,3), mean(3), ex(3)
    integer(int64) :: s
    integer  :: n, k, j, m
    real(rp) :: th
    ! A rotated principal frame and unequal half-extents.
    th = 0.6_rp
    frame = [cos(th), sin(th), 0.0_rp, -sin(th), cos(th), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]
    sig   = [0.30_rp, 0.21_rp, 0.12_rp]
    ampl  = [0.8_rp, -1.3_rp, 0.5_rp]
    s = 123456789_int64
    h = 1.0e-6_rp
    div = 0.0_rp; gradmax = 0.0_rp
    do n = 1, 2000
      d = [urand(s) - 0.5_rp, urand(s) - 0.5_rp, urand(s) - 0.5_rp] * 0.6_rp
      do j = 1, 3
        ex = 0.0_rp; ex(j) = h
        up = sem_dfsem_eddy_velocity(d + ex, frame, sig, ampl)
        um = sem_dfsem_eddy_velocity(d - ex, frame, sig, ampl)
        g(:, j) = (up - um) / (2.0_rp * h)
      end do
      div     = max(div, abs(g(1,1) + g(2,2) + g(3,3)))
      gradmax = max(gradmax, maxval(abs(g)))
    end do
    call report('single eddy: max |div u| / max |grad u|', div / gradmax < 1.0e-6_rp, div / gradmax, 1.0e-6_rp)

    ! Volume integral of u on a fine grid covering the support: zero for the compactly supported eddy.
    m = 80; mean = 0.0_rp; flux_sum = 0.0_rp
    dv = (0.8_rp / m)**3
    do k = 1, m
      do j = 1, m
        do n = 1, m
          d = [(n - 0.5_rp) / m - 0.5_rp, (j - 0.5_rp) / m - 0.5_rp, (k - 0.5_rp) / m - 0.5_rp] * 0.8_rp
          mean = mean + sem_dfsem_eddy_velocity(d, frame, sig, ampl) * dv
          flux_sum = flux_sum + sum(abs(sem_dfsem_eddy_velocity(d, frame, sig, ampl))) * dv
        end do
      end do
    end do
    call report('single eddy: |int u dV| / int |u| dV', maxval(abs(mean)) / flux_sum < 1.0e-3_rp, &
                maxval(abs(mean)) / flux_sum, 1.0e-3_rp)
  end subroutine check_single_eddy

  !> Expected covariance at a point from N identical eddies with uniformly random centres and independent
  !> random signs: (N / V_box) * sign-average of the integral of u u^T over the eddy support (midpoint rule).
  function quad_cov(st, ne) result(cov)
    type(SEMState), intent(in) :: st
    integer(ip),    intent(in) :: ne
    real(rp) :: cov(3,3)

    integer(ip), parameter :: m = 80
    real(rp)    :: d(3), u(3), ampl(3), h, dv
    integer(ip) :: i, j, k, c

    h = 2.0_rp * st%rsup / m
    dv = h**3
    cov = 0.0_rp
    do c = 0, 7
      ampl = abs(st%ampl(:, 1)) * [merge(1.0_rp, -1.0_rp, btest(c, 0)), merge(1.0_rp, -1.0_rp, btest(c, 1)), &
                                   merge(1.0_rp, -1.0_rp, btest(c, 2))]
      do k = 1, m
        do j = 1, m
          do i = 1, m
            d = [(i - 0.5_rp) * h - st%rsup, (j - 0.5_rp) * h - st%rsup, (k - 0.5_rp) * h - st%rsup]
            u = sem_dfsem_eddy_velocity(d, st%frame(:, 1), st%sig(:, 1), ampl)
            cov = cov + spread(u, 2, 3) * spread(u, 1, 3) * dv / 8.0_rp
          end do
        end do
      end do
    end do
    cov = cov * ne / st%box_volume
  end function quad_cov

  !> Eddies fill a doubly periodic unit inlet plane; sample the fluctuation on that plane and compare its
  !> covariance with the target tensor (in global axes).
  subroutine check_statistics(label, wind_deg, shear)
    character(len=*), intent(in) :: label
    real(rp),         intent(in) :: wind_deg
    logical,          intent(in) :: shear

    type(SEMState)  :: st
    type(ABLConfig) :: cfg
    integer(int64) :: s
    integer(ip), parameter :: ne = 6000, ns = 20000
    integer(ip) :: i, k, a, b
    real(rp) :: fc(3), um(3), uf(3), cov(3,3), Rm(3,3), Rg(3,3), Wr(3,3), th, err, scale, c(3), quad_err

    cfg%wind_dir_deg = wind_deg
    cfg%u_ref = 10.0_rp; cfg%z_ref = 0.5_rp; cfg%z0_m = 0.01_rp
    cfg%sem_type = 'dfsem'
    cfg%sem_sigma = 0.1_rp
    cfg%n_eddies = ne

    st%n_in = [1.0_rp, 0.0_rp, 0.0_rp]; st%t1 = [0.0_rp, 1.0_rp, 0.0_rp]; st%t2 = [0.0_rp, 0.0_rp, 1.0_rp]
    st%x_inlet = 0.0_rp
    st%per1 = .true.; st%per2 = .true.
    st%period1 = 1.0_rp; st%period2 = 1.0_rp; st%per1_lo = 0.0_rp; st%per2_lo = 0.0_rp
    st%t1_min = 0.0_rp; st%t1_max = 1.0_rp; st%t2_min = 0.0_rp; st%t2_max = 1.0_rp
    st%sigma = cfg%sem_sigma
    st%dfsem = .true.
    if (shear) then
      ! Height-independent full tensor with a strong shear stress (global axes): off the principal frame.
      st%has_profile = .true.; st%prof_global_axes = .true.; st%n_levels = 2
      allocate(st%prof_z(2), st%prof_U(3, 2), st%prof_R(6, 2))
      st%prof_z = [0.0_rp, 1.0_rp]
      st%prof_U = 0.0_rp; st%prof_U(1, :) = 10.0_rp
      st%prof_R(:, 1) = [4.0_rp, 2.5_rp, 1.2_rp, 0.3_rp, -1.0_rp, 0.1_rp]
      st%prof_R(:, 2) = st%prof_R(:, 1)
    end if
    call sem_reference_scales(st, cfg)
    st%box_volume = 2.0_rp * st%rsup * 1.0_rp * 1.0_rp
    st%n_eddies = ne
    allocate(st%pos(3, ne), st%amp(3, ne), st%frame(9, ne), st%sig(3, ne), st%ampl(3, ne))
    s = 88172645463325252_int64
    do i = 1, ne
      st%pos(:, i) = [(2.0_rp * urand(s) - 1.0_rp) * st%rsup, urand(s), urand(s)]
      do k = 1, 3
        st%amp(k, i) = sign(1.0_rp, urand(s) - 0.5_rp)
      end do
    end do
    st%initialized = .true.
    call sem_prepare(st, cfg)

    cov = 0.0_rp; c = 0.0_rp
    do i = 1, ns
      fc = [0.0_rp, urand(s), urand(s)]
      call sem_velocity_at_face(st, cfg, fc, um, uf)
      do a = 1, 3
        do b = 1, 3
          cov(a, b) = cov(a, b) + uf(a) * uf(b) / ns
        end do
      end do
      c = c + uf / ns
    end do

    ! Target tensor in global axes.
    call sem_stress_at(st, cfg, 0.5_rp, Rm, Wr)
    Rg = matmul(Wr, matmul(Rm, transpose(Wr)))
    scale = maxval(abs(Rg))
    err = maxval(abs(cov - Rg)) / scale
    quad_err = maxval(abs(quad_cov(st, ne) - Rg)) / scale
    th = maxval(abs(c)) / sqrt(scale)
    call report(label//': eddy second moment x density vs R (quadrature)', quad_err < 1.0e-2_rp, quad_err, 1.0e-2_rp)
    call report(label//': plane sampling |cov - R| / max|R|', err < 0.15_rp, err, 0.15_rp)
    call report(label//': |mean u| / sqrt(max R)', th < 0.05_rp, th, 0.05_rp)
    write(*,'(A,6ES11.3)') '        cov  uu vv ww uv uw vw: ', cov(1,1), cov(2,2), cov(3,3), cov(1,2), cov(1,3), cov(2,3)
    write(*,'(A,6ES11.3)') '        R    uu vv ww uv uw vw: ', Rg(1,1), Rg(2,2), Rg(3,3), Rg(1,2), Rg(1,3), Rg(2,3)
  end subroutine check_statistics

  !> Binned evaluation against the brute-force sum over all eddies on the same state.
  subroutine check_bins(label, kernel, per1, per2, width, sigma)
    character(len=*), intent(in) :: label, kernel
    logical,          intent(in) :: per1, per2
    real(rp),         intent(in) :: width, sigma

    type(SEMState)  :: st, ref
    type(ABLConfig) :: cfg
    integer(int64) :: s
    integer(ip), parameter :: nq = 1500
    integer(ip) :: iq, istep
    real(rp) :: fc(3), um(3), uf(3), ur(3), err, scale, tau, ext

    cfg%wind_dir_deg = 270.0_rp
    cfg%u_ref = 10.0_rp; cfg%z_ref = 0.5_rp; cfg%z0_m = 0.01_rp
    cfg%sem_type = kernel
    cfg%sem_sigma = sigma
    cfg%n_eddies = 4000

    st%n_in = [1.0_rp, 0.0_rp, 0.0_rp]; st%t1 = [0.0_rp, 1.0_rp, 0.0_rp]; st%t2 = [0.0_rp, 0.0_rp, 1.0_rp]
    st%x_inlet = 0.0_rp
    st%per1 = per1; st%per2 = per2
    st%period1 = width; st%period2 = width; st%per1_lo = 0.0_rp; st%per2_lo = 0.0_rp
    st%t1_min = 0.0_rp; st%t1_max = width; st%t2_min = 0.0_rp; st%t2_max = width
    call sem_setup(st, cfg)
    ext = st%rsup

    s = 1234567_int64
    err = 0.0_rp; scale = 0.0_rp
    do istep = 0, 3
      if (istep > 0) then
        call sem_advance(st, 0.37_rp * st%rsup / st%u_convect)
        call sem_prepare(st, cfg)
      end if
      ref = st
      ref%bins_ready = .false.
      do iq = 1, nq
        ! Points inside the box, beyond its edges and (open inlets) well outside it.
        fc = [0.0_rp, (1.6_rp * urand(s) - 0.3_rp) * width, (1.6_rp * urand(s) - 0.3_rp) * width]
        tau = 0.0_rp
        if (mod(iq, 3) == 0) tau = 1.3_rp * ext / st%u_convect * (2.0_rp * urand(s) - 1.0_rp)
        call sem_velocity_at_face(st, cfg, fc, um, uf, tau)
        call sem_velocity_at_face(ref, cfg, fc, um, ur, tau)
        err = max(err, maxval(abs(uf - ur)))
        scale = max(scale, maxval(abs(ur)))
      end do
    end do
    call report(label//': binned == all-eddy sum (max |diff| / max |u|)', &
                st%bins_ready .and. err <= 1.0e-12_rp * max(scale, 1.0_rp), err / max(scale, 1.0e-300_rp), 1.0e-12_rp)
  end subroutine check_bins

end program test_sem_dfsem
