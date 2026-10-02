! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Verification of the divergence-free SEM inflow generator against its own targets:
!>   1. div(u') = 0 at arbitrary points (finite differences, second-order convergence)
!>   2. mean, turbulence intensities and u'w' versus the prescribed profile, at several heights
!>   3. height-dependent integral length scale from the two-point (Taylor) autocorrelation
!>   4. the stress/eddy-strength map: exact round trip, and a loud error when unrealisable
!>   5. wind-tunnel style inputs: turbulence-intensity tables and namelist laws
!>   6. the analytic Panofsky-Dutton mode carries u'w' = -u*^2
program test_sem_inflow
  use mod_kinds
  use mod_abl_types
  use mod_abl_profiles
  use mod_sem_inflow
  implicit none

  integer :: n_fail = 0
  real(rp), parameter :: SMALL_T = 1.0e-300_rp

  call test_coeff_roundtrip()
  call test_divergence_free()
  call test_bins_and_periodic()
  call test_stats_and_scales()
  call test_ti_table()
  call test_ti_law()

  if (n_fail > 0) then
    write(*,'(A,I0,A)') 'test_sem_inflow: ', n_fail, ' FAILED'
    error stop 1
  end if
  write(*,'(A)') 'test_sem_inflow: all passed'

contains

  subroutine check(cond, msg)
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (cond) then
      write(*,'(A,A)') '  ok    ', msg
    else
      write(*,'(A,A)') '  FAIL  ', msg
      n_fail = n_fail + 1
    end if
  end subroutine check

  !> Westerly wind (+x) over a 60 m x 60 m inlet at x = 0, with scales growing like z^0.3.
  subroutine make_state(state, cfg, n_eddies)
    type(SEMState),  intent(out) :: state
    type(ABLConfig), intent(out) :: cfg
    integer,         intent(in)  :: n_eddies

    cfg%u_ref = 10.0_rp; cfg%z_ref = 50.0_rp; cfg%z0_m = 0.05_rp
    cfg%wind_dir_deg = 270.0_rp
    cfg%n_eddies = n_eddies
    cfg%sem_sigma = 6.0_rp
    cfg%sem_L_exponent = 0.3_rp
    state%x_inlet = 0.0_rp
    state%n_in = [1.0_rp, 0.0_rp, 0.0_rp]
    state%t1   = [0.0_rp, 1.0_rp, 0.0_rp]
    state%t2   = [0.0_rp, 0.0_rp, 1.0_rp]
    state%t1_min = -30.0_rp; state%t1_max = 30.0_rp
    state%t2_min = 0.0_rp;   state%t2_max = 60.0_rp
  end subroutine make_state

  subroutine test_coeff_roundtrip()
    real(rp) :: R(6), sig(3), C(3,3), err
    write(*,'(A)') 'stress <-> eddy strength map'

    ! Elongated eddies (3:1.5:1) realise Panofsky-like anisotropy with a shear stress.
    R = [5.76_rp, 3.61_rp, 1.56_rp, 0.0_rp, -0.4_rp, 0.0_rp]
    sig = [3.0_rp, 1.5_rp, 1.0_rp]
    call sem_eddy_coeffs(R, sig, 3.0_rp, C, err)
    call check(err < 1.0e-12_rp, 'attainable stresses are reproduced to rounding')

    ! Isotropic eddies cannot give uu > vv + ww: must be reported, not hidden.
    sig = [1.0_rp, 1.0_rp, 1.0_rp]
    call sem_eddy_coeffs(R, sig, 3.0_rp, C, err)
    call check(err > 0.05_rp, 'unrealisable stresses are flagged (err > 5%)')
  end subroutine test_coeff_roundtrip

  subroutine test_divergence_free()
    type(SEMState)  :: state
    type(ABLConfig) :: cfg
    real(rp) :: ratio(3), fac(3), smin
    integer  :: m

    write(*,'(A)') 'divergence of the fluctuation field'
    call make_state(state, cfg, 20000)
    call sem_setup(state, cfg)

    smin = minval(state%sig)
    fac = [0.04_rp, 0.02_rp, 0.01_rp]
    do m = 1, 3
      ratio(m) = div_over_grad(state, fac(m) * smin)
    end do
    write(*,'(A,3ES10.2)') '    rms(div)/rms(grad) at h = 0.04, 0.02, 0.01 sigma_min: ', ratio
    call check(ratio(3) < 1.0e-4_rp, 'rms(div u)/rms(grad u) < 1e-4 at the finest step')
    ! The kernel is C1, so u'' jumps at the support edge: the FD error converges at ~h^1.5, not h^2.
    call check(ratio(2) < 0.6_rp * ratio(1) .and. ratio(3) < 0.6_rp * ratio(2), &
               'divergence decreases with h (discretisation error only)')
  end subroutine test_divergence_free

  subroutine test_bins_and_periodic()
    type(SEMState)  :: state, sp
    type(ABLConfig) :: cfg
    real(rp) :: x(3), ub(3), ur(3), r(3), emax, scale, W, xs(3), u1(3), u2(3), t_bin, t_ref, h, gsum
    integer  :: ipt, it, c0, c1, crate

    write(*,'(A)') 'spatial bins and spanwise periodicity'

    ! Bins must reproduce the full loop to rounding, also after eddies have moved and been recycled.
    call make_state(state, cfg, 20000)
    call sem_setup(state, cfg)
    do it = 1, 5
      call sem_advance(state, cfg, 0.3_rp)
    end do
    emax = 0.0_rp; scale = 0.0_rp
    do ipt = 1, 400
      call random_number(r)
      x = [-8.0_rp + 16.0_rp * r(1), -34.0_rp + 68.0_rp * r(2), -3.0_rp + 70.0_rp * r(3)]
      call sem_fluct_at_point(state, x, ub)
      call sem_fluct_at_point_ref(state, x, ur)
      emax = max(emax, maxval(abs(ub - ur)))
      scale = max(scale, maxval(abs(ur)))
    end do
    call check(emax < 1.0e-12_rp * max(scale, 1.0_rp), 'binned evaluation equals the reference loop (1e-12)')

    call system_clock(c0, crate)
    gsum = 0.0_rp
    do ipt = 1, 2000
      call sem_fluct_at_point(state, [0.0_rp, -20.0_rp + 0.02_rp * real(ipt, rp), 20.0_rp], ub)
      gsum = gsum + ub(1)
    end do
    call system_clock(c1); t_bin = real(c1 - c0, rp) / real(crate, rp)
    do ipt = 1, 2000
      call sem_fluct_at_point_ref(state, [0.0_rp, -20.0_rp + 0.02_rp * real(ipt, rp), 20.0_rp], ur)
      gsum = gsum - ur(1)
    end do
    if (abs(gsum) > 1.0e-6_rp) write(*,'(A)') '    (timing loops disagree)'
    call system_clock(c0); t_ref = real(c0 - c1, rp) / real(crate, rp)
    write(*,'(A,ES9.2,A,ES9.2,A)') '    speedup of the binned evaluator: ', &
      t_ref / max(t_bin, 1.0e-9_rp), 'x  (binned time for 2000 points, s): ', t_bin
    call check(t_ref > 5.0_rp * t_bin, 'binning is at least 5x faster at 20000 eddies')

    ! Spanwise-periodic field: exact translation invariance by W, also divergence-free across the seam.
    call make_state(sp, cfg, 20000)
    W = 60.0_rp
    sp%periodic_t1 = .true.; sp%t1_period = W; sp%t1_min = 0.0_rp; sp%t1_max = W
    call sem_setup(sp, cfg)
    call sem_advance(sp, cfg, 0.2_rp)
    emax = 0.0_rp; scale = 0.0_rp
    do ipt = 1, 200
      call random_number(r)
      x = [-5.0_rp + 10.0_rp * r(1), W * r(2), 2.0_rp + 55.0_rp * r(3)]
      xs = x + [0.0_rp, W, 0.0_rp]
      call sem_fluct_at_point(sp, x, u1)
      call sem_fluct_at_point(sp, xs, u2)
      emax = max(emax, maxval(abs(u1 - u2))); scale = max(scale, maxval(abs(u1)))
    end do
    call check(emax < 1.0e-12_rp * scale, 'u(y + W) = u(y) to rounding')
    ! Divergence by central differences straddling the seam y = 0.
    h = 0.01_rp * minval(sp%sig)
    gsum = 0.0_rp; scale = 0.0_rp
    do ipt = 1, 100
      call random_number(r)
      x = [-5.0_rp + 10.0_rp * r(1), -1.0_rp + 2.0_rp * r(2), 2.0_rp + 55.0_rp * r(3)]
      call div_and_grad(sp, x, h, emax, scale)
      gsum = gsum + emax**2
      u1(1) = u1(1) + scale**2
    end do
    call check(sqrt(gsum / max(u1(1), SMALL_T)) < 1.0e-4_rp, 'divergence-free across the periodic seam')
  end subroutine test_bins_and_periodic

  !> div and |grad u|^2 by central differences at x.
  subroutine div_and_grad(state, x, h, dv, gn)
    type(SEMState), intent(in)  :: state
    real(rp),       intent(in)  :: x(3), h
    real(rp),       intent(out) :: dv, gn
    real(rp) :: xx(3), up(3), um(3)
    integer :: k, j
    dv = 0.0_rp; gn = 0.0_rp
    do k = 1, 3
      xx = x; xx(k) = x(k) + h; call sem_fluct_at_point(state, xx, up)
      xx = x; xx(k) = x(k) - h; call sem_fluct_at_point(state, xx, um)
      dv = dv + (up(k) - um(k)) / (2.0_rp * h)
      do j = 1, 3
        gn = gn + ((up(j) - um(j)) / (2.0_rp * h))**2
      end do
    end do
  end subroutine div_and_grad

  !> rms of the central-difference divergence over rms of the velocity-gradient tensor, same points every call.
  real(rp) function div_over_grad(state, h) result(q)
    type(SEMState), intent(in) :: state
    real(rp),       intent(in) :: h
    real(rp) :: x(3), up(3), um(3), dv, gn, sd, sg
    integer  :: ipt, k, j
    real(rp) :: r(3)
    integer  :: seed_n
    integer, allocatable :: seed(:)

    call random_seed(size=seed_n)
    allocate(seed(seed_n)); seed = 12345
    call random_seed(put=seed)
    sd = 0.0_rp; sg = 0.0_rp
    do ipt = 1, 300
      call random_number(r)
      x = [-5.0_rp + 10.0_rp * r(1), -25.0_rp + 50.0_rp * r(2), 2.0_rp + 55.0_rp * r(3)]
      dv = 0.0_rp; gn = 0.0_rp
      do k = 1, 3
        x(k) = x(k) + h;          call sem_fluct_at_point(state, x, up)
        x(k) = x(k) - 2.0_rp * h; call sem_fluct_at_point(state, x, um)
        x(k) = x(k) + h
        dv = dv + (up(k) - um(k)) / (2.0_rp * h)
        do j = 1, 3
          gn = gn + ((up(j) - um(j)) / (2.0_rp * h))**2
        end do
      end do
      sd = sd + dv**2
      sg = sg + gn
    end do
    q = sqrt(sd / sg)
  end function div_over_grad

  subroutine test_stats_and_scales()
    type(SEMState)  :: state
    type(ABLConfig) :: cfg
    real(rp), parameter :: zl(4) = [3.0_rp, 10.0_rp, 25.0_rp, 50.0_rp]
    real(rp), parameter :: yl(5) = [-20.0_rp, -10.0_rp, 0.0_rp, 10.0_rp, 20.0_rp]
    integer,  parameter :: nstep = 4000, nlag = 400
    real(rp), parameter :: zl_s(3) = [3.0_rp, 10.0_rp, 50.0_rp]
    real(rp) :: s1(3,4), s2(6,4), x(3), u(3), um(3), R(6), L(3), cnt, dt
    real(rp) :: series(nstep, 3), rr(0:nlag), lux, rho(0:nlag), mean_s, var_s, sx, uc, ratio
    integer  :: it, iz, iy, lag

    write(*,'(A)') 'mean, intensities, shear stress and length scale vs the prescribed profile'
    call make_state(state, cfg, 20000)
    call sem_setup(state, cfg)
    write(*,'(A,ES9.2,A,ES9.2)') '    n_sigma = ', state%n_sigma, ', realisability error = ', state%realiz_err
    call check(state%realiz_err < 1.0e-9_rp, 'default targets are exactly realisable at every height')

    s1 = 0.0_rp; s2 = 0.0_rp; cnt = 0.0_rp
    dt = 0.3_rp
    do it = 1, 800
      call sem_advance(state, cfg, dt)
      cnt = cnt + 1.0_rp
      do iz = 1, 4
        do iy = 1, 5
          x = [0.0_rp, yl(iy), zl(iz)]
          call sem_fluct_at_point(state, x, u)
          s1(:, iz) = s1(:, iz) + u
          s2(1, iz) = s2(1, iz) + u(1)**2; s2(2, iz) = s2(2, iz) + u(2)**2; s2(3, iz) = s2(3, iz) + u(3)**2
          s2(4, iz) = s2(4, iz) + u(1)*u(2); s2(5, iz) = s2(5, iz) + u(1)*u(3); s2(6, iz) = s2(6, iz) + u(2)*u(3)
        end do
      end do
    end do
    cnt = cnt * 5.0_rp
    do iz = 1, 4
      call sem_target_at(state, cfg, zl(iz), um, R, L)
      s1(:, iz) = s1(:, iz) / cnt
      s2(:, iz) = s2(:, iz) / cnt
      write(*,'(A,F5.1,A,3F7.3,A,3F7.3,A,2F7.3)') '    z=', zl(iz), '  sigma/target (u,v,w) =', &
        sqrt(s2(1:3, iz) / R(1:3)), '  mean/sigma_u =', s1(:, iz) / sqrt(R(1)), '  uw, target =', s2(5, iz), R(5)
      call check(all(abs(sqrt(s2(1:3, iz) / R(1:3)) - 1.0_rp) < 0.12_rp), 'rms of u,v,w within 12% of target')
      call check(abs(s2(5, iz) - R(5)) < 0.15_rp * sqrt(R(1) * R(3)), 'u''w'' within 0.15 sigma_u sigma_w of target')
      call check(maxval(abs(s1(:, iz))) < 0.1_rp * sqrt(R(1)), 'fluctuation mean below 0.1 sigma_u')
    end do

    ! Streamwise integral scale from the time series at a point (frozen turbulence at U(z)),
    ! integrating the autocorrelation to its first zero crossing.
    do iz = 1, 3
      call sem_target_at(state, cfg, zl_s(iz), um, R, L)
      uc = um(1)
      dt = 0.02_rp
      do it = 1, nstep
        call sem_advance(state, cfg, dt)
        do iy = 1, 3
          call sem_fluct_at_point(state, [0.0_rp, -15.0_rp + 15.0_rp * real(iy - 1, rp), zl_s(iz)], u)
          series(it, iy) = u(1)
        end do
      end do
      rho = 0.0_rp
      do iy = 1, 3  ! average the autocorrelation over lateral positions to beat sampling noise
        series(:, iy) = series(:, iy) - sum(series(:, iy)) / real(nstep, rp)
        var_s = sum(series(:, iy)**2) / real(nstep, rp)
        do lag = 0, nlag
          rr(lag) = sum(series(1:nstep-lag, iy) * series(1+lag:nstep, iy)) / real(nstep - lag, rp) / var_s
        end do
        rho = rho + rr / 3.0_rp
      end do
      lux = 0.5_rp * rho(0)
      do lag = 1, nlag
        if (rho(lag) <= 0.0_rp) exit
        lux = lux + rho(lag)
      end do
      lux = lux * uc * dt
      ratio = lux / L(1)
      write(*,'(A,F5.1,A,F7.2,A,F7.2,A,F6.2)') '    z=', zl_s(iz), '  Lux measured =', lux, ' m, target =', L(1), &
        ' m, ratio', ratio
      call check(ratio > 0.6_rp .and. ratio < 1.7_rp, 'streamwise integral scale within ~40% of target')
    end do
    ! The scale must actually grow with height (the point of the power-law model).
    call sem_target_at(state, cfg, 3.0_rp, um, R, L);  sx = L(1)
    call sem_target_at(state, cfg, 50.0_rp, um, R, L)
    call check(L(1) > 2.0_rp * sx, 'prescribed Lx(z=50) > 2 Lx(z=3)')
  end subroutine test_stats_and_scales

  subroutine test_ti_table()
    type(SEMState)  :: state
    type(ABLConfig) :: cfg
    integer :: iu
    real(rp) :: um(3), R(6), L(3), su, sw
    character(len=*), parameter :: fn = 'sem_ti_test_profile.csv'

    write(*,'(A)') 'wind-tunnel table: z, U, Iu, Iv, Iw, rho_uw, Lx'
    open(newunit=iu, file=fn, status='replace', action='write')
    write(iu,'(A)') '# z  U   Iu   Iv   Iw   rho   Lx'
    write(iu,'(A)') '0.5  4.0 0.30 0.22 0.15 -0.30 2.0'
    write(iu,'(A)') '60.0 12.0 0.10 0.075 0.05 -0.30 20.0'
    close(iu)

    call make_state(state, cfg, 5000)
    cfg%sem_Ly_ratio = 0.5_rp
    call sem_read_profile(fn, state, fmt='ti', rho_uw=-0.3_rp)
    call sem_setup(state, cfg)
    call sem_target_at(state, cfg, 0.5_rp, um, R, L)
    su = 0.30_rp * 4.0_rp; sw = 0.15_rp * 4.0_rp
    call check(abs(R(1) - su**2) < 1.0e-12_rp .and. abs(R(3) - sw**2) < 1.0e-12_rp, 'sigma = I * U at the first level')
    call check(abs(R(5) - (-0.3_rp) * su * sw) < 1.0e-12_rp, 'u''w'' = rho sigma_u sigma_w')
    call check(abs(L(1) - 2.0_rp) < 1.0e-12_rp .and. abs(L(2) - 1.0_rp) < 1.0e-12_rp, &
               'Lx from the table, Ly from the ratio')
    call sem_target_at(state, cfg, 30.0_rp, um, R, L)
    call check(abs(um(1) - (4.0_rp + 8.0_rp * 29.5_rp / 59.5_rp)) < 1.0e-9_rp, &
               'mean interpolates the table (westerly: +x)')
    call check(state%override_mean, 'solver is told to take the inlet mean from the SEM')
    call execute_command_line('rm -f '//fn)
  end subroutine test_ti_table

  subroutine test_ti_law()
    type(SEMState)  :: state
    type(ABLConfig) :: cfg
    real(rp) :: um(3), R(6), L(3), Iu_ref, umag, sn, sd

    write(*,'(A)') 'namelist laws: power-law mean, Iu(z), Lx(z)'
    call make_state(state, cfg, 5000)
    cfg%sem_mean_law = 'power'; cfg%sem_alpha = 0.25_rp
    cfg%sem_Iu_ref = 0.15_rp; cfg%sem_rho_uw = -0.25_rp
    call sem_setup(state, cfg)
    call sem_target_at(state, cfg, 50.0_rp, um, R, L)
    call check(abs(um(1) - 10.0_rp) < 1.0e-9_rp, 'U(z_ref) = u_ref')
    call check(abs(sqrt(R(1)) / 10.0_rp - 0.15_rp) < 1.0e-9_rp, 'Iu(z_ref) = sem_Iu_ref')
    call check(abs(sqrt(R(2) / R(1)) - 0.75_rp) < 1.0e-9_rp .and. abs(sqrt(R(3) / R(1)) - 0.5_rp) < 1.0e-9_rp, &
               'Iv/Iu = 0.75 and Iw/Iu = 0.5')
    call check(abs(R(5) / sqrt(R(1) * R(3)) + 0.25_rp) < 1.0e-9_rp, 'rho_uw from the namelist')
    call check(state%override_mean, 'power-law mean is handed to the solver')
    call sem_target_at(state, cfg, 5.0_rp, um, R, L)
    Iu_ref = sqrt(R(1)) / um(1)
    call check(Iu_ref > 0.15_rp, 'intensity increases towards the wall')
    umag = um(1)
    sn = L(1); sd = 6.0_rp * (5.0_rp / 50.0_rp) ** 0.3_rp
    call check(abs(sn - sd) < 1.0e-9_rp, 'Lx(z) = L_ref (z/z_ref)^0.3')
  end subroutine test_ti_law

end program test_sem_inflow
