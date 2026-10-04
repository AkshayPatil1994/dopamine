! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> The successor-side replay of an inflow-plane database:
!>   1. schedule: deterministic, no immediate repeats, shifts well spread
!>   2. interpolation order: exact for band-limited (span) x cubic (z) x quadratic (time) fields
!>   3. mesh-matched faces are an exact copy (no filtering)
!>   4. cross-fade keeps the variance; time interpolation loses < 1%
!>   5. no spectral line at the replay period with several segments (and a line with one)
!>   6. the lateral-periodicity stop check
program test_inflow_planes
  use mod_kinds
  use mod_parameters, only: BC_PERIODIC, BC_SLIP
  use mod_abl_types
  use mod_mesh
  use mod_inflow_db
  use mod_inflow_gen
  use mod_inflow_planes
  use mpi_f08
  implicit none

  integer :: n_fail = 0
  integer :: ierr
  real(rp), parameter :: PI = 3.14159265358979323846_rp

  call MPI_Init(ierr)
  call test_schedule()
  call test_interpolation_orders()
  call test_crossfade_and_time_interp()
  call test_spectral_line()
  call test_periodic_check()
  call MPI_Finalize(ierr)

  if (n_fail > 0) then
    write(*,'(A,I0,A)') 'test_inflow_planes: ', n_fail, ' FAILED'
    error stop 1
  end if
  write(*,'(A)') 'test_inflow_planes: all passed'

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

  !> Small SEM database on a 16 x 16 plane with `nseg` segments.
  subroutine make_sem_db(fn, nseg, turn, decorrelate)
    character(len=*), intent(in) :: fn
    integer(ip),      intent(in) :: nseg
    real(rp),         intent(in) :: turn
    logical,          intent(in), optional :: decorrelate
    type(ABLConfig) :: cfg
    type(InflowGenParams) :: gp
    logical :: ok
    cfg%u_ref = 1.0_rp; cfg%z_ref = 1.0_rp; cfg%z0_m = 0.001_rp; cfg%wind_dir_deg = 270.0_rp
    cfg%sem_sigma = 0.15_rp
    cfg%sem_mean_law = 'power'; cfg%sem_alpha = 0.15_rp; cfg%sem_Iu_ref = 0.1_rp
    gp%file = fn; gp%W = 1.0_rp; gp%y0 = 0.03125_rp; gp%ny = 16; gp%nz = 16
    gp%z_min = 0.03125_rp; gp%z_max = 0.96875_rp; gp%n_segments = nseg; gp%n_turnovers = turn
    gp%verify = .false.; gp%seed = 3
    if (present(decorrelate)) gp%decorrelate = decorrelate
    call inflow_gen_generate(cfg, gp, ok)
  end subroutine make_sem_db

  subroutine open_planes(fp, fn, seed, xcf)
    type(InflowPlanes), intent(inout) :: fp
    character(len=*),   intent(in)    :: fn
    integer(ip),        intent(in)    :: seed
    real(rp),           intent(in), optional :: xcf
    type(InflowConfig) :: ic
    type(ABLConfig) :: cfg
    ic%file = fn; ic%seed = seed
    if (present(xcf)) ic%crossfade_time = xcf
    cfg%wind_dir_deg = 270.0_rp
    call inflow_planes_open(fp, ic, cfg, 1_ip)
  end subroutine open_planes

  ! ----------------------------------------------------------------------------------------------
  subroutine test_schedule()
    type(InflowPlanes) :: a, b
    integer(ip) :: m, sa, sb, ha, hb, hprev, cnt(4), prev, d, mind, shifts_hist(8), bad_rep
    logical :: same

    write(*,'(A)') 'replay schedule'
    call make_sem_db('sched.nc', 4_ip, 1.0_rp)
    call open_planes(a, 'sched.nc', 11_ip)
    call open_planes(b, 'sched.nc', 11_ip)
    same = .true.; cnt = 0; prev = -1; bad_rep = 0; mind = 1000; shifts_hist = 0
    do m = 0, 1999
      call inflow_stream(a, m, sa, ha)
      call inflow_stream(b, m, sb, hb)
      same = same .and. (sa == sb) .and. (ha == hb)
      cnt(sa) = cnt(sa) + 1
      if (sa == prev) bad_rep = bad_rep + 1
      if (m > 0) then
        d = abs(ha - hprev); d = min(d, 16 - d); mind = min(mind, d)
      end if
      shifts_hist(ha / 2 + 1) = shifts_hist(ha / 2 + 1) + 1
      prev = sa; hprev = ha
    end do
    call check(same, 'same seed gives the same schedule')
    call check(bad_rep == 0, 'a segment is never replayed twice in a row')
    call check(minval(cnt) > 400 .and. maxval(cnt) < 600, 'all four segments are used about equally (2000 periods)')
    call check(mind >= 5, 'consecutive spanwise shifts differ by >= 5 of 16 nodes (golden-ratio spacing)')
    call check(minval(shifts_hist) > 150, 'shifts cover the whole period uniformly')
    call inflow_planes_close(a); call inflow_planes_close(b)
    call execute_command_line('rm -f sched.nc')
  end subroutine test_schedule

  ! ----------------------------------------------------------------------------------------------
  !> Synthetic field: sin(2 pi k (y-y0)/W) in span, cubic in z, quadratic in time.
  real(rp) function fexact(c, y, z, tau, W, y0)
    integer(ip), intent(in) :: c
    real(rp),    intent(in) :: y, z, tau, W, y0
    fexact = (0.5_rp + 0.7_rp * real(c, rp)) * sin(2.0_rp * PI * 2.0_rp * (y - y0) / W + 0.3_rp * c) &
           * (1.0_rp + 0.8_rp * z - 1.3_rp * z**2 + 0.6_rp * z**3) * (1.0_rp + 0.4_rp * tau + 0.1_rp * tau**2)
  end function fexact

  subroutine test_interpolation_orders()
    type(InflowDB) :: db
    type(InflowDBHeader) :: h
    type(InflowPlanes) :: fp
    real(rp), allocatable :: plane(:,:,:), xf(:,:), fl(:,:), zg(:)
    integer(ip), allocatable :: fb(:)
    integer(ip) :: is, it, j, k, c, nf, i, seg, shift
    real(rp) :: W, y0, dy, err, scale, t, tau, ex, yq, zq
    real(rp), parameter :: dt = 0.1_rp
    integer(ip), parameter :: ny = 12, nz = 9, nt = 40, nseg = 2

    write(*,'(A)') 'interpolation: exact for band-limited span, cubic z, quadratic time'
    W = 2.0_rp; y0 = 0.1_rp; dy = W / ny
    h%ny = ny; h%nz = nz; h%nt = nt; h%nseg = nseg; h%W = W; h%y0 = y0; h%dt_p = dt
    h%frame = reshape([1.0_rp,0.0_rp,0.0_rp, 0.0_rp,1.0_rp,0.0_rp, 0.0_rp,0.0_rp,1.0_rp], [3,3])
    h%n_in = [1.0_rp, 0.0_rp, 0.0_rp]; h%t1 = [0.0_rp, 1.0_rp, 0.0_rp]; h%wind_dir_deg = 270.0_rp
    allocate(h%z(nz), zg(nz))
    do k = 1, nz
      h%z(k) = 0.05_rp + 0.9_rp * (real(k - 1, rp) / (nz - 1))**1.5_rp   ! graded
    end do
    zg = h%z
    allocate(plane(ny, nz, 3))
    call inflow_db_create(db, 'synth.nc', h)
    do is = 1, nseg
      do it = 1, nt
        do c = 1, 3
          do k = 1, nz
            do j = 1, ny
              plane(j, k, c) = fexact(c, y0 + (j - 1) * dy, zg(k), (it - 1) * dt, W, y0)
            end do
          end do
        end do
        call inflow_db_write_frame(db, is, it, plane)
      end do
    end do
    call inflow_db_close(db)

    call open_planes(fp, 'synth.nc', 5_ip)
    ! faces at arbitrary positions (off the nodes in y and z)
    nf = 40
    allocate(xf(3, nf), fb(nf), fl(3, nf))
    do i = 1, nf
      xf(:, i) = [0.0_rp, y0 + W * (0.013_rp + 0.0237_rp * i) + 0.37_rp * dy, 0.08_rp + 0.83_rp * real(i, rp) / nf]
      fb(i) = i
    end do
    call inflow_planes_bind(fp, xf, fb, [1.0_rp, 0.0_rp, 0.0_rp], [0.0_rp, 1.0_rp, 0.0_rp], W, 1_ip)
    call inflow_stream(fp, 0_ip, seg, shift)
    err = 0.0_rp; scale = 0.0_rp
    do it = 1, 6
      t = 0.45_rp + 0.31_rp * it       ! inside period 0: no cross-fade, tau = t
      call inflow_planes_eval(fp, t, fl)
      tau = t
      do i = 1, nf
        yq = xf(2, i) + shift * dy
        zq = xf(3, i)
        do c = 1, 3
          ex = fexact(c, yq, zq, tau, W, y0)
          err = max(err, abs(fl(c, i) - ex)); scale = max(scale, abs(ex))
        end do
      end do
    end do
    write(*,'(A,ES10.2)') '    max error relative to the field amplitude: ', err / scale
    call check(err < 2.0e-5_rp * scale, 'off-node faces: error at float32 level (span exact, z cubic exact, time quadratic exact)')

    ! faces exactly on the nodes: an exact copy of the stored data
    call inflow_planes_close(fp)
    call open_planes(fp, 'synth.nc', 5_ip)
    nf = ny * nz
    deallocate(xf, fb, fl); allocate(xf(3, nf), fb(nf), fl(3, nf))
    i = 0
    do k = 1, nz
      do j = 1, ny
        i = i + 1
        xf(:, i) = [0.0_rp, y0 + (j - 1) * dy, zg(k)]
        fb(i) = i
      end do
    end do
    call inflow_planes_bind(fp, xf, fb, [1.0_rp, 0.0_rp, 0.0_rp], [0.0_rp, 1.0_rp, 0.0_rp], W, 1_ip)
    call check(fp%n_match_y == ny .and. fp%n_match_z == nz, 'every face recognised as coinciding with a database node')
    call inflow_stream(fp, 0_ip, seg, shift)
    t = 0.3_rp * dt * 7      ! between frames: time interpolation only
    call inflow_planes_eval(fp, t, fl)
    err = 0.0_rp
    i = 0
    do k = 1, nz
      do j = 1, ny
        i = i + 1
        do c = 1, 3
          ex = fexact(c, y0 + (j - 1 + shift) * dy, zg(k), t, W, y0)
          err = max(err, abs(fl(c, i) - ex))
        end do
      end do
    end do
    call check(err < 2.0e-5_rp * scale, 'matched faces: no spatial filtering, only the (exact) time interpolation')
    call inflow_planes_close(fp)
    call execute_command_line('rm -f synth.nc')
  end subroutine test_interpolation_orders

  ! ----------------------------------------------------------------------------------------------
  subroutine test_crossfade_and_time_interp()
    type(InflowPlanes) :: fp
    type(InflowDB) :: db
    real(rp), allocatable :: xf(:,:), fl(:,:), pl(:,:,:), fr(:,:,:)
    integer(ip), allocatable :: fb(:)
    integer(ip) :: nf, i, j, k, it, nper, n_in_fade, n_out
    real(rp) :: v_fade, v_out, t, P, Tx, vmid, vnode, tau

    write(*,'(A)') 'cross-fade and time interpolation'
    call make_sem_db('xf.nc', 4_ip, 2.0_rp)
    call open_planes(fp, 'xf.nc', 7_ip)
    P = fp%P; Tx = fp%T_x
    nf = 16 * 8
    allocate(xf(3, nf), fb(nf), fl(3, nf))
    call inflow_db_open(db, 'xf.nc')
    i = 0
    do k = 2, 16, 2
      do j = 1, 16
        i = i + 1
        xf(:, i) = [0.0_rp, 0.03125_rp + (j - 1) / 16.0_rp, db%h%z(k)]
        fb(i) = i
      end do
    end do
    call inflow_db_close(db)
    call inflow_planes_bind(fp, xf, fb, [1.0_rp, 0.0_rp, 0.0_rp], [0.0_rp, 1.0_rp, 0.0_rp], 1.0_rp, 1_ip)

    v_fade = 0.0_rp; v_out = 0.0_rp; n_in_fade = 0; n_out = 0
    do nper = 1, 60
      do it = 0, 39
        t = (nper + (it + 0.5_rp) / 40.0_rp) * P
        call inflow_planes_eval(fp, t, fl)
        tau = t - nper * P
        if (tau < Tx) then
          v_fade = v_fade + sum(fl(1, :)**2) / nf; n_in_fade = n_in_fade + 1
        else
          v_out = v_out + sum(fl(1, :)**2) / nf; n_out = n_out + 1
        end if
      end do
    end do
    v_fade = v_fade / max(n_in_fade, 1); v_out = v_out / max(n_out, 1)
    write(*,'(A,F6.3)') '    variance inside cross-fades / outside: ', v_fade / v_out
    call check(abs(v_fade / v_out - 1.0_rp) < 0.15_rp, 'cross-fade preserves the variance (within 15%, sampling noise)')

    ! time interpolation: mid-frame values must retain the variance of the frames themselves
    allocate(pl(16, 16, 3), fr(16, 16, 3))
    vmid = 0.0_rp; vnode = 0.0_rp
    do it = 3, fp%db%h%nt - 3
      call inflow_planes_time_field(fp, 2_ip, (it - 1 + 0.5_rp) * fp%db%h%dt_p, pl)
      call inflow_db_read_frame(fp%db, 2_ip, it, fr)
      vmid = vmid + sum(pl(:, :, 1)**2); vnode = vnode + sum(fr(:, :, 1)**2)
    end do
    write(*,'(A,F7.4)') '    variance at half-frame times / at frames: ', vmid / vnode
    call check(vmid / vnode > 0.99_rp, 'cubic time interpolation retains > 99% of the variance at half-frame times')
    call inflow_planes_close(fp)
    call execute_command_line('rm -f xf.nc')
  end subroutine test_crossfade_and_time_interp

  ! ----------------------------------------------------------------------------------------------
  !> Power of the spanwise-mean u(t) at the replay frequency relative to its neighbours.
  real(rp) function line_ratio(fn, seed) result(r)
    character(len=*), intent(in) :: fn
    integer(ip),      intent(in) :: seed
    type(InflowPlanes) :: fp
    type(InflowDB) :: db
    real(rp), allocatable :: xf(:,:), fl(:,:), series(:), pw(:)
    integer(ip), allocatable :: fb(:)
    integer(ip) :: j, n, Np, S, kf, kk, nn
    real(rp) :: t, re, im, P, neigh

    call open_planes(fp, fn, seed)
    P = fp%P
    call inflow_db_open(db, fn)
    allocate(xf(3, 16), fb(16), fl(3, 16))
    do j = 1, 16
      xf(:, j) = [0.0_rp, 0.03125_rp + (j - 1) / 16.0_rp, db%h%z(8)]; fb(j) = j
    end do
    call inflow_db_close(db)
    call inflow_planes_bind(fp, xf, fb, [1.0_rp, 0.0_rp, 0.0_rp], [0.0_rp, 1.0_rp, 0.0_rp], 1.0_rp, 1_ip)
    Np = 120; S = 32; n = Np * S
    allocate(series(n), pw(n / 2))
    do nn = 1, n
      t = (nn - 1) * P / S
      call inflow_planes_eval(fp, t, fl)
      series(nn) = sum(fl(1, :)) / 16.0_rp
    end do
    series = series - sum(series) / n
    do kk = 1, n / 2
      re = 0.0_rp; im = 0.0_rp
      do nn = 1, n
        re = re + series(nn) * cos(2.0_rp * PI * kk * (nn - 1) / n)
        im = im + series(nn) * sin(2.0_rp * PI * kk * (nn - 1) / n)
      end do
      pw(kk) = re**2 + im**2
    end do
    kf = Np     ! frequency 1/P
    neigh = 0.0_rp
    do kk = kf - 40, kf + 40
      if (abs(kk - kf) >= 4) neigh = neigh + pw(kk)
    end do
    neigh = neigh / count([(abs(kk - kf) >= 4, kk = kf - 40, kf + 40)])
    r = pw(kf) / neigh
    call inflow_planes_close(fp)
  end function line_ratio

  !> Variance of the spanwise-mean u over all frames and segments at mid-height.
  real(rp) function spanmean_var(fn) result(v)
    character(len=*), intent(in) :: fn
    type(InflowDB) :: db
    real(rp), allocatable :: pl(:,:,:)
    integer(ip) :: is, it
    real(rp) :: s1, s2, y
    call inflow_db_open(db, fn)
    allocate(pl(db%h%ny, db%h%nz, 3))
    s1 = 0.0_rp; s2 = 0.0_rp
    do is = 1, db%h%nseg
      do it = 1, db%h%nt
        call inflow_db_read_frame(db, is, it, pl)
        y = sum(pl(:, 8, 1)) / db%h%ny
        s1 = s1 + y; s2 = s2 + y**2
      end do
    end do
    v = s2 / (db%h%nseg * db%h%nt) - (s1 / (db%h%nseg * db%h%nt))**2
    call inflow_db_close(db)
  end function spanmean_var

  subroutine test_spectral_line()
    real(rp) :: r1, r4, r4d, frac, vplain, vdec
    write(*,'(A)') 'periodic-injection line in the spanwise-mean signal (the part shifts cannot move)'
    call make_sem_db('line1.nc', 1_ip, 2.0_rp, .false.)
    call make_sem_db('line4.nc', 4_ip, 2.0_rp, .false.)
    r1 = line_ratio('line1.nc', 3_ip)
    r4 = line_ratio('line4.nc', 3_ip)
    vplain = spanmean_var('line4.nc')
    call inflow_db_decorrelate_span_mean('line4.nc', frac)
    r4d = line_ratio('line4.nc', 3_ip)
    vdec = spanmean_var('line4.nc')
    write(*,'(A,ES9.2)') '    line / neighbours, one segment looped:                ', r1
    write(*,'(A,ES9.2)') '    line / neighbours, four random segments:              ', r4
    write(*,'(A,ES9.2)') '    line / neighbours, four segments, span-mean decorrelated: ', r4d
    write(*,'(A,F7.3,A,F7.3)') '    spanwise-mean variance before / after decorrelation: ', vplain, ' /', vdec
    call check(r1 > 20.0_rp, 'detector sees the line when a single segment is looped')
    call check(r4 > 20.0_rp, 'random segment order alone does NOT remove it (line carries ~1/K of the spanwise-mean power)')
    call check(r4d < 10.0_rp, 'ensemble-mean removal leaves no line at 1/P')
    call check(abs(vdec / vplain - 1.0_rp) < 0.2_rp, 'spanwise-mean variance is preserved (within sampling noise)')
    call check(frac < 0.2_rp, 'the spanwise-mean part is a small share of the plane variance')
    call execute_command_line('rm -f line1.nc line4.nc')
  end subroutine test_spectral_line

  ! ----------------------------------------------------------------------------------------------
  subroutine test_periodic_check()
    type(Mesh) :: m
    real(rp) :: W
    logical :: ok

    write(*,'(A)') 'lateral-periodicity stop check'
    m%comm = MPI_COMM_WORLD
    m%n_int_faces = 0
    allocate(m%faces(2), m%patches(1))
    m%faces(1)%normal = [0.0_rp, -1.0_rp, 0.0_rp]; m%faces(1)%centre = [0.5_rp, 0.0_rp, 0.5_rp]; m%faces(1)%area = 1.0_rp
    m%faces(2)%normal = [0.0_rp,  1.0_rp, 0.0_rp]; m%faces(2)%centre = [0.5_rp, 2.5_rp, 0.5_rp]; m%faces(2)%area = 1.0_rp
    m%n_patches = 1
    m%patches(1)%bc_id = BC_PERIODIC; m%patches(1)%start_face = 1; m%patches(1)%n_faces = 2
    call inflow_check_lateral_periodic(m, [0.0_rp, 1.0_rp, 0.0_rp], W, ok)
    call check(ok .and. abs(W - 2.5_rp) < 1.0e-12_rp, 'periodic y planes found, width = 2.5')
    m%patches(1)%bc_id = BC_SLIP
    call inflow_check_lateral_periodic(m, [0.0_rp, 1.0_rp, 0.0_rp], W, ok)
    call check(.not. ok, 'a domain with no periodic y planes is reported as not periodic')
    m%patches(1)%bc_id = BC_PERIODIC
    m%faces(2)%area = 0.5_rp
    call inflow_check_lateral_periodic(m, [0.0_rp, 1.0_rp, 0.0_rp], W, ok)
    call check(.not. ok, 'unequal periodic areas (unpaired planes) are rejected')
  end subroutine test_periodic_check

end program test_inflow_planes
