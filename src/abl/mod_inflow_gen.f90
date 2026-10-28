! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Offline generation of an inflow-plane database from the SEM.
module mod_inflow_gen
  use mod_kinds
  use mod_abl_types
  use mod_abl_profiles, only: abl_wind_vector
  use mod_sem_inflow
  use mod_inflow_db
  implicit none
  private

  real(rp), parameter :: PI = 3.14159265358979323846_rp
  real(rp), parameter :: SMALL_GEN = 1.0e-16_rp

  public :: inflow_gen_generate, inflow_gen_verify
  public :: inflow_write_grid_file, inflow_read_grid_file

contains

  pure function cross3(a, b) result(c)
    real(rp), intent(in) :: a(3), b(3)
    real(rp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  !> Plane grid description shared between the solver (export) and the generator (import), so a
  !> database can be generated on exactly the successor's inlet face positions.
  subroutine inflow_write_grid_file(filename, W, y0, ny, n_in, z)
    character(len=*), intent(in) :: filename
    real(rp),         intent(in) :: W, y0, n_in(3), z(:)
    integer(ip),      intent(in) :: ny

    integer :: iu, k
    open(newunit=iu, file=trim(filename), status='replace', action='write')
    write(iu,'(A)') '# DOPAMINE inflow plane grid (written by the solver for dopamine-inflow-gen)'
    write(iu,'(A,ES24.16)') 'W  ', W
    write(iu,'(A,ES24.16)') 'y0 ', y0
    write(iu,'(A,I0)')      'ny ', ny
    write(iu,'(A,3ES24.16)') 'n_in', n_in
    write(iu,'(A)') 'z'
    do k = 1, size(z)
      write(iu,'(ES24.16)') z(k)
    end do
    close(iu)
  end subroutine inflow_write_grid_file

  subroutine inflow_read_grid_file(filename, W, y0, ny, n_in, z)
    character(len=*),      intent(in)  :: filename
    real(rp),              intent(out) :: W, y0, n_in(3)
    integer(ip),           intent(out) :: ny
    real(rp), allocatable, intent(out) :: z(:)

    integer :: iu, ios, nz
    character(len=256) :: line, key
    real(rp), allocatable :: tmp(:)
    real(rp) :: v
    logical :: in_z

    W = 0.0_rp; y0 = 0.0_rp; ny = 0; n_in = 0.0_rp
    allocate(tmp(0))
    nz = 0
    in_z = .false.
    open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'inflow grid file: cannot open '//trim(filename)
    do
      read(iu, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      if (in_z) then
        read(line, *, iostat=ios) v
        if (ios == 0) then
          tmp = [tmp, v]; nz = nz + 1
        end if
        cycle
      end if
      read(line, *, iostat=ios) key
      select case (trim(key))
      case ('W');    read(line(len_trim(key)+1:), *) W
      case ('y0');   read(line(len_trim(key)+1:), *) y0
      case ('ny');   read(line(len_trim(key)+1:), *) ny
      case ('n_in'); read(line(len_trim(key)+1:), *) n_in
      case ('z');    in_z = .true.
      end select
    end do
    close(iu)
    if (W <= 0.0_rp .or. ny < 2 .or. nz < 2) error stop 'inflow grid file: needs W, ny and at least two z levels'
    allocate(z(nz)); z = tmp
  end subroutine inflow_read_grid_file

  !> Generate the database described by `gp` for the inflow described by the &dopamine_abl block.
  subroutine inflow_gen_generate(abl_cfg, gp, ok)
    type(ABLConfig),       intent(in)  :: abl_cfg
    type(InflowGenParams), intent(in)  :: gp
    logical,               intent(out), optional :: ok

    type(SEMState)  :: base, pilot, st
    type(ABLConfig) :: cfg
    type(InflowDB)  :: db
    type(InflowDBHeader) :: h
    real(rp), allocatable :: z(:), plane(:,:,:), gz(:)
    real(rp) :: W, y0, n_in(3), t1(3), um(3), R(6), Te, tseg, dt, zk, x(3), ug(3), z_ref, Rm(3,3), Wr(3,3)
    integer(ip) :: ny, nz, nt, k, iseg, it, j, n_eddies, n_pilot
    logical :: verified
    real(rp) :: spanvar

    if (present(ok)) ok = .true.
    cfg = abl_cfg

    ! ---- grid ----
    if (len_trim(gp%grid_file) > 0) then
      call inflow_read_grid_file(gp%grid_file, W, y0, ny, n_in, gz)
      nz = size(gz); allocate(z(nz)); z = gz
      n_in = n_in / max(sqrt(sum(n_in**2)), 1.0e-300_rp)
      if (maxval(abs(n_in - abs_wind(cfg))) > 1.0e-3_rp) then
        write(*,'(A)') 'ERROR: the grid file inlet normal differs from the wind direction in &dopamine_abl.'
        error stop 'inflow_gen: wind_dir_deg is not aligned with the inlet normal of the grid file'
      end if
    else
      W = gp%W; y0 = gp%y0; ny = gp%ny; nz = gp%nz
      if (W <= 0.0_rp .or. ny < 2 .or. nz < 2 .or. gp%z_max <= gp%z_min) &
        error stop 'inflow_gen: set gen_W/gen_ny/gen_nz/gen_z_* or gen_grid_file'
      allocate(z(nz))
      do k = 1, nz
        z(k) = gp%z_min + (gp%z_max - gp%z_min) * real(k - 1, rp) / real(nz - 1, rp)
      end do
      n_in = abs_wind(cfg)
    end if
    if (gp%n_segments < 1 .or. gp%n_turnovers <= 0.0_rp .or. gp%n_sigma <= 0.0_rp) &
      error stop 'inflow_gen: need gen_n_segments >= 1, gen_n_turnovers > 0, gen_n_sigma > 0'
    if (abs(n_in(3)) > 1.0e-6_rp) error stop 'inflow_gen: the flow direction through the plane must be horizontal'
    t1 = cross3([0.0_rp, 0.0_rp, 1.0_rp], n_in)
    z_ref = cfg%z_ref

    ! ---- template: geometry plus (optional) profile table ----
    if (len_trim(cfg%sem_profile_file) > 0) &
      call sem_read_profile(cfg%sem_profile_file, base, verbose=.true., fmt=cfg%sem_profile_format, rho_uw=cfg%sem_rho_uw)
    base%x_inlet = 0.0_rp
    base%n_in = n_in
    base%t1 = t1
    base%t2 = cross3(n_in, t1)
    base%per1 = .true.; base%period1 = W; base%per1_lo = y0
    base%per2 = .false.; base%period2 = 0.0_rp
    base%t1_min = y0; base%t1_max = y0 + W
    base%t2_min = z(1); base%t2_max = z(nz)

    ! ---- pilot setup: eddy scale, box, turnover time, frame spacing ----
    pilot = base
    n_pilot = 400
    cfg%n_eddies = n_pilot
    cfg%sem_seed = 1
    call sem_setup(pilot, cfg, verbose=.false.)
    call sem_target_at(pilot, cfg, z_ref, um, R)
    Te = gp%T_turnover
    if (Te <= 0.0_rp) then
      if (R(1) <= 0.0_rp) error stop 'inflow_gen: zero streamwise stress at z_ref, set gen_T_turnover'
      Te = pilot%sigma / sqrt(R(1))
    end if
    tseg = gp%n_turnovers * Te

    ! Eddies are convected rigidly at u_convect: at least ~8 frames across the shortest eddy passage.
    dt = gp%dt
    if (dt <= 0.0_rp) dt = max(pilot%sig_min, SMALL_GEN) / max(pilot%u_convect, SMALL_GEN) / 4.0_rp
    nt = max(8_ip, int(ceiling(tseg / dt), ip) + 1_ip)
    ! gen_n_sigma eddies per sigma^3 of box volume.
    n_eddies = max(n_pilot, int(ceiling(gp%n_sigma * pilot%box_volume / pilot%sigma**3), ip))

    write(*,'(A,ES10.3,A,ES10.3,A,I0,A,ES9.2,A)') ' inflow-gen: turnover time ', Te, ' s, segment ', &
      (nt - 1) * dt, ' s (', nt, ' frames of dt = ', dt, ' s)'
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,F7.1,A)') ' inflow-gen: ', ny, ' x ', nz, ' plane, ', gp%n_segments, &
      ' segments, ', n_eddies, ' eddies, ', &
      real(ny, rp) * nz * 3 * 4 * nt * gp%n_segments / 1.0e6_rp, ' MB'

    ! How many plane nodes span the smallest eddy (diameter 2 sigma)? Only matters if the database is
    ! later interpolated to a different mesh; on the mesh's own lattice nothing is interpolated.
    call report_eddy_resolution(pilot, z, W / real(ny, rp))

    ! ---- header ----
    h%ny = ny; h%nz = nz; h%nt = nt; h%nseg = gp%n_segments
    h%W = W; h%y0 = y0; h%dt_p = dt
    call sem_stress_at(pilot, cfg, z_ref, Rm, Wr)
    h%frame = Wr
    h%n_in = n_in; h%t1 = t1
    h%wind_dir_deg = cfg%wind_dir_deg
    h%u_ref = cfg%u_ref; h%z_ref = z_ref; h%T_turnover = Te
    h%seed = gp%seed
    h%generator = 'dopamine-inflow-gen sem'
    allocate(h%z(nz), h%U_mean(nz), h%R_target(6, nz))
    h%z = z
    do k = 1, nz
      call sem_target_at(pilot, cfg, z(k), um, R)
      h%U_mean(k) = sqrt(sum(um**2))
      h%R_target(:, k) = R
    end do
    call inflow_db_create(db, trim(gp%file), h)

    ! ---- segments ----
    allocate(plane(ny, nz, 3))
    cfg%n_eddies = n_eddies
    do iseg = 1, gp%n_segments
      st = base
      cfg%sem_seed = 7919 * gp%seed + 104729 * int(iseg, ip)
      call sem_setup(st, cfg, verbose=(iseg == 1))
      do it = 1, nt
        if (it > 1) then
          call sem_advance(st, dt)
          call sem_prepare(st, cfg)  ! re-seeded eddies sit at new heights
        end if
        !$omp parallel do collapse(2) private(j, k, x, ug, um, zk) schedule(static)
        do k = 1, nz
          do j = 1, ny
            zk = z(k)
            x = (y0 + real(j - 1, rp) * W / real(ny, rp)) * t1 + zk * st%t2
            call sem_velocity_at_face(st, cfg, x, um, ug)
            plane(j, k, 1) = dot_product(Wr(:, 1), ug)
            plane(j, k, 2) = dot_product(Wr(:, 2), ug)
            plane(j, k, 3) = dot_product(Wr(:, 3), ug)
          end do
        end do
        !$omp end parallel do
        call inflow_db_write_frame(db, iseg, it, plane)
      end do
      write(*,'(A,I0,A,I0)') ' inflow-gen: segment ', iseg, ' of ', gp%n_segments
    end do
    call inflow_db_close(db)
    write(*,'(A)') ' inflow-gen: wrote '//trim(gp%file)
    if (gp%decorrelate) then
      call inflow_db_decorrelate_span_mean(trim(gp%file), spanvar)
      if (gp%n_segments >= 2) write(*,'(A,F6.2,A)') ' inflow-gen: removed the segment-mean spanwise-mean waveform (', &
        100.0_rp * spanvar, '% of the plane variance), so replay has no line at 1/P'
      if (gp%n_segments < 4) write(*,'(A)') ' WARNING: fewer than 4 segments; with so few the spanwise-mean content ' // &
        'of replayed segments is strongly repeated. Use gen_n_segments >= 4.'
    end if

    if (gp%verify) then
      call inflow_gen_verify(trim(gp%file), 0.15_rp, verified)
      if (present(ok)) ok = verified
    end if
  end subroutine inflow_gen_generate

  subroutine report_eddy_resolution(st, z, dy)
    type(SEMState),  intent(in) :: st
    real(rp),        intent(in) :: z(:), dy
    integer(ip) :: k, nz
    real(rp) :: dz, dzmin, ny_min, nz_min
    nz = size(z)
    dzmin = huge(1.0_rp)
    do k = 1, nz - 1
      dzmin = min(dzmin, z(k + 1) - z(k))
    end do
    dz = max(dzmin, SMALL_GEN)
    ! The smallest eddy has diameter 2 sig_min.
    ny_min = 2.0_rp * st%sig_min / dy
    nz_min = 2.0_rp * st%sig_min / dz
    write(*,'(A,F6.1,A,F6.1,A)') ' inflow-gen: smallest eddy spans ', ny_min, ' nodes across the span and ', nz_min, &
      ' in height'
    if (min(ny_min, nz_min) < 4.0_rp) write(*,'(A)') &
      ' inflow-gen: NOTE: fewer than 4 nodes per eddy; fine if the successor mesh uses this very lattice ' // &
      '(gen_grid_file), but interpolating to another mesh would filter the small scales.'
  end subroutine report_eddy_resolution

  !> Flow direction unit vector of the wind (horizontal), the direction the eddies are convected.
  function abs_wind(cfg) result(n)
    type(ABLConfig), intent(in) :: cfg
    real(rp) :: n(3)
    n = abl_wind_vector(1.0_rp, cfg%wind_dir_deg)
  end function abs_wind

  !> Re-read a database and compare plane statistics (over y, time and segments) with the
  !> targets stored in its header: rms of each component and the shear stress.
  subroutine inflow_gen_verify(filename, tol, ok)
    character(len=*), intent(in)  :: filename
    real(rp),         intent(in)  :: tol
    logical,          intent(out) :: ok

    type(InflowDB) :: db
    real(rp), allocatable :: plane(:,:,:), s1(:,:), s2(:,:)
    real(rp) :: n, rms(3), mean(3), uw, tu, tw, worst
    integer(ip) :: is, it, k, j, c

    call inflow_db_open(db, filename)
    allocate(plane(db%h%ny, db%h%nz, 3), s1(3, db%h%nz), s2(4, db%h%nz))
    s1 = 0.0_rp; s2 = 0.0_rp
    do is = 1, db%h%nseg
      do it = 1, db%h%nt
        call inflow_db_read_frame(db, is, it, plane)
        do k = 1, db%h%nz
          do j = 1, db%h%ny
            do c = 1, 3
              s1(c, k) = s1(c, k) + plane(j, k, c)
              s2(c, k) = s2(c, k) + plane(j, k, c)**2
            end do
            s2(4, k) = s2(4, k) + plane(j, k, 1) * plane(j, k, 3)
          end do
        end do
      end do
    end do
    n = real(db%h%ny, rp) * db%h%nt * db%h%nseg
    ok = .true.
    worst = 0.0_rp
    write(*,'(A)') ' inflow-gen verify:     z   rms(u,v,w)/target            u''w''/(su sw target)   mean/su'
    do k = 1, db%h%nz
      mean = s1(:, k) / n
      rms  = sqrt(max(s2(1:3, k) / n - mean**2, 0.0_rp))
      uw   = s2(4, k) / n - mean(1) * mean(3)
      tu = sqrt(max(db%h%R_target(1, k), 0.0_rp)); tw = sqrt(max(db%h%R_target(3, k), 0.0_rp))
      if (tu <= 1.0e-12_rp .or. tw <= 1.0e-12_rp) cycle
      write(*,'(A,F9.3,3F8.3,A,2F9.3,A,F7.3)') '                      ', db%h%z(k), &
        rms(1) / tu, rms(2) / sqrt(db%h%R_target(2, k)), rms(3) / tw, '   ', &
        uw / (tu * tw), db%h%R_target(5, k) / (tu * tw), '  ', mean(1) / tu
      worst = max(worst, abs(rms(1) / tu - 1.0_rp), abs(rms(2) / sqrt(db%h%R_target(2, k)) - 1.0_rp), &
                  abs(rms(3) / tw - 1.0_rp), abs(uw - db%h%R_target(5, k)) / (tu * tw))
    end do
    ok = (worst < tol)
    write(*,'(A,F6.3,A,F5.2,A)') ' inflow-gen verify: worst deviation ', worst, '  (tolerance ', tol, ')'
    if (.not. ok) write(*,'(A)') ' inflow-gen verify: FAILED'
    call inflow_db_close(db)
  end subroutine inflow_gen_verify

end module mod_inflow_gen
