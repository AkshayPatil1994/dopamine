! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Inflow-plane database: NetCDF-4 file of velocity-fluctuation planes u'(y, z, comp) sampled at a
!> fixed time step, stored as several independent segments.
!>
!> Producers (dopamine-inflow-gen, the precursor recorder) write this format; the successor run
!> reads it. Conventions:
!>   - the plane is spanwise periodic: y_j = y0 + j*W/ny, j = 0..ny-1 (uniform), y along `t1`
!>   - z is height above ground, ascending, any spacing (graded is fine)
!>   - u holds the FLUCTUATION in the profile frame (columns of `frame`; 1 = streamwise)
!>   - frame i of a segment is at t = (i-1)*dt_p; a segment spans (nt-1)*dt_p
!>   - plane position along the normal is implicit: the database is the inlet
module mod_inflow_db
  use mod_kinds
  use netcdf
  implicit none
  private


  !> Everything about the planes except the planes themselves.
  type, public :: InflowDBHeader
    integer(ip) :: ny = 0, nz = 0, nt = 0, nseg = 0
    real(rp)    :: W = 0.0_rp          ! spanwise period [m]
    real(rp)    :: y0 = 0.0_rp         ! t1 coordinate of the first node (absolute) [m]
    real(rp)    :: dt_p = 0.0_rp       ! frame spacing [s]
    real(rp)    :: frame(3,3) = 0.0_rp ! columns: streamwise, lateral, vertical axes of the components
    real(rp)    :: n_in(3) = 0.0_rp    ! unit flow direction through the plane
    real(rp)    :: t1(3) = 0.0_rp      ! unit spanwise direction (y axis of the plane)
    real(rp)    :: wind_dir_deg = 0.0_rp
    real(rp)    :: u_ref = 0.0_rp, z_ref = 0.0_rp, T_turnover = 0.0_rp
    integer(ip) :: seed = 0
    character(len=64)  :: generator = ''
    real(rp), allocatable :: z(:)          ! (nz) heights [m]
    real(rp), allocatable :: U_mean(:)     ! (nz) mean streamwise speed [m/s] (informational)
    real(rp), allocatable :: R_target(:,:) ! (6, nz) target uu,vv,ww,uv,uw,vw in the profile frame (QA)
  end type InflowDBHeader

  type, public :: InflowDB
    integer :: ncid = -1
    integer :: varid_u = -1
    integer :: varid_ub = -1, varid_rs = -1
    logical :: writable = .false.
    type(InflowDBHeader) :: h
  end type InflowDB

  public :: inflow_db_create, inflow_db_open, inflow_db_close
  public :: inflow_db_write_frame, inflow_db_read_frame, inflow_db_node_y, inflow_db_segment_time
  public :: inflow_db_decorrelate_span_mean, inflow_db_write_profiles

contains

  subroutine chk(status, what)
    integer,          intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) then
      write(*,'(A)') 'ERROR: inflow database: '//trim(what)//': '//trim(nf90_strerror(status))
      error stop 'inflow database I/O error'
    end if
  end subroutine chk

  !> y coordinate (absolute, along t1) of node j (1-based).
  pure real(rp) function inflow_db_node_y(h, j) result(y)
    type(InflowDBHeader), intent(in) :: h
    integer(ip),          intent(in) :: j
    y = h%y0 + real(j - 1, rp) * h%W / real(h%ny, rp)
  end function inflow_db_node_y

  !> Duration of one segment [s]: first to last frame.
  pure real(rp) function inflow_db_segment_time(h) result(T)
    type(InflowDBHeader), intent(in) :: h
    T = real(h%nt - 1, rp) * h%dt_p
  end function inflow_db_segment_time

  !> Create a new database file with the given header (planes are written afterwards).
  subroutine inflow_db_create(db, filename, h)
    type(InflowDB),       intent(inout) :: db
    character(len=*),     intent(in)    :: filename
    type(InflowDBHeader), intent(in)    :: h

    integer :: dy, dz, dc, dt, ds, dsix, vz, vy, vU, vR, ny, nz, k
    real(rp), allocatable :: yy(:)

    if (h%ny < 2 .or. h%nz < 2 .or. h%nt < 4 .or. h%nseg < 1) &
      error stop 'inflow_db_create: need ny>=2, nz>=2, nt>=4, nseg>=1'
    if (h%W <= 0.0_rp .or. h%dt_p <= 0.0_rp) error stop 'inflow_db_create: W and dt_p must be positive'
    if (.not. allocated(h%z)) error stop 'inflow_db_create: header has no z levels'
    ny = h%ny; nz = h%nz
    do k = 2, nz
      if (h%z(k) <= h%z(k-1)) error stop 'inflow_db_create: z must be strictly ascending'
    end do

    db%h = h
    db%writable = .true.
    call chk(nf90_create(trim(filename), ior(nf90_clobber, nf90_netcdf4), db%ncid), 'create '//trim(filename))
    call chk(nf90_def_dim(db%ncid, 'y', ny, dy), 'def y')
    call chk(nf90_def_dim(db%ncid, 'z', nz, dz), 'def z')
    call chk(nf90_def_dim(db%ncid, 'comp', 3, dc), 'def comp')
    call chk(nf90_def_dim(db%ncid, 'time', h%nt, dt), 'def time')
    call chk(nf90_def_dim(db%ncid, 'seg', h%nseg, ds), 'def seg')
    call chk(nf90_def_dim(db%ncid, 'rcomp', 6, dsix), 'def rcomp')

    call chk(nf90_def_var(db%ncid, 'y', nf90_double, [dy], vy), 'def var y')
    call chk(nf90_def_var(db%ncid, 'z', nf90_double, [dz], vz), 'def var z')
    call chk(nf90_def_var(db%ncid, 'U_mean', nf90_double, [dz], vU), 'def var U_mean')
    call chk(nf90_def_var(db%ncid, 'R_target', nf90_double, [dsix, dz], vR), 'def var R_target')
    db%varid_ub = vU; db%varid_rs = vR
    call chk(nf90_def_var(db%ncid, 'u', nf90_float, [dy, dz, dc, dt, ds], db%varid_u), 'def var u')
    ! One frame per chunk: the consumer reads whole frames, in time order.
    call chk(nf90_def_var_chunking(db%ncid, db%varid_u, nf90_chunked, [ny, nz, 3, 1, 1]), 'chunking u')

    call chk(nf90_put_att(db%ncid, nf90_global, 'format', 'dopamine-inflow-planes-1'), 'att format')
    call chk(nf90_put_att(db%ncid, nf90_global, 'W', h%W), 'att W')
    call chk(nf90_put_att(db%ncid, nf90_global, 'y0', h%y0), 'att y0')
    call chk(nf90_put_att(db%ncid, nf90_global, 'dt_p', h%dt_p), 'att dt_p')
    call chk(nf90_put_att(db%ncid, nf90_global, 'frame', reshape(h%frame, [9])), 'att frame')
    call chk(nf90_put_att(db%ncid, nf90_global, 'n_in', h%n_in), 'att n_in')
    call chk(nf90_put_att(db%ncid, nf90_global, 't1', h%t1), 'att t1')
    call chk(nf90_put_att(db%ncid, nf90_global, 'wind_dir_deg', h%wind_dir_deg), 'att wind_dir_deg')
    call chk(nf90_put_att(db%ncid, nf90_global, 'u_ref', h%u_ref), 'att u_ref')
    call chk(nf90_put_att(db%ncid, nf90_global, 'z_ref', h%z_ref), 'att z_ref')
    call chk(nf90_put_att(db%ncid, nf90_global, 'T_turnover', h%T_turnover), 'att T_turnover')
    call chk(nf90_put_att(db%ncid, nf90_global, 'seed', int(h%seed)), 'att seed')
    call chk(nf90_put_att(db%ncid, nf90_global, 'generator', trim(h%generator)), 'att generator')
    call chk(nf90_enddef(db%ncid), 'enddef')

    allocate(yy(ny))
    do k = 1, ny
      yy(k) = inflow_db_node_y(h, int(k, ip))
    end do
    call chk(nf90_put_var(db%ncid, vy, yy), 'put y')
    call chk(nf90_put_var(db%ncid, vz, h%z), 'put z')
    if (allocated(h%U_mean)) then
      call chk(nf90_put_var(db%ncid, vU, h%U_mean), 'put U_mean')
    else
      call chk(nf90_put_var(db%ncid, vU, spread(0.0_rp, 1, nz)), 'put U_mean')
    end if
    if (allocated(h%R_target)) then
      call chk(nf90_put_var(db%ncid, vR, h%R_target), 'put R_target')
    else
      call chk(nf90_put_var(db%ncid, vR, spread(spread(0.0_rp, 1, 6), 2, nz)), 'put R_target')
    end if
  end subroutine inflow_db_create

  !> Write one frame (segment iseg, frame it, both 1-based). plane(y, z, comp) in the profile frame.
  subroutine inflow_db_write_frame(db, iseg, it, plane)
    type(InflowDB), intent(in) :: db
    integer(ip),    intent(in) :: iseg, it
    real(rp),       intent(in) :: plane(:,:,:)

    if (.not. db%writable) error stop 'inflow_db_write_frame: database opened read-only'
    if (size(plane,1) /= db%h%ny .or. size(plane,2) /= db%h%nz .or. size(plane,3) /= 3) &
      error stop 'inflow_db_write_frame: plane has the wrong shape'
    if (iseg < 1 .or. iseg > db%h%nseg .or. it < 1 .or. it > db%h%nt) &
      error stop 'inflow_db_write_frame: segment or frame index out of range'
    call chk(nf90_put_var(db%ncid, db%varid_u, real(plane, sp), &
             start=[1, 1, 1, int(it), int(iseg)], count=[db%h%ny, db%h%nz, 3, 1, 1]), 'put frame')
  end subroutine inflow_db_write_frame

  !> Set the informational mean-speed and stress profiles after the planes are known (precursor recording).
  subroutine inflow_db_write_profiles(db, U_mean, R_target)
    type(InflowDB), intent(inout) :: db
    real(rp),       intent(in)    :: U_mean(:), R_target(:,:)
    if (.not. db%writable) error stop 'inflow_db_write_profiles: database opened read-only'
    call chk(nf90_put_var(db%ncid, db%varid_ub, U_mean), 'put U_mean')
    call chk(nf90_put_var(db%ncid, db%varid_rs, R_target), 'put R_target')
  end subroutine inflow_db_write_profiles

  !> Open an existing database read-only and validate its header.
  subroutine inflow_db_open(db, filename, writable)
    type(InflowDB),   intent(inout) :: db
    character(len=*), intent(in)    :: filename
    logical,          intent(in), optional :: writable

    integer :: dimid, ierr, vz, vU, vR, k, ival
    integer :: nyy, nzz, ntt, nss
    real(rp) :: f9(9), dy
    character(len=64) :: fmt
    character(len=len(db%h%generator)) :: gen
    real(rp), allocatable :: yy(:)

    db%writable = .false.
    if (present(writable)) db%writable = writable
    ierr = nf90_open(trim(filename), merge(nf90_write, nf90_nowrite, db%writable), db%ncid)
    if (ierr /= nf90_noerr) then
      write(*,'(A)') 'ERROR: cannot open inflow database '//trim(filename)//': '//trim(nf90_strerror(ierr))
      error stop 'inflow_db_open: open failed'
    end if
    fmt = ''
    call chk(nf90_get_att(db%ncid, nf90_global, 'format', fmt), 'read format attribute')
    if (trim(fmt) /= 'dopamine-inflow-planes-1') then
      write(*,'(A)') 'ERROR: '//trim(filename)//' is not a dopamine inflow database (format = "'//trim(fmt)//'")'
      error stop 'inflow_db_open: wrong file format'
    end if

    call chk(nf90_inq_dimid(db%ncid, 'y', dimid), 'dim y');    call chk(nf90_inquire_dimension(db%ncid, dimid, len=nyy), 'len y')
    call chk(nf90_inq_dimid(db%ncid, 'z', dimid), 'dim z');    call chk(nf90_inquire_dimension(db%ncid, dimid, len=nzz), 'len z')
    call chk(nf90_inq_dimid(db%ncid, 'time', dimid), 'dim t'); call chk(nf90_inquire_dimension(db%ncid, dimid, len=ntt), 'len t')
    call chk(nf90_inq_dimid(db%ncid, 'seg', dimid), 'dim s');  call chk(nf90_inquire_dimension(db%ncid, dimid, len=nss), 'len s')
    db%h%ny = nyy; db%h%nz = nzz; db%h%nt = ntt; db%h%nseg = nss

    call chk(nf90_get_att(db%ncid, nf90_global, 'W', db%h%W), 'att W')
    call chk(nf90_get_att(db%ncid, nf90_global, 'y0', db%h%y0), 'att y0')
    call chk(nf90_get_att(db%ncid, nf90_global, 'dt_p', db%h%dt_p), 'att dt_p')
    call chk(nf90_get_att(db%ncid, nf90_global, 'frame', f9), 'att frame')
    db%h%frame = reshape(f9, [3, 3])
    call chk(nf90_get_att(db%ncid, nf90_global, 'n_in', db%h%n_in), 'att n_in')
    call chk(nf90_get_att(db%ncid, nf90_global, 't1', db%h%t1), 'att t1')
    call chk(nf90_get_att(db%ncid, nf90_global, 'wind_dir_deg', db%h%wind_dir_deg), 'att wind_dir_deg')
    call chk(nf90_get_att(db%ncid, nf90_global, 'u_ref', db%h%u_ref), 'att u_ref')
    call chk(nf90_get_att(db%ncid, nf90_global, 'z_ref', db%h%z_ref), 'att z_ref')
    call chk(nf90_get_att(db%ncid, nf90_global, 'T_turnover', db%h%T_turnover), 'att T_turnover')
    call chk(nf90_get_att(db%ncid, nf90_global, 'seed', ival), 'att seed')
    db%h%seed = ival
    gen = ''
    ierr = nf90_get_att(db%ncid, nf90_global, 'generator', gen)
    db%h%generator = gen

    if (allocated(db%h%z)) deallocate(db%h%z)
    if (allocated(db%h%U_mean)) deallocate(db%h%U_mean)
    if (allocated(db%h%R_target)) deallocate(db%h%R_target)
    allocate(db%h%z(nzz), db%h%U_mean(nzz), db%h%R_target(6, nzz), yy(nyy))
    call chk(nf90_inq_varid(db%ncid, 'z', vz), 'var z');            call chk(nf90_get_var(db%ncid, vz, db%h%z), 'get z')
    call chk(nf90_inq_varid(db%ncid, 'U_mean', vU), 'var U_mean');  call chk(nf90_get_var(db%ncid, vU, db%h%U_mean), 'get U_mean')
    call chk(nf90_inq_varid(db%ncid, 'R_target', vR), 'var R');     call chk(nf90_get_var(db%ncid, vR, db%h%R_target), 'get R')
    call chk(nf90_inq_varid(db%ncid, 'u', db%varid_u), 'var u')

    ! Header sanity: the consumer relies on all of these.
    if (db%h%ny < 2 .or. db%h%nz < 2 .or. db%h%nt < 4) error stop 'inflow database: too few planes or frames'
    if (db%h%W <= 0.0_rp .or. db%h%dt_p <= 0.0_rp) error stop 'inflow database: non-positive W or dt_p'
    do k = 2, nzz
      if (db%h%z(k) <= db%h%z(k-1)) error stop 'inflow database: z levels not strictly ascending'
    end do
    call chk(nf90_inq_varid(db%ncid, 'y', vz), 'var y'); call chk(nf90_get_var(db%ncid, vz, yy), 'get y')
    dy = db%h%W / real(nyy, rp)
    do k = 1, nyy
      if (abs(yy(k) - inflow_db_node_y(db%h, int(k, ip))) > 1.0e-6_rp * dy) &
        error stop 'inflow database: y nodes are not uniform over the period (y0 + j W/ny)'
    end do
  end subroutine inflow_db_open

  !> Read one frame into plane(y, z, comp) (profile frame, double precision).
  subroutine inflow_db_read_frame(db, iseg, it, plane)
    type(InflowDB), intent(in)  :: db
    integer(ip),    intent(in)  :: iseg, it
    real(rp),       intent(out) :: plane(:,:,:)

    real(sp), allocatable :: buf(:,:,:,:,:)

    if (iseg < 1 .or. iseg > db%h%nseg .or. it < 1 .or. it > db%h%nt) &
      error stop 'inflow_db_read_frame: segment or frame index out of range'
    allocate(buf(db%h%ny, db%h%nz, 3, 1, 1))
    call chk(nf90_get_var(db%ncid, db%varid_u, buf, start=[1, 1, 1, int(it), int(iseg)], &
             count=[db%h%ny, db%h%nz, 3, 1, 1]), 'get frame')
    plane = real(buf(:, :, :, 1, 1), rp)
  end subroutine inflow_db_read_frame


  !> Remove the spanwise-mean waveform shifts cannot hide: subtract the per-frame mean over K >= 2 segments (rescaled).
  subroutine inflow_db_decorrelate_span_mean(filename, span_var_fraction)
    character(len=*), intent(in)  :: filename
    real(rp),         intent(out) :: span_var_fraction  ! spanwise-mean share of the plane variance (before)

    type(InflowDB) :: db
    real(rp), allocatable :: ybar(:,:,:,:), mbar(:,:,:), plane(:,:,:)
    real(rp) :: scale, tot, spn
    integer(ip) :: is, it, ny, nz, nt, ns, c, kz

    call inflow_db_open(db, filename, writable=.true.)
    ny = db%h%ny; nz = db%h%nz; nt = db%h%nt; ns = db%h%nseg
    span_var_fraction = 0.0_rp
    if (ns < 2) then
      write(*,'(A)') ' WARNING: one segment only; the spanwise-mean line at the replay period cannot be removed.'
      call inflow_db_close(db)
      return
    end if
    allocate(ybar(nz, 3, nt, ns), mbar(nz, 3, nt), plane(ny, nz, 3))
    tot = 0.0_rp
    do is = 1, ns
      do it = 1, nt
        call inflow_db_read_frame(db, is, it, plane)
        do c = 1, 3
          ybar(:, c, it, is) = sum(plane(:, :, c), dim=1) / real(ny, rp)
        end do
        tot = tot + sum(plane**2)
      end do
    end do
    mbar = sum(ybar, dim=4) / real(ns, rp)
    scale = sqrt(real(ns, rp) / real(ns - 1, rp))
    spn = real(ny, rp) * sum(ybar**2)
    span_var_fraction = spn / max(tot, 1.0e-300_rp)

    do is = 1, ns
      do it = 1, nt
        call inflow_db_read_frame(db, is, it, plane)
        do c = 1, 3
          do kz = 1, nz
            plane(:, kz, c) = plane(:, kz, c) + (scale * (ybar(kz, c, it, is) - mbar(kz, c, it)) - ybar(kz, c, it, is))
          end do
        end do
        call inflow_db_write_frame(db, is, it, plane)
      end do
    end do
    call inflow_db_close(db)
  end subroutine inflow_db_decorrelate_span_mean

  subroutine inflow_db_close(db)
    type(InflowDB), intent(inout) :: db
    if (db%ncid >= 0) call chk(nf90_close(db%ncid), 'close')
    db%ncid = -1
  end subroutine inflow_db_close

end module mod_inflow_db
