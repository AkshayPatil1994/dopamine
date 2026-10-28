! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Successor-side reader of an inflow-plane database: replays the stored segments in a
!> pseudo-random order with spanwise shifts and cross-fades between them, and interpolates the
!> planes to the inlet faces.
!>
!> Replay (a pure function of time and seed, so a restart needs no extra state):
!>   period P = T_seg - T_x, with T_x the cross-fade length. Period m uses segment k_m (random,
!>   never the same twice in a row) shifted spanwise by s_m nodes (golden-ratio sequence). During
!>   the first T_x of a period the previous period's stream fades out: u = sin(th) u_new + cos(th) u_old.
!>   The weights are uniform in space, so a divergence-free plane stays so; sin^2 + cos^2 = 1
!>   keeps the variance for decorrelated segments.
!>
!> Interpolation:
!>   time    - Catmull-Rom cubic over four frames
!>   span    - exact trigonometric (Fourier) interpolation on the periodic uniform nodes; a shift is
!>             a whole number of nodes, so it costs no accuracy and a matching face is untouched
!>   height  - six-point Lagrange on the (possibly graded) z nodes; a matching level is untouched
module mod_inflow_planes
  use mod_kinds
  use mod_parameters, only: BC_PERIODIC
  use mod_abl_types
  use mod_inflow_db
  use mod_mesh
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  integer(ip), parameter :: NCACHE = 16
  integer(int64), parameter :: HASH_SEED = 88172645463325252_int64
  real(rp), parameter :: PI = 3.14159265358979323846_rp
  real(rp), parameter :: PHI_GOLD = 0.6180339887498949_rp

  type, public :: InflowPlanes
    logical :: ready = .false.
    type(InflowDB) :: db
    ! replay schedule
    integer(ip) :: seed = 1
    real(rp) :: t_offset = 0.0_rp, T_seg = 0.0_rp, T_x = 0.0_rp, P = 0.0_rp, frac0 = 0.0_rp
    integer(ip), allocatable :: sched_seg(:)   ! segment of period m, m = 0..n_sched-1 (grown on demand)
    integer(ip) :: n_sched = 0
    ! mapped faces
    integer(ip) :: nf = 0
    integer(ip), allocatable :: f_bnd(:)       ! solver boundary-face index of each mapped face
    integer(ip), allocatable :: uy(:), uz(:)   ! unique-y / unique-z id of each face
    integer(ip) :: nuy = 0, nuz = 0
    integer(ip), allocatable :: jm(:)          ! unique y: matching node (0-based) or -1
    real(rp),    allocatable :: wy(:,:)        ! (ny, nuy) trigonometric weights
    integer(ip), allocatable :: zs0(:), zn(:)  ! unique z: first stencil node and stencil size
    real(rp),    allocatable :: wz(:,:)        ! (6, nuz)
    integer(ip) :: n_match_y = 0, n_match_z = 0
    ! frame cache
    integer(ip) :: ckey_seg(NCACHE) = -1, ckey_it(NCACHE) = -1, clast(NCACHE) = 0
    integer(ip) :: tick = 0
    real(rp), allocatable :: cdata(:,:,:,:)    ! (ny, nz, 3, NCACHE)
  end type InflowPlanes

  public :: inflow_planes_open, inflow_planes_bind, inflow_planes_eval, inflow_planes_close
  public :: inflow_stream, inflow_check_lateral_periodic, inflow_planes_time_field, inflow_unique_sorted_global

contains

  ! ------------------------------------------------------------------ hashing / schedule

  integer(int64) function mix64(a, b) result(h)
    integer(int64), intent(in) :: a, b
    integer :: r
    h = ieor(ieor(a, ishft(b, 21)), HASH_SEED)
    do r = 1, 4
      h = ieor(h, ishft(h, 13)); h = ieor(h, ishft(h, -7)); h = ieor(h, ishft(h, 17))
    end do
  end function mix64

  real(rp) function u01(h) result(u)
    integer(int64), intent(in) :: h
    u = real(iand(h, 4503599627370495_int64), rp) / 4503599627370496.0_rp
  end function u01

  !> Segment and spanwise shift (whole nodes) of replay period m >= 0.
  subroutine inflow_stream(fp, m, seg, shift)
    type(InflowPlanes), intent(inout) :: fp
    integer(ip),        intent(in)    :: m
    integer(ip),        intent(out)   :: seg, shift

    integer(ip) :: mm, k, K_n, prev
    integer(int64) :: h
    integer(ip), allocatable :: tmp(:)
    real(rp) :: f

    K_n = fp%db%h%nseg
    if (.not. allocated(fp%sched_seg)) allocate(fp%sched_seg(0:63))
    do while (m >= fp%n_sched)
      mm = fp%n_sched
      if (mm > ubound(fp%sched_seg, 1)) then
        allocate(tmp(0:2 * ubound(fp%sched_seg, 1) + 1)); tmp(0:ubound(fp%sched_seg, 1)) = fp%sched_seg
        call move_alloc(tmp, fp%sched_seg)
      end if
      h = mix64(int(fp%seed, int64), int(mm, int64))
      k = int(u01(h) * real(K_n, rp), ip) + 1_ip
      k = min(k, K_n)
      if (mm > 0 .and. K_n > 1) then
        prev = fp%sched_seg(mm - 1)
        if (k == prev) then
          h = ieor(h, ishft(h, 13)); h = ieor(h, ishft(h, -7)); h = ieor(h, ishft(h, 17))
          k = modulo(k - 1 + 1 + int(u01(h) * real(K_n - 1, rp), ip), K_n) + 1_ip
        end if
      end if
      fp%sched_seg(mm) = k
      fp%n_sched = mm + 1
    end do
    seg = fp%sched_seg(m)
    f = fp%frac0 + real(m, rp) * PHI_GOLD
    f = f - floor(f)
    shift = min(int(f * real(fp%db%h%ny, rp), ip), fp%db%h%ny - 1_ip)
  end subroutine inflow_stream

  ! ------------------------------------------------------------------ open / close

  !> Open the database named in `icfg` and set up the replay schedule.
  subroutine inflow_planes_open(fp, icfg, abl_cfg, rank)
    type(InflowPlanes), intent(inout) :: fp
    type(InflowConfig), intent(in)    :: icfg
    type(ABLConfig),    intent(in)    :: abl_cfg
    integer(ip),        intent(in)    :: rank

    integer(int64) :: h
    real(rp) :: dtp

    call inflow_db_open(fp%db, trim(icfg%file))
    associate (hd => fp%db%h)
      if (abs(hd%wind_dir_deg - abl_cfg%wind_dir_deg) > 0.5_rp) then
        write(*,'(A,F8.2,A,F8.2,A)') 'ERROR: the inflow database was generated for wind_dir_deg = ', hd%wind_dir_deg, &
          ' but the case uses ', abl_cfg%wind_dir_deg, '.'
        error stop 'inflow planes: wind direction differs from the database'
      end if
      dtp = hd%dt_p
      fp%T_seg = inflow_db_segment_time(hd)
      fp%T_x = icfg%crossfade_time
      if (fp%T_x <= 0.0_rp) fp%T_x = 0.1_rp * fp%T_seg
      fp%T_x = min(max(fp%T_x, 2.0_rp * dtp), 0.5_rp * fp%T_seg)
      fp%P = fp%T_seg - fp%T_x
      fp%seed = icfg%seed
      fp%t_offset = icfg%t_offset
      h = mix64(int(fp%seed, int64), -1_int64)
      fp%frac0 = u01(h)
      if (allocated(fp%cdata)) deallocate(fp%cdata)
      if (allocated(fp%sched_seg)) deallocate(fp%sched_seg)
      allocate(fp%cdata(hd%ny, hd%nz, 3, NCACHE))
    end associate
    fp%ckey_seg = -1; fp%ckey_it = -1; fp%clast = 0; fp%tick = 0
    fp%n_sched = 0

    if (rank == 0) then
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' Inflow planes: ', fp%db%h%ny, ' x ', fp%db%h%nz, ' plane, ', fp%db%h%nseg, &
        ' segments of ', fp%db%h%nt, ' frames'
      write(*,'(A,ES9.2,A,ES9.2,A,ES9.2,A)') ' Inflow planes: dt_p = ', fp%db%h%dt_p, ' s, segment ', fp%T_seg, &
        ' s, cross-fade ', fp%T_x, ' s'
      if (fp%db%h%nseg < 4) &
        write(*,'(A)') ' WARNING: fewer than 4 inflow segments. The spanwise-mean content is not moved by shifts and ' // &
          'repeats every replay period unless the database was decorrelated (gen_decorrelate) with enough segments.'
    end if
  end subroutine inflow_planes_open

  subroutine inflow_planes_close(fp)
    type(InflowPlanes), intent(inout) :: fp
    call inflow_db_close(fp%db)
    fp%ready = .false.
    if (allocated(fp%cdata)) deallocate(fp%cdata)
    if (allocated(fp%sched_seg)) deallocate(fp%sched_seg)
    fp%n_sched = 0
  end subroutine inflow_planes_close

  ! ------------------------------------------------------------------ mesh checks

  !> Is the domain periodic in the spanwise direction t1? Looks for periodic boundary faces with
  !> normals +-t1 and returns the period W (distance between the two periodic planes).
  subroutine inflow_check_lateral_periodic(msh, t1, W, ok)
    type(Mesh), intent(in)  :: msh
    real(rp),   intent(in)  :: t1(3)
    real(rp),   intent(out) :: W
    logical,    intent(out) :: ok

    integer(ip) :: p, f, fg, ierr
    real(rp) :: lbuf(4), gbuf(4), dn, c

    ! 1 area on +t1 faces, 2 area on -t1 faces, 3 max centre on +t1 faces, 4 -(min centre on -t1 faces)
    lbuf = [0.0_rp, 0.0_rp, -huge(1.0_rp), -huge(1.0_rp)]
    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_PERIODIC) cycle
      do f = 1, msh%patches(p)%n_faces
        fg = msh%n_int_faces + msh%patches(p)%start_face + f - 1
        if (fg < 1 .or. fg > size(msh%faces)) cycle
        dn = dot_product(msh%faces(fg)%normal, t1)
        c  = dot_product(msh%faces(fg)%centre, t1)
        if (dn > 0.99_rp) then
          lbuf(1) = lbuf(1) + msh%faces(fg)%area; lbuf(3) = max(lbuf(3), c)
        else if (dn < -0.99_rp) then
          lbuf(2) = lbuf(2) + msh%faces(fg)%area; lbuf(4) = max(lbuf(4), -c)
        end if
      end do
    end do
    call MPI_Allreduce(lbuf(1:2), gbuf(1:2), 2, MPI_DOUBLE_PRECISION, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(lbuf(3:4), gbuf(3:4), 2, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)
    W = gbuf(3) + gbuf(4)
    ok = (gbuf(1) > 0.0_rp .and. gbuf(2) > 0.0_rp .and. abs(gbuf(1) - gbuf(2)) <= 1.0e-6_rp * gbuf(1) .and. W > 0.0_rp)
  end subroutine inflow_check_lateral_periodic


  !> Sorted distinct values (to within tol) of the per-rank arrays vloc over all ranks; the result is
  !> returned on every rank. Used to recover the regular (y, z) lattice of a plane of cells or faces.
  subroutine inflow_unique_sorted_global(comm, vloc, tol, uniq)
    type(MPI_Comm),        intent(in)  :: comm
    real(rp),              intent(in)  :: vloc(:)
    real(rp),              intent(in)  :: tol
    real(rp), allocatable, intent(out) :: uniq(:)

    integer :: rank, nr, ierr, nl, i, ntot, nu
    integer, allocatable :: counts(:), displs(:)
    real(rp), allocatable :: all(:)
    real(rp) :: prev

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nr, ierr)
    nl = size(vloc)
    allocate(counts(nr), displs(nr))
    call MPI_Allgather(nl, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, comm, ierr)
    displs(1) = 0
    do i = 2, nr
      displs(i) = displs(i - 1) + counts(i - 1)
    end do
    ntot = sum(counts)
    allocate(all(max(ntot, 1)))
    call MPI_Allgatherv(vloc, nl, MPI_DOUBLE_PRECISION, all, counts, displs, MPI_DOUBLE_PRECISION, comm, ierr)
    call heap_sort(all, ntot)
    allocate(uniq(max(ntot, 1)))
    nu = 0
    prev = 0.0_rp
    do i = 1, ntot
      if (nu == 0) then
        nu = 1; uniq(1) = all(i); prev = all(i)
      else if (all(i) - prev > tol) then
        nu = nu + 1; uniq(nu) = all(i); prev = all(i)
      end if
    end do
    uniq = uniq(1:nu)
  end subroutine inflow_unique_sorted_global

  subroutine heap_sort(a, n)
    real(rp),    intent(inout) :: a(:)
    integer,     intent(in)    :: n
    integer :: i, j, k, m
    real(rp) :: t
    if (n < 2) return
    do i = n / 2, 1, -1
      call sift(i, n)
    end do
    do m = n, 2, -1
      t = a(1); a(1) = a(m); a(m) = t
      call sift(1, m - 1)
    end do
  contains
    subroutine sift(root, last)
      integer, intent(in) :: root, last
      integer :: r, c
      real(rp) :: v
      r = root; v = a(r)
      do
        c = 2 * r
        if (c > last) exit
        if (c < last) then
          if (a(c + 1) > a(c)) c = c + 1
        end if
        if (a(c) <= v) exit
        a(r) = a(c); r = c
      end do
      a(r) = v
    end subroutine sift
  end subroutine heap_sort

  ! ------------------------------------------------------------------ face binding

  !> Lagrange weights of the nodes xs(1:n) at x.
  pure subroutine lagrange_weights(xs, n, x, w)
    real(rp),    intent(in)  :: xs(:), x
    integer(ip), intent(in)  :: n
    real(rp),    intent(out) :: w(:)
    integer(ip) :: i, j
    do i = 1, n
      w(i) = 1.0_rp
      do j = 1, n
        if (j /= i) w(i) = w(i) * (x - xs(j)) / (xs(i) - xs(j))
      end do
    end do
  end subroutine lagrange_weights

  !> Trigonometric interpolation weight of periodic node j (0-based) at fractional position xi,
  !> N nodes. Even N uses the symmetric (Nyquist-split) kernel.
  pure real(rp) function trig_weight(N, xi, j) result(w)
    integer(ip), intent(in) :: N, j
    real(rp),    intent(in) :: xi
    real(rp) :: d, s
    d = xi - real(j, rp)
    s = sin(PI * d)
    if (abs(s) < 1.0e-300_rp) then
      w = merge(1.0_rp, 0.0_rp, abs(d - anint(d / real(N, rp)) * real(N, rp)) < 1.0e-9_rp)
      return
    end if
    if (mod(N, 2_ip) == 1_ip) then
      w = s / (real(N, rp) * sin(PI * d / real(N, rp)))
    else
      w = s / (real(N, rp) * tan(PI * d / real(N, rp)))
    end if
  end function trig_weight

  !> Map the inlet faces onto the database grid. xf(:, i) is the centre of face i, f_bnd(i) its
  !> boundary index in the solver. n_in / t1 / W_dom describe the successor inlet and domain.
  subroutine inflow_planes_bind(fp, xf, f_bnd, n_in, t1, W_dom, rank)
    type(InflowPlanes), intent(inout) :: fp
    real(rp),           intent(in)    :: xf(:,:), n_in(3), t1(3), W_dom
    integer(ip),        intent(in)    :: f_bnd(:), rank

    integer(ip) :: i, u, j, k, ny, nz, nf, n, s0
    real(rp), allocatable :: yq(:), zq(:), xs(:)
    real(rp) :: xi, dy, zr, zf, tolz, wk(6), frac_y, frac_z

    associate (hd => fp%db%h)
      ny = hd%ny; nz = hd%nz; nf = size(f_bnd)
      if (sum(hd%n_in * n_in) < 0.999_rp) then
        write(*,'(A,3F8.4,A,3F8.4)') 'ERROR: database flow direction ', hd%n_in, ' differs from the inlet normal ', n_in
        error stop 'inflow planes: the database was generated for a different inlet normal'
      end if
      if (sum(hd%t1 * t1) < 0.999_rp) error stop 'inflow planes: spanwise direction differs from the database'
      if (abs(hd%W - W_dom) > 1.0e-6_rp * W_dom) then
        write(*,'(A,ES11.4,A,ES11.4,A)') 'ERROR: database spanwise period W = ', hd%W, ' m but the domain is ', W_dom, ' m wide.'
        error stop 'inflow planes: database period does not match the domain width'
      end if
      dy = hd%W / real(ny, rp)
      zr = hd%z(nz) - hd%z(1)
      tolz = 1.0e-9_rp * zr

      fp%nf = nf
      if (allocated(fp%f_bnd)) deallocate(fp%f_bnd, fp%uy, fp%uz)
      allocate(fp%f_bnd(nf), fp%uy(nf), fp%uz(nf))
      fp%f_bnd = f_bnd

      ! unique spanwise offsets and heights
      allocate(yq(nf), zq(nf))
      fp%nuy = 0; fp%nuz = 0
      do i = 1, nf
        xi = (dot_product(xf(:, i), t1) - hd%y0) / dy
        xi = xi - real(ny, rp) * floor(xi / real(ny, rp))
        if (abs(xi - anint(xi)) < 1.0e-6_rp) xi = anint(xi)
        if (xi >= real(ny, rp)) xi = 0.0_rp
        zf = xf(3, i)
        if (zf < hd%z(1) - 1.0e-6_rp * zr .or. zf > hd%z(nz) + 1.0e-6_rp * zr) then
          write(*,'(A,ES11.4,A,ES11.4,A,ES11.4,A)') 'ERROR: inlet face at z = ', zf, ' m lies outside the database heights [', &
            hd%z(1), ', ', hd%z(nz), '] m.'
          error stop 'inflow planes: database does not cover the inlet heights'
        end if
        u = 0
        do j = 1, fp%nuy
          if (abs(yq(j) - xi) < 1.0e-9_rp) then
            u = j; exit
          end if
        end do
        if (u == 0) then
          fp%nuy = fp%nuy + 1; yq(fp%nuy) = xi; u = fp%nuy
        end if
        fp%uy(i) = u
        u = 0
        do j = 1, fp%nuz
          if (abs(zq(j) - zf) < tolz + 1.0e-12_rp) then
            u = j; exit
          end if
        end do
        if (u == 0) then
          fp%nuz = fp%nuz + 1; zq(fp%nuz) = zf; u = fp%nuz
        end if
        fp%uz(i) = u
      end do

      if (allocated(fp%jm)) deallocate(fp%jm, fp%wy, fp%zs0, fp%zn, fp%wz)
      allocate(fp%jm(fp%nuy), fp%wy(ny, fp%nuy), fp%zs0(fp%nuz), fp%zn(fp%nuz), fp%wz(6, fp%nuz))
      fp%wy = 0.0_rp; fp%wz = 0.0_rp
      fp%n_match_y = 0; fp%n_match_z = 0

      do u = 1, fp%nuy
        xi = yq(u)
        if (abs(xi - anint(xi)) < 1.0e-9_rp) then
          fp%jm(u) = modulo(int(anint(xi), ip), ny)
          fp%n_match_y = fp%n_match_y + 1
        else
          fp%jm(u) = -1
          do j = 0, ny - 1
            fp%wy(j + 1, u) = trig_weight(ny, xi, j)
          end do
        end if
      end do

      allocate(xs(6))
      do u = 1, fp%nuz
        zf = zq(u)
        ! interval k with z_k <= zf <= z_(k+1)
        k = 1
        do while (k < nz - 1 .and. hd%z(k + 1) <= zf)
          k = k + 1
        end do
        if (abs(zf - hd%z(k)) <= tolz) then
          fp%zs0(u) = k; fp%zn(u) = 1; fp%wz(1, u) = 1.0_rp; fp%n_match_z = fp%n_match_z + 1
        else if (abs(zf - hd%z(k + 1)) <= tolz) then
          fp%zs0(u) = k + 1; fp%zn(u) = 1; fp%wz(1, u) = 1.0_rp; fp%n_match_z = fp%n_match_z + 1
        else
          n = min(6_ip, nz)
          s0 = min(max(k - 2_ip, 1_ip), nz - n + 1_ip)
          fp%zs0(u) = s0; fp%zn(u) = n
          xs(1:n) = hd%z(s0:s0 + n - 1)
          call lagrange_weights(xs, n, zf, wk)
          fp%wz(1:n, u) = wk(1:n)
        end if
      end do

      frac_y = real(count([(fp%jm(fp%uy(i)) >= 0, i = 1, nf)]), rp) / real(max(nf, 1_ip), rp)
      frac_z = real(count([(fp%zn(fp%uz(i)) == 1, i = 1, nf)]), rp) / real(max(nf, 1_ip), rp)
    end associate
    fp%ready = .true.

    if (rank == 0) then
      write(*,'(A,I0,A,I0,A,I0,A)') ' Inflow planes: ', nf, ' inlet faces (rank 0 share) on ', fp%nuy, ' spanwise offsets, ', &
        fp%nuz, ' heights'
      write(*,'(A,F6.1,A,F6.1,A)') ' Inflow planes: faces coinciding with database nodes - spanwise ', 100.0_rp * frac_y, &
        '%, vertical ', 100.0_rp * frac_z, '% (copied, not filtered)'
      if (frac_y < 1.0_rp .or. frac_z < 1.0_rp) &
        write(*,'(A)') ' Inflow planes: other faces are interpolated (span: exact Fourier, height: 6-point Lagrange).'
      call inflow_planes_diagnose(fp, frac_z < 1.0_rp)
    end if
  end subroutine inflow_planes_bind


  !> Resolution report for the stored data: how well the frame spacing and the (possibly graded)
  !> height spacing resolve the signal. Rank 0 only; reads a few frames of the first segment.
  !>   time   - correlation of consecutive frames (the cubic interpolation needs it near 1)
  !>   height - predict every interior level from its six nearest neighbours (omitting it, so a
  !>            gap twice the node spacing, i.e. worse than any real offset) and report the
  !>            smallest variance retained
  subroutine inflow_planes_diagnose(fp, check_z)
    type(InflowPlanes), intent(inout) :: fp
    logical,            intent(in)    :: check_z   ! only heights that are interpolated are worth reporting

    integer(ip) :: nt, nz, ny, ns, i0, i, k, c, j
    real(rp), allocatable :: pa(:,:,:), pb(:,:,:)
    real(rp) :: sab, saa, sbb, xs(6), w(6), pred, ret, ret_min, k_min, z_at
    real(rp), allocatable :: sp(:), st(:)

    nt = fp%db%h%nt; nz = fp%db%h%nz; ny = fp%db%h%ny
    ns = min(nt - 1_ip, 24_ip)
    i0 = max(1_ip, nt / 2 - ns / 2)
    allocate(pa(ny, nz, 3), pb(ny, nz, 3), sp(nz), st(nz))
    sab = 0.0_rp; saa = 0.0_rp; sbb = 0.0_rp; sp = 0.0_rp; st = 0.0_rp
    do i = i0, min(i0 + ns - 1_ip, nt - 1_ip)
      call inflow_db_read_frame(fp%db, 1_ip, i, pa)
      call inflow_db_read_frame(fp%db, 1_ip, i + 1_ip, pb)
      sab = sab + sum(pa * pb); saa = saa + sum(pa * pa); sbb = sbb + sum(pb * pb)
      do k = 4, nz - 3
        xs(1:3) = fp%db%h%z(k - 3:k - 1); xs(4:6) = fp%db%h%z(k + 1:k + 3)
        call lagrange_weights(xs, 6_ip, fp%db%h%z(k), w)
        do c = 1, 3
          do j = 1, ny
            pred = dot_product(w(1:3), pa(j, k - 3:k - 1, c)) + dot_product(w(4:6), pa(j, k + 1:k + 3, c))
            sp(k) = sp(k) + pred**2
            st(k) = st(k) + pa(j, k, c)**2
          end do
        end do
      end do
    end do
    write(*,'(A,F7.4,A)') ' Inflow planes: frame-to-frame correlation ', sab / sqrt(max(saa * sbb, 1.0e-300_rp)), &
      ' (time interpolation is accurate when this is close to 1)'
    ret_min = 1.0e30_rp; k_min = 0.0_rp
    do k = 4, nz - 3
      if (st(k) <= 0.0_rp) cycle
      ret = sp(k) / st(k)
      if (ret < ret_min) then
        ret_min = ret; k_min = fp%db%h%z(k)
      end if
    end do
    z_at = k_min
    if (check_z .and. nz >= 8 .and. ret_min < 1.0e29_rp) then
      write(*,'(A,F6.3,A,ES9.2,A)') ' Inflow planes: vertical resolution - smallest variance ratio of the 6-point stencil ', &
        ret_min, ' (at z = ', z_at, ' m; omit-one test, pessimistic)'
      if (ret_min < 0.97_rp .or. ret_min > 1.03_rp) write(*,'(A)') ' WARNING: the database heights under-resolve the ' // &
        'vertical structure; interpolation to off-node faces will damp or distort it. Regenerate with more levels, or on ' // &
        'the inlet face heights (dopamine-inflow-gen with gen_grid_file).'
    end if
    if (sab / sqrt(max(saa * sbb, 1.0e-300_rp)) < 0.95_rp) &
      write(*,'(A)') ' WARNING: consecutive frames are poorly correlated; the frame spacing under-resolves the signal.'
  end subroutine inflow_planes_diagnose

  ! ------------------------------------------------------------------ evaluation

  !> Catmull-Rom weights for the frames at offsets -1, 0, 1, 2 and fraction th in [0, 1).
  pure subroutine inflow_catmull_rom(th, w)
    real(rp), intent(in)  :: th
    real(rp), intent(out) :: w(4)
    w(1) = -0.5_rp * th**3 + th**2 - 0.5_rp * th
    w(2) =  1.5_rp * th**3 - 2.5_rp * th**2 + 1.0_rp
    w(3) = -1.5_rp * th**3 + 2.0_rp * th**2 + 0.5_rp * th
    w(4) =  0.5_rp * th**3 - 0.5_rp * th**2
  end subroutine inflow_catmull_rom

  !> Cache slot holding frame (seg, it), reading it on a miss (least recently used slot is replaced).
  integer(ip) function frame_slot(fp, seg, it) result(slot)
    type(InflowPlanes), intent(inout) :: fp
    integer(ip),        intent(in)    :: seg, it
    integer(ip) :: c, oldest
    real(rp), allocatable :: buf(:,:,:)

    fp%tick = fp%tick + 1
    do c = 1, NCACHE
      if (fp%ckey_seg(c) == seg .and. fp%ckey_it(c) == it) then
        slot = c; fp%clast(c) = fp%tick; return
      end if
    end do
    oldest = 1
    do c = 2, NCACHE
      if (fp%clast(c) < fp%clast(oldest)) oldest = c
    end do
    allocate(buf(fp%db%h%ny, fp%db%h%nz, 3))
    call inflow_db_read_frame(fp%db, seg, it, buf)
    fp%cdata(:, :, :, oldest) = buf
    fp%ckey_seg(oldest) = seg; fp%ckey_it(oldest) = it; fp%clast(oldest) = fp%tick
    slot = oldest
  end function frame_slot

  !> The cubic-in-time plane of stream (seg) at local time tau, shifted by `shift` nodes, before
  !> spatial interpolation: out(ny, nz, 3) in the profile frame.
  subroutine inflow_planes_time_field(fp, seg, tau, out)
    type(InflowPlanes), intent(inout) :: fp
    integer(ip),        intent(in)    :: seg
    real(rp),           intent(in)    :: tau
    real(rp),           intent(out)   :: out(:,:,:)

    integer(ip) :: i0, c, it, slot, nt
    real(rp) :: x, th, w(4)

    nt = fp%db%h%nt
    x = min(max(tau, 0.0_rp), fp%T_seg) / fp%db%h%dt_p
    i0 = int(floor(x), ip) + 1_ip
    th = x - floor(x)
    if (i0 >= nt) then
      i0 = nt - 1_ip; th = 1.0_rp
    end if
    call inflow_catmull_rom(th, w)
    out = 0.0_rp
    do c = 1, 4
      it = min(max(i0 + c - 2_ip, 1_ip), nt)
      slot = frame_slot(fp, seg, it)
      out = out + w(c) * fp%cdata(:, :, :, slot)
    end do
  end subroutine inflow_planes_time_field

  !> Spatially interpolate a time-field plane to the faces, applying a spanwise shift of `shift`
  !> nodes; result in the profile frame.
  subroutine plane_to_faces(fp, plane, shift, ul)
    type(InflowPlanes), intent(in)  :: fp
    real(rp),           intent(in)  :: plane(:,:,:)
    integer(ip),        intent(in)  :: shift
    real(rp),           intent(out) :: ul(:,:)

    integer(ip) :: u, j, k, c, ny, nz, i, n, s0, jj
    real(rp), allocatable :: g(:,:,:)
    real(rp) :: acc

    ny = fp%db%h%ny; nz = fp%db%h%nz
    allocate(g(fp%nuy, nz, 3))
    do c = 1, 3
      do k = 1, nz
        do u = 1, fp%nuy
          if (fp%jm(u) >= 0) then
            g(u, k, c) = plane(modulo(fp%jm(u) + shift, ny) + 1, k, c)
          else
            acc = 0.0_rp
            do j = 0, ny - 1
              jj = modulo(j + shift, ny) + 1
              acc = acc + fp%wy(j + 1, u) * plane(jj, k, c)
            end do
            g(u, k, c) = acc
          end if
        end do
      end do
    end do
    do i = 1, fp%nf
      u = fp%uz(i); n = fp%zn(u); s0 = fp%zs0(u)
      do c = 1, 3
        acc = 0.0_rp
        do k = 1, n
          acc = acc + fp%wz(k, u) * g(fp%uy(i), s0 + k - 1, c)
        end do
        ul(c, i) = acc
      end do
    end do
  end subroutine plane_to_faces

  !> Fluctuation (global axes) at the mapped faces at physical time t.
  subroutine inflow_planes_eval(fp, t, fluc)
    type(InflowPlanes), intent(inout) :: fp
    real(rp),           intent(in)    :: t
    real(rp),           intent(out)   :: fluc(:,:)   ! (3, nf)

    integer(ip) :: m, seg, shift, i, c
    real(rp) :: te, tau, th, wn, wo
    real(rp), allocatable :: pl(:,:,:), ua(:,:), ub(:,:)

    if (.not. fp%ready) error stop 'inflow_planes_eval: faces are not bound'
    te = max(t + fp%t_offset, 0.0_rp)
    m = int(floor(te / fp%P), ip)
    tau = te - real(m, rp) * fp%P

    allocate(pl(fp%db%h%ny, fp%db%h%nz, 3), ua(3, fp%nf))
    call inflow_stream(fp, m, seg, shift)
    call inflow_planes_time_field(fp, seg, tau, pl)
    call plane_to_faces(fp, pl, shift, ua)

    if (m >= 1 .and. tau < fp%T_x) then
      allocate(ub(3, fp%nf))
      call inflow_stream(fp, m - 1_ip, seg, shift)
      call inflow_planes_time_field(fp, seg, tau + fp%P, pl)
      call plane_to_faces(fp, pl, shift, ub)
      th = 0.5_rp * PI * tau / fp%T_x
      wn = sin(th); wo = cos(th)
      ua = wn * ua + wo * ub
    end if

    ! profile frame -> global axes
    do i = 1, fp%nf
      do c = 1, 3
        fluc(c, i) = fp%db%h%frame(c, 1) * ua(1, i) + fp%db%h%frame(c, 2) * ua(2, i) + fp%db%h%frame(c, 3) * ua(3, i)
      end do
    end do
  end subroutine inflow_planes_eval

end module mod_inflow_planes
