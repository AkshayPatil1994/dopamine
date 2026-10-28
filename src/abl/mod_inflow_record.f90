! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Precursor-side recorder: samples a plane of the running simulation into an inflow-plane database.
!>
!> The plane is a layer of cell centres (all cells whose centroid lies on the plane through
!> `record_point` with normal `record_normal`), so the stored values are the solver's own cell
!> values with no interpolation. The layer must form a regular (y, z) lattice, uniform across a
!> spanwise-periodic domain; otherwise the run stops. Frames are taken at fixed times
!> t_start + n*dt (linear interpolation between solver steps), so adaptive time stepping is fine.
!> When the requested segments are complete the mean profile is removed (the database holds
!> fluctuations only), the measured stresses are stored, and the spanwise-mean waveform shared
!> by the segments is removed (see inflow_db_decorrelate_span_mean).
module mod_inflow_record
  use mod_kinds
  use mod_mesh
  use mod_abl_types
  use mod_inflow_db
  use mod_inflow_planes, only: inflow_check_lateral_periodic, inflow_unique_sorted_global
  use mod_fv_solver
  use mpi_f08
  implicit none
  private

  real(rp), parameter :: PI = 3.14159265358979323846_rp

  type, public :: InflowRecorder
    logical :: active = .false., done = .false., warned_dt = .false.
    type(InflowRecordConfig) :: cfg
    type(InflowDB) :: db                       ! rank 0 only
    integer(ip) :: ny = 0, nz = 0, nt = 0, nseg = 0, rank = 0
    integer(ip) :: n_cell = 0
    integer(ip), allocatable :: cell(:), iy(:), iz(:)
    real(rp) :: n_in(3) = 0.0_rp, t1(3) = 0.0_rp, frame(3,3) = 0.0_rp
    real(rp), allocatable :: P_prev(:,:,:), P_now(:,:,:)
    real(rp) :: t_prev = 0.0_rp
    integer(ip) :: n_next = 0, n_total = 0
    real(rp), allocatable :: S1(:,:), S2(:,:)  ! rank 0: sums of u and of the six products, per height
    real(rp) :: n_samp = 0.0_rp
    character(len=256) :: file = ''
  end type InflowRecorder

  public :: inflow_record_init, inflow_record_step, inflow_record_finalize

contains

  pure function cross3(a, b) result(c)
    real(rp), intent(in) :: a(3), b(3)
    real(rp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  !> Find the cell layer, build the lattice and open the database. Collective.
  subroutine inflow_record_init(rec, solver, cfg)
    type(InflowRecorder), intent(inout) :: rec
    type(FVSolver),       intent(inout) :: solver
    type(InflowRecordConfig), intent(in) :: cfg

    type(Mesh), pointer :: msh
    integer(ip) :: c, i, n_loc, j, ierr, jy, kz
    real(rp) :: n_in(3), t1(3), s, sabs_loc, sabs, h_loc, h, tol, s_layer, s_loc, W
    real(rp), allocatable :: al(:), zl(:), ua(:), uz(:)
    integer(ip), allocatable :: hits(:,:), hits_g(:,:)
    logical :: ok
    type(InflowDBHeader) :: hd
    real(rp) :: dy, y0, zscale

    msh => solver%mesh
    rec%cfg = cfg
    rec%rank = msh%rank
    rec%file = cfg%file
    n_in = cfg%normal / sqrt(sum(cfg%normal**2))
    if (abs(n_in(3)) > 1.0e-6_rp) error stop 'inflow record: the plane normal (flow direction) must be horizontal'
    t1 = cross3([0.0_rp, 0.0_rp, 1.0_rp], n_in)
    rec%n_in = n_in; rec%t1 = t1
    rec%frame(:, 1) = n_in; rec%frame(:, 2) = t1; rec%frame(:, 3) = [0.0_rp, 0.0_rp, 1.0_rp]

    call inflow_check_lateral_periodic(msh, t1, W, ok)
    if (.not. ok) then
      if (rec%rank == 0) write(*,'(A)') 'ERROR: inflow recording needs a domain that is periodic in the spanwise direction.'
      error stop 'inflow record: spanwise periodicity required'
    end if

    ! the layer: cells whose centroid is on the plane closest to record_point (positive side on ties)
    sabs_loc = huge(1.0_rp); h_loc = huge(1.0_rp)
    do c = 1, msh%n_cells
      s = dot_product(msh%cells(c)%centroid - cfg%point, n_in)
      if (abs(s) < sabs_loc) then
        sabs_loc = abs(s); h_loc = msh%cells(c)%volume**(1.0_rp / 3.0_rp)
      end if
    end do
    call MPI_Allreduce(sabs_loc, sabs, 1, MPI_DOUBLE_PRECISION, MPI_MIN, msh%comm, ierr)
    call MPI_Allreduce(h_loc, h, 1, MPI_DOUBLE_PRECISION, MPI_MIN, msh%comm, ierr)
    tol = 1.0e-3_rp * h
    s_loc = -huge(1.0_rp)
    do c = 1, msh%n_cells
      s = dot_product(msh%cells(c)%centroid - cfg%point, n_in)
      if (abs(abs(s) - sabs) <= tol) s_loc = max(s_loc, s)
    end do
    call MPI_Allreduce(s_loc, s_layer, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)

    n_loc = 0
    do c = 1, msh%n_cells
      s = dot_product(msh%cells(c)%centroid - cfg%point, n_in)
      if (abs(s - s_layer) <= tol) n_loc = n_loc + 1
    end do
    allocate(rec%cell(n_loc), al(n_loc), zl(n_loc))
    i = 0
    do c = 1, msh%n_cells
      s = dot_product(msh%cells(c)%centroid - cfg%point, n_in)
      if (abs(s - s_layer) > tol) cycle
      i = i + 1
      rec%cell(i) = c
      al(i) = dot_product(msh%cells(c)%centroid, t1)
      zl(i) = msh%cells(c)%centroid(3)
    end do
    rec%n_cell = n_loc

    zscale = 1.0_rp
    if (n_loc > 0) zscale = max(maxval(abs(zl)), 1.0_rp)
    call MPI_Allreduce(MPI_IN_PLACE, zscale, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)
    call inflow_unique_sorted_global(msh%comm, al, 1.0e-9_rp * W, ua)
    call inflow_unique_sorted_global(msh%comm, zl, 1.0e-9_rp * zscale, uz)
    rec%ny = size(ua); rec%nz = size(uz)
    if (rec%ny < 2 .or. rec%nz < 2) error stop 'inflow record: the plane has fewer than 2 nodes in y or z'
    y0 = ua(1); dy = W / real(rec%ny, rp)
    do j = 1, rec%ny
      if (abs(ua(j) - (y0 + real(j - 1, rp) * dy)) > 1.0e-6_rp * dy) then
        if (rec%rank == 0) write(*,'(A)') 'ERROR: the cells on the recording plane are not uniformly spaced across the span.'
        error stop 'inflow record: spanwise spacing of the recording plane must be uniform'
      end if
    end do

    ! every (y, z) node must be hit exactly once
    allocate(rec%iy(n_loc), rec%iz(n_loc), hits(rec%ny, rec%nz), hits_g(rec%ny, rec%nz))
    hits = 0
    do i = 1, n_loc
      jy = minloc(abs(ua - al(i)), 1); kz = minloc(abs(uz - zl(i)), 1)
      rec%iy(i) = jy; rec%iz(i) = kz
      hits(jy, kz) = hits(jy, kz) + 1
    end do
    call MPI_Allreduce(hits, hits_g, int(rec%ny * rec%nz), MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    if (any(hits_g /= 1)) then
      if (rec%rank == 0) write(*,'(A,I0,A,I0,A)') 'ERROR: the recording plane is not a regular ', rec%ny, ' x ', rec%nz, &
        ' lattice of cells (a node is missing or doubly covered).'
      error stop 'inflow record: plane is not a regular y-z lattice of cells'
    end if

    rec%nt = int(nint(cfg%segment_time / cfg%dt), ip) + 1_ip
    rec%nseg = cfg%n_segments
    if (rec%nt < 8) error stop 'inflow record: record_segment_time must span at least 8 frames'
    rec%n_total = rec%nt * rec%nseg
    allocate(rec%P_prev(rec%ny, rec%nz, 3), rec%P_now(rec%ny, rec%nz, 3))

    if (rec%rank == 0) then
      hd%ny = rec%ny; hd%nz = rec%nz; hd%nt = rec%nt; hd%nseg = rec%nseg
      hd%W = W; hd%y0 = y0; hd%dt_p = cfg%dt
      hd%frame = rec%frame; hd%n_in = n_in; hd%t1 = t1
      hd%wind_dir_deg = modulo(atan2(-n_in(1), -n_in(2)) * 180.0_rp / PI, 360.0_rp)
      hd%generator = 'dopamine precursor recorder'
      allocate(hd%z(rec%nz), hd%U_mean(rec%nz), hd%R_target(6, rec%nz))
      hd%z = uz; hd%U_mean = 0.0_rp; hd%R_target = 0.0_rp
      call inflow_db_create(rec%db, trim(cfg%file), hd)
      allocate(rec%S1(rec%nz, 3), rec%S2(rec%nz, 6))
      rec%S1 = 0.0_rp; rec%S2 = 0.0_rp
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' Inflow record: ', rec%ny, ' x ', rec%nz, ' cell layer, ', rec%nseg, &
        ' segments of ', rec%nt, ' frames'
      write(*,'(A,ES9.2,A,ES9.2,A,F7.2,A)') ' Inflow record: dt = ', cfg%dt, ' s, from t = ', cfg%t_start, &
        ' s, wind_dir_deg of the database = ', hd%wind_dir_deg, ' (set the successor to match)'
    end if

    rec%n_next = max(0_ip, int(ceiling((solver%time - cfg%t_start) / cfg%dt - 1.0e-9_rp), ip))
    call sample_plane(rec, solver, rec%P_prev)
    rec%t_prev = solver%time
    rec%active = .true.
    rec%done = .false.
  end subroutine inflow_record_init

  !> Gather the layer's velocity onto the full lattice (profile frame), identical on all ranks.
  subroutine sample_plane(rec, solver, P)
    type(InflowRecorder), intent(in)    :: rec
    type(FVSolver),       intent(inout) :: solver
    real(rp),             intent(out)   :: P(:,:,:)

    integer(ip) :: i, ierr
    real(rp) :: u(3)
    real(rp), allocatable :: loc(:,:,:)

    !$acc update self(solver%U) if(solver%on_device)
    allocate(loc(rec%ny, rec%nz, 3))
    loc = 0.0_rp
    do i = 1, rec%n_cell
      u = solver%U(:, rec%cell(i))
      loc(rec%iy(i), rec%iz(i), 1) = dot_product(rec%frame(:, 1), u)
      loc(rec%iy(i), rec%iz(i), 2) = dot_product(rec%frame(:, 2), u)
      loc(rec%iy(i), rec%iz(i), 3) = dot_product(rec%frame(:, 3), u)
    end do
    call MPI_Allreduce(loc, P, int(rec%ny * rec%nz * 3), MPI_DOUBLE_PRECISION, MPI_SUM, solver%mesh%comm, ierr)
  end subroutine sample_plane

  !> Call after every completed step. Writes the frames that fall in the step just taken.
  subroutine inflow_record_step(rec, solver, dt_step)
    type(InflowRecorder), intent(inout) :: rec
    type(FVSolver),       intent(inout) :: solver
    real(rp),             intent(in)    :: dt_step

    real(rp) :: t_now, t_f, th
    real(rp), allocatable :: P(:,:,:)
    integer(ip) :: iseg, it, k, c

    if (.not. rec%active .or. rec%done) return
    if (dt_step > rec%cfg%dt .and. .not. rec%warned_dt) then
      rec%warned_dt = .true.
      if (rec%rank == 0) write(*,'(A,ES9.2,A,ES9.2,A)') ' WARNING: inflow record: solver step ', dt_step, &
        ' s exceeds record_dt = ', rec%cfg%dt, ' s; frames are linearly interpolated across a step.'
    end if
    t_now = solver%time
    call sample_plane(rec, solver, rec%P_now)

    allocate(P(rec%ny, rec%nz, 3))
    do while (rec%n_next < rec%n_total)
      t_f = rec%cfg%t_start + real(rec%n_next, rp) * rec%cfg%dt
      if (t_f > t_now + 1.0e-12_rp * max(1.0_rp, abs(t_now))) exit
      th = 0.0_rp
      if (t_now > rec%t_prev) th = min(max((t_f - rec%t_prev) / (t_now - rec%t_prev), 0.0_rp), 1.0_rp)
      P = (1.0_rp - th) * rec%P_prev + th * rec%P_now
      if (rec%rank == 0) then
        iseg = rec%n_next / rec%nt + 1_ip
        it   = mod(rec%n_next, rec%nt) + 1_ip
        call inflow_db_write_frame(rec%db, iseg, it, P)
        do c = 1, 3
          do k = 1, rec%nz
            rec%S1(k, c) = rec%S1(k, c) + sum(P(:, k, c))
          end do
        end do
        do k = 1, rec%nz
          rec%S2(k, 1) = rec%S2(k, 1) + sum(P(:, k, 1) * P(:, k, 1))
          rec%S2(k, 2) = rec%S2(k, 2) + sum(P(:, k, 2) * P(:, k, 2))
          rec%S2(k, 3) = rec%S2(k, 3) + sum(P(:, k, 3) * P(:, k, 3))
          rec%S2(k, 4) = rec%S2(k, 4) + sum(P(:, k, 1) * P(:, k, 2))
          rec%S2(k, 5) = rec%S2(k, 5) + sum(P(:, k, 1) * P(:, k, 3))
          rec%S2(k, 6) = rec%S2(k, 6) + sum(P(:, k, 2) * P(:, k, 3))
        end do
        rec%n_samp = rec%n_samp + real(rec%ny, rp)
      end if
      rec%n_next = rec%n_next + 1
    end do
    rec%P_prev = rec%P_now
    rec%t_prev = t_now

    if (rec%n_next >= rec%n_total) call inflow_record_finalize(rec)
  end subroutine inflow_record_step

  !> Remove the mean, store the measured profiles, decorrelate the spanwise mean, close the file.
  !> Called automatically when all frames are written; calling it earlier leaves an incomplete file.
  subroutine inflow_record_finalize(rec)
    type(InflowRecorder), intent(inout) :: rec

    real(rp), allocatable :: Ub(:,:), R(:,:), Um(:), plane(:,:,:)
    integer(ip) :: is, it, c, k
    real(rp) :: frac

    if (.not. rec%active .or. rec%done) return
    rec%done = .true.
    if (rec%rank /= 0) return

    if (rec%n_next < rec%n_total) then
      write(*,'(A,I0,A,I0,A)') ' WARNING: inflow record incomplete (', rec%n_next, ' of ', rec%n_total, &
        ' frames); the file is not usable. Extend the run or shorten record_segment_time.'
      call inflow_db_close(rec%db)
      return
    end if

    allocate(Ub(rec%nz, 3), R(6, rec%nz), Um(rec%nz), plane(rec%ny, rec%nz, 3))
    Ub = rec%S1 / rec%n_samp
    do k = 1, rec%nz
      R(1, k) = rec%S2(k, 1) / rec%n_samp - Ub(k, 1)**2
      R(2, k) = rec%S2(k, 2) / rec%n_samp - Ub(k, 2)**2
      R(3, k) = rec%S2(k, 3) / rec%n_samp - Ub(k, 3)**2
      R(4, k) = rec%S2(k, 4) / rec%n_samp - Ub(k, 1) * Ub(k, 2)
      R(5, k) = rec%S2(k, 5) / rec%n_samp - Ub(k, 1) * Ub(k, 3)
      R(6, k) = rec%S2(k, 6) / rec%n_samp - Ub(k, 2) * Ub(k, 3)
    end do
    Um = Ub(:, 1)

    do is = 1, rec%nseg
      do it = 1, rec%nt
        call inflow_db_read_frame(rec%db, is, it, plane)
        do c = 1, 3
          do k = 1, rec%nz
            plane(:, k, c) = plane(:, k, c) - Ub(k, c)
          end do
        end do
        call inflow_db_write_frame(rec%db, is, it, plane)
      end do
    end do
    call inflow_db_write_profiles(rec%db, Um, R)
    call inflow_db_close(rec%db)
    write(*,'(A)') ' Inflow record: mean removed (fluctuation-only database), measured stresses stored'
    if (rec%nseg >= 2) then
      call inflow_db_decorrelate_span_mean(trim(rec%file), frac)
      write(*,'(A,F6.2,A)') ' Inflow record: removed the segment-mean spanwise-mean waveform (', 100.0_rp * frac, &
        '% of the plane variance)'
      if (rec%nseg < 4) write(*,'(A)') ' WARNING: fewer than 4 recorded segments; use record_n_segments >= 4.'
    else
      write(*,'(A)') ' WARNING: a single recorded segment cannot be decorrelated; replay will repeat its spanwise-mean content.'
    end if
    write(*,'(A)') ' Inflow record: wrote '//trim(rec%file)
  end subroutine inflow_record_finalize

end module mod_inflow_record
