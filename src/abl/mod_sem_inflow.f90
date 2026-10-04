! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Synthetic Eddy Method (SEM) for correlated turbulent inflow generation
module mod_sem_inflow
  use mod_kinds
  use mod_parameters, only: BC_ABL_INLET, BC_PERIODIC
  use mod_abl_types
  use mod_abl_profiles
  use mod_mesh
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  real(rp), parameter :: PI = 3.14159265358979323846_rp
  real(rp), parameter :: SMALL = 1.0e-16_rp

  ! Tent kernel max(0,1-r) needs this factor so that integral f^2 d3y = 1.
  real(rp), parameter :: TENT_NORM = 1.5450968080927583_rp  ! sqrt(15/(2*pi))

  ! Divergence-free eddy of Poletto et al. (2013), Eq. 8: u_b = sigma_b (1 - d^2) (r x alpha)_b with
  ! r_j = (x_j - x_j^k) / sigma_j, d = |r| < 1. DF_I = int (1-d^2)^2 r_2^2 d3r over the unit ball.
  real(rp), parameter :: DF_I = 32.0_rp * PI / 945.0_rp
  ! Smallest eddy half-extent relative to the largest (guards a vanishing principal stress).
  real(rp), parameter :: DF_FMIN = 0.2_rp

  integer(int64), parameter :: RNG_SEED = 88172645463325252_int64
  integer(int64), save      :: rng_state = RNG_SEED

  public :: sem_init, sem_advance, sem_velocity_at_face
  public :: sem_get_rng_state, sem_set_rng_state
  public :: sem_read_profile, sem_profile_at, sem_cholesky
  public :: sem_prepare, sem_eigen3, sem_dfsem_eddy_velocity, sem_reference_scales, sem_stress_at
  public :: sem_init_targets, sem_setup, sem_target_at, sem_mean_at

contains

  !> Random-number state, for checkpoints: the eddy field is a pure function of (pos, amp, rng_state).
  integer(int64) function sem_get_rng_state() result(st)
    st = rng_state
  end function sem_get_rng_state

  subroutine sem_set_rng_state(st)
    integer(int64), intent(in) :: st
    rng_state = st
  end subroutine sem_set_rng_state

  !> Portable xorshift64 uniform in [0,1). Bit ops only, so no reliance on integer overflow.
  !> Every rank runs the identical call sequence, keeping the eddy field consistent.
  real(rp) function sem_rand() result(u)
    rng_state = ieor(rng_state, ishft(rng_state,  13))
    rng_state = ieor(rng_state, ishft(rng_state,  -7))
    rng_state = ieor(rng_state, ishft(rng_state,  17))
    u = real(iand(rng_state, 4503599627370495_int64), rp) / 4503599627370496.0_rp
  end function sem_rand

  !> Build the inlet-plane frame from all BC_ABL_INLET faces, reduced across ranks.
  !> Returns found = .false. when no rank owns an ABL inlet face.
  subroutine sem_inlet_frame(msh, state, found)
    type(Mesh),     intent(in)    :: msh
    type(SEMState), intent(inout) :: state
    logical,        intent(out)   :: found

    integer(ip) :: p, f, f_global, ierr
    real(rp)    :: nsum(3), xsum(3), lbuf(7), gbuf(7), area, area_sum
    real(rp)    :: nlen, d(3), s1, s2, lo(2), hi(2), glo(2), ghi(2)

    nsum = 0.0_rp; xsum = 0.0_rp; area_sum = 0.0_rp
    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, msh%patches(p)%n_faces
        f_global = msh%n_int_faces + msh%patches(p)%start_face + f - 1
        if (f_global < 1 .or. f_global > size(msh%faces)) cycle
        area     = msh%faces(f_global)%area
        nsum     = nsum + area * msh%faces(f_global)%normal
        xsum     = xsum + area * msh%faces(f_global)%centre
        area_sum = area_sum + area
      end do
    end do

    lbuf(1:3) = nsum; lbuf(4:6) = xsum; lbuf(7) = area_sum
    call MPI_Allreduce(lbuf, gbuf, 7, MPI_DOUBLE_PRECISION, MPI_SUM, msh%comm, ierr)

    found = (gbuf(7) > SMALL)
    if (.not. found) return

    ! Inward normal is the convection direction; faces store the outward normal.
    nlen = sqrt(sum(gbuf(1:3)**2))
    if (nlen <= SMALL) then
      found = .false.
      return
    end if
    state%n_in    = -gbuf(1:3) / nlen
    state%x_inlet =  gbuf(4:6) / gbuf(7)

    ! In-plane basis; e_z is a safe seed because ABL inlets have a horizontal normal.
    if (abs(state%n_in(3)) < 0.9_rp) then
      state%t1 = cross3([0.0_rp, 0.0_rp, 1.0_rp], state%n_in)
    else
      state%t1 = cross3([1.0_rp, 0.0_rp, 0.0_rp], state%n_in)
    end if
    state%t1 = state%t1 / max(sqrt(sum(state%t1**2)), SMALL)
    state%t2 = cross3(state%n_in, state%t1)

    ! In-plane extent of the inlet, again reduced across ranks.
    lo =  huge(1.0_rp)
    hi = -huge(1.0_rp)
    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, msh%patches(p)%n_faces
        f_global = msh%n_int_faces + msh%patches(p)%start_face + f - 1
        if (f_global < 1 .or. f_global > size(msh%faces)) cycle
        d  = msh%faces(f_global)%centre - state%x_inlet
        s1 = dot_product(d, state%t1)
        s2 = dot_product(d, state%t2)
        lo(1) = min(lo(1), s1); hi(1) = max(hi(1), s1)
        lo(2) = min(lo(2), s2); hi(2) = max(hi(2), s2)
      end do
    end do
    call MPI_Allreduce(lo, glo, 2, MPI_DOUBLE_PRECISION, MPI_MIN, msh%comm, ierr)
    call MPI_Allreduce(hi, ghi, 2, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)

    state%t1_min = glo(1); state%t1_max = ghi(1)
    state%t2_min = glo(2); state%t2_max = ghi(2)

    call sem_detect_periodic(msh, state)
  end subroutine sem_inlet_frame

  !> A periodic patch pair whose faces are normal to t1 (or t2) makes the inlet plane periodic along it.
  !> The period is the distance between the two planes of periodic faces, reduced across ranks.
  subroutine sem_detect_periodic(msh, state)
    type(Mesh),     intent(in)    :: msh
    type(SEMState), intent(inout) :: state

    integer(ip) :: p, f, f_global, ierr, k
    real(rp)    :: s, tk(3), lo(2), hi(2), glo(2), ghi(2)

    state%per1 = .false.; state%per2 = .false.
    state%period1 = 0.0_rp; state%period2 = 0.0_rp
    lo =  huge(1.0_rp)
    hi = -huge(1.0_rp)
    do k = 1, 2
      tk = merge(state%t1, state%t2, k == 1)
      do p = 1, msh%n_patches
        if (msh%patches(p)%bc_id /= BC_PERIODIC) cycle
        do f = 1, msh%patches(p)%n_faces
          f_global = msh%n_int_faces + msh%patches(p)%start_face + f - 1
          if (f_global < 1 .or. f_global > size(msh%faces)) cycle
          if (abs(dot_product(msh%faces(f_global)%normal, tk)) < 0.99_rp) cycle
          s = dot_product(msh%faces(f_global)%centre - state%x_inlet, tk)
          lo(k) = min(lo(k), s); hi(k) = max(hi(k), s)
        end do
      end do
    end do
    call MPI_Allreduce(lo, glo, 2, MPI_DOUBLE_PRECISION, MPI_MIN, msh%comm, ierr)
    call MPI_Allreduce(hi, ghi, 2, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)

    ! Both planes must exist and enclose the inlet; a lone periodic patch is not a pair.
    if (ghi(1) - glo(1) > SMALL .and. glo(1) <= state%t1_min + SMALL .and. ghi(1) >= state%t1_max - SMALL) then
      state%per1 = .true.; state%period1 = ghi(1) - glo(1); state%per1_lo = glo(1)
    end if
    if (ghi(2) - glo(2) > SMALL .and. glo(2) <= state%t2_min + SMALL .and. ghi(2) >= state%t2_max - SMALL) then
      state%per2 = .true.; state%period2 = ghi(2) - glo(2); state%per2_lo = glo(2)
    end if
  end subroutine sem_detect_periodic

  pure function cross3(a, b) result(c)
    real(rp), intent(in) :: a(3), b(3)
    real(rp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  !> Initialize SEM: seed eddies in a recycle box one sigma upstream of the inlet plane.
  subroutine sem_init(state, abl_cfg, msh)
    type(SEMState),   intent(inout) :: state
    type(ABLConfig),  intent(in)    :: abl_cfg
    type(Mesh),       intent(in)    :: msh

    logical :: found

    call sem_init_targets(state, abl_cfg, msh, found)
    if (.not. found) then
      state%initialized = .false.
      return
    end if

    call sem_setup(state, abl_cfg, verbose=(msh%rank == 0))
    if (state%dfsem) call sem_check_resolution(state, msh)
  end subroutine sem_init

  !> Profile table (if any) and inlet frame, no eddies: for a source that supplies its own fluctuation
  !> (the inflow-plane database) but still takes the mean profile from the SEM laws and tables.
  subroutine sem_init_targets(state, abl_cfg, msh, found)
    type(SEMState),   intent(inout) :: state
    type(ABLConfig),  intent(in)    :: abl_cfg
    type(Mesh),       intent(in)    :: msh
    logical,          intent(out)   :: found

    if (len_trim(abl_cfg%sem_profile_file) > 0 .and. .not. state%has_profile) &
      call sem_read_profile(abl_cfg%sem_profile_file, state, verbose=(msh%rank == 0), &
                            fmt=abl_cfg%sem_profile_format, rho_uw=abl_cfg%sem_rho_uw)
    call sem_inlet_frame(msh, state, found)
    if (.not. found) return
    state%override_mean = state%has_profile .or. trim(abl_cfg%sem_mean_law) == 'power'
  end subroutine sem_init_targets

  !> Everything after the inlet frame is known: eddy scale, box, seeding. Needs state%{x_inlet,n_in,t1,t2,
  !> t1_min..t2_max, per1/per2 with period and lower plane} and, for tabulated inflow, the profile.
  subroutine sem_setup(state, abl_cfg, verbose)
    type(SEMState),   intent(inout) :: state
    type(ABLConfig),  intent(in)    :: abl_cfg
    logical, intent(in), optional   :: verbose

    integer(ip) :: i, n_e
    real(rp)    :: span1, span2, r
    logical     :: talk

    talk = .false.
    if (present(verbose)) talk = verbose
    state%override_mean = state%has_profile .or. trim(abl_cfg%sem_mean_law) == 'power'

    state%n_eddies = max(1_ip, abl_cfg%n_eddies)
    n_e = state%n_eddies
    if (allocated(state%pos))   deallocate(state%pos)
    if (allocated(state%amp))   deallocate(state%amp)
    if (allocated(state%frame)) deallocate(state%frame)
    if (allocated(state%sig))   deallocate(state%sig)
    if (allocated(state%ampl))  deallocate(state%ampl)
    allocate(state%pos(3, n_e), source=0.0_rp)
    allocate(state%amp(3, n_e), source=0.0_rp)

    span1 = state%t1_max - state%t1_min
    span2 = state%t2_max - state%t2_min

    ! Eddy scale defaults to a tenth of the smaller in-plane extent.
    if (abl_cfg%sem_sigma <= SMALL) then
      state%sigma = 0.1_rp * min(span1, span2)
    else
      state%sigma = abl_cfg%sem_sigma
    end if
    if (state%sigma <= SMALL) state%sigma = 0.1_rp * abl_cfg%z_ref

    ! Divergence-free eddies are ellipsoidal (half-extents follow the principal stresses, the largest
    ! being sigma); rsup bounds the support half-extent along any direction and replaces sigma for the box.
    state%dfsem = (trim(abl_cfg%sem_type) /= 'jarrin')
    state%rsup  = state%sigma
    if (state%dfsem) then
      allocate(state%frame(9, n_e), source=0.0_rp)
      allocate(state%sig(3, n_e),   source=0.0_rp)
      allocate(state%ampl(3, n_e),  source=0.0_rp)
      call sem_reference_scales(state, abl_cfg)
    end if

    ! A periodic direction needs exactly one period of eddies (they are wrapped when evaluated, so
    ! the fluctuation is periodic across the seam) and a kernel that fits inside half of it.
    if ((state%per1 .and. state%rsup >= 0.5_rp * state%period1) .or. &
        (state%per2 .and. state%rsup >= 0.5_rp * state%period2)) then
      if (talk) write(*,'(A,ES10.3,A)') 'ERROR: the SEM eddy support (', state%rsup, &
        ' m) must be below half the periodic width of the inlet plane; reduce sem_sigma'
      error stop 'sem_setup: eddy scale too large for the periodic inlet'
    end if
    ! Otherwise pad the box by the support so eddies straddling the edges are represented.
    if (state%per1) then
      state%t1_min = state%per1_lo; state%t1_max = state%per1_lo + state%period1
    else
      state%t1_min = state%t1_min - state%rsup; state%t1_max = state%t1_max + state%rsup
    end if
    if (state%per2) then
      state%t2_min = state%per2_lo; state%t2_max = state%per2_lo + state%period2
    else
      state%t2_min = state%t2_min - state%rsup; state%t2_max = state%t2_max + state%rsup
    end if
    ! Box spans -rsup to +rsup about the plane so a face there sees every eddy that can touch it.
    state%box_volume = 2.0_rp * state%rsup &
                     * (state%t1_max - state%t1_min) * (state%t2_max - state%t2_min)

    state%u_convect = max(abl_cfg%u_ref, SMALL)

    rng_state = RNG_SEED
    if (abl_cfg%sem_seed /= 0) then
      rng_state = ieor(RNG_SEED, int(abl_cfg%sem_seed, int64) * 2654435761_int64)
      if (rng_state == 0_int64) rng_state = RNG_SEED
      do i = 1, 8
        r = sem_rand()  ! warm up so nearby seeds decorrelate
      end do
    end if
    do i = 1, n_e
      call sem_place_eddy(state, i, state%rsup * (2.0_rp * sem_rand() - 1.0_rp), .true.)
    end do
    call sem_balance_amplitudes(state)
    if (state%dfsem) call sem_prepare(state, abl_cfg)

    state%initialized = .true.
  end subroutine sem_setup

  !> Force each amplitude component to sum to zero so the eddies add no mean velocity.
  subroutine sem_balance_amplitudes(state)
    type(SEMState), intent(inout) :: state
    integer(ip) :: k, i
    real(rp)    :: s_sum

    do k = 1, 3
      s_sum = sum(state%amp(k, 1:state%n_eddies))
      do i = 1, state%n_eddies
        if (abs(s_sum) < 1.0_rp) exit
        ! Flip members of the majority sign until the running sum is balanced.
        if (state%amp(k, i) * s_sum > 0.0_rp) then
          state%amp(k, i) = -state%amp(k, i)
          s_sum = s_sum + 2.0_rp * state%amp(k, i)
        end if
      end do
    end do
  end subroutine sem_balance_amplitudes

  !> Place eddy i at signed distance s along the inward normal from the inlet plane.
  !> Amplitude signs are redrawn only at seeding, so the zero-mean balance survives recycling.
  subroutine sem_place_eddy(state, i, s, new_amp)
    type(SEMState), intent(inout) :: state
    integer(ip),    intent(in)    :: i
    real(rp),       intent(in)    :: s
    logical,        intent(in)    :: new_amp

    real(rp) :: a, b
    integer(ip) :: k

    a = state%t1_min + (state%t1_max - state%t1_min) * sem_rand()
    b = state%t2_min + (state%t2_max - state%t2_min) * sem_rand()

    state%pos(:, i) = state%x_inlet + s * state%n_in + a * state%t1 + b * state%t2
    if (new_amp) then
      do k = 1, 3
        state%amp(k, i) = sign(1.0_rp, sem_rand() - 0.5_rp)
      end do
    end if
  end subroutine sem_place_eddy

  !> Advance SEM: convect eddies along the inward normal and recycle those leaving the box.
  subroutine sem_advance(state, dt)
    type(SEMState), intent(inout) :: state
    real(rp),       intent(in)    :: dt

    integer(ip) :: i
    real(rp)    :: displacement, s

    if (.not. state%initialized) return

    displacement = state%u_convect * dt

    do i = 1, state%n_eddies
      state%pos(:, i) = state%pos(:, i) + displacement * state%n_in
      s = dot_product(state%pos(:, i) - state%x_inlet, state%n_in)
      if (s > state%rsup) call sem_place_eddy(state, i, -state%rsup, .false.)
    end do
  end subroutine sem_advance

  !> Count whitespace/comma separated tokens on a line.
  integer(ip) function count_tokens(line) result(n)
    character(len=*), intent(in) :: line
    integer(ip) :: i
    logical :: in_tok
    character :: ch

    n = 0
    in_tok = .false.
    do i = 1, len_trim(line)
      ch = line(i:i)
      if (ch == ' ' .or. ch == ',' .or. ch == achar(9) .or. ch == ';') then
        in_tok = .false.
      else if (.not. in_tok) then
        in_tok = .true.
        n = n + 1
      end if
    end do
  end function count_tokens

  !> Read prescribed inflow profiles (heights need not be sorted; '#' comments allowed).
  !>
  !> fmt = 'stress' (default):  z U uu vv ww uw            6 columns, wind-aligned
  !>                            z U V W uu vv ww uv uw vw  10 columns, global axes
  !> fmt = 'ti':                z U Iu Iv Iw [rho_uw]      5-6 columns, wind-aligned; I = sigma/U
  !> Length-scale columns are not used: the eddy sizes follow from sem_sigma and the principal stresses.
  subroutine sem_read_profile(filename, state, verbose, fmt, rho_uw)
    character(len=*), intent(in)    :: filename
    type(SEMState),   intent(inout) :: state
    logical,          intent(in), optional :: verbose  ! rank-0 only in parallel runs
    character(len=*), intent(in), optional :: fmt
    real(rp),         intent(in), optional :: rho_uw   ! correlation used when a 'ti' table has none

    integer(ip) :: iunit, ios, k, n, ncol
    character(len=512) :: line
    real(rp) :: v(10), tmp, rho_def, su, sv, sw, rho
    logical  :: is_ti, ok

    is_ti = .false.
    if (present(fmt)) is_ti = (trim(fmt) == 'ti')
    rho_def = -0.3_rp
    if (present(rho_uw)) rho_def = rho_uw

    open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'ERROR: cannot open SEM profile file: '//trim(filename)
      error stop 'sem_read_profile: open failed'
    end if

    ! First pass: settle the column count from the first numeric row, then count valid rows.
    n = 0
    ncol = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      if (ncol == 0) then
        read(line, *, iostat=ios) tmp
        if (ios /= 0) cycle  ! header row
        ncol = count_tokens(line)
        if (is_ti) then
          ok = (ncol == 5 .or. ncol == 6)
        else
          ok = (ncol == 6 .or. ncol == 10)
        end if
        if (.not. ok) then
          write(*,'(A,I0,A)') 'ERROR: SEM profile has ', ncol, ' columns, not valid for sem_profile_format = ' &
            //merge("'ti'    ", "'stress'", is_ti)//' in '//trim(filename)
          error stop 'sem_read_profile: bad column count'
        end if
      end if
      read(line, *, iostat=ios) v(1:ncol)
      if (ios == 0) n = n + 1
    end do

    if (n < 2) then
      write(*,'(A)') 'ERROR: SEM profile needs at least 2 levels: '//trim(filename)
      error stop 'sem_read_profile: too few levels'
    end if

    state%n_levels = n
    state%prof_global_axes = (.not. is_ti) .and. (ncol == 10)
    if (allocated(state%prof_z)) deallocate(state%prof_z, state%prof_U, state%prof_R)
    allocate(state%prof_z(n), source=0.0_rp)
    allocate(state%prof_U(3, n), source=0.0_rp)
    allocate(state%prof_R(6, n), source=0.0_rp)

    rewind(iunit)
    k = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) v(1:ncol)
      if (ios /= 0) cycle
      k = k + 1
      if (k > n) exit
      state%prof_z(k) = v(1)
      if (is_ti) then
        state%prof_U(1, k) = v(2)
        su = v(3) * v(2); sv = v(4) * v(2); sw = v(5) * v(2)
        rho = rho_def
        if (ncol >= 6) rho = v(6)
        state%prof_R(1, k) = su**2; state%prof_R(2, k) = sv**2; state%prof_R(3, k) = sw**2
        state%prof_R(5, k) = rho * su * sw
      else if (ncol == 6) then
        state%prof_U(1, k)   = v(2)  ! streamwise magnitude, rotated later
        state%prof_R(1:3, k) = v(3:5)
        state%prof_R(5, k)   = v(6)  ! uw
      else
        state%prof_U(:, k) = v(2:4)
        state%prof_R(:, k) = v(5:10)
      end if
    end do
    close(iunit)

    call sem_sort_profile(state)
    state%has_profile = .true.
    if (present(verbose)) then
      if (verbose) write(*,'(A,I0,A,A,A,I0,A,ES10.3,A,ES10.3,A)') &
        ' SEM inflow profiles: ', n, ' levels, format ', merge('ti    ', 'stress', is_ti), ', ', ncol, &
        ' columns, z = ', state%prof_z(1), ' to ', state%prof_z(n), ' m'
    end if
  end subroutine sem_read_profile

  !> Insertion sort by height; tables are short and often already ascending.
  subroutine sem_sort_profile(state)
    type(SEMState), intent(inout) :: state

    integer(ip) :: i, j
    real(rp) :: zt, ut(3), rt(6)

    do i = 2, state%n_levels
      zt = state%prof_z(i); ut = state%prof_U(:, i); rt = state%prof_R(:, i)
      j = i - 1
      do while (j >= 1)
        if (state%prof_z(j) <= zt) exit
        state%prof_z(j+1)    = state%prof_z(j)
        state%prof_U(:, j+1) = state%prof_U(:, j)
        state%prof_R(:, j+1) = state%prof_R(:, j)
        j = j - 1
      end do
      state%prof_z(j+1) = zt; state%prof_U(:, j+1) = ut; state%prof_R(:, j+1) = rt
    end do
  end subroutine sem_sort_profile

  !> Linear interpolation in z, clamped to the table ends.
  subroutine sem_profile_at(state, z, u_mean, R)
    type(SEMState), intent(in)  :: state
    real(rp),       intent(in)  :: z
    real(rp),       intent(out) :: u_mean(3), R(6)

    integer(ip) :: k, n
    real(rp) :: w

    n = state%n_levels
    if (z <= state%prof_z(1)) then
      u_mean = state%prof_U(:, 1); R = state%prof_R(:, 1); return
    else if (z >= state%prof_z(n)) then
      u_mean = state%prof_U(:, n); R = state%prof_R(:, n); return
    end if

    k = 1
    do while (k < n - 1 .and. state%prof_z(k+1) < z)
      k = k + 1
    end do
    w = (z - state%prof_z(k)) / max(state%prof_z(k+1) - state%prof_z(k), SMALL)
    u_mean = (1.0_rp - w) * state%prof_U(:, k) + w * state%prof_U(:, k+1)
    R      = (1.0_rp - w) * state%prof_R(:, k) + w * state%prof_R(:, k+1)
  end subroutine sem_profile_at

  !> Cholesky factor of the full symmetric Reynolds tensor R = [uu,vv,ww,uv,uw,vw].
  !> Non-positive-definite input is clipped rather than failing the run.
  pure function sem_cholesky(R) result(L)
    real(rp), intent(in) :: R(6)
    real(rp) :: L(3,3)

    L = 0.0_rp
    L(1,1) = sqrt(max(R(1), SMALL))
    L(2,1) = R(4) / L(1,1)
    L(3,1) = R(5) / L(1,1)
    L(2,2) = sqrt(max(R(2) - L(2,1)**2, SMALL))
    L(3,2) = (R(6) - L(3,1) * L(2,1)) / L(2,2)
    L(3,3) = sqrt(max(R(3) - L(3,1)**2 - L(3,2)**2, SMALL))
  end function sem_cholesky

  !> Target mean velocity (global axes) and stresses R = [uu,vv,ww,uv,uw,vw] at height z. R is in the
  !> mean-wind frame (e_s, e_n, e_z), except for 10-column tables, which are already global.
  subroutine sem_target_at(state, cfg, z, u_mean, R)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: cfg
    real(rp),        intent(in)  :: z
    real(rp),        intent(out) :: u_mean(3), R(6)

    real(rp) :: u_star, umag, zz, zi, Iu, su, sv, sw, z_floor

    if (state%has_profile) then
      call sem_profile_at(state, z, u_mean, R)
      ! 6/5-column tables give a streamwise speed; 10-column ones are already global.
      if (.not. state%prof_global_axes) u_mean = abl_wind_vector(u_mean(1), cfg%wind_dir_deg)
      return
    end if

    z_floor = max(cfg%z0_m, 1.0e-6_rp)
    zz      = max(z, z_floor)
    u_star  = abl_ustar_from_ref(cfg%u_ref, cfg%z_ref, cfg%z0_m)
    if (trim(cfg%sem_mean_law) == 'power') then
      umag = cfg%u_ref * (zz / max(cfg%z_ref, SMALL)) ** cfg%sem_alpha
    else
      umag = abl_u_profile(z, u_star, cfg%z0_m, 0.0_rp)
    end if
    u_mean = abl_wind_vector(umag, cfg%wind_dir_deg)

    R = 0.0_rp
    if (cfg%sem_Iu_ref > 0.0_rp) then
      ! Wind-tunnel style: intensity decays like the log law, Iu = Iu_ref ln(z_ref/z0)/ln(z/z0),
      ! evaluated no lower than 10 z0 where the log law stops being meaningful.
      zi = max(z, 10.0_rp * cfg%z0_m)
      Iu = cfg%sem_Iu_ref * log(max(cfg%z_ref, 1.0001_rp * cfg%z0_m) / cfg%z0_m) / log(zi / cfg%z0_m)
      if (trim(cfg%sem_mean_law) == 'power') then
        umag = cfg%u_ref * (max(zi, z_floor) / max(cfg%z_ref, SMALL)) ** cfg%sem_alpha
      else
        umag = abl_u_profile(zi, u_star, cfg%z0_m, 0.0_rp)
      end if
      su = Iu * umag
      sv = cfg%sem_Iv_ratio * su
      sw = cfg%sem_Iw_ratio * su
      R(1) = su**2; R(2) = sv**2; R(3) = sw**2
      R(5) = cfg%sem_rho_uw * su * sw
    else
      ! Panofsky & Dutton neutral surface layer: sigma_u,v,w / u* = 2.4, 1.9, 1.25 and a constant
      ! momentum flux <u'w'> = -u*^2.
      R(1) = (2.4_rp * u_star)**2
      R(2) = (1.9_rp * u_star)**2
      R(3) = (1.25_rp * u_star)**2
      R(5) = -u_star**2
    end if
  end subroutine sem_target_at

  !> Mean velocity at height z (global axes).
  subroutine sem_mean_at(state, abl_cfg, z, u_mean)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: z
    real(rp),        intent(out) :: u_mean(3)

    real(rp) :: R6(6)

    call sem_target_at(state, abl_cfg, z, u_mean, R6)
  end subroutine sem_mean_at

  !> Target Reynolds tensor at height z as a symmetric matrix Rm in the stress frame, and the matrix Wr
  !> whose columns are the stress-frame axes in global coordinates (R_global = Wr Rm Wr^T). The stress
  !> frame is the mean-wind frame (e_s, e_n, e_z) except for 10-column tables, which are already global.
  subroutine sem_stress_at(state, abl_cfg, z, Rm, Wr)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: z
    real(rp),        intent(out) :: Rm(3,3), Wr(3,3)

    real(rp) :: R6(6), u_mean(3), wind_rad

    call sem_target_at(state, abl_cfg, z, u_mean, R6)
    Rm(1,1) = R6(1); Rm(2,2) = R6(2); Rm(3,3) = R6(3)
    Rm(1,2) = R6(4); Rm(2,1) = R6(4)
    Rm(1,3) = R6(5); Rm(3,1) = R6(5)
    Rm(2,3) = R6(6); Rm(3,2) = R6(6)

    Wr = 0.0_rp
    if (state%has_profile .and. state%prof_global_axes) then
      Wr(1,1) = 1.0_rp; Wr(2,2) = 1.0_rp; Wr(3,3) = 1.0_rp
    else
      wind_rad = abl_cfg%wind_dir_deg * PI / 180.0_rp
      Wr(:,1) = [-sin(wind_rad), -cos(wind_rad), 0.0_rp]  ! e_s along the mean wind
      Wr(:,2) = [ cos(wind_rad), -sin(wind_rad), 0.0_rp]  ! e_n horizontal, normal to it
      Wr(:,3) = [0.0_rp, 0.0_rp, 1.0_rp]
    end if
  end subroutine sem_stress_at

  !> Eigen-decomposition of a symmetric 3x3 matrix by cyclic Jacobi: A = V diag(lam) V^T, columns of V
  !> orthonormal. Eigenvalues are returned in the order they converge (not sorted).
  pure subroutine sem_eigen3(A_in, lam, V)
    real(rp), intent(in)  :: A_in(3,3)
    real(rp), intent(out) :: lam(3), V(3,3)

    real(rp)    :: A(3,3), off, th, t, c, sn, tau, apq
    integer(ip) :: sweep, p, q, k
    real(rp)    :: akp, akq, vkp, vkq

    A = A_in
    V = 0.0_rp
    V(1,1) = 1.0_rp; V(2,2) = 1.0_rp; V(3,3) = 1.0_rp
    do sweep = 1, 30
      off = abs(A(1,2)) + abs(A(1,3)) + abs(A(2,3))
      if (off <= 1.0e-15_rp * (abs(A(1,1)) + abs(A(2,2)) + abs(A(3,3)) + SMALL)) exit
      do p = 1, 2
        do q = p + 1, 3
          apq = A(p,q)
          if (abs(apq) <= SMALL) cycle
          th  = (A(q,q) - A(p,p)) / (2.0_rp * apq)
          t   = sign(1.0_rp, th) / (abs(th) + sqrt(th * th + 1.0_rp))
          c   = 1.0_rp / sqrt(t * t + 1.0_rp)
          sn  = t * c
          tau = sn / (1.0_rp + c)
          A(p,p) = A(p,p) - t * apq
          A(q,q) = A(q,q) + t * apq
          A(p,q) = 0.0_rp; A(q,p) = 0.0_rp
          do k = 1, 3
            if (k /= p .and. k /= q) then
              akp = A(k,p); akq = A(k,q)
              A(k,p) = akp - sn * (akq + tau * akp); A(p,k) = A(k,p)
              A(k,q) = akq + sn * (akp - tau * akq); A(q,k) = A(k,q)
            end if
            vkp = V(k,p); vkq = V(k,q)
            V(k,p) = vkp - sn * (vkq + tau * vkp)
            V(k,q) = vkq + sn * (vkp - tau * vkq)
          end do
        end do
      end do
    end do
    lam = [A(1,1), A(2,2), A(3,3)]
  end subroutine sem_eigen3

  !> Largest principal stress over the inflow (sets the eddy sizes). Eddy half-extents are
  !> sigma_i = sigma * max(DF_FMIN, sqrt(lam_i / lam_ref)): the longest eddy is along the strongest principal
  !> stress, which keeps the eddy intensities real for any positive-definite tensor. The ellipsoidal support
  !> reaches at most sigma along any direction, so rsup = sigma.
  subroutine sem_reference_scales(state, abl_cfg)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: abl_cfg

    integer(ip) :: k, nl
    real(rp)    :: Rm(3,3), Wr(3,3), lam(3), V(3,3), lmax, fmin, zlo, zhi
    real(rp), allocatable :: zl(:)

    ! Heights to scan: the table levels, or (analytic laws) the inlet height range, where the intensity law can
    ! vary by a large factor; the largest principal stress over all of them must set sigma_i <= sigma.
    if (state%has_profile) then
      nl = state%n_levels
      allocate(zl(nl))
      zl = state%prof_z(1:nl)
    else
      nl = 65
      allocate(zl(nl))
      zlo = state%x_inlet(3) + state%t2_min
      zhi = state%x_inlet(3) + state%t2_max
      do k = 1, nl - 1
        zl(k) = zlo + (zhi - zlo) * real(k - 1, rp) / real(nl - 2, rp)
      end do
      zl(nl) = abl_cfg%z_ref
    end if

    lmax = SMALL
    do k = 1, nl
      call sem_stress_at(state, abl_cfg, zl(k), Rm, Wr)
      call sem_eigen3(Rm, lam, V)
      lmax = max(lmax, maxval(lam))
    end do
    state%lam_ref = lmax
    state%rsup    = state%sigma
    ! Smallest half-extent where the inflow is energetic (levels with a principal stress above a quarter of the
    ! largest); weak-stress levels sit at the DF_FMIN floor and carry little pressure.
    fmin = 1.0_rp
    do k = 1, nl
      call sem_stress_at(state, abl_cfg, zl(k), Rm, Wr)
      call sem_eigen3(Rm, lam, V)
      if (maxval(lam) >= 0.25_rp * lmax) fmin = min(fmin, max(DF_FMIN, sqrt(max(minval(lam), 0.0_rp) / lmax)))
    end do
    state%sig_min = state%sigma * fmin
    deallocate(zl)
  end subroutine sem_reference_scales

  !> Warn when the smallest eddy is poorly resolved by the inlet cells. The first cell's convective term has an
  !> O(dx) closure error at the inlet face (the interior's two face errors cancel, the exact inlet value does
  !> not), which the pressure absorbs: roughly p ~ ((k dx)^2 / 8) U u', i.e. (dx/sigma)^1.5-2.
  subroutine sem_check_resolution(state, msh)
    type(SEMState), intent(in) :: state
    type(Mesh),     intent(in) :: msh

    integer(ip) :: p, f, f_global, o, ierr
    real(rp)    :: h_loc, h_glob, area, h

    h_loc = huge(1.0_rp)
    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_ABL_INLET) cycle
      do f = 1, msh%patches(p)%n_faces
        f_global = msh%n_int_faces + msh%patches(p)%start_face + f - 1
        if (f_global < 1 .or. f_global > size(msh%faces)) cycle
        area = msh%faces(f_global)%area
        o    = msh%faces(f_global)%owner
        h    = sqrt(area)  ! in-plane spacing
        if (o >= 1 .and. o <= msh%n_cells) h = max(h, msh%cells(o)%volume / max(area, SMALL))  ! and the normal one
        h_loc = min(h_loc, h)
      end do
    end do
    call MPI_Allreduce(h_loc, h_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, msh%comm, ierr)
    if (msh%rank /= 0 .or. h_glob >= 0.5_rp * huge(1.0_rp)) return
    write(*,'(A,ES9.2,A,F5.1,A)') ' SEM: smallest eddy half-extent ', state%sig_min, ' m = ', &
      state%sig_min / max(h_glob, SMALL), ' inlet cells'
    if (state%sig_min < 5.0_rp * h_glob) &
      write(*,'(A)') ' WARNING: SEM eddies are under-resolved (< 5 cells); the inlet pressure fluctuation grows ' // &
        'roughly as (dx/sigma)^2 - raise sem_sigma or refine the inlet (>= 8 cells keeps it near the nonlinear floor)'
  end subroutine sem_check_resolution

  !> Derive the per-eddy principal frame, half-extents and intensities from the target stress at each eddy
  !> centre. The independent random signs make the cross-stresses vanish in the principal frame and give
  !> lam_1 = K sigma_1^2 (A_2 + A_3) and cyclic, with A_i = alpha_i^2 and
  !> K = N sigma_1 sigma_2 sigma_3 DF_I / V_box, so A_1 = (b_2 + b_3 - b_1)/2 with b_i = lam_i/(K sigma_i^2).
  !> Cheap (O(N_eddies)); call after the eddies move or are restored.
  subroutine sem_prepare(state, abl_cfg)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: abl_cfg

    integer(ip) :: i, j
    real(rp)    :: Rm(3,3), Wr(3,3), lam(3), V(3,3), Vg(3,3), s(3), K, b(3), A(3)

    if (.not. state%dfsem) return
    if (.not. allocated(state%ampl)) return
    do i = 1, state%n_eddies
      call sem_stress_at(state, abl_cfg, state%pos(3, i), Rm, Wr)
      call sem_eigen3(Rm, lam, V)
      Vg = matmul(Wr, V)  ! principal axes in global coordinates
      do j = 1, 3
        state%frame(3*j-2:3*j, i) = Vg(:, j)
        s(j) = state%sigma * max(DF_FMIN, sqrt(max(lam(j), 0.0_rp) / state%lam_ref))
      end do
      state%sig(:, i) = s
      K = dble(state%n_eddies) * s(1) * s(2) * s(3) * DF_I / state%box_volume
      b = max(lam, 0.0_rp) / (K * s**2)
      A(1) = 0.5_rp * (b(2) + b(3) - b(1))
      A(2) = 0.5_rp * (b(1) + b(3) - b(2))
      A(3) = 0.5_rp * (b(1) + b(2) - b(3))
      state%ampl(:, i) = sqrt(max(A, 0.0_rp)) * state%amp(:, i)
    end do
  end subroutine sem_prepare

  !> Velocity of one divergence-free eddy at offset d from its centre (global axes), Poletto et al. Eq. 8 in
  !> the principal frame, then rotated back: u_b = sigma_b (1 - d^2) (r x alpha)_b, r_j = d.e_j / sigma_j.
  !> Zero outside the ellipsoid d < 1. Exactly solenoidal: du_b/dx_b = -2 r_b (r x alpha)_b (no sum),
  !> and the three terms cancel.
  pure function sem_dfsem_eddy_velocity(d, frame, sig, ampl) result(u)
    real(rp), intent(in) :: d(3), frame(9), sig(3), ampl(3)
    real(rp) :: u(3)

    real(rp) :: r(3), d2, g, up(3)
    integer(ip) :: j

    u = 0.0_rp
    d2 = 0.0_rp
    do j = 1, 3
      r(j) = dot_product(d, frame(3*j-2:3*j)) / sig(j)
      d2 = d2 + r(j)**2
    end do
    if (d2 >= 1.0_rp) return
    g = 1.0_rp - d2
    up(1) = sig(1) * g * (r(2) * ampl(3) - r(3) * ampl(2))
    up(2) = sig(2) * g * (r(3) * ampl(1) - r(1) * ampl(3))
    up(3) = sig(3) * g * (r(1) * ampl(2) - r(2) * ampl(1))
    u = up(1) * frame(1:3) + up(2) * frame(4:6) + up(3) * frame(7:9)
  end function sem_dfsem_eddy_velocity

  !> Mean velocity and synthetic fluctuation at an inlet face centre.
  !> tau (default 0) is a time look-ahead from the eddy state: the eddies are convected rigidly, so the field
  !> at time t + tau is the field at t sampled upstream by u_convect * tau along the inward normal.
  subroutine sem_velocity_at_face(state, abl_cfg, face_centre, u_mean, u_fluc_out, tau)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: face_centre(3)
    real(rp),        intent(out) :: u_mean(3)
    real(rp),        intent(out) :: u_fluc_out(3)
    real(rp), optional, intent(in) :: tau

    integer(ip) :: i
    real(rp)    :: d_vec(3), rs2, xs(3)

    xs = face_centre
    if (present(tau)) xs = face_centre - state%u_convect * tau * state%n_in

    if (.not. state%dfsem) then
      call sem_velocity_jarrin(state, abl_cfg, xs, u_mean, u_fluc_out, face_centre)
      return
    end if

    call sem_mean_at(state, abl_cfg, face_centre(3), u_mean)
    u_fluc_out = 0.0_rp
    rs2 = state%rsup**2
    do i = 1, state%n_eddies
      d_vec = xs - state%pos(:, i)
      if (state%per1) d_vec = d_vec - state%period1 * anint(dot_product(d_vec, state%t1) / state%period1) * state%t1
      if (state%per2) d_vec = d_vec - state%period2 * anint(dot_product(d_vec, state%t2) / state%period2) * state%t2
      if (sum(d_vec**2) >= rs2) cycle
      u_fluc_out = u_fluc_out + sem_dfsem_eddy_velocity(d_vec, state%frame(:, i), state%sig(:, i), state%ampl(:, i))
    end do
  end subroutine sem_velocity_at_face

  subroutine sem_velocity_jarrin(state, abl_cfg, face_centre, u_mean, u_fluc_out, face_orig)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: face_centre(3)  ! sample point (look-ahead shifted)
    real(rp),        intent(in)  :: face_orig(3)    ! the face centre itself (heights)
    real(rp),        intent(out) :: u_mean(3)
    real(rp),        intent(out) :: u_fluc_out(3)

    real(rp)    :: wind_rad, u_fluc_iso(3), L_chol(3,3), u_fluc_aniso(3), R6(6)
    integer(ip) :: i
    real(rp)    :: r_dist, tent_func, eddy_contrib(3), d_vec(3), sig, norm

    ! Mean and target stresses at the face height, shared with the divergence-free path.
    call sem_target_at(state, abl_cfg, face_orig(3), u_mean, R6)
    wind_rad = abl_cfg%wind_dir_deg * PI / 180.0_rp

    ! Isotropic unit-variance fluctuation from the eddy superposition (Jarrin et al. 2006).
    u_fluc_iso = 0.0_rp
    sig  = max(state%sigma, SMALL)
    norm = TENT_NORM * sqrt(state%box_volume / sig**3) / sqrt(dble(state%n_eddies))

    do i = 1, state%n_eddies
      ! Compact support: reject on the bounding box before taking the square root.
      d_vec = face_centre - state%pos(:, i)
      if (state%per1) d_vec = d_vec - state%period1 * anint(dot_product(d_vec, state%t1) / state%period1) * state%t1
      if (state%per2) d_vec = d_vec - state%period2 * anint(dot_product(d_vec, state%t2) / state%period2) * state%t2
      if (abs(d_vec(1)) >= sig .or. abs(d_vec(2)) >= sig .or. abs(d_vec(3)) >= sig) cycle

      r_dist = sqrt(sum(d_vec**2)) / sig
      tent_func = max(0.0_rp, 1.0_rp - r_dist)
      if (tent_func <= 0.0_rp) cycle

      eddy_contrib = norm * tent_func * state%amp(:, i)
      u_fluc_iso = u_fluc_iso + eddy_contrib
    end do

    ! Full Cholesky of the target tensor; 10-column tables are already in global axes.
    L_chol = sem_cholesky(R6)
    u_fluc_aniso(1) = L_chol(1,1) * u_fluc_iso(1)
    u_fluc_aniso(2) = L_chol(2,1) * u_fluc_iso(1) + L_chol(2,2) * u_fluc_iso(2)
    u_fluc_aniso(3) = L_chol(3,1) * u_fluc_iso(1) + L_chol(3,2) * u_fluc_iso(2) &
                    + L_chol(3,3) * u_fluc_iso(3)
    if (state%has_profile .and. state%prof_global_axes) then
      u_fluc_out = u_fluc_aniso
    else
      ! Rotate into the mean-wind frame: e_s along the mean, e_n horizontal and normal to it.
      u_fluc_out(1) = -u_fluc_aniso(1) * sin(wind_rad) + u_fluc_aniso(2) * cos(wind_rad)
      u_fluc_out(2) = -u_fluc_aniso(1) * cos(wind_rad) - u_fluc_aniso(2) * sin(wind_rad)
      u_fluc_out(3) =  u_fluc_aniso(3)
    end if
  end subroutine sem_velocity_jarrin

end module mod_sem_inflow
