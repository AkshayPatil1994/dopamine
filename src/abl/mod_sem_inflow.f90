! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Synthetic Eddy Method (SEM) for correlated turbulent inflow generation
module mod_sem_inflow
  use mod_kinds
  use mod_parameters, only: BC_ABL_INLET
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

  integer(int64), parameter :: RNG_SEED = 88172645463325252_int64
  integer(int64), save      :: rng_state = RNG_SEED

  public :: sem_init, sem_advance, sem_velocity_at_face
  public :: sem_read_profile, sem_profile_at, sem_cholesky

contains

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
  end subroutine sem_inlet_frame

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

    integer(ip) :: i, n_e
    real(rp)    :: span1, span2
    logical     :: found

    if (len_trim(abl_cfg%sem_profile_file) > 0) &
      call sem_read_profile(abl_cfg%sem_profile_file, state, verbose=(msh%rank == 0))

    call sem_inlet_frame(msh, state, found)
    if (.not. found) then
      state%initialized = .false.
      return
    end if

    state%n_eddies = max(1_ip, abl_cfg%n_eddies)
    n_e = state%n_eddies
    if (allocated(state%pos)) deallocate(state%pos)
    if (allocated(state%amp)) deallocate(state%amp)
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

    ! Pad the box by sigma so eddies straddling the edges are represented.
    state%t1_min = state%t1_min - state%sigma; state%t1_max = state%t1_max + state%sigma
    state%t2_min = state%t2_min - state%sigma; state%t2_max = state%t2_max + state%sigma
    ! Box spans -sigma to +sigma about the plane so a face there sees full kernel support.
    state%box_volume = 2.0_rp * state%sigma &
                     * (state%t1_max - state%t1_min) * (state%t2_max - state%t2_min)

    state%u_convect = max(abl_cfg%u_ref, SMALL)

    rng_state = RNG_SEED
    do i = 1, n_e
      call sem_place_eddy(state, i, state%sigma * (2.0_rp * sem_rand() - 1.0_rp), .true.)
    end do
    call sem_balance_amplitudes(state)

    state%initialized = .true.
  end subroutine sem_init

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
  !> Amplitudes are redrawn only at seeding, so the zero-mean balance survives recycling.
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
      if (s > state%sigma) call sem_place_eddy(state, i, -state%sigma, .false.)
    end do
  end subroutine sem_advance

  !> Evaluate SEM synthetic velocity fluctuation at a point
  !> Returns fluctuation vector u'(x) computed as a sum of eddy contributions:

  !> Read prescribed inflow profiles: 6 columns z,U,uu,vv,ww,uw (wind-aligned)
  !> or 10 columns z,U,V,W,uu,vv,ww,uv,uw,vw (global axes).
  subroutine sem_read_profile(filename, state, verbose)
    character(len=*), intent(in)    :: filename
    type(SEMState),   intent(inout) :: state
    logical, intent(in), optional   :: verbose  ! rank-0 only in parallel runs

    integer(ip) :: iunit, ios, k, n, ncol
    character(len=512) :: line
    real(rp) :: v10(10), v6(6)

    open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'ERROR: cannot open SEM profile file: '//trim(filename)
      error stop 'sem_read_profile: open failed'
    end if

    ! First pass: count valid rows and settle the column count from the first one.
    n = 0
    ncol = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      if (ncol == 0) then
        read(line, *, iostat=ios) v10
        if (ios == 0) then
          ncol = 10
        else
          read(line, *, iostat=ios) v6
          if (ios /= 0) cycle
          ncol = 6
        end if
      end if
      if (ncol == 10) then
        read(line, *, iostat=ios) v10
      else
        read(line, *, iostat=ios) v6
      end if
      if (ios == 0) n = n + 1
    end do

    if (n < 2) then
      write(*,'(A)') 'ERROR: SEM profile needs at least 2 levels: '//trim(filename)
      error stop 'sem_read_profile: too few levels'
    end if

    state%n_levels = n
    state%prof_global_axes = (ncol == 10)
    allocate(state%prof_z(n), source=0.0_rp)
    allocate(state%prof_U(3, n), source=0.0_rp)
    allocate(state%prof_R(6, n), source=0.0_rp)

    rewind(iunit)
    k = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      if (ncol == 10) then
        read(line, *, iostat=ios) v10
        if (ios /= 0) cycle
        k = k + 1
        if (k > n) exit
        state%prof_z(k)    = v10(1)
        state%prof_U(:, k) = v10(2:4)
        state%prof_R(:, k) = v10(5:10)
      else
        read(line, *, iostat=ios) v6
        if (ios /= 0) cycle
        k = k + 1
        if (k > n) exit
        state%prof_z(k)      = v6(1)
        state%prof_U(1, k)   = v6(2)  ! streamwise magnitude, rotated later
        state%prof_R(1:3, k) = v6(3:5)
        state%prof_R(5, k)   = v6(6)  ! uw
      end if
    end do
    close(iunit)

    call sem_sort_profile(state)
    state%has_profile = .true.
    if (present(verbose)) then
      if (verbose) write(*,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A)') &
        ' SEM inflow profiles: ', n, ' levels, ', ncol, ' columns, z = ', &
        state%prof_z(1), ' to ', state%prof_z(n)
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

  subroutine sem_velocity_at_face(state, abl_cfg, face_centre, u_mean, u_fluc_out)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: face_centre(3)
    real(rp),        intent(out) :: u_mean(3)
    real(rp),        intent(out) :: u_fluc_out(3)

    real(rp)    :: z, wind_rad, u_mag
    real(rp)    :: u_star, u_fluc_iso(3), R_iso(3,3), L_chol(3,3)
    real(rp)    :: u_fluc_aniso(3), R_prof(6)
    integer(ip) :: i
    real(rp)    :: r_dist, tent_func, eddy_contrib(3), d_vec(3), sig, norm

    ! Height above ground
    z = face_centre(3)

    wind_rad = abl_cfg%wind_dir_deg * PI / 180.0_rp
    if (state%has_profile) then
      call sem_profile_at(state, z, u_mean, R_prof)
      ! 6-column tables give a streamwise magnitude; 10-column ones are already global.
      if (.not. state%prof_global_axes) u_mean = abl_wind_vector(u_mean(1), abl_cfg%wind_dir_deg)
    else
      ! Mean velocity: neutral MOST log-law, same helpers as the ABL inlet BC.
      u_star = abl_ustar_from_ref(abl_cfg%u_ref, abl_cfg%z_ref, abl_cfg%z0_m)
      u_mag  = abl_u_profile(z, u_star, abl_cfg%z0_m, 0.0_rp)
      u_mean = abl_wind_vector(u_mag, abl_cfg%wind_dir_deg)
    end if

    ! Isotropic unit-variance fluctuation from the eddy superposition (Jarrin et al. 2006).
    u_fluc_iso = 0.0_rp
    sig  = max(state%sigma, SMALL)
    norm = TENT_NORM * sqrt(state%box_volume / sig**3) / sqrt(dble(state%n_eddies))

    do i = 1, state%n_eddies
      ! Compact support: reject on the bounding box before taking the square root.
      d_vec = face_centre - state%pos(:, i)
      if (abs(d_vec(1)) >= sig .or. abs(d_vec(2)) >= sig .or. abs(d_vec(3)) >= sig) cycle

      r_dist = sqrt(sum(d_vec**2)) / sig
      tent_func = max(0.0_rp, 1.0_rp - r_dist)
      if (tent_func <= 0.0_rp) cycle

      eddy_contrib = norm * tent_func * state%amp(:, i)
      u_fluc_iso = u_fluc_iso + eddy_contrib
    end do

    ! Prescribed tensor path: full Cholesky, and 10-column tables need no rotation.
    if (state%has_profile) then
      L_chol = sem_cholesky(R_prof)
      u_fluc_aniso(1) = L_chol(1,1) * u_fluc_iso(1)
      u_fluc_aniso(2) = L_chol(2,1) * u_fluc_iso(1) + L_chol(2,2) * u_fluc_iso(2)
      u_fluc_aniso(3) = L_chol(3,1) * u_fluc_iso(1) + L_chol(3,2) * u_fluc_iso(2) &
                      + L_chol(3,3) * u_fluc_iso(3)
      if (state%prof_global_axes) then
        u_fluc_out = u_fluc_aniso
      else
        u_fluc_out(1) = -u_fluc_aniso(1) * sin(wind_rad) + u_fluc_aniso(2) * cos(wind_rad)
        u_fluc_out(2) = -u_fluc_aniso(1) * cos(wind_rad) - u_fluc_aniso(2) * sin(wind_rad)
        u_fluc_out(3) =  u_fluc_aniso(3)
      end if
      return
    end if

    ! Target Reynolds stresses (Panofsky & Dutton profiles for surface layer)
    ! sigma_u/u_star ~= 2.4, sigma_v/u_star ~= 1.9, sigma_w/u_star ~= 1.25
    R_iso(1, 1) = (2.4_rp * u_star) ** 2
    R_iso(2, 2) = (1.9_rp * u_star) ** 2
    R_iso(3, 3) = (1.25_rp * u_star) ** 2
    R_iso(1, 2) = -0.1_rp * u_star**2  ! small negative correlation
    R_iso(2, 1) = R_iso(1, 2)
    R_iso(1, 3) = 0.0_rp
    R_iso(3, 1) = 0.0_rp
    R_iso(2, 3) = 0.0_rp
    R_iso(3, 2) = 0.0_rp

    ! Cholesky decomposition of R: R = L * L^T
    ! (simplified: here we just apply diagonal scaling)
    L_chol = 0.0_rp
    L_chol(1, 1) = sqrt(max(R_iso(1, 1), SMALL))
    L_chol(2, 1) = R_iso(2, 1) / max(L_chol(1, 1), SMALL)
    L_chol(2, 2) = sqrt(max(R_iso(2, 2) - L_chol(2, 1)**2, SMALL))
    L_chol(3, 3) = sqrt(max(R_iso(3, 3), SMALL))

    ! Apply to isotropic fluctuations
    u_fluc_aniso(1) = L_chol(1, 1) * u_fluc_iso(1)
    u_fluc_aniso(2) = L_chol(2, 1) * u_fluc_iso(1) + L_chol(2, 2) * u_fluc_iso(2)
    u_fluc_aniso(3) = L_chol(3, 3) * u_fluc_iso(3)

    ! Rotate into the mean-wind frame: e_s along the mean, e_n horizontal and normal to it.
    u_fluc_out(1) = -u_fluc_aniso(1) * sin(wind_rad) + u_fluc_aniso(2) * cos(wind_rad)
    u_fluc_out(2) = -u_fluc_aniso(1) * cos(wind_rad) - u_fluc_aniso(2) * sin(wind_rad)
    u_fluc_out(3) =  u_fluc_aniso(3)

  end subroutine sem_velocity_at_face

end module mod_sem_inflow
