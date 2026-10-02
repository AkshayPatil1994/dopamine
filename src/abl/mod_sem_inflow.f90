! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Divergence-free Synthetic Eddy Method (SEM) for correlated turbulent inflow generation.
!>
!> Each eddy is the curl of a compact vector potential A = phi(s) * a, with
!> s^2 = sum_i (d_i / sig_i)^2 the ellipsoidal distance from the eddy centre:
!>
!>     u(x) = grad(phi) x a,     grad(phi)_i = -6 (1 - s^2)^2 d_i / sig_i^2.
!>
!> This is solenoidal for any sig_i and any a, which is what lets every eddy carry its own
!> height-dependent scales and stress-fitted strength without breaking div(u) = 0. The
!> kernel is C1, so vorticity is continuous. The strength a = C xi (xi random +-1) is chosen
!> from the target Reynolds stresses R and the eddy scales by `sem_eddy_coeffs`; the stress
!> map is closed-form because the kernel is even:
!>
!>     uu = Qy Szz + Qz Syy      uv = -Qz Sxy        Q_i = N_sigma * QJ / sig_i^2
!>     vv = Qz Sxx + Qx Szz      uw = -Qy Sxz        S = C C^T
!>     ww = Qy Sxx + Qx Syy      vw = -Qx Syz
!>
!> Not every (R, sig) pair has S positive semi-definite (e.g. isotropic eddies cannot give
!> uu > vv + ww). Those cases are clipped to the nearest realisable S and reported, never silent.
!>
!> References: Jarrin et al. (2006) for the SEM; Poletto et al. (2013) for the
!> divergence-free construction.
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

  ! Int s^4 (1-s^2)^4 ds over [0,1] = 384/45045, so Int |grad phi|^2 dV / sx sy sz = 144 pi * that.
  ! QJ = that / 3: the per-axis second moment of grad(phi) per unit eddy density (sig_i = 1).
  real(rp), parameter :: QJ = 48.0_rp * PI * 384.0_rp / 45045.0_rp
  ! Kernel semi-axis per streamwise integral length scale of u. Calibrated on the isotropic
  ! kernel (R_uu(r e_x) integrated numerically: L_ux = 0.596 sigma).
  real(rp), parameter :: SIGMA_PER_L = 1.678_rp

  integer(int64), parameter :: RNG_SEED = 88172645463325252_int64
  integer(int64), save      :: rng_state = RNG_SEED

  public :: sem_init, sem_setup, sem_advance, sem_velocity_at_face, sem_fluct_at_point
  public :: sem_fluct_at_point_ref, sem_bins_fill, sem_setup_targets, sem_init_targets
  public :: sem_read_profile, sem_profile_at, sem_cholesky, sem_target_at, sem_eddy_coeffs
  public :: SIGMA_PER_L

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

  !> Initialize SEM from the mesh: read any profile table, find the inlet frame, then set up eddies.
  subroutine sem_init(state, abl_cfg, msh)
    type(SEMState),   intent(inout) :: state
    type(ABLConfig),  intent(in)    :: abl_cfg
    type(Mesh),       intent(in)    :: msh

    logical :: found

    if (len_trim(abl_cfg%sem_profile_file) > 0) &
      call sem_read_profile(abl_cfg%sem_profile_file, state, verbose=(msh%rank == 0), &
                            fmt=abl_cfg%sem_profile_format, rho_uw=abl_cfg%sem_rho_uw)

    call sem_inlet_frame(msh, state, found)
    if (.not. found) then
      state%initialized = .false.
      return
    end if

    call sem_setup(state, abl_cfg, verbose=(msh%rank == 0))
  end subroutine sem_init

  !> Targets and inlet frame only, no eddies: for a source that supplies its own fluctuation
  !> (the inflow-plane database) but still takes the mean profile from the SEM laws and tables.
  subroutine sem_init_targets(state, abl_cfg, msh, found, verbose)
    type(SEMState),   intent(inout) :: state
    type(ABLConfig),  intent(in)    :: abl_cfg
    type(Mesh),       intent(in)    :: msh
    logical,          intent(out)   :: found
    logical,          intent(in), optional :: verbose

    if (len_trim(abl_cfg%sem_profile_file) > 0 .and. .not. state%has_profile) &
      call sem_read_profile(abl_cfg%sem_profile_file, state, verbose=(msh%rank == 0), &
                            fmt=abl_cfg%sem_profile_format, rho_uw=abl_cfg%sem_rho_uw)
    call sem_inlet_frame(msh, state, found)
    if (.not. found) return
    call sem_setup_targets(state, abl_cfg)
  end subroutine sem_init_targets

  !> Everything after the inlet frame is known: scales, density, eddy seeding, realisability scan.
  !> Needs state%{x_inlet,n_in,t1,t2,t1_min..t2_max} and, for tabulated inflow, the profile.
  subroutine sem_setup(state, abl_cfg, verbose)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: abl_cfg
    logical,         intent(in), optional :: verbose

    integer(ip) :: i, n_e
    real(rp)    :: span1, span2, rad
    logical     :: talk

    talk = .false.
    if (present(verbose)) talk = verbose

    rng_state = RNG_SEED
    if (abl_cfg%sem_seed /= 0) rng_state = ieor(RNG_SEED, int(abl_cfg%sem_seed, int64) * 2654435761_int64)
    if (rng_state == 0_int64) rng_state = RNG_SEED
    do i = 1, 8
      rad = sem_rand()  ! warm up so nearby seeds decorrelate
    end do

    n_e = max(1_ip, abl_cfg%n_eddies)
    state%n_eddies = n_e
    if (allocated(state%pos))   deallocate(state%pos)
    if (allocated(state%amp))   deallocate(state%amp)
    if (allocated(state%sig))   deallocate(state%sig)
    if (allocated(state%uconv)) deallocate(state%uconv)
    allocate(state%pos(3, n_e), source=0.0_rp)
    allocate(state%amp(3, n_e), source=0.0_rp)
    allocate(state%sig(3, n_e), source=1.0_rp)
    allocate(state%uconv(n_e),  source=0.0_rp)

    ! Profile frame: wind-aligned unless a table is given in global axes.
    if (state%has_profile .and. state%prof_global_axes) then
      state%frame = reshape([1.0_rp,0.0_rp,0.0_rp, 0.0_rp,1.0_rp,0.0_rp, 0.0_rp,0.0_rp,1.0_rp], [3,3])
    else
      rad = abl_cfg%wind_dir_deg * PI / 180.0_rp
      state%frame(:, 1) = [-sin(rad), -cos(rad), 0.0_rp]
      state%frame(:, 2) = [ cos(rad), -sin(rad), 0.0_rp]
      state%frame(:, 3) = [0.0_rp, 0.0_rp, 1.0_rp]
    end if

    span1 = state%t1_max - state%t1_min
    span2 = state%t2_max - state%t2_min
    state%L_ref = abl_cfg%sem_sigma
    if (state%L_ref <= SMALL) state%L_ref = 0.05_rp * min(span1, span2)
    if (state%L_ref <= SMALL) state%L_ref = 0.05_rp * abl_cfg%z_ref
    state%z_floor    = max(abl_cfg%z0_m, 1.0e-6_rp)
    state%u_conv_min = 0.05_rp * max(abl_cfg%u_ref, SMALL)
    state%override_mean = state%has_profile .or. trim(abl_cfg%sem_mean_law) == 'power'

    if (state%has_profile) call sem_fill_table_lengths(state, abl_cfg)
    call sem_build_density(state, abl_cfg)

    do i = 1, n_e
      call sem_place_eddy(state, abl_cfg, i, state%s_half * (2.0_rp * sem_rand() - 1.0_rp), .false.)
    end do
    call sem_bins_fill(state)

    state%initialized = .true.
    call sem_report(state, abl_cfg, talk)
  end subroutine sem_setup

  !> Targets only (no eddies): lets a prescribed-fluctuation source such as the inflow-plane
  !> database reuse the mean-profile laws and tables for its mean velocity.
  subroutine sem_setup_targets(state, abl_cfg)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: abl_cfg

    state%L_ref      = max(abl_cfg%sem_sigma, SMALL)
    state%z_floor    = max(abl_cfg%z0_m, 1.0e-6_rp)
    state%u_conv_min = 0.05_rp * max(abl_cfg%u_ref, SMALL)
    state%override_mean = state%has_profile .or. trim(abl_cfg%sem_mean_law) == 'power'
    if (state%has_profile) call sem_fill_table_lengths(state, abl_cfg)
  end subroutine sem_setup_targets

  !> Table rows without length scales take them from the power-law / ratio model.
  subroutine sem_fill_table_lengths(state, cfg)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: cfg
    integer(ip) :: k
    real(rp)    :: Llaw(3)

    do k = 1, state%n_levels
      Llaw = sem_length_law(state, cfg, state%prof_z(k))
      if (state%prof_L(1, k) <= 0.0_rp) state%prof_L(1, k) = Llaw(1)
      if (state%prof_L(2, k) <= 0.0_rp) state%prof_L(2, k) = cfg%sem_Ly_ratio * state%prof_L(1, k)
      if (state%prof_L(3, k) <= 0.0_rp) state%prof_L(3, k) = cfg%sem_Lz_ratio * state%prof_L(1, k)
    end do
  end subroutine sem_fill_table_lengths

  !> Integral length scales Lx,Ly,Lz [m] from the namelist model: Lx = L_ref (z/z_ref)^m.
  function sem_length_law(state, cfg, z) result(L)
    type(SEMState),  intent(in) :: state
    type(ABLConfig), intent(in) :: cfg
    real(rp),        intent(in) :: z
    real(rp) :: L(3), zz

    L(1) = state%L_ref
    if (abs(cfg%sem_L_exponent) > SMALL) then
      zz   = max(z, state%z_floor)
      L(1) = state%L_ref * (zz / max(cfg%z_ref, SMALL)) ** cfg%sem_L_exponent
    end if
    L(2) = cfg%sem_Ly_ratio * L(1)
    L(3) = cfg%sem_Lz_ratio * L(1)
  end function sem_length_law

  !> Target mean velocity (global axes), stresses R = [uu,vv,ww,uv,uw,vw] and integral scales
  !> (profile frame) at height z.
  subroutine sem_target_at(state, cfg, z, u_mean, R, L)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: cfg
    real(rp),        intent(in)  :: z
    real(rp),        intent(out) :: u_mean(3), R(6), L(3)

    real(rp) :: u_star, umag, zz, zi, Iu, su, sv, sw

    if (state%has_profile) then
      call sem_profile_at(state, z, u_mean, R, L)
      ! 6/5-column tables give a streamwise speed; 10-column ones are already global.
      if (.not. state%prof_global_axes) u_mean = abl_wind_vector(u_mean(1), cfg%wind_dir_deg)
      L = max(L, cfg%sem_L_min)
      return
    end if

    zz     = max(z, state%z_floor)
    u_star = abl_ustar_from_ref(cfg%u_ref, cfg%z_ref, cfg%z0_m)
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
        umag = cfg%u_ref * (max(zi, state%z_floor) / max(cfg%z_ref, SMALL)) ** cfg%sem_alpha
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

    L = max(sem_length_law(state, cfg, z), cfg%sem_L_min)
  end subroutine sem_target_at

  !> Stress-fitted potential strength for one eddy. Given the target stresses R (profile frame)
  !> and kernel semi-axes sig, returns C with S = C C^T so that a = C xi (xi unit-variance,
  !> independent) reproduces R. err is the worst relative stress error after the S
  !> realisability clip: 0 to rounding when R is attainable with these scales.
  pure subroutine sem_eddy_coeffs(R, sig, n_sigma, C, err)
    real(rp), intent(in)  :: R(6), sig(3), n_sigma
    real(rp), intent(out) :: C(3,3), err

    real(rp) :: a, b, cz, S(3,3), Ra(6), dmax, e
    integer(ip) :: k

    ! Q_i = n_sigma * QJ / sig_i^2
    a  = n_sigma * QJ / max(sig(1), SMALL)**2
    b  = n_sigma * QJ / max(sig(2), SMALL)**2
    cz = n_sigma * QJ / max(sig(3), SMALL)**2

    S(1,1) = (-a * R(1) + b * R(2) + cz * R(3)) / (2.0_rp * b * cz)
    S(2,2) = ( a * R(1) - b * R(2) + cz * R(3)) / (2.0_rp * a * cz)
    S(3,3) = ( a * R(1) + b * R(2) - cz * R(3)) / (2.0_rp * a * b)
    S(2,1) = -R(4) / cz
    S(3,1) = -R(5) / b
    S(3,2) = -R(6) / a

    ! Clip to a PSD matrix with a Cholesky that floors its pivots.
    dmax = max(abs(S(1,1)), abs(S(2,2)), abs(S(3,3)), SMALL)
    C = 0.0_rp
    C(1,1) = sqrt(max(S(1,1), 1.0e-12_rp * dmax))
    C(2,1) = S(2,1) / C(1,1)
    C(3,1) = S(3,1) / C(1,1)
    C(2,2) = sqrt(max(S(2,2) - C(2,1)**2, 1.0e-12_rp * dmax))
    C(3,2) = (S(3,2) - C(3,1) * C(2,1)) / C(2,2)
    C(3,3) = sqrt(max(S(3,3) - C(3,1)**2 - C(3,2)**2, 1.0e-12_rp * dmax))

    ! Stresses actually delivered by the clipped S.
    S = matmul(C, transpose(C))
    Ra(1) = b * S(3,3) + cz * S(2,2)
    Ra(2) = cz * S(1,1) + a * S(3,3)
    Ra(3) = b * S(1,1) + a * S(2,2)
    Ra(4) = -cz * S(2,1)
    Ra(5) = -b * S(3,1)
    Ra(6) = -a * S(3,2)

    err = 0.0_rp
    do k = 1, 3
      e = abs(Ra(k) - R(k)) / max(abs(R(k)), SMALL)
      err = max(err, e)
    end do
    do k = 4, 6
      e = abs(Ra(k) - R(k)) / max(sqrt(abs(R(1) * R(3))), SMALL)
      err = max(err, e)
    end do
  end subroutine sem_eddy_coeffs

  !> Extent of the kernel ellipsoid along the unit vector v (global axes), given semi-axes in the frame.
  pure real(rp) function kernel_extent(state, sg, v) result(ext)
    type(SEMState), intent(in) :: state
    real(rp),       intent(in) :: sg(3), v(3)
    real(rp) :: vl(3)
    vl  = matmul(transpose(state%frame), v)
    ext = sqrt(sum((sg * vl)**2))
  end function kernel_extent

  !> Pad the eddy box so eddies straddling its faces are represented, and tabulate the height
  !> density along t2: eddies per volume ~ 1/(sx sy sz) so every kernel holds ~ n_sigma eddies.
  subroutine sem_build_density(state, cfg)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: cfg

    integer(ip), parameter :: NCDF = 512
    integer(ip) :: it, k
    real(rp) :: b0, b1, db, b, x(3), um(3), R(6), L(3), sg(3)
    real(rp) :: pad1, pad2, s_h, tlo1, thi1, tlo2, thi2, w(NCDF), wr(NCDF), vol_w, uc

    tlo1 = state%t1_min; thi1 = state%t1_max
    tlo2 = state%t2_min; thi2 = state%t2_max
    pad1 = 0.0_rp; pad2 = 0.0_rp; s_h = 0.0_rp

    ! Kernel sizes depend on height, which depends on the padded box: iterate to a fixed point.
    do it = 1, 3
      b0 = tlo2 - pad2; b1 = thi2 + pad2
      db = (b1 - b0) / real(NCDF - 1, rp)
      pad1 = 0.0_rp; pad2 = 0.0_rp; s_h = 0.0_rp
      do k = 1, NCDF
        b = b0 + real(k - 1, rp) * db
        x = state%x_inlet + b * state%t2
        call sem_target_at(state, cfg, x(3), um, R, L)
        sg = SIGMA_PER_L * L
        s_h  = max(s_h,  kernel_extent(state, sg, state%n_in))
        pad1 = max(pad1, kernel_extent(state, sg, state%t1))
        pad2 = max(pad2, kernel_extent(state, sg, state%t2))
      end do
    end do

    state%ext = 1.02_rp * [pad1, pad2, max(s_h, SMALL)]
    if (state%periodic_t1) then
      ! The box is exactly one period; images are handled by minimum-image distances, which are
      ! only correct while a kernel cannot reach across more than half the period.
      if (state%t1_period < 2.0_rp * pad1) then
        write(*,'(A,ES10.3,A,ES10.3,A)') 'ERROR: spanwise period ', state%t1_period, &
          ' m is shorter than twice the largest eddy half-extent (', pad1, ' m).'
        error stop 'sem_setup: period too short for the requested length scales (reduce Lx or enlarge the period)'
      end if
      state%t1_min = tlo1; state%t1_max = tlo1 + state%t1_period
    else
      state%t1_min = tlo1 - pad1; state%t1_max = thi1 + pad1
    end if
    state%t2_min = tlo2 - pad2; state%t2_max = thi2 + pad2
    state%s_half = max(s_h, SMALL)

    b0 = state%t2_min
    db = (state%t2_max - state%t2_min) / real(NCDF - 1, rp)
    do k = 1, NCDF
      b = b0 + real(k - 1, rp) * db
      x = state%x_inlet + b * state%t2
      call sem_target_at(state, cfg, x(3), um, R, L)
      sg = SIGMA_PER_L * L
      w(k)  = 1.0_rp / max(sg(1) * sg(2) * sg(3), SMALL)
      uc    = max(dot_product(um, state%n_in), state%u_conv_min)
      wr(k) = w(k) * uc
    end do

    if (allocated(state%cdf_seed)) deallocate(state%cdf_seed)
    if (allocated(state%cdf_rec))  deallocate(state%cdf_rec)
    allocate(state%cdf_seed(NCDF), state%cdf_rec(NCDF))
    state%n_cdf  = NCDF
    state%cdf_b0 = b0
    state%cdf_db = db
    state%cdf_seed(1) = 0.0_rp
    state%cdf_rec(1)  = 0.0_rp
    do k = 2, NCDF
      state%cdf_seed(k) = state%cdf_seed(k-1) + 0.5_rp * (w(k-1)  + w(k))  * db
      state%cdf_rec(k)  = state%cdf_rec(k-1)  + 0.5_rp * (wr(k-1) + wr(k)) * db
    end do
    vol_w = state%cdf_seed(NCDF)
    state%cdf_seed = state%cdf_seed / max(state%cdf_seed(NCDF), SMALL)
    state%cdf_rec  = state%cdf_rec  / max(state%cdf_rec(NCDF),  SMALL)

    ! N = n_sigma * (2 s_half) * (t1 span) * Int dz / (sx sy sz)
    state%n_sigma = real(state%n_eddies, rp) &
                  / max(2.0_rp * state%s_half * (state%t1_max - state%t1_min) * vol_w, SMALL)
  end subroutine sem_build_density

  !> Draw a height coordinate b along t2 from the tabulated density by inverse CDF.
  real(rp) function sem_sample_b(state, recycle) result(b)
    type(SEMState), intent(in) :: state
    logical,        intent(in) :: recycle

    integer(ip) :: lo, hi, mid
    real(rp)    :: u, frac, c0, c1

    u  = sem_rand()
    lo = 1; hi = state%n_cdf
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (recycle) then
        c0 = state%cdf_rec(mid)
      else
        c0 = state%cdf_seed(mid)
      end if
      if (c0 <= u) then
        lo = mid
      else
        hi = mid
      end if
    end do
    if (recycle) then
      c0 = state%cdf_rec(lo);  c1 = state%cdf_rec(hi)
    else
      c0 = state%cdf_seed(lo); c1 = state%cdf_seed(hi)
    end if
    frac = (u - c0) / max(c1 - c0, SMALL)
    b = state%cdf_b0 + (real(lo - 1, rp) + frac) * state%cdf_db
  end function sem_sample_b

  !> Place eddy i at signed distance s along the inward normal and draw its height, lateral
  !> position, scales and strength. Recycled eddies draw heights flux-weighted by their
  !> convection speed, which keeps the height density stationary when slow eddies dwell longer.
  subroutine sem_place_eddy(state, cfg, i, s, recycle)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: cfg
    integer(ip),     intent(in)    :: i
    real(rp),        intent(in)    :: s
    logical,         intent(in)    :: recycle

    real(rp)    :: a, b, um(3), R(6), L(3), C(3,3), xi(3), e
    integer(ip) :: k

    b = sem_sample_b(state, recycle)
    a = state%t1_min + (state%t1_max - state%t1_min) * sem_rand()
    state%pos(:, i) = state%x_inlet + s * state%n_in + a * state%t1 + b * state%t2

    call sem_target_at(state, cfg, state%pos(3, i), um, R, L)
    state%sig(:, i) = SIGMA_PER_L * L
    call sem_eddy_coeffs(R, state%sig(:, i), state%n_sigma, C, e)
    do k = 1, 3
      xi(k) = sign(1.0_rp, sem_rand() - 0.5_rp)
    end do
    state%amp(:, i)  = matmul(C, xi)
    state%uconv(i)   = max(dot_product(um, state%n_in), state%u_conv_min)
  end subroutine sem_place_eddy

  !> Advance SEM: convect each eddy along the inward normal at its own mean speed (Taylor's
  !> hypothesis with the local U(z)) and recycle those leaving the box.
  subroutine sem_advance(state, abl_cfg, dt)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: abl_cfg
    real(rp),        intent(in)    :: dt

    integer(ip) :: i
    real(rp)    :: s

    if (.not. state%initialized) return

    do i = 1, state%n_eddies
      state%pos(:, i) = state%pos(:, i) + state%uconv(i) * dt * state%n_in
      s = dot_product(state%pos(:, i) - state%x_inlet, state%n_in)
      if (s > state%s_half) then
        ! Wrap rather than reset so a long step keeps the uniform spacing along the normal.
        s = modulo(s + state%s_half, 2.0_rp * state%s_half) - state%s_half
        call sem_place_eddy(state, abl_cfg, i, s, .true.)
      end if
    end do
    call sem_bins_fill(state)
  end subroutine sem_advance

  !> Report setup diagnostics: kernel counts and the realisability of the targets. Rank 0 only.
  subroutine sem_report(state, cfg, talk)
    type(SEMState),  intent(inout) :: state
    type(ABLConfig), intent(in)    :: cfg
    logical,         intent(in)    :: talk

    integer(ip), parameter :: NS = 33
    integer(ip) :: k
    real(rp)    :: b, x(3), um(3), R(6), L(3), C(3,3), e, worst, zworst, overlap

    worst = 0.0_rp; zworst = 0.0_rp
    do k = 1, NS
      b = state%t2_min + (state%t2_max - state%t2_min) * real(k - 1, rp) / real(NS - 1, rp)
      x = state%x_inlet + b * state%t2
      call sem_target_at(state, cfg, x(3), um, R, L)
      call sem_eddy_coeffs(R, SIGMA_PER_L * L, state%n_sigma, C, e)
      if (e > worst) then
        worst = e; zworst = x(3)
      end if
    end do
    state%realiz_err = worst

    if (.not. talk) return
    overlap = state%n_sigma * 4.0_rp * PI / 3.0_rp
    write(*,'(A,I0,A,ES9.2,A,ES9.2,A)') ' SEM inflow: ', state%n_eddies, ' divergence-free eddies, ', &
      state%n_sigma, ' per kernel volume (~', overlap, ' overlapping at a point)'
    write(*,'(A,ES9.2,A,ES9.2,A,ES9.2,A)') ' SEM inflow: Lx(z_ref) = ', state%L_ref, ' m, box depth ', &
      2.0_rp * state%s_half, ' m, exponent ', cfg%sem_L_exponent, ''
    if (state%n_sigma < 1.0_rp) &
      write(*,'(A,I0,A)') ' WARNING: SEM has few eddies per kernel; the field will be intermittent. n_eddies >= ', &
        ceiling(real(state%n_eddies, rp) * 2.0_rp / max(state%n_sigma, SMALL)), ' gives n_sigma = 2.'
    if (worst > 0.02_rp) then
      write(*,'(A,F6.1,A,ES9.2,A)') ' WARNING: SEM stresses are not realisable by divergence-free eddies with these ' // &
        'length scales (worst error ', 100.0_rp * worst, '% at z = ', zworst, ' m).'
      write(*,'(A)') '          Make Ly/Lx and Lz/Lx smaller (elongate the eddies) or lower uu relative to vv+ww.'
    end if
  end subroutine sem_report

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
  !> Any of these may be followed by Lx Ly Lz (integral length scales [m], 0 = use the model),
  !> or for 'ti' just Lx or Lx Ly.
  subroutine sem_read_profile(filename, state, verbose, fmt, rho_uw)
    character(len=*), intent(in)    :: filename
    type(SEMState),   intent(inout) :: state
    logical,          intent(in), optional :: verbose  ! rank-0 only in parallel runs
    character(len=*), intent(in), optional :: fmt
    real(rp),         intent(in), optional :: rho_uw   ! correlation used when a 'ti' table has none

    integer(ip) :: iunit, ios, k, n, ncol
    character(len=512) :: line
    real(rp) :: v(13), tmp, rho_def
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
        ok = .false.
        if (is_ti) then
          ok = (ncol >= 5 .and. ncol <= 9)
        else
          ok = (ncol == 6 .or. ncol == 9 .or. ncol == 10 .or. ncol == 13)
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
    state%prof_global_axes = (.not. is_ti) .and. (ncol == 10 .or. ncol == 13)
    if (allocated(state%prof_z)) deallocate(state%prof_z, state%prof_U, state%prof_R, state%prof_L)
    allocate(state%prof_z(n), source=0.0_rp)
    allocate(state%prof_U(3, n), source=0.0_rp)
    allocate(state%prof_R(6, n), source=0.0_rp)
    allocate(state%prof_L(3, n), source=0.0_rp)

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
      call sem_parse_row(v, ncol, is_ti, rho_def, state%prof_z(k), state%prof_U(:, k), &
                         state%prof_R(:, k), state%prof_L(:, k))
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

  !> One table row to height, mean, stresses and scales. For wind-aligned tables U(1) is a speed.
  pure subroutine sem_parse_row(v, ncol, is_ti, rho_def, z, U, R, L)
    real(rp),    intent(in)  :: v(13), rho_def
    integer(ip), intent(in)  :: ncol
    logical,     intent(in)  :: is_ti
    real(rp),    intent(out) :: z, U(3), R(6), L(3)

    real(rp) :: su, sv, sw, rho
    integer(ip) :: nb

    U = 0.0_rp; R = 0.0_rp; L = 0.0_rp
    z = v(1)
    if (is_ti) then
      U(1) = v(2)
      su = v(3) * v(2); sv = v(4) * v(2); sw = v(5) * v(2)
      rho = rho_def
      if (ncol >= 6) rho = v(6)
      R(1) = su**2; R(2) = sv**2; R(3) = sw**2
      R(5) = rho * su * sw
      if (ncol >= 7) L(1) = v(7)
      if (ncol >= 8) L(2) = v(8)
      if (ncol >= 9) L(3) = v(9)
    else if (ncol == 6 .or. ncol == 9) then
      U(1)   = v(2)  ! streamwise magnitude, rotated later
      R(1:3) = v(3:5)
      R(5)   = v(6)  ! uw
      nb = 6
      if (ncol == 9) L = v(nb+1:nb+3)
    else
      U    = v(2:4)
      R    = v(5:10)
      nb = 10
      if (ncol == 13) L = v(nb+1:nb+3)
    end if
  end subroutine sem_parse_row

  !> Insertion sort by height; tables are short and often already ascending.
  subroutine sem_sort_profile(state)
    type(SEMState), intent(inout) :: state

    integer(ip) :: i, j
    real(rp) :: zt, ut(3), rt(6), lt(3)

    do i = 2, state%n_levels
      zt = state%prof_z(i); ut = state%prof_U(:, i); rt = state%prof_R(:, i); lt = state%prof_L(:, i)
      j = i - 1
      do while (j >= 1)
        if (state%prof_z(j) <= zt) exit
        state%prof_z(j+1)    = state%prof_z(j)
        state%prof_U(:, j+1) = state%prof_U(:, j)
        state%prof_R(:, j+1) = state%prof_R(:, j)
        state%prof_L(:, j+1) = state%prof_L(:, j)
        j = j - 1
      end do
      state%prof_z(j+1) = zt; state%prof_U(:, j+1) = ut; state%prof_R(:, j+1) = rt
      state%prof_L(:, j+1) = lt
    end do
  end subroutine sem_sort_profile

  !> Linear interpolation in z, clamped to the table ends.
  subroutine sem_profile_at(state, z, u_mean, R, L)
    type(SEMState), intent(in)  :: state
    real(rp),       intent(in)  :: z
    real(rp),       intent(out) :: u_mean(3), R(6), L(3)

    integer(ip) :: k, n
    real(rp) :: w

    n = state%n_levels
    if (z <= state%prof_z(1)) then
      u_mean = state%prof_U(:, 1); R = state%prof_R(:, 1); L = state%prof_L(:, 1); return
    else if (z >= state%prof_z(n)) then
      u_mean = state%prof_U(:, n); R = state%prof_R(:, n); L = state%prof_L(:, n); return
    end if

    k = 1
    do while (k < n - 1 .and. state%prof_z(k+1) < z)
      k = k + 1
    end do
    w = (z - state%prof_z(k)) / max(state%prof_z(k+1) - state%prof_z(k), SMALL)
    u_mean = (1.0_rp - w) * state%prof_U(:, k) + w * state%prof_U(:, k+1)
    R      = (1.0_rp - w) * state%prof_R(:, k) + w * state%prof_R(:, k+1)
    L      = (1.0_rp - w) * state%prof_L(:, k) + w * state%prof_L(:, k+1)
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

  !> Add the curl of eddy i at point x to the profile-frame velocity ul.
  subroutine sem_eddy_add(state, i, x, ul)
    type(SEMState), intent(in)    :: state
    integer(ip),    intent(in)    :: i
    real(rp),       intent(in)    :: x(3)
    real(rp),       intent(inout) :: ul(3)

    real(rp) :: d(3), dl(3), q(3), s2, w, g(3), shift

    d = x - state%pos(:, i)
    if (state%periodic_t1) then
      shift = state%t1_period * anint(dot_product(d, state%t1) / state%t1_period)
      d = d - shift * state%t1
    end if
    dl(1) = dot_product(state%frame(:, 1), d)
    dl(2) = dot_product(state%frame(:, 2), d)
    dl(3) = dot_product(state%frame(:, 3), d)
    q = dl / state%sig(:, i)
    ! Compact support: reject on the bounding box before the square.
    if (abs(q(1)) >= 1.0_rp .or. abs(q(2)) >= 1.0_rp .or. abs(q(3)) >= 1.0_rp) return
    s2 = q(1)**2 + q(2)**2 + q(3)**2
    if (s2 >= 1.0_rp) return

    w = -6.0_rp * (1.0_rp - s2)**2
    g = w * q / state%sig(:, i)
    ul = ul + cross3(g, state%amp(:, i))
  end subroutine sem_eddy_add

  !> Bin geometry (once per setup) and counting sort of the eddies into (t1, t2, n_in) bins whose
  !> size is at least the largest kernel half-extent, so a point only ever needs its 3x3x3 neighbourhood.
  subroutine sem_bins_fill(state)
    type(SEMState), intent(inout) :: state

    integer(ip) :: i, ib, nbt, ia, ic, is, nb(3)
    real(rp)    :: lo(3), span(3), c(3), xr(3)
    integer(ip), allocatable :: cnt(:), fillpos(:), bid(:)

    if (.not. allocated(state%pos)) return
    lo   = [state%t1_min, state%t2_min, -state%s_half]
    span = [state%t1_max - state%t1_min, state%t2_max - state%t2_min, 2.0_rp * state%s_half]
    do i = 1, 3
      nb(i) = max(1_ip, int(span(i) / max(state%ext(i), SMALL), ip))
      nb(i) = min(nb(i), 256_ip)  ! keep the table small; larger bins only cost a few extra tests
    end do
    state%nbin  = nb
    state%bin_lo = lo
    state%bin_h  = span / real(nb, rp)
    nbt = nb(1) * nb(2) * nb(3)

    allocate(cnt(nbt), bid(state%n_eddies), source=0_ip)
    do i = 1, state%n_eddies
      xr = state%pos(:, i) - state%x_inlet
      c = [dot_product(xr, state%t1), dot_product(xr, state%t2), dot_product(xr, state%n_in)]
      ia = min(nb(1), max(1_ip, int(floor((c(1) - lo(1)) / state%bin_h(1)), ip) + 1_ip))
      ic = min(nb(2), max(1_ip, int(floor((c(2) - lo(2)) / state%bin_h(2)), ip) + 1_ip))
      is = min(nb(3), max(1_ip, int(floor((c(3) - lo(3)) / state%bin_h(3)), ip) + 1_ip))
      bid(i) = ia + nb(1) * ((ic - 1) + nb(2) * (is - 1))
      cnt(bid(i)) = cnt(bid(i)) + 1
    end do

    if (allocated(state%bin_start)) deallocate(state%bin_start)
    if (allocated(state%bin_list))  deallocate(state%bin_list)
    allocate(state%bin_start(nbt + 1), state%bin_list(state%n_eddies))
    state%bin_start(1) = 1
    do ib = 1, nbt
      state%bin_start(ib + 1) = state%bin_start(ib) + cnt(ib)
    end do
    allocate(fillpos(nbt))
    fillpos = state%bin_start(1:nbt)
    do i = 1, state%n_eddies
      state%bin_list(fillpos(bid(i))) = i
      fillpos(bid(i)) = fillpos(bid(i)) + 1
    end do
    state%bins_ready = .true.
  end subroutine sem_bins_fill

  !> Bin indices to visit along one axis: the point's bin and its neighbours, without repeats.
  subroutine bin_range(k, nb, periodic, idx, n)
    integer(ip), intent(in)  :: k, nb
    logical,     intent(in)  :: periodic
    integer(ip), intent(out) :: idx(3), n

    integer(ip) :: j
    n = 0
    if (periodic) then
      if (nb <= 3) then
        do j = 1, nb
          n = n + 1; idx(n) = j
        end do
      else
        do j = -1, 1
          n = n + 1; idx(n) = modulo(k - 1 + j, nb) + 1
        end do
      end if
    else
      do j = max(1_ip, k - 1), min(nb, k + 1)
        n = n + 1; idx(n) = j
      end do
    end if
  end subroutine bin_range

  !> Synthetic velocity fluctuation (global axes) at any point: the sum of the eddy curls.
  !> Exactly divergence-free in 3D; the stresses it carries are those of the eddies near x.
  !> Uses the spatial bins when they are current; otherwise (or for points well outside the
  !> eddy box) falls back to the full loop, which is the reference.
  subroutine sem_fluct_at_point(state, x, u_fluc)
    type(SEMState), intent(in)  :: state
    real(rp),       intent(in)  :: x(3)
    real(rp),       intent(out) :: u_fluc(3)

    integer(ip) :: i, k(3), ia(3), ib(3), ic(3), na, nbn, ns, j1, j2, j3, bin, m
    real(rp)    :: xr(3), c(3), ul(3), a1
    logical     :: brute

    if (.not. state%bins_ready) then
      call sem_fluct_at_point_ref(state, x, u_fluc)
      return
    end if

    xr = x - state%x_inlet
    c = [dot_product(xr, state%t1), dot_product(xr, state%t2), dot_product(xr, state%n_in)]
    if (state%periodic_t1) then
      a1 = c(1) - state%bin_lo(1)
      a1 = a1 - state%t1_period * floor(a1 / state%t1_period)
      c(1) = state%bin_lo(1) + a1
    end if
    k = int(floor((c - state%bin_lo) / state%bin_h), ip) + 1_ip
    ! Beyond one bin outside the box no eddy can reach the point only if the box really is the
    ! support region; to stay safe in that corner case use the reference loop.
    brute = (k(2) < 0 .or. k(2) > state%nbin(2) + 1 .or. k(3) < 0 .or. k(3) > state%nbin(3) + 1)
    if (.not. state%periodic_t1) brute = brute .or. k(1) < 0 .or. k(1) > state%nbin(1) + 1
    if (brute) then
      call sem_fluct_at_point_ref(state, x, u_fluc)
      return
    end if

    call bin_range(k(1), state%nbin(1), state%periodic_t1, ia, na)
    call bin_range(k(2), state%nbin(2), .false., ib, nbn)
    call bin_range(k(3), state%nbin(3), .false., ic, ns)
    ul = 0.0_rp
    do j3 = 1, ns
      do j2 = 1, nbn
        do j1 = 1, na
          bin = ia(j1) + state%nbin(1) * ((ib(j2) - 1) + state%nbin(2) * (ic(j3) - 1))
          do m = state%bin_start(bin), state%bin_start(bin + 1) - 1
            call sem_eddy_add(state, state%bin_list(m), x, ul)
          end do
        end do
      end do
    end do
    u_fluc = matmul(state%frame, ul)
  end subroutine sem_fluct_at_point

  !> Reference evaluation: loop over every eddy.
  subroutine sem_fluct_at_point_ref(state, x, u_fluc)
    type(SEMState), intent(in)  :: state
    real(rp),       intent(in)  :: x(3)
    real(rp),       intent(out) :: u_fluc(3)

    integer(ip) :: i
    real(rp) :: ul(3)

    ul = 0.0_rp
    do i = 1, state%n_eddies
      call sem_eddy_add(state, i, x, ul)
    end do
    u_fluc = matmul(state%frame, ul)
  end subroutine sem_fluct_at_point_ref

  !> Mean velocity and synthetic fluctuation at an inlet face.
  subroutine sem_velocity_at_face(state, abl_cfg, face_centre, u_mean, u_fluc_out)
    type(SEMState),  intent(in)  :: state
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(in)  :: face_centre(3)
    real(rp),        intent(out) :: u_mean(3)
    real(rp),        intent(out) :: u_fluc_out(3)

    real(rp) :: R(6), L(3)

    call sem_target_at(state, abl_cfg, face_centre(3), u_mean, R, L)
    call sem_fluct_at_point(state, face_centre, u_fluc_out)
  end subroutine sem_velocity_at_face

end module mod_sem_inflow
