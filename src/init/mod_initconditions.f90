! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Initial-condition setup routines.
module mod_initconditions
  use mod_kinds
  use mod_parameters
  use mod_fv_solver, only: FVSolver
  use mpi_f08
  implicit none
  private

  public :: apply_initial_conditions

contains

  subroutine report_u_bounds(slv, stage)
    type(FVSolver),   intent(in) :: slv
    character(len=*), intent(in) :: stage

    integer(ip) :: ierr_local
    real(rp) :: umin_l, umax_l, umin_g, umax_g

    if (slv%mesh%n_cells > 0) then
      umin_l = minval(slv%U(1, 1:slv%mesh%n_cells))
      umax_l = maxval(slv%U(1, 1:slv%mesh%n_cells))
    else
      umin_l = huge(1.0_rp)
      umax_l = -huge(1.0_rp)
    end if

    call MPI_Allreduce(umin_l, umin_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, slv%mesh%comm, ierr_local)
    call MPI_Allreduce(umax_l, umax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, slv%mesh%comm, ierr_local)

    if (slv%mesh%rank == 0) then
      write(*, '(A,A,A,2ES12.4)') ' IC: ', trim(stage), ' Umin/Umax = ', umin_g, umax_g
    end if
  end subroutine report_u_bounds

  subroutine apply_initial_conditions(slv, ic_type, ic_uref, ic_lref, noise_level, wall_normal_dir, noise_seed)
    type(FVSolver),    intent(inout) :: slv
    character(len=*),  intent(in)    :: ic_type
    real(rp),          intent(in)    :: ic_uref
    real(rp),          intent(in)    :: ic_lref
    real(rp),          intent(in)    :: noise_level
    integer(ip),       intent(in)    :: wall_normal_dir
    integer(ip), optional, intent(in) :: noise_seed  ! > 0: reproducible noise; 0/absent: clock

    real(rp) :: uref, lref_eff
    integer(ip) :: wnd
    real(rp) :: xmin_g(3), xmax_g(3)
    logical  :: need_bounds

    uref = ic_uref
    lref_eff = ic_lref
    wnd = max(1, min(3, wall_normal_dir))
    need_bounds = .false.

    select case (trim(ic_type))
    case ('linear', 'channel_linear', 'parabolic', 'channel_parabolic', 'poiseuille', &
          'loglaw', 'log_law', 'channel_loglaw', 'boundary_layer', 'channel_boundary_layer', 'bl')
      need_bounds = .true.
    end select
    if (need_bounds) call compute_domain_bounds(slv, xmin_g, xmax_g)

    select case (trim(ic_type))
    case ('uniform', '')
      ! Uniform streamwise velocity uref everywhere.
      slv%U(1, :) = uref
      slv%U(2, :) = 0.0_rp
      slv%U(3, :) = 0.0_rp
      if (slv%mesh%rank == 0) then
        write(*, '(A,ES12.4)') ' IC: uniform Uref = ', uref
      end if
    case ('taylor_green', 'TaylorGreen', 'TAYLOR_GREEN')
      call initialize_taylor_green_ic(slv, uref)
    case ('scalar_gauss', 'scalar_pulse')
      call initialize_scalar_advection_ic(slv, uref, trim(ic_type))
    case ('linear', 'channel_linear')
      call initialize_channel_profile_ic(slv, uref, lref_eff, 'linear', wnd, xmin_g, xmax_g)
    case ('parabolic', 'channel_parabolic', 'poiseuille')
      call initialize_channel_profile_ic(slv, uref, lref_eff, 'parabolic', wnd, xmin_g, xmax_g)
    case ('loglaw', 'log_law', 'channel_loglaw')
      call initialize_channel_profile_ic(slv, uref, lref_eff, 'loglaw', wnd, xmin_g, xmax_g)
    case ('boundary_layer', 'channel_boundary_layer', 'bl')
      call initialize_channel_profile_ic(slv, uref, lref_eff, 'boundary_layer', wnd, xmin_g, xmax_g)
    case default
      ! Unknown type - leave field as-is (solver BC-based init already ran).
    end select

    call report_u_bounds(slv, 'after profile')

    if (noise_level > 0.0_rp) then
      call add_random_velocity_noise(slv, noise_level, noise_seed)
      call report_u_bounds(slv, 'after noise')
      ! Extra projection pass: multiplicative noise perturbs dU_1/dx, making
      ! the Poisson RHS O(noisexuref/dx). With rtol=1e-8 the first pass leaves
      call slv%project_once(dt_for_projection(slv))
      call report_u_bounds(slv, 'after pre-projection')
    end if

    ! Always project (div U=0) before time loop; for RK3W use stage-1 dt so initial p is at
    ! stage-1 scale, avoiding a startup grad(p) residual that can destabilise transient flows.
    call slv%project_once(dt_for_projection(slv))
    call report_u_bounds(slv, 'after final projection')
  end subroutine apply_initial_conditions

  !> Returns dt for project_once scaled to match the first RK3W stage's incremental pre-correction.
  !> Uses GAMMA(1)*dt=(8/15)*dt so stage-1 correction exactly cancels the initial projection.
  pure real(rp) function dt_for_projection(slv) result(dt_proj)
    type(FVSolver), intent(in) :: slv
    dt_proj = RK3W_GAMMA(1) * slv%params%dt
  end function dt_for_projection

  subroutine initialize_channel_profile_ic(slv, uref, lref_in, profile_kind, wall_normal_dir, xmin_in, xmax_in)
    type(FVSolver),   intent(inout) :: slv
    real(rp),         intent(in)    :: uref
    real(rp),         intent(in)    :: lref_in
    character(len=*), intent(in)    :: profile_kind
    integer(ip),      intent(in)    :: wall_normal_dir
    real(rp),         intent(in), optional :: xmin_in(3), xmax_in(3)

    integer(ip) :: c, wnd, stream_dir
    real(rp) :: xmin(3), xmax(3), H, zeta, zc, ywall
    real(rp) :: lref, delta_bl, uval

    wnd = max(1, min(3, wall_normal_dir))
    stream_dir = 1
    if (stream_dir == wnd) stream_dir = 2

    if (present(xmin_in) .and. present(xmax_in)) then
      xmin = xmin_in
      xmax = xmax_in
    else
      call compute_domain_bounds(slv, xmin, xmax)
    end if
    H = max(xmax(wnd) - xmin(wnd), SMALL)
    if (lref_in > 0.0_rp) then
      lref = lref_in
    else
      lref = H
    end if
    lref = max(lref, SMALL)

    do c = 1, slv%mesh%n_cells
      zc = slv%mesh%cells(c)%centroid(wnd)
      zeta = (zc - xmin(wnd)) / lref
      zeta = max(0.0_rp, min(1.0_rp, zeta))

      select case (trim(profile_kind))
      case ('linear')
        uval = uref * zeta

      case ('parabolic')
        zeta = (zc - xmin(wnd)) / max(H, SMALL)
        zeta = max(0.0_rp, min(1.0_rp, zeta))
        uval = 4.0_rp * uref * zeta * (1.0_rp - zeta)

      case ('loglaw')
        ywall = min(zc - xmin(wnd), xmax(wnd) - zc)
        ywall = max(ywall, 0.0_rp)
        uval = uref * log(1.0_rp + 9.0_rp * min(1.0_rp, max(0.0_rp, ywall) / lref)) / log(10.0_rp)

      case ('boundary_layer')
        ywall = min(zc - xmin(wnd), xmax(wnd) - zc)
        ywall = max(ywall, 0.0_rp)
        delta_bl = 0.15_rp * lref
        uval = uref * (1.0_rp - exp(-max(0.0_rp, ywall) / (delta_bl + SMALL)))

      case default
        uval = 0.0_rp
      end select

      slv%U(:, c) = 0.0_rp
      slv%U(stream_dir, c) = uval
      slv%p(c)    = 0.0_rp
      slv%T(c)    = slv%params%T_ref
    end do

    if (slv%mesh%n_ghost > 0) then
      slv%U(:, slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%p(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%T(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = slv%params%T_ref
    end if

    call slv%mesh%exchange_halos_state(slv%U, slv%p, slv%T)

    if (slv%mesh%rank == 0) then
      write(*, '(A,A)') ' IC: channel profile set to ', trim(profile_kind)
      write(*, '(A,ES12.4)') ' IC: Uref = ', uref
      write(*, '(A,ES12.4)') ' IC: Lref = ', lref
      write(*, '(A,I0)') ' IC: wall_normal_direction = ', wnd
      write(*, '(A,3ES12.4)') ' IC: domain bounds xmin = ', xmin
      write(*, '(A,3ES12.4)') ' IC: domain bounds xmax = ', xmax
      write(*, '(A,ES12.4)') ' IC: wall half-height = ', (xmax(wnd) - xmin(wnd)) / 2.0_rp
    end if
  end subroutine initialize_channel_profile_ic

  ! Multiplicative relative noise: each velocity component is perturbed by
  ! +/-noise_frac of its own local value. Components that are zero (e.g.
  subroutine add_random_velocity_noise(slv, noise_frac, noise_seed)
    type(FVSolver), intent(inout) :: slv
    real(rp),       intent(in)    :: noise_frac  ! fraction of local |U_k| to perturb
    integer(ip), optional, intent(in) :: noise_seed

    integer(ip) :: c, k, seed_size, ierr_local
    integer, allocatable :: seed(:)
    integer :: seed_values(8), base_seed, max_seed
    real(rp) :: r

    call random_seed(size=seed_size)
    allocate(seed(seed_size))

    if (slv%mesh%rank == 0) then
      call date_and_time(values=seed_values)
      base_seed = seed_values(8)
      do k = 1, 8
        base_seed = modulo(1664525 * base_seed + seed_values(k) + 1013904223, huge(1) - 1)
      end do
      if (base_seed <= 0) base_seed = 13579
      if (present(noise_seed)) then
        if (noise_seed > 0) base_seed = noise_seed
      end if
    end if
    call MPI_Bcast(base_seed, 1, MPI_INTEGER, 0, slv%mesh%comm, ierr_local)

    max_seed = huge(1) - 1
    do k = 1, seed_size
      seed(k) = modulo(base_seed + 7919 * slv%mesh%rank + 104729 * k, max_seed)
      if (seed(k) <= 0) seed(k) = k
    end do
    call random_seed(put=seed)
    deallocate(seed)

    do c = 1, slv%mesh%n_cells
      do k = 1, 3
        call random_number(r)
        slv%U(k, c) = slv%U(k, c) * (1.0_rp + noise_frac * (2.0_rp * r - 1.0_rp))
      end do
    end do

    if (slv%mesh%n_ghost > 0) then
      slv%U(:, slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
    end if

    call slv%mesh%exchange_halos_vec(slv%U)

    if (slv%mesh%rank == 0) then
      write(*, '(A,ES12.4)') ' IC: added multiplicative noise fraction = ', noise_frac
    end if
  end subroutine add_random_velocity_noise

  subroutine compute_domain_bounds(slv, xmin_g, xmax_g)
    type(FVSolver), intent(in) :: slv
    real(rp),       intent(out) :: xmin_g(3), xmax_g(3)

    integer(ip) :: c, ierr_local
    real(rp) :: xmin_l(3), xmax_l(3)

    xmin_l = huge(1.0_rp)
    xmax_l = -huge(1.0_rp)

    if (allocated(slv%mesh%vertices) .and. slv%mesh%n_vertices > 0) then
      xmin_l = minval(slv%mesh%vertices, dim=2)
      xmax_l = maxval(slv%mesh%vertices, dim=2)
    else
      do c = 1, slv%mesh%n_cells
        xmin_l = min(xmin_l, slv%mesh%cells(c)%centroid)
        xmax_l = max(xmax_l, slv%mesh%cells(c)%centroid)
      end do
    end if

    call MPI_Allreduce(xmin_l, xmin_g, 3, MPI_DOUBLE_PRECISION, MPI_MIN, slv%mesh%comm, ierr_local)
    call MPI_Allreduce(xmax_l, xmax_g, 3, MPI_DOUBLE_PRECISION, MPI_MAX, slv%mesh%comm, ierr_local)
  end subroutine compute_domain_bounds

  subroutine initialize_taylor_green_ic(slv, amp)
    type(FVSolver), intent(inout) :: slv
    real(rp),       intent(in)    :: amp

    integer(ip) :: c, ierr_local
    real(rp)    :: xmin_l(3), xmax_l(3), xmin_g(3), xmax_g(3)
    real(rp)    :: L(3), xloc(3), k(3), u0
    real(rp)    :: theta_x, theta_y, theta_z

    xmin_l = huge(1.0_rp)
    xmax_l = -huge(1.0_rp)

    if (allocated(slv%mesh%vertices) .and. slv%mesh%n_vertices > 0) then
      xmin_l = minval(slv%mesh%vertices, dim=2)
      xmax_l = maxval(slv%mesh%vertices, dim=2)
    else
      do c = 1, slv%mesh%n_cells
        xmin_l = min(xmin_l, slv%mesh%cells(c)%centroid)
        xmax_l = max(xmax_l, slv%mesh%cells(c)%centroid)
      end do
    end if

    call MPI_Allreduce(xmin_l, xmin_g, 3, MPI_DOUBLE_PRECISION, MPI_MIN, slv%mesh%comm, ierr_local)
    call MPI_Allreduce(xmax_l, xmax_g, 3, MPI_DOUBLE_PRECISION, MPI_MAX, slv%mesh%comm, ierr_local)

    L = xmax_g - xmin_g
    where (L <= 0.0_rp) L = 1.0_rp
    k = TWOPI / L
    u0 = amp

    do c = 1, slv%mesh%n_cells
      xloc(1) = modulo(slv%mesh%cells(c)%centroid(1) - xmin_g(1), L(1))
      xloc(2) = modulo(slv%mesh%cells(c)%centroid(2) - xmin_g(2), L(2))
      xloc(3) = modulo(slv%mesh%cells(c)%centroid(3) - xmin_g(3), L(3))

      theta_x = k(1) * xloc(1)
      theta_y = k(2) * xloc(2)
      theta_z = k(3) * xloc(3)

      slv%U(1,c) = u0 * sin(theta_x) * cos(theta_y) * cos(theta_z)
      slv%U(2,c) = -u0 * cos(theta_x) * sin(theta_y) * cos(theta_z)
      slv%U(3,c) = 0.0_rp
      slv%p(c)   = 0.0_rp
      slv%T(c)   = slv%params%T_ref
    end do

    if (slv%mesh%n_ghost > 0) then
      slv%U(:, slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%p(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%T(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = slv%params%T_ref
    end if

    call slv%mesh%exchange_halos_state(slv%U, slv%p, slv%T)

    if (slv%mesh%rank == 0) then
      write(*, '(A,3F10.5)') ' TGV IC domain lengths:', L
      write(*, '(A,F10.5)')  ' TGV IC amplitude     :', u0
    end if
  end subroutine initialize_taylor_green_ic

  !> Passive-scalar advection test: uniform streamwise velocity uref and a temperature profile
  !> T = T_ref + f(x) along x on the periodic domain.
  !>   'scalar_gauss' : Gaussian, centre 0.25 L, sigma 0.05 L (smooth; for order-of-accuracy tests)
  !>   'scalar_pulse' : the same Gaussian plus a square pulse on [0.55, 0.75] L (for scheme comparison)
  subroutine initialize_scalar_advection_ic(slv, uref, kind)
    type(FVSolver),   intent(inout) :: slv
    real(rp),         intent(in)    :: uref
    character(len=*), intent(in)    :: kind

    integer(ip) :: c
    real(rp)    :: xmin_g(3), xmax_g(3), lx, s, d, f

    call compute_domain_bounds(slv, xmin_g, xmax_g)
    lx = max(xmax_g(1) - xmin_g(1), SMALL)

    do c = 1, slv%mesh%n_cells
      s = modulo(slv%mesh%cells(c)%centroid(1) - xmin_g(1), lx) / lx
      d = s - 0.25_rp
      d = d - anint(d)
      f = exp(-0.5_rp * (d / 0.05_rp)**2)
      if (trim(kind) == 'scalar_pulse') then
        if (s >= 0.55_rp .and. s < 0.75_rp) f = f + 1.0_rp
      end if
      slv%U(1, c) = uref
      slv%U(2, c) = 0.0_rp
      slv%U(3, c) = 0.0_rp
      slv%p(c)    = 0.0_rp
      slv%T(c)    = slv%params%T_ref + f
    end do

    if (slv%mesh%n_ghost > 0) then
      slv%U(:, slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%p(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = 0.0_rp
      slv%T(slv%mesh%n_cells+1:slv%mesh%n_cells+slv%mesh%n_ghost) = slv%params%T_ref
    end if
    call slv%mesh%exchange_halos_state(slv%U, slv%p, slv%T)

    if (slv%mesh%rank == 0) write(*, '(A,A,A,ES12.4)') ' IC: ', trim(kind), ', Uref = ', uref
  end subroutine initialize_scalar_advection_ic

end module mod_initconditions
