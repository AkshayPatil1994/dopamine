! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Equilibrium log-law wall model (Newton iteration for u_tau).
module mod_wall_model
  use mod_kinds
  use mod_parameters
  use mod_mesh
  implicit none
  private

  public :: compute_wall_stress
  public :: apply_wall_model_bc
  public :: wall_model_free
  public :: wall_model_reset_cap_count
  public :: wall_model_get_cap_count
  public :: reichardt_uplus
  public :: C_BUF

  ! Reichardt (1951) buffer-layer constants; A2 sets the near-wall exponential decay.
  real(rp), parameter :: REICH_A1 = 11.0_rp
  real(rp), parameter :: REICH_A2 = 0.33_rp
  real(rp), parameter :: REICH_A1_INV = 1.0_rp / REICH_A1

  ! Buffer constant fixed by the log-law asymptote, so Reichardt matches KAPPA and BLOG.
  real(rp), parameter :: C_BUF = BLOG - log(KAPPA) / KAPPA

  !> Counts matching points below y+ = 2; non-zero means the mesh is too fine for a wall model.
  integer(ip), save :: utau_cap_count = 0_ip

  !> Per-boundary-face sampling cache; built once on first call.
  real(rp),    allocatable, save :: y_w_arr(:)       ! wall-normal distance to the sampling cell
  integer(ip), allocatable, save :: sample_cell(:)   ! cell index sampled by the wall model
  real(rp),    allocatable, save :: u_par_filt(:,:)  ! (3, n_bnd) time-filtered sampled velocity

  !> False until the time filter has been seeded with its first instantaneous sample.
  logical, save :: filt_ready = .false.

contains

  !> Reichardt (1951) law of the wall: u+ as a function of y+, valid from the wall through the log layer.
  pure subroutine reichardt_uplus(yplus, c_buf, uplus, duplus)
    real(rp), intent(in)  :: yplus, c_buf
    real(rp), intent(out) :: uplus, duplus
    real(rp) :: e1, e2

    e1  = exp(-yplus * REICH_A1_INV)
    e2  = exp(-REICH_A2 * yplus)

    uplus  = log(1.0_rp + KAPPA * yplus) / KAPPA &
           + c_buf * (1.0_rp - e1 - yplus * REICH_A1_INV * e2)
    duplus = 1.0_rp / (1.0_rp + KAPPA * yplus) &
           + c_buf * (REICH_A1_INV * e1 - REICH_A1_INV * e2 + yplus * REICH_A1_INV * REICH_A2 * e2)
  end subroutine reichardt_uplus

  !> Compute friction velocity and shear stress magnitude from wall-parallel velocity.
  subroutine compute_wall_stress(u_par, y_w, nu, u_tau, tau_w_mag)
    real(rp), intent(in)  :: u_par, y_w, nu
    real(rp), intent(out) :: u_tau, tau_w_mag

    real(rp) :: u_tau_new, yplus, f, df, uplus, duplus, inv_nu, y_w_inv_nu
    integer(ip) :: iter
    integer(ip), parameter :: MAX_ITER = 50
    real(rp), parameter :: TOL = 1.0e-8_rp

    if (u_par <= SMALL .or. y_w <= SMALL) then
      u_tau = 0.0_rp
      tau_w_mag = 0.0_rp
      return
    end if

    ! Initial guess via Blasius
    u_tau = max(0.01_rp * u_par, sqrt(nu * u_par / (y_w + SMALL)))

    inv_nu     = 1.0_rp / (nu + SMALL)
    y_w_inv_nu = y_w * inv_nu  ! loop-invariant: nu and y_w are fixed across the Newton iteration

    do iter = 1, MAX_ITER
      yplus = u_tau * y_w_inv_nu
      call reichardt_uplus(yplus, C_BUF, uplus, duplus)

      ! F(u_tau) = u_par/u_tau - u+(y+); dF/du_tau adds the chain term y_w/nu.
      f  =  u_par / (u_tau + SMALL) - uplus
      df = -u_par / (u_tau**2 + SMALL) - duplus * y_w_inv_nu

      u_tau_new = u_tau - f / (df + sign(SMALL, df))
      ! Bound the step: F is monotone, so a large jump means an overshoot.
      u_tau_new = min(max(u_tau_new, 0.1_rp * u_tau), 10.0_rp * u_tau)

      if (abs(u_tau_new - u_tau) < TOL * u_tau_new) then
        u_tau = u_tau_new
        exit
      end if
      u_tau = u_tau_new
    end do

    ! Diagnostic only: counts faces where the matching point falls in the viscous sublayer
    ! (mesh too fine for a wall model there); u_tau itself is not clamped.
    if (u_tau * y_w_inv_nu < 1.0_rp) utau_cap_count = utau_cap_count + 1_ip

    tau_w_mag = u_tau**2
  end subroutine compute_wall_stress

  !> Walk inward from the wall face to the first cell at least h_wm from the wall.
  !> Returns the wall-adjacent owner when h_wm <= 0 or the walk cannot proceed.
  subroutine find_sample_cell(msh, f_global, h_wm, cell_out, y_out)
    type(Mesh),  intent(in)  :: msh
    integer(ip), intent(in)  :: f_global
    real(rp),    intent(in)  :: h_wm
    integer(ip), intent(out) :: cell_out
    real(rp),    intent(out) :: y_out

    integer(ip), parameter :: MAX_WALK = 20
    integer(ip) :: c, step, fi, fj, other, best
    real(rp)    :: n_in(3), x_w(3), dist, best_dist

    n_in     = -msh%faces(f_global)%normal  ! inward
    x_w      = msh%faces(f_global)%centre
    cell_out = msh%faces(f_global)%owner
    y_out    = SMALL
    if (cell_out < 1 .or. cell_out > size(msh%cells)) return

    y_out = abs(dot_product(msh%cells(cell_out)%centroid - x_w, n_in)) + SMALL
    if (h_wm <= 0.0_rp .or. .not. allocated(msh%cell_faces)) return

    c = cell_out
    do step = 1, MAX_WALK
      if (y_out >= h_wm) exit
      ! Ghost cells carry valid halo data but no face list, so the walk stops there.
      if (c > msh%n_cells) exit
      best = -1_ip
      best_dist = y_out
      do fi = msh%cells(c)%face_offset, msh%cells(c)%face_offset + msh%cells(c)%n_faces - 1
        fj = msh%cell_faces(fi)
        if (fj < 1 .or. fj > size(msh%faces)) cycle
        if (msh%faces(fj)%neighbour < 1) cycle
        other = msh%faces(fj)%owner
        if (other == c) other = msh%faces(fj)%neighbour
        if (other < 1 .or. other > size(msh%cells)) cycle
        dist = dot_product(msh%cells(other)%centroid - x_w, n_in)
        if (dist > best_dist) then
          best_dist = dist
          best = other
        end if
      end do
      if (best < 1) exit
      c        = best
      cell_out = best
      y_out    = best_dist
    end do
  end subroutine find_sample_cell

  !> Apply the Reichardt wall model to all BC_WALL_MODEL patches.
  !> h_wm > 0 samples at that wall distance; t_filt > 0 exponentially filters the sampled velocity.
  subroutine apply_wall_model_bc(U, nu, msh, tau_w_bnd, h_wm, t_filt, dt)
    real(rp),   intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    real(rp),   intent(in)    :: nu
    type(Mesh), intent(in)    :: msh
    real(rp),   intent(inout) :: tau_w_bnd(:,:)  ! (3, n_bnd_faces) output
    real(rp),   intent(in), optional :: h_wm    ! sampling height; <= 0 uses the wall-adjacent cell
    real(rp),   intent(in), optional :: t_filt  ! filter time constant; <= 0 disables filtering
    real(rp),   intent(in), optional :: dt      ! time increment for this filter update

    integer(ip) :: p, f, f_bnd, f_global, o, c_s
    real(rp)    :: n_hat(3), u_par_vec(3), u_par_mag
    real(rp)    :: y_w, u_tau, tau_mag, h_use, t_use, dt_use, w_filt
    integer(ip) :: n_bnd
    logical     :: needs_rebuild

    h_use  = 0.0_rp; if (present(h_wm))   h_use  = h_wm
    t_use  = 0.0_rp; if (present(t_filt)) t_use  = t_filt
    dt_use = 0.0_rp; if (present(dt))     dt_use = dt

    ! Rebuild the sampling cache if uninitialized or the mesh size has changed.
    ! .or. does not guarantee short-circuit in Fortran, so size(y_w_arr) must not be
    ! evaluated in the same expression as the allocated() check.
    n_bnd = msh%n_faces - msh%n_int_faces
    needs_rebuild = .not. allocated(y_w_arr)
    if (.not. needs_rebuild) needs_rebuild = (size(y_w_arr) /= n_bnd)
    if (needs_rebuild) then
      if (allocated(y_w_arr))      deallocate(y_w_arr)
      if (allocated(sample_cell))  deallocate(sample_cell)
      if (allocated(u_par_filt))   deallocate(u_par_filt)
      if (n_bnd > 0) then
        allocate(y_w_arr(n_bnd), source=0.0_rp)
        allocate(sample_cell(n_bnd), source=-1_ip)
        allocate(u_par_filt(3, n_bnd), source=0.0_rp)
        do p = 1, msh%n_patches
          if (msh%patches(p)%bc_id /= BC_WALL_MODEL) cycle
          do f = 1, msh%patches(p)%n_faces
            f_bnd    = msh%patches(p)%start_face + f - 1
            f_global = msh%n_int_faces + f_bnd
            if (f_bnd < 1 .or. f_bnd > n_bnd) cycle
            call find_sample_cell(msh, f_global, h_use, sample_cell(f_bnd), y_w_arr(f_bnd))
          end do
        end do
      end if
    end if

    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_WALL_MODEL) cycle

      do f = 1, msh%patches(p)%n_faces
        f_bnd    = msh%patches(p)%start_face + f - 1
        f_global = msh%n_int_faces + f_bnd
        o        = msh%faces(f_global)%owner
        n_hat    = msh%faces(f_global)%normal  ! outward unit normal

        c_s = o
        y_w = SMALL
        ! Split into nested ifs: .and. does not short-circuit in Fortran, so size(y_w_arr)/
        ! sample_cell(f_bnd) must never be evaluated in the same expression as allocated().
        if (allocated(sample_cell)) then
          if (f_bnd >= 1 .and. f_bnd <= size(y_w_arr)) then
            if (sample_cell(f_bnd) >= 1) c_s = sample_cell(f_bnd)
            y_w = y_w_arr(f_bnd)
          end if
        end if

        u_par_vec = U(:, c_s) - dot_product(U(:, c_s), n_hat) * n_hat

        ! Exponential time filter (Yang et al. 2017) damps the instantaneous sampled velocity.
        if (t_use > SMALL .and. allocated(u_par_filt)) then
          if (.not. filt_ready) then
            u_par_filt(:, f_bnd) = u_par_vec  ! seed, else tau_w ramps from zero
          else
            w_filt = dt_use / (t_use + dt_use)
            u_par_filt(:, f_bnd) = u_par_filt(:, f_bnd) + w_filt * (u_par_vec - u_par_filt(:, f_bnd))
          end if
          u_par_vec = u_par_filt(:, f_bnd)
        end if

        u_par_mag = sqrt(sum(u_par_vec**2))

        call compute_wall_stress(u_par_mag, y_w, nu, u_tau, tau_mag)

        ! Wall shear stress vector opposes the wall-parallel velocity.
        tau_w_bnd(:, f_bnd) = (-tau_mag / (u_par_mag + SMALL)) * u_par_vec

      end do
    end do

    if (t_use > SMALL) filt_ready = .true.
  end subroutine apply_wall_model_bc

  subroutine wall_model_reset_cap_count()
    utau_cap_count = 0_ip
  end subroutine wall_model_reset_cap_count

  integer(ip) function wall_model_get_cap_count() result(n)
    n = utau_cap_count
  end function wall_model_get_cap_count

  !> Free module-level persistent allocations (call on solver teardown).
  subroutine wall_model_free()
    if (allocated(y_w_arr))     deallocate(y_w_arr)
    if (allocated(sample_cell)) deallocate(sample_cell)
    if (allocated(u_par_filt))  deallocate(u_par_filt)
    filt_ready = .false.
  end subroutine wall_model_free

end module mod_wall_model
