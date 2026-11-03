! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> MOST surface layer wall model with Businger-Dyer stability functions
module mod_abl_wall
  use mod_kinds
  use mod_parameters, only: BC_WALL_MODEL, SMALL, KAPPA, G_STD
  use mod_abl_types
  use mod_abl_profiles
  use mod_bc, only: BCConfig
  use mod_mesh
  use mpi_f08
  implicit none
  private

  real(rp), parameter :: L_CLIP = 5000.0_rp
  ! Inverse Stanton number kB^-1 = ln(z0_m/z0_h); ~2 over vegetated land.
  real(rp), parameter :: KB_INV = 2.0_rp

  public :: abl_wall_model

  integer(ip), allocatable, save :: wm_fbnd(:)  ! boundary-face index of each wall-model face
  real(rp),    allocatable, save :: wm_z0m(:), wm_z0h(:)  ! per-face roughness lengths
  logical, save :: wm_staged = .false.

contains

  !> Cache the wall-model boundary faces and their roughness; the mesh is static, so this runs once.
  subroutine wm_build_faces(msh, abl_cfg, bc_configs, dev)
    type(Mesh),      intent(in) :: msh
    type(ABLConfig), intent(in) :: abl_cfg
    type(BCConfig),  intent(in) :: bc_configs(:)
    logical,         intent(in) :: dev
    integer(ip) :: p, f, n
    real(rp)    :: z0_m_p, z0_h_p

    if (.not. allocated(wm_fbnd)) then
      n = 0
      do p = 1, msh%n_patches
        if (msh%patches(p)%bc_id == BC_WALL_MODEL) n = n + msh%patches(p)%n_faces
      end do
      allocate(wm_fbnd(n), wm_z0m(n), wm_z0h(n))
      n = 0
      do p = 1, msh%n_patches
        if (msh%patches(p)%bc_id /= BC_WALL_MODEL) cycle
        ! Patch roughness overrides the global ABL value when set
        call patch_roughness(bc_configs(p), abl_cfg, z0_m_p, z0_h_p)
        do f = 1, msh%patches(p)%n_faces
          n = n + 1
          wm_fbnd(n) = msh%patches(p)%start_face + f - 1
          wm_z0m(n)  = z0_m_p
          wm_z0h(n)  = z0_h_p
        end do
      end do
    end if
    if (dev .and. .not. wm_staged) then
      !$acc enter data copyin(wm_fbnd, wm_z0m, wm_z0h)
      wm_staged = .true.
    end if
  end subroutine wm_build_faces

  !> Compute MOST surface fluxes (u_star, T_star) at all wall-model faces; on_device: U, T and the outputs are resident.
  subroutine abl_wall_model(U, T, msh, abl_cfg, bc_configs, tau_w_bnd, q_w_bnd, L_bnd, update_L, on_device)
    real(rp),         intent(in)    :: U(:,:)  ! (3, nc+ng)
    real(rp),         intent(in)    :: T(:)  ! (nc+ng)
    type(Mesh),       intent(in), target :: msh
    type(ABLConfig),  intent(inout) :: abl_cfg
    type(BCConfig),   intent(in)    :: bc_configs(:)  ! per-patch roughness
    real(rp),         intent(inout) :: tau_w_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),         intent(inout) :: q_w_bnd(:)  ! (n_bnd_faces)
    real(rp),         intent(inout) :: L_bnd(:)  ! (n_bnd_faces) per-face Obukhov length
    logical,          intent(in), optional :: update_L  ! reduce L once per step, not per RK stage
    logical,          intent(in), optional :: on_device

    integer(ip) :: i, f_bnd, f_global, o, n_int, nwf
    real(rp)    :: n1, n2, n3, y_w, un, up1, up2, up3, U_par_mag, T_o
    real(rp)    :: u_star, T_star, L_new, L_guess, T_surf, L_cfg
    real(rp)    :: L_sum, rbuf(2), gbuf(2)
    integer(ip) :: ierr
    logical     :: do_update_L, dev
    integer(ip), pointer :: fown(:)
    real(rp),    pointer :: fnorm(:,:), fcen(:,:), ccen(:,:)

    do_update_L = .true.
    if (present(update_L)) do_update_L = update_L
    dev = .false.
    if (present(on_device)) dev = on_device
    call wm_build_faces(msh, abl_cfg, bc_configs, dev)
    nwf = size(wm_fbnd)
    n_int = msh%n_int_faces
    T_surf = abl_cfg%T_surface
    L_cfg = abl_cfg%L_obukhov
    fown  => msh%face_owner_soa
    fnorm => msh%face_normal_soa
    fcen  => msh%face_centre_soa
    ccen  => msh%cell_centroid_soa
    L_sum = 0.0_rp

    !$acc parallel loop if(dev) default(present) reduction(+:L_sum) &
    !$acc      private(f_bnd, f_global, o, n1, n2, n3, y_w, un, up1, up2, up3, U_par_mag, T_o, &
    !$acc              u_star, T_star, L_new, L_guess)
    do i = 1, nwf
      f_bnd    = wm_fbnd(i)
      f_global = n_int + f_bnd
      o        = fown(f_global)
      n1 = fnorm(1, f_global);  n2 = fnorm(2, f_global);  n3 = fnorm(3, f_global)

      ! Wall-normal distance from owner cell centroid to face
      y_w = abs((fcen(1, f_global) - ccen(1, o)) * n1 + (fcen(2, f_global) - ccen(2, o)) * n2 &
              + (fcen(3, f_global) - ccen(3, o)) * n3) + SMALL

      T_o = T(o)

      ! Tangential velocity (remove normal component)
      un = U(1, o) * n1 + U(2, o) * n2 + U(3, o) * n3
      up1 = U(1, o) - un * n1;  up2 = U(2, o) - un * n2;  up3 = U(3, o) - un * n3
      U_par_mag = sqrt(up1 * up1 + up2 * up2 + up3 * up3)

      ! Normalized tangential direction
      if (U_par_mag > SMALL) then
        up1 = up1 / U_par_mag;  up2 = up2 / U_par_mag;  up3 = up3 / U_par_mag
      else
        up1 = 0.0_rp;  up2 = 0.0_rp;  up3 = 0.0_rp
      end if

      ! Start from this face's own converged L; fall back to the config value on the first call.
      L_guess = L_bnd(f_bnd)
      if (abs(L_guess) < SMALL) L_guess = L_cfg

      ! Solve coupled MOST system for u_star and T_star
      call most_iterate(U_par_mag, T_o, T_surf, y_w, wm_z0m(i), wm_z0h(i), L_guess, u_star, T_star, L_new)

      ! Keep L local: a single domain mean is meaningless over heterogeneous surfaces.
      L_bnd(f_bnd) = L_new
      L_sum = L_sum + L_new

      ! Wall stress: tau_w = -u_star^2 * (U_par / |U_par|); the sign opposes the flow
      tau_w_bnd(1, f_bnd) = -u_star * u_star * up1
      tau_w_bnd(2, f_bnd) = -u_star * u_star * up2
      tau_w_bnd(3, f_bnd) = -u_star * u_star * up3

      ! Kinematic heat flux [K m/s]
      q_w_bnd(f_bnd) = -u_star * T_star
    end do

    ! Domain-mean Obukhov length: diagnostic and first-call seed only, the solve uses L_bnd.
    if (abl_cfg%compute_L .and. do_update_L) then
      rbuf(1) = L_sum
      rbuf(2) = real(nwf, rp)
      call MPI_Allreduce(rbuf, gbuf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, msh%comm, ierr)
      if (gbuf(2) > 0.0_rp) abl_cfg%L_obukhov = gbuf(1) / gbuf(2)
    end if

  end subroutine abl_wall_model

  !> Coupled Newton iteration for u_star and T_star (Monin-Obukhov).
  subroutine most_iterate(U_par_mag, T_o, T_surface, y_w, z0_m, z0_h, L_init, &
                           u_star_out, T_star_out, L_out)
    !$acc routine seq
    real(rp), intent(in)  :: U_par_mag, T_o, T_surface, y_w, z0_m, z0_h, L_init
    real(rp), intent(out) :: u_star_out, T_star_out, L_out

    real(rp) :: u_star, T_star, L, dFu_du, dFu_dT, dFT_du, dFT_dT
    real(rp) :: Fu, FT, Fu_p, FT_p, denom, du_star, dT_star
    integer(ip) :: iter, max_iter
    real(rp) :: T_ref, du_eps, dT_eps

    max_iter = 20
    T_ref = T_o

    ! Initial guess
    u_star = max(SMALL, KAPPA * U_par_mag / max(SMALL, log(y_w / z0_m)))
    T_star = 0.1_rp * (T_o - T_surface)
    ! L_obukhov = 0 means neutral, so start from the clip rather than a near-zero (very stable) L.
    if (abs(L_init) < SMALL) then
      L = L_CLIP
    else
      L = L_init
    end if

    ! Newton iteration
    do iter = 1, max_iter

      ! Clip L to avoid singularity
      L = sign(max(abs(L), 0.001_rp), L)
      if (abs(L) > L_CLIP) L = sign(L_CLIP, L)

      ! Residuals at the current state
      call most_residuals(u_star, T_star, L, U_par_mag, T_o, T_surface, &
                          y_w, z0_m, z0_h, Fu, FT)

      ! Check convergence
      if (abs(Fu) < 1.0e-8_rp .and. abs(FT) < 1.0e-8_rp) exit

      ! Finite-difference Jacobian (full coupling through L = -u*^3 T_ref/(kappa g T*)).
      du_eps = max(1.0e-6_rp, 1.0e-6_rp * abs(u_star))
      dT_eps = max(1.0e-6_rp, 1.0e-6_rp * abs(T_star))

      call most_residuals(u_star + du_eps, T_star, L_from(u_star + du_eps, T_star, T_ref), &
                          U_par_mag, T_o, T_surface, y_w, z0_m, z0_h, Fu_p, FT_p)
      dFu_du = (Fu_p - Fu) / du_eps
      dFT_du = (FT_p - FT) / du_eps

      call most_residuals(u_star, T_star + dT_eps, L_from(u_star, T_star + dT_eps, T_ref), &
                          U_par_mag, T_o, T_surface, y_w, z0_m, z0_h, Fu_p, FT_p)
      dFu_dT = (Fu_p - Fu) / dT_eps
      dFT_dT = (FT_p - FT) / dT_eps

      ! Solve linear system: J * [du, dT]^T = -[Fu, FT]^T
      denom = dFu_du * dFT_dT - dFu_dT * dFT_du
      if (abs(denom) < SMALL) exit

      du_star = -(dFT_dT * Fu - dFu_dT * FT) / denom
      dT_star = -(dFu_du * FT - dFT_du * Fu) / denom

      ! Update with relaxation for stability
      u_star = u_star + 0.5_rp * du_star
      T_star = T_star + 0.5_rp * dT_star

      ! Prevent negative u_star and cap it
      u_star = max(SMALL, min(u_star, 0.5_rp * U_par_mag))

      ! Update Obukhov length
      L = L_from(u_star, T_star, T_ref)

    end do

    u_star_out = u_star
    T_star_out = T_star
    L_out = L

  end subroutine most_iterate

  !> L = -u*^3 T/(kappa g <w'T'>) with <w'T'> = -u* T*, so the u* powers cancel to u*^2.
  pure real(rp) function L_from(u_star, T_star, T_ref) result(L)
    !$acc routine seq
    real(rp), intent(in) :: u_star, T_star, T_ref
    if (abs(T_star) > SMALL) then
      L = u_star**2 * T_ref / (KAPPA * G_STD * T_star)
      ! Keep L away from zero to avoid singular psi functions.
      L = sign(max(abs(L), 0.001_rp), L)
      if (abs(L) > L_CLIP) L = sign(L_CLIP, L)
    else
      ! Neutral: L -> inf; use finite proxy so psi ~= 0.
      L = L_CLIP
    end if
  end function L_from

  !> Resolve patch roughness; non-positive values inherit the global ABL setting.
  pure subroutine patch_roughness(bc, abl_cfg, z0_m, z0_h)
    type(BCConfig),  intent(in)  :: bc
    type(ABLConfig), intent(in)  :: abl_cfg
    real(rp),        intent(out) :: z0_m, z0_h

    z0_m = bc%abl_z0_m
    if (z0_m <= 0.0_rp) z0_m = abl_cfg%z0_m
    z0_m = max(z0_m, SMALL)

    z0_h = bc%abl_z0_h
    ! A non-positive thermal roughness means "derive from z0_m via kB^-1".
    if (z0_h <= 0.0_rp) z0_h = z0_m * exp(-KB_INV)
    z0_h = max(z0_h, SMALL)
  end subroutine patch_roughness

  subroutine most_residuals(u_star, T_star, L, U_par_mag, T_o, T_surface, &
                            y_w, z0_m, z0_h, Fu, FT)
    !$acc routine seq
    real(rp), intent(in)  :: u_star, T_star, L, U_par_mag, T_o, T_surface
    real(rp), intent(in)  :: y_w, z0_m, z0_h
    real(rp), intent(out) :: Fu, FT
    real(rp) :: psi_m_z, psi_m_z0m, psi_h_z, psi_h_z0h, ln_ratio_m, ln_ratio_h
    real(rp) :: denom_m, denom_h
    psi_m_z   = abl_psi_m(y_w / L)
    psi_m_z0m = abl_psi_m(z0_m / L)
    psi_h_z   = abl_psi_h(y_w / L)
    psi_h_z0h = abl_psi_h(z0_h / L)
    ln_ratio_m = log(y_w / z0_m)
    ln_ratio_h = log(y_w / z0_h)
    denom_m = ln_ratio_m - psi_m_z + psi_m_z0m
    denom_h = ln_ratio_h - psi_h_z + psi_h_z0h
    if (abs(denom_m) < SMALL) denom_m = sign(SMALL, denom_m)
    if (abs(denom_h) < SMALL) denom_h = sign(SMALL, denom_h)
    Fu = u_star - KAPPA * U_par_mag        / denom_m
    FT = T_star - KAPPA * (T_o - T_surface)/ denom_h
  end subroutine most_residuals

end module mod_abl_wall
