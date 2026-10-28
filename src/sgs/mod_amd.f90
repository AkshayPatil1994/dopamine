! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Anisotropic minimum-dissipation (AMD) subgrid-scale eddy-viscosity model (Rozema et al. 2015).
module mod_amd
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  use mod_vreman, only: ensure_filter_tensor, g_filt_arr
  use mod_wall_model, only: wall_stress_newton
  implicit none
  private

  public :: amd_sgs_viscosity
  public :: amd_point_viscosity

  ! Floor on alpha_ij alpha_ij below which the flow is treated as irrotational and nu_t is zero.
  real(rp), parameter :: AMD_ALPHA_SQ_MIN = 1.0e-20_rp

  ! The numerator cancels to round-off in near-irrotational flow: values below this fraction of the
  ! Smagorinsky-type scale tr(G) |alpha|^3 are dropped.
  real(rp), parameter :: AMD_ROUNDOFF_FRAC = 1.0e-12_rp

  ! Cap on nu_sgs as a fraction of the explicit diffusive stability limit, as in mod_vreman.
  real(rp), parameter :: NU_SGS_STAB_FRAC = 0.4_rp

contains

  !> AMD closure at one point: nu_t = C max(0, -beta_mn S_mn) / (alpha_ij alpha_ij), with
  !> beta_mn = G_lk alpha_lm alpha_kn and S the symmetric part of alpha_ij = du_j/dx_i.
  pure subroutine amd_point_viscosity(alpha, G, c_amd, nu_t)
    !$acc routine seq
    real(rp), intent(in)  :: alpha(3,3), G(3,3), c_amd
    real(rp), intent(out) :: nu_t

    integer(ip) :: i, j, m, n
    real(rp)    :: beta, num, alpha_sq, tr_g

    alpha_sq = 0.0_rp
    do j = 1, 3
      do i = 1, 3
        alpha_sq = alpha_sq + alpha(i, j)**2
      end do
    end do

    ! -beta_mn S_mn = -beta_mn alpha_mn, beta being symmetric.
    num = 0.0_rp
    do n = 1, 3
      do m = 1, 3
        beta = 0.0_rp
        do j = 1, 3
          do i = 1, 3
            beta = beta + G(i, j) * alpha(i, m) * alpha(j, n)
          end do
        end do
        num = num - beta * alpha(m, n)
      end do
    end do

    tr_g = G(1,1) + G(2,2) + G(3,3)
    if (alpha_sq > AMD_ALPHA_SQ_MIN .and. num > AMD_ROUNDOFF_FRAC * tr_g * alpha_sq * sqrt(alpha_sq)) then
      nu_t = c_amd * num / alpha_sq
    else
      nu_t = 0.0_rp
    end if
  end subroutine amd_point_viscosity

  !> Compute the AMD SGS viscosity at owned cells; wall_gradient uses the log-law gradient u_tau/(kappa y) in wall-model cells.
  subroutine amd_sgs_viscosity(U, msh, nu_sgs, nu, c_amd, wall_gradient, U_bnd, grad_U_ext, dt_stab)
    real(rp),   intent(in)    :: U(:,:)  ! (3, n_cells+n_ghost)
    type(Mesh), intent(inout), target :: msh  ! inout: the shared filter-tensor build halo-exchanges G
    real(rp),   intent(inout) :: nu_sgs(:)  ! (n_cells)
    real(rp),   intent(in)    :: nu         ! molecular viscosity, for the wall-law u_tau
    real(rp),   intent(in)    :: c_amd
    logical,    intent(in)    :: wall_gradient
    real(rp),   intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),   intent(in), optional, target :: grad_U_ext(:,:,:)  ! (3,3,n_cells[+n_ghost])
    real(rp),   intent(in), optional :: dt_stab  ! caps nu_sgs at NU_SGS_STAB_FRAC * V^(2/3) / dt_stab

    real(rp), allocatable, target :: grad_U_local(:,:,:)
    real(rp), pointer             :: grad_U(:,:,:)
    integer(ip) :: c, nc, f_bnd, f_global, o, n_int, n_bnd, n_patches, patch_id
    real(rp)    :: alpha(3,3), G(3,3), nu_c, nu_cap, dt_stab_local
    real(rp)    :: nrm(3), un, ut(3), ut_mag, y_w, u_tau, tau_mag, g_wall, d_old(3), d_t(3), dn
    logical     :: has_dt_stab, do_wall, capped
    integer(ip) :: i, j
    real(rp), pointer :: cvol_local(:), fnrm(:,:), fcen(:,:), ccen(:,:)
    integer(ip), pointer :: fown(:), fpid(:), pbc(:)

    if (present(grad_U_ext)) then
      grad_U => grad_U_ext(:, :, 1:msh%n_cells)
    else
      allocate(grad_U_local(3, 3, msh%n_cells))
      if (present(U_bnd)) then
        call green_gauss_gradient_vector(U, msh, grad_U_local, U_bnd)
      else
        call ls_gradient_vector(U, msh, grad_U_local)
      end if
      grad_U => grad_U_local
    end if

    call ensure_filter_tensor(msh)

    nc = msh%n_cells
    has_dt_stab = present(dt_stab)
    dt_stab_local = 1.0_rp
    if (has_dt_stab) dt_stab_local = dt_stab
    cvol_local => msh%cell_volume_soa

    n_int = msh%n_int_faces
    n_bnd = msh%n_faces - msh%n_int_faces
    n_patches = msh%n_patches
    fown => msh%face_owner_soa
    fpid => msh%face_patch_id_soa
    pbc  => msh%patch_bc_id_soa
    fnrm => msh%face_normal_soa
    fcen => msh%face_centre_soa
    ccen => msh%cell_centroid_soa
    do_wall = wall_gradient .and. n_bnd > 0

#ifdef _OPENACC
    !$acc data copyin(grad_U(:,:,1:nc), g_filt_arr(:,:,1:nc), cvol_local(1:nc)) &
    !$acc      copyout(nu_sgs(1:nc))
    !$acc parallel loop private(alpha, G, nu_c, nu_cap)
#else
    !$omp parallel do private(c, alpha, G, nu_c, nu_cap) schedule(static)
#endif
    do c = 1, nc
      alpha = grad_U(:, :, c)
      G     = g_filt_arr(:, :, c)
      call amd_point_viscosity(alpha, G, c_amd, nu_c)
      if (has_dt_stab) then
        nu_cap = NU_SGS_STAB_FRAC * cvol_local(c)**(2.0_rp/3.0_rp) / max(dt_stab_local, SMALL)
        nu_c = min(nu_c, nu_cap)
      end if
      nu_sgs(c) = nu_c
    end do
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end parallel do
#endif

    ! Wall-model cells, recomputed from the base gradient with the log-law wall-normal derivative. A cell
    ! owning two wall faces (a one-cell-thick channel) is overwritten by whichever face comes last.
    if (do_wall) then
#ifdef _OPENACC
      !$acc parallel loop default(present) copyin(U) &
      !$acc      private(f_global, patch_id, o, alpha, G, nrm, un, ut, ut_mag, y_w, u_tau, tau_mag, capped, &
      !$acc              g_wall, d_old, d_t, dn, nu_c, nu_cap, i, j)
#else
      !$omp parallel do private(f_bnd, f_global, patch_id, o, alpha, G, nrm, un, ut, ut_mag, y_w, u_tau, &
      !$omp&                    tau_mag, capped, g_wall, d_old, d_t, dn, nu_c, nu_cap, i, j) schedule(static)
#endif
      do f_bnd = 1, n_bnd
        f_global = n_int + f_bnd
        patch_id = fpid(f_global)
        if (patch_id < 1 .or. patch_id > n_patches) cycle
        if (pbc(patch_id) /= BC_WALL_MODEL) cycle
        o = fown(f_global)
        if (o < 1 .or. o > nc) cycle

        nrm = -fnrm(:, f_global)  ! unit normal pointing into the fluid
        y_w = abs(dot_product(fcen(:, f_global) - ccen(:, o), nrm))
        un  = dot_product(U(:, o), nrm)
        ut  = U(:, o) - un * nrm
        ut_mag = sqrt(ut(1)**2 + ut(2)**2 + ut(3)**2)
        if (ut_mag <= 1.0e-14_rp .or. y_w <= SMALL) cycle

        call wall_stress_newton(ut_mag, y_w, nu, u_tau, tau_mag, capped)
        g_wall = u_tau / (KAPPA * y_w)

        ! d_j = n_i alpha_ij is the wall-normal derivative of u_j; swap its wall-parallel part.
        alpha = grad_U(:, :, o)
        G     = g_filt_arr(:, :, o)
        do j = 1, 3
          d_old(j) = nrm(1)*alpha(1, j) + nrm(2)*alpha(2, j) + nrm(3)*alpha(3, j)
        end do
        dn  = dot_product(d_old, nrm)
        d_t = g_wall * ut / ut_mag - (d_old - dn * nrm)
        do j = 1, 3
          do i = 1, 3
            alpha(i, j) = alpha(i, j) + nrm(i) * d_t(j)
          end do
        end do

        call amd_point_viscosity(alpha, G, c_amd, nu_c)
        if (has_dt_stab) then
          nu_cap = NU_SGS_STAB_FRAC * cvol_local(o)**(2.0_rp/3.0_rp) / max(dt_stab_local, SMALL)
          nu_c = min(nu_c, nu_cap)
        end if
        nu_sgs(o) = nu_c
      end do
#ifdef _OPENACC
      !$acc end parallel loop
#else
      !$omp end parallel do
#endif
    end if
#ifdef _OPENACC
    !$acc end data
#endif

    nullify(cvol_local, fnrm, fcen, ccen, fown, fpid, pbc)
    if (allocated(grad_U_local)) deallocate(grad_U_local)
    nullify(grad_U)
  end subroutine amd_sgs_viscosity

end module mod_amd
