! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Viscous flux assembly for momentum and scalar transport.
module mod_viscous
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  implicit none
  private

  public :: viscous_rhs_vector
  public :: viscous_rhs_scalar

  ! Persistent scratch buffers, same rationale as mod_gradient.f90: mesh topology is static
  ! for a run, so grow these once instead of allocate/deallocate every RK substage.
  integer(ip), allocatable, target, save :: fcp_identity_buf(:)  ! shared by both routines below
  real(rp),    allocatable, save :: vec_U_bnd_buf(:,:)           ! viscous_rhs_vector, (3, n)
  real(rp),    allocatable, save :: sca_phi_bnd_buf(:)           ! viscous_rhs_scalar
  real(rp),    allocatable, target, save :: visc_inv_vol_buf(:)  ! fallback 1/volume when inv_vol not passed in

contains

  !> Jasak min-correction orthogonal coefficient A*cos(theta); limiting the
  !> dimensionless cosine (not the m^2 projection) keeps the clamp unit-safe.
  pure real(rp) function ortho_weight(s_dot_dhat, area) result(w)
    !$acc routine seq
    real(rp), intent(in) :: s_dot_dhat, area
    w = area * max(0.0_rp, min(1.0_rp, s_dot_dhat / (area + SMALL)))
  end function ortho_weight

  !> Viscous RHS for momentum: rhs(k,c) = (1/V_c) sum_f F_visc_k_f
  subroutine viscous_rhs_vector(U, nu, nu_sgs, grad_U, msh, rhs, U_bnd, inv_vol)
    real(rp),   intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    real(rp),   intent(in)    :: nu  ! molecular kinematic viscosity
    real(rp),   intent(in)    :: nu_sgs(:)  ! (n_cells + n_ghost) SGS viscosity
    real(rp),   intent(in)    :: grad_U(:,:,:)  ! (3, 3, n_cells + n_ghost)
    type(Mesh), intent(in), target :: msh
    real(rp),   intent(inout) :: rhs(:,:)  ! (3, n_cells)
    real(rp),   intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),   intent(in), optional, target :: inv_vol(:)  ! (n_cells) cached 1/volume

    integer(ip) :: f, o, nb, k, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd
    real(rp)    :: nu_eff_f, dU_n(3), stress_f(3), sf(3)
    real(rp)    :: n_hat(3), grad_Uo(3,3), grad_Unb(3,3)
    real(rp)    :: grad_Uf(3,3), d_on(3), d_on_len, gc
    real(rp)    :: s_dot_dhat, dn_perp(3), d_per, gc_per, weight_ortho
    integer(ip) :: col, idx, n_loops, lo, hi, nc, ng, nf, n_int, nbnd, n_patches
    integer(ip) :: n_bnd_buf
    logical     :: use_color, has_bnd
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:)
    real(rp),    pointer :: sf_soa(:,:), gc_soa(:), fnorm_soa(:,:), farea_soa(:)
    real(rp),    pointer :: fcentre_soa(:,:), fdcf_soa(:), ppair_dcf_soa(:)
    real(rp),    pointer :: ccent_soa(:,:), cvol_soa(:), vinv(:)

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'viscous_rhs_vector: msh%soa_ready is false; call mesh%build_soa_mirror first'

    rhs = 0.0_rp
    nc = msh%n_cells
    ng = msh%n_ghost
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches
    fo_soa      => msh%face_owner_soa
    fn_soa      => msh%face_neighbour_soa
    sf_soa      => msh%face_sf_soa
    gc_soa      => msh%face_gc_soa
    fnorm_soa   => msh%face_normal_soa
    farea_soa   => msh%face_area_soa
    fcentre_soa => msh%face_centre_soa
    fdcf_soa    => msh%face_dcf_soa
    fpatch_soa  => msh%face_patch_id_soa
    pbc_soa     => msh%patch_bc_id_soa
    pneigh_soa  => msh%periodic_neighbor_soa
    ppair_soa   => msh%periodic_pair_soa
    ppair_dcf_soa => msh%periodic_pair_dcf_soa
    ccent_soa   => msh%cell_centroid_soa
    cvol_soa    => msh%cell_volume_soa

    if (present(inv_vol)) then
      vinv => inv_vol
    else
      call ensure_buf_1d(visc_inv_vol_buf, nc)
      visc_inv_vol_buf(1:nc) = 1.0_rp / cvol_soa(1:nc)
      vinv => visc_inv_vol_buf(1:nc)
    end if

    ! U_bnd copied to a plain (always-present) local array: an absent optional
    ! can't be referenced inside a statically-mapped !$acc region.
    has_bnd = present(U_bnd)
    n_bnd_buf = max(merge(size(U_bnd, 2), nbnd, has_bnd), 1_ip)
    call ensure_buf_2d3(vec_U_bnd_buf, n_bnd_buf)
    if (has_bnd) then
      vec_U_bnd_buf(:, 1:n_bnd_buf) = U_bnd
    else
      vec_U_bnd_buf(:, 1:n_bnd_buf) = 0.0_rp
    end if

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets) .and. msh%soa_ready)
    n_loops = merge(msh%n_colors, 1_ip, use_color)

    if (use_color) then
      fcp => msh%face_color_perm
    else
      call ensure_fcp_identity(fcp_identity_buf, nf)
      fcp => fcp_identity_buf(1:nf)
    end if

#ifdef _OPENACC
    !$acc data copyin(U(:,:), nu_sgs(:), grad_U(:,:,:), vec_U_bnd_buf(:,1:n_bnd_buf), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), sf_soa(:,1:nf), gc_soa(1:nf), &
    !$acc      fnorm_soa(:,1:nf), farea_soa(1:nf), fcentre_soa(:,1:nf), fdcf_soa(1:nf), &
    !$acc      fpatch_soa(1:nf), pbc_soa(1:max(n_patches,1_ip)), &
    !$acc      pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      ppair_dcf_soa(1:max(nbnd,1_ip)), ccent_soa(:,1:nc+ng), vinv(1:nc), &
    !$acc      fcp(1:nf)) &
    !$acc      copy(rhs(:,1:nc)) if(use_color)
#endif
    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) private(f, o, nb, k, f_bnd, patch_id, bc_id, &
      !$acc      pair_face, nb_per, pair_bnd, nu_eff_f, dU_n, stress_f, sf, n_hat, gc, &
      !$acc      grad_Uo, grad_Unb, grad_Uf, d_on, d_on_len, s_dot_dhat, dn_perp, &
      !$acc      d_per, gc_per, weight_ortho)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, k, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd, &
      !$omp&         nu_eff_f, dU_n, stress_f, sf, n_hat, gc, grad_Uo, grad_Unb, &
      !$omp&         grad_Uf, d_on, d_on_len, s_dot_dhat, dn_perp, d_per, gc_per, weight_ortho)
#endif
      do idx = lo, hi
        if (use_color) then
          f = fcp(idx)
        else
          f = idx
        end if
        o  = fo_soa(f)
        nb = fn_soa(f)
        sf = sf_soa(:, f)
        n_hat = fnorm_soa(:, f)

        if (nb >= 1 .and. nb <= nc + ng) then
          ! Owned and ghost neighbours both valid (nu_sgs and grad_U halo-exchanged).
          gc = gc_soa(f)
          grad_Uo  = grad_U(:, :, o)
          grad_Unb = grad_U(:, :, nb)
          ! gc-weighted mean SGS viscosity; nu_sgs is a continuous field, not a material property.
          nu_eff_f = nu + gc * nu_sgs(o) + (1.0_rp - gc) * nu_sgs(nb)

          ! gc-weighted mean gradient at face (Jasak 1996 Sect. 3.3.1)
          grad_Uf = gc * grad_Uo + (1.0_rp - gc) * grad_Unb

          ! Jasak min-correction split; smooth limiter avoids hard discontinuity at Sf.d_on=0.
          d_on       = ccent_soa(:, nb) - ccent_soa(:, o)
          d_on_len   = sqrt(sum(d_on**2)) + SMALL
          s_dot_dhat = dot_product(sf, d_on) / d_on_len
          dn_perp    = sf - (s_dot_dhat / d_on_len) * d_on
          weight_ortho = ortho_weight(s_dot_dhat, farea_soa(f))

          do k = 1, 3
            stress_f(k) = nu_eff_f * ( &
                 (U(k, nb) - U(k, o)) * weight_ortho / d_on_len &
               + dot_product(grad_Uf(:, k), dn_perp) &
               + dot_product(grad_Uf(k, :), sf) )
          end do
        else
          ! Boundary face
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          bc_id = -1
          if (patch_id >= 1 .and. patch_id <= n_patches) bc_id = pbc_soa(patch_id)

          if (bc_id == BC_PERIODIC) then
            nb_per = -1
            pair_face = -1
            if (f_bnd >= 1 .and. f_bnd <= nbnd) then
              nb_per = pneigh_soa(f_bnd)
              pair_bnd = ppair_soa(f_bnd)
              if (pair_bnd > 0) then
                pair_face = n_int + pair_bnd
                if (pair_face >= 1 .and. pair_face <= nf) then
                  if (nb_per <= 0) nb_per = fo_soa(pair_face)
                else
                  pair_face = -1
                end if
              end if
              if (nb_per <= 0 .and. n_int + f_bnd <= nf) nb_per = fn_soa(n_int + f_bnd)
              if (nb_per == o) nb_per = -1
            end if

            if (nb_per >= 1 .and. nb_per <= nc + ng) then
              ! Compute gc_per from d_cf distances (mesh%faces%gc = 1 for boundary/periodic faces).
              if (pair_face >= 1 .and. pair_face <= nf) then
                d_per  = fdcf_soa(f) + fdcf_soa(pair_face)
                gc_per = fdcf_soa(pair_face) / (d_per + SMALL)
              else if (f_bnd >= 1 .and. f_bnd <= nbnd .and. ppair_dcf_soa(f_bnd) > SMALL) then
                d_per  = fdcf_soa(f) + ppair_dcf_soa(f_bnd)
                gc_per = ppair_dcf_soa(f_bnd) / (d_per + SMALL)
              else
                d_per  = 2.0_rp * fdcf_soa(f)
                gc_per = 0.5_rp
              end if
              ! gc-weighted mean SGS viscosity on periodic faces.
              nu_eff_f = nu + gc_per * nu_sgs(o) + (1.0_rp - gc_per) * nu_sgs(nb_per)
              ! Jasak min-correction on periodic faces: d_on = d_per * n_hat.
              d_on       = d_per * n_hat
              d_on_len   = d_per + SMALL
              s_dot_dhat = dot_product(sf, d_on) / d_on_len
              dn_perp    = sf - (s_dot_dhat / d_on_len) * d_on
              grad_Uo    = grad_U(:, :, o)
              grad_Unb   = grad_U(:, :, nb_per)
              grad_Uf    = gc_per * grad_Uo + (1.0_rp - gc_per) * grad_Unb
              weight_ortho = ortho_weight(s_dot_dhat, farea_soa(f))
              do k = 1, 3
                stress_f(k) = nu_eff_f * ( &
                     (U(k, nb_per) - U(k, o)) * weight_ortho / d_on_len &
                   + dot_product(grad_Uf(:, k), dn_perp) &
                   + dot_product(grad_Uf(k, :), sf) )
              end do
            else
              nu_eff_f = nu + nu_sgs(o)
              if (has_bnd) then
                dU_n = (vec_U_bnd_buf(:, max(f_bnd,1_ip)) - U(:, o)) / (fdcf_soa(f) + SMALL)
                stress_f = nu_eff_f * dU_n * farea_soa(f)
              else
                stress_f = 0.0_rp
              end if
            end if
          else if (bc_id == BC_WALL_MODEL) then
            ! All wall stress is applied by apply_wall_stress_to_predictor.
            stress_f = 0.0_rp
          else
            ! Covers BC_WALL_NOSLIP/VELOCITY/SLIP/SYMMETRY/INLET/OUTLET/ABL_INLET/
            ! FIXEDTEMP/ZEROGRADIENT and the invalid-patch fallback alike (Jasak split).
            nu_eff_f = nu + nu_sgs(o)
            if (has_bnd) then
              d_on       = fcentre_soa(:, f) - ccent_soa(:, o)
              d_on_len   = sqrt(sum(d_on**2)) + SMALL
              s_dot_dhat = dot_product(sf, d_on) / d_on_len
              dn_perp    = sf - (s_dot_dhat / d_on_len) * d_on
              grad_Uo    = grad_U(:, :, o)
              weight_ortho = ortho_weight(s_dot_dhat, farea_soa(f))
              do k = 1, 3
                stress_f(k) = nu_eff_f * ( &
                     (vec_U_bnd_buf(k, max(f_bnd,1_ip)) - U(k, o)) * weight_ortho / d_on_len &
                   + dot_product(grad_Uo(:, k), dn_perp) &
                   + dot_product(grad_Uo(k, :), sf) )
              end do
            else
              stress_f = 0.0_rp
            end if
          end if
        end if

        if (o  >= 1 .and. o  <= nc) rhs(:, o)  = rhs(:, o)  + stress_f * vinv(o)
        if (nb >= 1 .and. nb <= nc) rhs(:, nb) = rhs(:, nb) - stress_f * vinv(nb)
      end do
#ifdef _OPENACC
      !$acc end parallel loop
#else
      !$omp end parallel do
#endif
    end do
#ifdef _OPENACC
    !$acc end data
#endif

    nullify(fo_soa, fn_soa, fcp, sf_soa, gc_soa, fnorm_soa, farea_soa, fcentre_soa, &
            fdcf_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, ppair_dcf_soa, &
            ccent_soa, cvol_soa, vinv)
  end subroutine viscous_rhs_vector

  !> Viscous RHS for a scalar (temperature): rhs(c) = (1/V_c) sum_f F_T_f
  subroutine viscous_rhs_scalar(phi, alpha_mol, nu_sgs, Pr_t, grad_phi, msh, rhs, phi_bnd, inv_vol)
    real(rp),   intent(in)    :: phi(:)  ! (n_cells + n_ghost)
    real(rp),   intent(in)    :: alpha_mol  ! molecular diffusivity
    real(rp),   intent(in)    :: nu_sgs(:)  ! (n_cells + n_ghost) SGS viscosity
    real(rp),   intent(in)    :: Pr_t  ! turbulent Prandtl number (~0.85)
    real(rp),   intent(in)    :: grad_phi(:,:)  ! (3, n_cells + n_ghost) cell gradients
    type(Mesh), intent(in), target :: msh
    real(rp),   intent(inout) :: rhs(:)  ! (n_cells)
    real(rp),   intent(in), optional :: phi_bnd(:)  ! (n_bnd_faces)
    real(rp),   intent(in), optional, target :: inv_vol(:)  ! (n_cells) cached 1/volume

    integer(ip) :: f, o, nb, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd
    real(rp)    :: alpha_eff, dphi_n, flux, sf(3), d_on(3), gc
    real(rp)    :: dn_len, s_dot_dhat, dn_perp(3), d_per, gc_per
    real(rp)    :: grad_o(3), grad_nb(3), grad_f(3)
    integer(ip) :: col, idx, n_loops, lo, hi, nc, ng, nf, n_int, nbnd, n_patches
    integer(ip) :: n_bnd_buf
    logical     :: use_color, has_bnd
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:)
    real(rp),    pointer :: sf_soa(:,:), gc_soa(:), fnorm_soa(:,:), farea_soa(:)
    real(rp),    pointer :: fcentre_soa(:,:), fdcf_soa(:), ppair_dcf_soa(:)
    real(rp),    pointer :: ccent_soa(:,:), cvol_soa(:), vinv(:)

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'viscous_rhs_scalar: msh%soa_ready is false; call mesh%build_soa_mirror first'

    rhs = 0.0_rp
    nc = msh%n_cells
    ng = msh%n_ghost
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches
    fo_soa      => msh%face_owner_soa
    fn_soa      => msh%face_neighbour_soa
    sf_soa      => msh%face_sf_soa
    gc_soa      => msh%face_gc_soa
    fnorm_soa   => msh%face_normal_soa
    farea_soa   => msh%face_area_soa
    fcentre_soa => msh%face_centre_soa
    fdcf_soa    => msh%face_dcf_soa
    fpatch_soa  => msh%face_patch_id_soa
    pbc_soa     => msh%patch_bc_id_soa
    pneigh_soa  => msh%periodic_neighbor_soa
    ppair_soa   => msh%periodic_pair_soa
    ppair_dcf_soa => msh%periodic_pair_dcf_soa
    ccent_soa   => msh%cell_centroid_soa
    cvol_soa    => msh%cell_volume_soa

    if (present(inv_vol)) then
      vinv => inv_vol
    else
      call ensure_buf_1d(visc_inv_vol_buf, nc)
      visc_inv_vol_buf(1:nc) = 1.0_rp / cvol_soa(1:nc)
      vinv => visc_inv_vol_buf(1:nc)
    end if

    has_bnd = present(phi_bnd)
    n_bnd_buf = max(merge(size(phi_bnd), nbnd, has_bnd), 1_ip)
    call ensure_buf_1d(sca_phi_bnd_buf, n_bnd_buf)
    if (has_bnd) then
      sca_phi_bnd_buf(1:n_bnd_buf) = phi_bnd
    else
      sca_phi_bnd_buf(1:n_bnd_buf) = 0.0_rp
    end if

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets) .and. msh%soa_ready)
    n_loops = merge(msh%n_colors, 1_ip, use_color)

    if (use_color) then
      fcp => msh%face_color_perm
    else
      call ensure_fcp_identity(fcp_identity_buf, nf)
      fcp => fcp_identity_buf(1:nf)
    end if

#ifdef _OPENACC
    !$acc data copyin(phi(:), nu_sgs(:), grad_phi(:,:), sca_phi_bnd_buf(1:n_bnd_buf), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), sf_soa(:,1:nf), gc_soa(1:nf), &
    !$acc      fnorm_soa(:,1:nf), farea_soa(1:nf), fcentre_soa(:,1:nf), fdcf_soa(1:nf), &
    !$acc      fpatch_soa(1:nf), pbc_soa(1:max(n_patches,1_ip)), &
    !$acc      pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      ppair_dcf_soa(1:max(nbnd,1_ip)), ccent_soa(:,1:nc+ng), vinv(1:nc), &
    !$acc      fcp(1:nf)) &
    !$acc      copy(rhs(1:nc)) if(use_color)
#endif
    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) private(f, o, nb, f_bnd, patch_id, bc_id, &
      !$acc      pair_face, nb_per, pair_bnd, alpha_eff, dphi_n, flux, sf, d_on, gc, &
      !$acc      dn_len, s_dot_dhat, dn_perp, d_per, gc_per, grad_o, grad_nb, grad_f)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd, &
      !$omp&         alpha_eff, dphi_n, flux, sf, d_on, gc, dn_len, s_dot_dhat, &
      !$omp&         dn_perp, d_per, gc_per, grad_o, grad_nb, grad_f)
#endif
      do idx = lo, hi
        if (use_color) then
          f = fcp(idx)
        else
          f = idx
        end if
        o  = fo_soa(f)
        nb = fn_soa(f)
        sf = sf_soa(:, f)

        if (nb >= 1 .and. nb <= nc + ng) then
          ! nu_sgs and grad_phi carry valid halo values (exchanged in solver).
          gc = gc_soa(f)
          ! gc-weighted mean SGS viscosity in thermal diffusivity.
          alpha_eff = alpha_mol + (gc * nu_sgs(o) + (1.0_rp - gc) * nu_sgs(nb)) / Pr_t
          grad_nb   = grad_phi(:, nb)
          grad_o    = grad_phi(:, o)
          d_on       = ccent_soa(:, nb) - ccent_soa(:, o)
          dn_len     = sqrt(sum(d_on**2)) + SMALL
          s_dot_dhat = dot_product(sf, d_on) / dn_len
          dn_perp    = sf - (s_dot_dhat / dn_len) * d_on
          ! gc-weighted mean gradient at face
          grad_f     = gc * grad_o + (1.0_rp - gc) * grad_nb
          dphi_n     = (phi(nb) - phi(o)) * ortho_weight(s_dot_dhat, farea_soa(f)) / dn_len &
                     + dot_product(grad_f, dn_perp)
          flux = alpha_eff * dphi_n
        else
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          bc_id = -1
          if (patch_id >= 1 .and. patch_id <= n_patches) bc_id = pbc_soa(patch_id)

          if (bc_id == BC_PERIODIC) then
            nb_per = -1
            pair_face = -1
            if (f_bnd >= 1 .and. f_bnd <= nbnd) then
              nb_per = pneigh_soa(f_bnd)
              pair_bnd = ppair_soa(f_bnd)
              if (pair_bnd > 0) then
                pair_face = n_int + pair_bnd
                if (pair_face >= 1 .and. pair_face <= nf) then
                  if (nb_per <= 0) nb_per = fo_soa(pair_face)
                else
                  pair_face = -1
                end if
              end if
              if (nb_per <= 0 .and. n_int + f_bnd <= nf) nb_per = fn_soa(n_int + f_bnd)
              if (nb_per == o) nb_per = -1
            end if

            if (nb_per >= 1 .and. nb_per <= nc + ng) then
              ! Compute gc_per from d_cf distances (mesh%faces%gc = 1 for boundary/periodic faces).
              if (pair_face >= 1 .and. pair_face <= nf) then
                d_per  = fdcf_soa(f) + fdcf_soa(pair_face)
                gc_per = fdcf_soa(pair_face) / (d_per + SMALL)
              else if (f_bnd >= 1 .and. f_bnd <= nbnd .and. ppair_dcf_soa(f_bnd) > SMALL) then
                d_per  = fdcf_soa(f) + ppair_dcf_soa(f_bnd)
                gc_per = ppair_dcf_soa(f_bnd) / (d_per + SMALL)
              else
                d_per  = 2.0_rp * fdcf_soa(f)
                gc_per = 0.5_rp
              end if
              ! gc-weighted mean SGS viscosity on periodic temperature face.
              alpha_eff = alpha_mol + (gc_per * nu_sgs(o) + (1.0_rp - gc_per) * nu_sgs(nb_per)) / Pr_t
              d_on       = d_per * fnorm_soa(:, f)
              dn_len     = d_per + SMALL
              s_dot_dhat = dot_product(sf, d_on) / dn_len
              dn_perp    = sf - (s_dot_dhat / dn_len) * d_on
              grad_o     = grad_phi(:, o)
              grad_nb    = grad_phi(:, nb_per)
              ! gc_per-weighted mean gradient at periodic face
              grad_f     = gc_per * grad_o + (1.0_rp - gc_per) * grad_nb
              dphi_n     = (phi(nb_per) - phi(o)) * ortho_weight(s_dot_dhat, farea_soa(f)) / dn_len &
                         + dot_product(grad_f, dn_perp)
              flux       = alpha_eff * dphi_n
            else
              alpha_eff = alpha_mol + nu_sgs(o)/Pr_t
              if (has_bnd) then
                ! Jasak split at boundary face
                d_on       = fcentre_soa(:, f) - ccent_soa(:, o)
                dn_len     = sqrt(sum(d_on**2)) + SMALL
                s_dot_dhat = dot_product(sf, d_on) / dn_len
                dn_perp    = sf - (s_dot_dhat / dn_len) * d_on
                grad_o     = grad_phi(:, o)
                dphi_n     = (sca_phi_bnd_buf(max(f_bnd,1_ip)) - phi(o)) * ortho_weight(s_dot_dhat, farea_soa(f)) / dn_len &
                           + dot_product(grad_o, dn_perp)
                flux       = alpha_eff * dphi_n
              else
                flux = 0.0_rp
              end if
            end if
          else
            ! Covers every non-periodic BC alike (no separate wall-stress carve-out) plus invalid-patch fallback.
            alpha_eff = alpha_mol + nu_sgs(o)/Pr_t
            if (has_bnd) then
              d_on       = fcentre_soa(:, f) - ccent_soa(:, o)
              dn_len     = sqrt(sum(d_on**2)) + SMALL
              s_dot_dhat = dot_product(sf, d_on) / dn_len
              dn_perp    = sf - (s_dot_dhat / dn_len) * d_on
              grad_o     = grad_phi(:, o)
              dphi_n     = (sca_phi_bnd_buf(max(f_bnd,1_ip)) - phi(o)) * ortho_weight(s_dot_dhat, farea_soa(f)) / dn_len &
                         + dot_product(grad_o, dn_perp)
              flux       = alpha_eff * dphi_n
            else
              flux = 0.0_rp
            end if
          end if
        end if

        if (o  >= 1 .and. o  <= nc) rhs(o)  = rhs(o)  + flux * vinv(o)
        if (nb >= 1 .and. nb <= nc) rhs(nb) = rhs(nb) - flux * vinv(nb)
      end do
#ifdef _OPENACC
      !$acc end parallel loop
#else
      !$omp end parallel do
#endif
    end do
#ifdef _OPENACC
    !$acc end data
#endif

    nullify(fo_soa, fn_soa, fcp, sf_soa, gc_soa, fnorm_soa, farea_soa, fcentre_soa, &
            fdcf_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, ppair_dcf_soa, &
            ccent_soa, cvol_soa, vinv)
  end subroutine viscous_rhs_scalar

  ! ensure_buf_1d/ensure_buf_2d3/ensure_fcp_identity now live in mod_gradient (used here via `use mod_gradient`).

end module mod_viscous
