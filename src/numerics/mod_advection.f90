! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Advection flux assembly: CENTRAL/KEP, QUICK, VANLEER_TVD, KOREN_TVD schemes.
module mod_advection
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: advection_rhs_scalar
  public :: advection_rhs_vector
  public :: advection_reset_tvd_fallback
  ! Exposed so the reconstruction linear-exactness test can call them directly.
  public :: quick_reconstruction
  public :: kep_reconstruction
  public :: advection_get_tvd_fallback_count
  public :: advection_reset_ke_budget
  public :: advection_get_ke_budget
  public :: advection_set_ke_budget

  !> KEP keeps the plain arithmetic mean; the midpoint correction trades conservation for order.
  logical, parameter, public :: KEP_PURE_MEAN = .true.

  !> Count of TVD faces per step that fell back to 1st-order upwind (missing far-upwind stencil).
  integer(ip), save :: tvd_fallback_count = 0_ip

  !> Per-face convective KE production -F_f*phi_f.(u_o-u_nb), split smooth/refinement; exactly zero for ideal KEP.
  real(rp), save :: ke_prod_smooth = 0.0_rp
  real(rp), save :: ke_prod_refine = 0.0_rp
  !> KE carried by the split form's advective half; equals the divergence error the projection left.
  real(rp), save :: ke_prod_split  = 0.0_rp
  !> Boundary-face KE flux and the exact rate it must equal; their difference is the conservation error.
  real(rp), save :: ke_prod_bnd    = 0.0_rp
  real(rp), save :: ke_bnd_exact   = 0.0_rp
  logical,  save :: ke_budget_on   = .false.

  !> Per-cell discrete divergence sum_f sigma*F_f, rebuilt each KEP call for the split form.
  real(rp), allocatable, save :: div_work(:)

  ! Persistent scratch buffers grown once instead of allocate/deallocate every RK substage, same rationale as mod_gradient.f90/mod_viscous.f90; calls into this module must stay serial, never concurrent, since these buffers are shared mutable state.
  integer(ip), allocatable, target, save :: fcp_identity_buf(:)
  real(rp),    allocatable, save :: cds_phi_bnd_buf(:)        ! advection_rhs_scalar_cds_soa
  real(rp),    allocatable, save :: recon_s_phi_bnd_buf(:)    ! advection_rhs_scalar_recon_soa
  real(rp),    allocatable, save :: recon_v_Ut_buf(:,:)       ! advection_rhs_vector_recon_soa, (n, 3)
  real(rp),    allocatable, save :: recon_v_gradt_buf(:,:,:)  ! advection_rhs_vector_recon_soa, (3, nc, 3)

contains

  !> Assemble advective RHS for a scalar field phi:
  !> rhs(c) = -1/V_c * sum_f F_f
  subroutine advection_rhs_scalar(phi, U, msh, scheme, rhs, phi_bnd, U_bnd, face_flux, inv_vol, grad_phi)
    real(rp),    intent(in)    :: phi(:)  ! (n_cells + n_ghost)
    real(rp),    intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    type(Mesh),  intent(in)    :: msh
    integer(ip), intent(in)    :: scheme  ! SCHEME_QUICK, SCHEME_VANLEER_TVD, SCHEME_KOREN_TVD
    real(rp),    intent(inout) :: rhs(:)  ! (n_cells)
    real(rp),    intent(in), optional :: phi_bnd(:)  ! (n_bnd_faces)
    real(rp),    intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),    intent(in), optional :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in), optional :: inv_vol(:)  ! (n_cells) cached 1/volume
    real(rp),    intent(in), optional :: grad_phi(:,:)  ! (3, n_cells) for O(h^2) CDS eps-correction on skewed meshes

    integer(ip) :: f, o, nb, f_bnd, patch_id, nb_per
    real(rp)    :: Fn, phi_f, phi_f_plain, w_blend
    integer(ip) :: col, idx, n_loops, lo, hi
    logical     :: use_color

    ! GPU/SoA fast path: CDS with grad_phi (eps-correction) and a projected face_flux, the
    ! production configuration for temperature/scalar transport, has a fixed owner/neighbour
    ! (+periodic) stencil and no reconstruction-function calls.
    if (scheme == SCHEME_CDS .and. present(grad_phi) .and. present(face_flux) .and. present(inv_vol) &
        .and. msh%soa_ready) then
      call advection_rhs_scalar_cds_soa(phi, msh, rhs, phi_bnd, face_flux, inv_vol, grad_phi)
      return
    end if

    ! GPU/SoA fast path: QUICK/VANLEER_TVD/KOREN_TVD, the same production configuration
    ! (grad_phi + projected face_flux) as CDS above, via the cell_face_offset_soa CSR mirror.
    if ((scheme == SCHEME_QUICK .or. scheme == SCHEME_VANLEER_TVD .or. scheme == SCHEME_KOREN_TVD) &
        .and. present(grad_phi) .and. present(face_flux) .and. present(inv_vol) .and. msh%soa_ready) then
      call advection_rhs_scalar_recon_soa(phi, msh, scheme, rhs, phi_bnd, face_flux, inv_vol, grad_phi)
      return
    end if

    rhs = 0.0_rp

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets))
    if (use_color) then
      n_loops = msh%n_colors
    else
      n_loops = 1
    end if

    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = msh%n_faces
      end if
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, f_bnd, patch_id, nb_per, Fn, phi_f, phi_f_plain, w_blend)
      do idx = lo, hi
        if (use_color) then
          f = msh%face_color_perm(idx)
        else
          f = idx
        end if
      f_bnd = 0
      patch_id = 0
      nb_per = -1
      o  = msh%faces(f)%owner
      nb = msh%faces(f)%neighbour
      if (nb < 1) f_bnd = f - msh%n_int_faces

      ! Face-normal volumetric flux
      ! Use projected face flux on interior and periodic faces.
      if (nb >= 1) then
        if (present(face_flux) .and. f >= 1 .and. f <= size(face_flux)) then
          Fn = face_flux(f)
        else
          Fn = dot_product(face_interpolate_vector(U, f, msh), msh%faces(f)%sf)
        end if
      else
        patch_id = msh%faces(f)%patch_id
        if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
          ! Reuse projected flux on periodic/outlet faces (single conservative mass flux).
          if (present(face_flux) .and. f >= 1 .and. f <= size(face_flux) .and. &
              (msh%patches(patch_id)%bc_id == BC_PERIODIC .or. &
               msh%patches(patch_id)%bc_id == BC_OUTLET)) then
            Fn = face_flux(f)
          else if (present(U_bnd)) then
            Fn = dot_product(U_bnd(:, f_bnd), msh%faces(f)%sf)
          else
            Fn = dot_product(U(:, o), msh%faces(f)%sf)
          end if
        else if (present(U_bnd)) then
          Fn = dot_product(U_bnd(:, f_bnd), msh%faces(f)%sf)
        else
          Fn = dot_product(U(:, o), msh%faces(f)%sf)
        end if
      end if
      if (.not. present(face_flux) .and. nb < 1 .and. f_bnd >= 1) then
        patch_id = msh%faces(f)%patch_id
        if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
          if (msh%patches(patch_id)%bc_id == BC_PERIODIC) then
            nb_per = resolve_periodic_nb(msh, f_bnd, o)
            if (nb_per >= 1 .and. nb_per <= size(U, 2)) then
              Fn = dot_product(0.5_rp*(U(:, o) + U(:, nb_per)), msh%faces(f)%sf)
            else if (present(U_bnd)) then
              Fn = dot_product(U_bnd(:, f_bnd), msh%faces(f)%sf)
            else
              Fn = dot_product(U(:, o), msh%faces(f)%sf)
            end if
          else
            select case (msh%patches(patch_id)%bc_id)
            case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
              ! Solid walls and slip boundaries: zero mass flux (no normal flow)
              Fn = 0.0_rp
            case (BC_INLET, BC_OUTLET, BC_ABL_INLET, BC_FIXEDTEMP, BC_ZEROGRADIENT)
              ! Open boundaries and inlet/outlet: advective mass flux from boundary velocity
              if (present(U_bnd)) then
                Fn = dot_product(U_bnd(:, f_bnd), msh%faces(f)%sf)
              else
                Fn = dot_product(U(:, o), msh%faces(f)%sf)
              end if
            case default
              ! Unknown BC type: fallback to boundary velocity flux
              if (present(U_bnd)) then
                Fn = dot_product(U_bnd(:, f_bnd), msh%faces(f)%sf)
              else
                Fn = dot_product(U(:, o), msh%faces(f)%sf)
              end if
            end select
          end if
        end if
      end if

      ! Reconstruct phi at face; CDS with grad_phi uses eps-correction (Kim & Choi 2000, Eq. 26);
      ! TVD/QUICK use directional upwind stencils (eps-correction non-standard).
      if (nb >= 1 .or. face_has_reconstruction_neighbor(f, msh)) then
        if (scheme == SCHEME_CDS .and. present(grad_phi)) then
          phi_f = face_interpolate_scalar(phi, f, msh, grad_phi=grad_phi)
          ! Fade the eps-correction out toward refinement interfaces; see refine_blend_weight.
          if (f <= msh%n_int_faces) then
            w_blend = refine_blend_weight(msh%faces(f)%skew)
            if (w_blend < 1.0_rp) then
              phi_f_plain = face_interpolate_scalar(phi, f, msh)
              phi_f = phi_f_plain + w_blend * (phi_f - phi_f_plain)
            end if
          end if
        else if (present(grad_phi)) then
          phi_f = reconstruct_face_value(phi, f, Fn, msh, scheme, grad_phi)
        else
          phi_f = reconstruct_face_value(phi, f, Fn, msh, scheme)
        end if
        if (.not. ieee_is_finite(phi_f)) then
          if (Fn >= 0.0_rp) then
            phi_f = phi(o)
          else if (nb >= 1 .and. nb <= size(phi)) then
            phi_f = phi(nb)
          else
            phi_f = phi(o)
          end if
        end if
      else
        patch_id = msh%faces(f)%patch_id
        if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
          if (msh%patches(patch_id)%bc_id == BC_PERIODIC) then
            ! F5: nb_per resolved earlier in this face iteration; reuse.
            if (nb_per >= 1 .and. nb_per <= size(phi)) then
              phi_f = 0.5_rp * (phi(o) + phi(nb_per))
            else if (present(phi_bnd)) then
              phi_f = phi_bnd(f_bnd)
            else
              phi_f = phi(o)
            end if
          else if (present(phi_bnd)) then
            phi_f = phi_bnd(f_bnd)
          else
            phi_f = phi(o)
          end if
        else if (present(phi_bnd)) then
          phi_f = phi_bnd(f_bnd)
        else
          phi_f = phi(o)
        end if
      end if

      if (.not. ieee_is_finite(Fn)) Fn = 0.0_rp
      if (.not. ieee_is_finite(phi_f)) phi_f = 0.0_rp

      ! Accumulate flux
      if (o >= 1 .and. o <= msh%n_cells) then
        if (present(inv_vol)) then
          rhs(o) = rhs(o) - Fn * phi_f * inv_vol(o)
        else
          rhs(o) = rhs(o) - Fn * phi_f / msh%cells(o)%volume
        end if
      end if

      if (nb >= 1 .and. nb <= msh%n_cells) then
        if (present(inv_vol)) then
          rhs(nb) = rhs(nb) + Fn * phi_f * inv_vol(nb)
        else
          rhs(nb) = rhs(nb) + Fn * phi_f / msh%cells(nb)%volume
        end if
      end if
    end do
      !$omp end parallel do
    end do
  end subroutine advection_rhs_scalar

  !> Device-safe finite check: ieee_is_finite itself is not device-linkable under nvfortran's OpenACC codegen.
  logical function soa_is_finite(x) result(ok)
    !$acc routine seq
    real(rp), intent(in) :: x
    ok = (x == x) .and. (abs(x) <= huge(x))
  end function soa_is_finite

  !> CDS fast path via the flat SoA mesh mirror, for scheme=SCHEME_CDS with grad_phi and a
  !> projected face_flux present (the production configuration): fixed owner/neighbour(+periodic)
  !> stencil, no reconstruction-function calls, matches face_interpolate_scalar's math exactly.
  subroutine advection_rhs_scalar_cds_soa(phi, msh, rhs, phi_bnd, face_flux, inv_vol, grad_phi)
    real(rp),    intent(in)    :: phi(:)  ! (n_cells + n_ghost)
    type(Mesh),  intent(in), target :: msh
    real(rp),    intent(inout) :: rhs(:)  ! (n_cells)
    real(rp),    intent(in), optional :: phi_bnd(:)  ! (n_bnd_faces)
    real(rp),    intent(in)    :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in)    :: inv_vol(:)  ! (n_cells) cached 1/volume
    real(rp),    intent(in)    :: grad_phi(:,:)  ! (3, n_cells) for O(h^2) eps-correction

    integer(ip) :: f, o, nb, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd, nb_simple
    real(rp)    :: Fn, phi_f, phi_f_plain, gc, epsilon(3), grad_mean(3), d_on(3), w_blend
    integer(ip) :: col, idx, n_loops, lo, hi, nc, nf, n_int, nbnd, n_patches
    integer(ip) :: n_bnd_buf
    logical     :: use_color, has_bnd, has_reconstruction
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:)
    real(rp),    pointer :: gc_soa(:), cvol_soa(:), fcentre_soa(:,:), ccent_soa(:,:)
    real(rp),    pointer :: eps_soa(:,:), skew_soa(:), pgc_soa(:)

    rhs = 0.0_rp
    nc = msh%n_cells
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches

    fo_soa      => msh%face_owner_soa
    fn_soa      => msh%face_neighbour_soa
    gc_soa      => msh%face_gc_soa
    pgc_soa     => msh%periodic_gc_soa
    cvol_soa    => msh%cell_volume_soa
    fpatch_soa  => msh%face_patch_id_soa
    pbc_soa     => msh%patch_bc_id_soa
    pneigh_soa  => msh%periodic_neighbor_soa
    ppair_soa   => msh%periodic_pair_soa
    fcentre_soa => msh%face_centre_soa
    ccent_soa   => msh%cell_centroid_soa
    eps_soa     => msh%face_eps_soa
    skew_soa    => msh%face_skew_soa

    has_bnd = present(phi_bnd)
    n_bnd_buf = max(merge(size(phi_bnd), nbnd, has_bnd), 1_ip)
    call ensure_buf_1d(cds_phi_bnd_buf, n_bnd_buf)
    if (has_bnd) then
      cds_phi_bnd_buf(1:n_bnd_buf) = phi_bnd
    else
      cds_phi_bnd_buf(1:n_bnd_buf) = 0.0_rp
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
    !$acc data copyin(phi(:), face_flux(1:nf), cds_phi_bnd_buf(1:n_bnd_buf), inv_vol(1:nc), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), gc_soa(1:nf), pgc_soa(1:max(nbnd,1_ip)), cvol_soa(1:nc), &
    !$acc      fpatch_soa(1:nf), pbc_soa(1:max(n_patches,1_ip)), &
    !$acc      pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      fcentre_soa(:,1:nf), ccent_soa(:,1:nc+msh%n_ghost), eps_soa(:,1:nf), &
    !$acc      grad_phi(:,1:nc+msh%n_ghost), &
    !$acc      skew_soa(1:nf), fcp(1:nf)) &
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
      !$acc      pair_face, nb_per, pair_bnd, nb_simple, Fn, phi_f, phi_f_plain, gc, epsilon, &
      !$acc      grad_mean, d_on, w_blend, has_reconstruction)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd, &
      !$omp&         nb_simple, Fn, phi_f, phi_f_plain, gc, epsilon, grad_mean, d_on, w_blend, &
      !$omp&         has_reconstruction)
#endif
      do idx = lo, hi
        if (use_color) then
          f = fcp(idx)
        else
          f = idx
        end if
        o  = fo_soa(f)
        nb = fn_soa(f)
        f_bnd = 0
        patch_id = 0
        bc_id = -1
        nb_per = -1

        ! Face-normal flux: projected face_flux on interior/periodic/outlet faces,
        ! boundary-velocity flux (zero on walls) elsewhere. Matches advection_rhs_scalar
        ! exactly for the present(face_flux) case (the only case this fast path handles).
        if (nb >= 1) then
          Fn = face_flux(f)
        else
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          if (patch_id >= 1 .and. patch_id <= n_patches) then
            bc_id = pbc_soa(patch_id)
            if (bc_id == BC_PERIODIC .or. bc_id == BC_OUTLET) then
              Fn = face_flux(f)
            else
              select case (bc_id)
              case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
                Fn = 0.0_rp
              case default
                Fn = face_flux(f)
              end select
            end if
          else
            Fn = face_flux(f)
          end if
        end if

        ! Reconstruction branch: rich periodic-partner resolution decides interior-style
        ! (face_interpolate_scalar) vs. boundary (phi_bnd) reconstruction, matching
        ! face_has_reconstruction_neighbor's own logic.
        has_reconstruction = (nb >= 1)
        if (.not. has_reconstruction .and. bc_id == BC_PERIODIC) then
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
          has_reconstruction = (nb_per >= 1)
        end if

        if (has_reconstruction) then
          ! face_interpolate_scalar's own (simpler) periodic resolution: periodic_neighbor
          ! only, no periodic_pair fallback -- replicated here for bit-identical behavior.
          if (nb >= 1) then
            nb_simple = nb
          else
            nb_simple = -1
            if (f_bnd >= 1 .and. f_bnd <= nbnd) nb_simple = pneigh_soa(f_bnd)
          end if

          if (nb_simple >= 1 .and. nb_simple <= nc + msh%n_ghost) then
            if (f > n_int) then
              gc = pgc_soa(f_bnd)
            else
              gc = gc_soa(f)
            end if
            phi_f_plain = gc * phi(o) + (1.0_rp - gc) * phi(nb_simple)
            if (f <= n_int) then
              epsilon = eps_soa(:, f)
            else
              d_on = ccent_soa(:, nb_simple) - ccent_soa(:, o)
              epsilon = fcentre_soa(:, f) - (ccent_soa(:, o) + (1.0_rp - gc) * d_on)
            end if
            grad_mean = gc * grad_phi(:, o) + (1.0_rp - gc) * grad_phi(:, nb_simple)
            phi_f = phi_f_plain + dot_product(grad_mean, epsilon)
            ! Fade the eps-correction out toward refinement interfaces; see refine_blend_weight.
            if (f <= n_int) then
              w_blend = refine_blend_weight(skew_soa(f))
              if (w_blend < 1.0_rp) phi_f = phi_f_plain + w_blend * (phi_f - phi_f_plain)
            end if
          else
            phi_f = phi(o)
          end if

          if (.not. soa_is_finite(phi_f)) then
            if (Fn >= 0.0_rp) then
              phi_f = phi(o)
            else if (nb >= 1 .and. nb <= size(phi)) then
              phi_f = phi(nb)
            else
              phi_f = phi(o)
            end if
          end if
        else
          if (has_bnd) then
            phi_f = cds_phi_bnd_buf(max(f_bnd, 1_ip))
          else
            phi_f = phi(o)
          end if
        end if

        if (.not. soa_is_finite(Fn)) Fn = 0.0_rp
        if (.not. soa_is_finite(phi_f)) phi_f = 0.0_rp

        if (o  >= 1 .and. o  <= nc) rhs(o)  = rhs(o)  - Fn * phi_f * inv_vol(o)
        if (nb >= 1 .and. nb <= nc) rhs(nb) = rhs(nb) + Fn * phi_f * inv_vol(nb)
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

    nullify(fo_soa, fn_soa, fcp, gc_soa, pgc_soa, cvol_soa, fpatch_soa, pbc_soa, pneigh_soa, &
            ppair_soa, fcentre_soa, ccent_soa, eps_soa, skew_soa)
  end subroutine advection_rhs_scalar_cds_soa

  !> QUICK/VANLEER_TVD/KOREN_TVD fast path via the flat SoA mesh mirror (including the
  !> cell_face_offset_soa CSR mirror for the upwind-cell stencil walk): matches
  !> quick_reconstruction/vanleer_tvd_reconstruction/koren_tvd_reconstruction's math exactly,
  !> minus the tvd_fallback_count diagnostic (dropped, not device-safe to increment).
  subroutine advection_rhs_scalar_recon_soa(phi, msh, scheme, rhs, phi_bnd, face_flux, inv_vol, grad_phi)
    real(rp),    intent(in)    :: phi(:)  ! (n_cells + n_ghost)
    type(Mesh),  intent(in), target :: msh
    integer(ip), intent(in)    :: scheme
    real(rp),    intent(inout) :: rhs(:)  ! (n_cells)
    real(rp),    intent(in), optional :: phi_bnd(:)  ! (n_bnd_faces)
    real(rp),    intent(in)    :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in)    :: inv_vol(:)  ! (n_cells) cached 1/volume
    real(rp),    intent(in)    :: grad_phi(:,:)  ! (3, n_cells), only used by QUICK

    integer(ip) :: f, o, nb, f_bnd, patch_id, bc_id, up, down, far_up
    real(rp)    :: Fn, phi_f, x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3), s_face, s_down
    integer(ip) :: col, idx, n_loops, lo, hi, nc, nf, n_int, nbnd, n_patches
    integer(ip) :: n_bnd_buf
    logical     :: use_color, has_bnd, ok
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:), cfo_soa(:), cell_faces_p(:)
    real(rp),    pointer :: cvol_soa(:), fcentre_soa(:,:), ccent_soa(:,:)
    real(rp),    pointer :: ppair_dcf_soa(:), fdcf_soa(:), fnorm_soa(:,:)

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'advection_rhs_scalar_recon_soa: msh%soa_ready is false; call mesh%build_soa_mirror first'

    rhs = 0.0_rp
    nc = msh%n_cells
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches

    fo_soa        => msh%face_owner_soa
    fn_soa        => msh%face_neighbour_soa
    cvol_soa      => msh%cell_volume_soa
    fpatch_soa    => msh%face_patch_id_soa
    pbc_soa       => msh%patch_bc_id_soa
    pneigh_soa    => msh%periodic_neighbor_soa
    ppair_soa     => msh%periodic_pair_soa
    ppair_dcf_soa => msh%periodic_pair_dcf_soa
    fdcf_soa      => msh%face_dcf_soa
    fnorm_soa     => msh%face_normal_soa
    fcentre_soa   => msh%face_centre_soa
    ccent_soa     => msh%cell_centroid_soa
    cfo_soa       => msh%cell_face_offset_soa
    cell_faces_p  => msh%cell_faces

    has_bnd = present(phi_bnd)
    n_bnd_buf = max(merge(size(phi_bnd), nbnd, has_bnd), 1_ip)
    call ensure_buf_1d(recon_s_phi_bnd_buf, n_bnd_buf)
    if (has_bnd) then
      recon_s_phi_bnd_buf(1:n_bnd_buf) = phi_bnd
    else
      recon_s_phi_bnd_buf(1:n_bnd_buf) = 0.0_rp
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
    !$acc data copyin(phi(:), face_flux(1:nf), recon_s_phi_bnd_buf(1:n_bnd_buf), inv_vol(1:nc), &
    !$acc      grad_phi(:,1:nc+msh%n_ghost), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), cvol_soa(1:nc), fpatch_soa(1:nf), &
    !$acc      pbc_soa(1:max(n_patches,1_ip)), pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      ppair_dcf_soa(1:max(nbnd,1_ip)), fdcf_soa(1:nf), fnorm_soa(:,1:nf), fcentre_soa(:,1:nf), &
    !$acc      ccent_soa(:,1:nc+msh%n_ghost), cfo_soa(1:nc+1), cell_faces_p(1:size(cell_faces_p)), &
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
      !$acc parallel loop if(use_color) private(f, o, nb, f_bnd, patch_id, bc_id, Fn, phi_f, &
      !$acc      up, down, far_up, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down, ok)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, f_bnd, patch_id, bc_id, Fn, phi_f, &
      !$omp&         up, down, far_up, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down, ok)
#endif
      do idx = lo, hi
        if (use_color) then
          f = fcp(idx)
        else
          f = idx
        end if
        o  = fo_soa(f)
        nb = fn_soa(f)
        f_bnd = 0
        patch_id = 0
        bc_id = -1

        ! Face-normal flux: identical to the CDS fast path (present(face_flux) case only).
        if (nb >= 1) then
          Fn = face_flux(f)
        else
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          if (patch_id >= 1 .and. patch_id <= n_patches) then
            bc_id = pbc_soa(patch_id)
            if (bc_id == BC_PERIODIC .or. bc_id == BC_OUTLET) then
              Fn = face_flux(f)
            else
              select case (bc_id)
              case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
                Fn = 0.0_rp
              case default
                Fn = face_flux(f)
              end select
            end if
          else
            Fn = face_flux(f)
          end if
        end if

        ok = soa_resolve_reconstruction_context_dev(f, Fn, n_int, nbnd, n_patches, nf, nc + msh%n_ghost, &
                  fo_soa, fn_soa, fcentre_soa, ccent_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, &
                  ppair_dcf_soa, fdcf_soa, fnorm_soa, up, down, x_up, x_down, x_face, up_shift, dir_hat, &
                  s_face, s_down)

        if (ok) then
          select case (scheme)
          case (SCHEME_QUICK)
            phi_f = soa_quick_reconstruction_ctx_dev(phi, up, down, x_up, x_down, up_shift, dir_hat, &
                      s_face, s_down, nc, cfo_soa, cell_faces_p, fo_soa, fn_soa, ccent_soa, grad_phi, &
                      nc + msh%n_ghost)
          case (SCHEME_VANLEER_TVD)
            far_up = soa_find_far_upwind_dev(up, down, x_up, up_shift, dir_hat, nc + msh%n_ghost, &
                       n_int, n_patches, nf, nc, cfo_soa, cell_faces_p, fo_soa, fn_soa, fpatch_soa, &
                       pbc_soa, pneigh_soa, ppair_soa, ccent_soa)
            phi_f = soa_vanleer_tvd_ctx_dev(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, &
                      far_up, ccent_soa)
          case (SCHEME_KOREN_TVD)
            far_up = soa_find_far_upwind_dev(up, down, x_up, up_shift, dir_hat, nc + msh%n_ghost, &
                       n_int, n_patches, nf, nc, cfo_soa, cell_faces_p, fo_soa, fn_soa, fpatch_soa, &
                       pbc_soa, pneigh_soa, ppair_soa, ccent_soa)
            phi_f = soa_koren_tvd_ctx_dev(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, &
                      far_up, ccent_soa)
          case default
            phi_f = 0.5_rp * (phi(up) + phi(down))
          end select

          if (.not. soa_is_finite(phi_f)) then
            if (Fn >= 0.0_rp) then
              phi_f = phi(o)
            else if (nb >= 1 .and. nb <= size(phi)) then
              phi_f = phi(nb)
            else
              phi_f = phi(o)
            end if
          end if
        else
          if (has_bnd) then
            phi_f = recon_s_phi_bnd_buf(max(f_bnd, 1_ip))
          else
            phi_f = phi(o)
          end if
        end if

        if (.not. soa_is_finite(Fn)) Fn = 0.0_rp
        if (.not. soa_is_finite(phi_f)) phi_f = 0.0_rp

        if (o  >= 1 .and. o  <= nc) rhs(o)  = rhs(o)  - Fn * phi_f * inv_vol(o)
        if (nb >= 1 .and. nb <= nc) rhs(nb) = rhs(nb) + Fn * phi_f * inv_vol(nb)
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

    nullify(fo_soa, fn_soa, fcp, cvol_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, &
            ppair_dcf_soa, fdcf_soa, fnorm_soa, fcentre_soa, ccent_soa, cfo_soa, cell_faces_p)
  end subroutine advection_rhs_scalar_recon_soa

  !> Assemble advective RHS for velocity (each component treated as scalar).
  !> Use scheme = KEP for the energy-conserving (skew-symmetric) form.
  subroutine advection_rhs_vector(U, msh, scheme, rhs, U_bnd, face_flux, inv_vol, grad_U)
    real(rp),    intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    type(Mesh),  intent(in)    :: msh
    integer(ip), intent(in)    :: scheme
    real(rp),    intent(inout) :: rhs(:,:)  ! (3, n_cells)
    real(rp),    intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),    intent(in), optional :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in), optional :: inv_vol(:)  ! (n_cells) cached 1/volume
    real(rp),    intent(in), optional :: grad_U(:,:,:)  ! (3, 3, n_cells) for O(h^2) correction

    integer(ip) :: f, o, nb, k, f_bnd, patch_id, nb_per
    real(rp)    :: Fn, Uf(3), Uf_plain(3), phi_f, phi_f_plain, ke_dot, inv_v, w_blend
    integer(ip) :: col, idx, n_loops, lo, hi
    logical     :: use_color, ke_count, split_form
    ! Reconstruction context for QUICK/VANLEER_TVD/KOREN_TVD: resolved once per face
    ! (component-independent) instead of once per velocity component inside the k loop.
    logical     :: ctx_ok
    integer(ip) :: ctx_up, ctx_down, ctx_far_up
    real(rp)    :: ctx_x_up(3), ctx_x_down(3), ctx_x_face(3), ctx_up_shift(3), ctx_dir_hat(3)
    real(rp)    :: ctx_s_face, ctx_s_down
    real(rp)    :: phi_f_quick(3)

    ! GPU/SoA fast path: KEP has a fixed owner/neighbour(+periodic) stencil, no reconstruction-function
    ! calls, so it maps onto !$acc parallel loop directly; other schemes stay on the AoS path below.
    if (scheme == SCHEME_KEP .and. KEP_PURE_MEAN .and. msh%soa_ready .and. .not. ke_budget_on &
        .and. present(U_bnd) .and. present(face_flux) .and. present(inv_vol)) then
      call advection_rhs_vector_kep_soa(U, msh, rhs, U_bnd, face_flux, inv_vol)
      return
    end if

    ! GPU/SoA fast path: QUICK/VANLEER_TVD/KOREN_TVD, mirrors the scalar recon SoA path.
    if ((scheme == SCHEME_QUICK .or. scheme == SCHEME_VANLEER_TVD .or. scheme == SCHEME_KOREN_TVD) &
        .and. present(U_bnd) .and. present(face_flux) .and. present(inv_vol) .and. present(grad_U) &
        .and. msh%soa_ready) then
      call advection_rhs_vector_recon_soa(U, msh, scheme, rhs, U_bnd, face_flux, inv_vol, grad_U)
      return
    end if

    rhs = 0.0_rp

    ! Skew-symmetric split: KE conservation then survives a non-zero discrete divergence.
    split_form = (scheme == SCHEME_KEP)
    if (split_form) then
      if (allocated(div_work)) then
        if (size(div_work) < msh%n_cells) deallocate(div_work)
      end if
      if (.not. allocated(div_work)) allocate(div_work(msh%n_cells))
      div_work(1:msh%n_cells) = 0.0_rp
    end if

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets))
    if (use_color) then
      n_loops = msh%n_colors
    else
      n_loops = 1
    end if

    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = msh%n_faces
      end if
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, k, f_bnd, patch_id, nb_per, Fn, Uf, Uf_plain, phi_f, &
      !$omp&         phi_f_plain, w_blend, ke_dot, ke_count, &
      !$omp&         ctx_ok, ctx_up, ctx_down, ctx_far_up, ctx_x_up, ctx_x_down, ctx_x_face, &
      !$omp&         ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down, phi_f_quick) &
      !$omp& reduction(+:ke_prod_smooth, ke_prod_refine, ke_prod_bnd, ke_bnd_exact)
      do idx = lo, hi
        if (use_color) then
          f = msh%face_color_perm(idx)
        else
          f = idx
        end if
      f_bnd = 0
      patch_id = 0
      nb_per = -1
      w_blend = 1.0_rp
      o   = msh%faces(f)%owner
      nb  = msh%faces(f)%neighbour
      if (nb >= 1) then
        if (scheme == SCHEME_KEP) then
          ! KEP: gc-weighted Uf preserves discrete KE on non-uniform grids (Verstappen 2003).
          Uf = msh%faces(f)%gc * U(:, o) + (1.0_rp - msh%faces(f)%gc) * U(:, nb)
        else
          if (present(grad_U)) then
            Uf = face_interpolate_vector(U, f, msh, grad_U=grad_U)
            ! Fade the eps-correction out toward refinement interfaces; see refine_blend_weight.
            w_blend = refine_blend_weight(msh%faces(f)%skew)
            if (w_blend < 1.0_rp) then
              Uf_plain = face_interpolate_vector(U, f, msh)
              Uf = Uf_plain + w_blend * (Uf - Uf_plain)
            end if
          else
            Uf = face_interpolate_vector(U, f, msh)
          end if
        end if
      else
        f_bnd = f - msh%n_int_faces
        patch_id = msh%faces(f)%patch_id
        if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
          if (msh%patches(patch_id)%bc_id == BC_PERIODIC) then
            nb_per = resolve_periodic_nb(msh, f_bnd, o)
            if (nb_per >= 1 .and. nb_per <= size(U, 2)) then
              Uf = msh%periodic_gc_soa(f_bnd) * U(:, o) + (1.0_rp - msh%periodic_gc_soa(f_bnd)) * U(:, nb_per)
            else if (present(U_bnd)) then
              Uf = U_bnd(:, f_bnd)
            else
              Uf = U(:, o)
            end if
          else if (present(U_bnd)) then
            Uf = U_bnd(:, f_bnd)
          else
            Uf = U(:, o)
          end if
        else if (present(U_bnd)) then
          Uf = U_bnd(:, f_bnd)
        else
          Uf = U(:, o)
        end if
      end if
      if (nb >= 1) then
        if (present(face_flux) .and. f >= 1 .and. f <= size(face_flux)) then
          Fn = face_flux(f)
        else
          Fn = dot_product(Uf, msh%faces(f)%sf)
        end if
      else
        patch_id = msh%faces(f)%patch_id
        if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
          ! Reuse projected flux on periodic/outlet faces (single conservative mass flux).
          if (present(face_flux) .and. f >= 1 .and. f <= size(face_flux) .and. &
              (msh%patches(patch_id)%bc_id == BC_PERIODIC .or. &
               msh%patches(patch_id)%bc_id == BC_OUTLET)) then
            Fn = face_flux(f)
          else
            select case (msh%patches(patch_id)%bc_id)
            case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
              ! Solid walls and slip boundaries: zero mass flux
              Fn = 0.0_rp
            case (BC_INLET, BC_OUTLET, BC_ABL_INLET, BC_FIXEDTEMP, BC_ZEROGRADIENT)
              ! Open boundaries: advective mass flux from boundary velocity
              Fn = dot_product(Uf, msh%faces(f)%sf)
            case default
              ! Unknown BC type: fallback to boundary velocity flux
              Fn = dot_product(Uf, msh%faces(f)%sf)
            end select
          end if
        else
          Fn = dot_product(Uf, msh%faces(f)%sf)
        end if
      end if

      if (split_form) then
        if (o  >= 1 .and. o  <= msh%n_cells) div_work(o)  = div_work(o)  + Fn
        if (nb >= 1 .and. nb <= msh%n_cells) div_work(nb) = div_work(nb) - Fn
      end if

      ! Count each interior face once: a processor face is visited by both ranks.
      ke_dot   = 0.0_rp
      ke_count = .false.
      if (ke_budget_on .and. nb >= 1 .and. o >= 1 .and. o <= msh%n_cells) then
        if (nb <= msh%n_cells) then
          ke_count = .true.
        else if (allocated(msh%global_id)) then
          ke_count = (msh%global_id(o) < msh%global_id(nb))
        end if
      end if

      ! Reconstruction geometry is component-independent (depends only on f, Fn, msh);
      ! resolve once here instead of once per component inside the k loop below.
      ctx_ok = .false.
      if (scheme == SCHEME_QUICK .or. scheme == SCHEME_VANLEER_TVD .or. scheme == SCHEME_KOREN_TVD) then
        if (nb >= 1 .or. face_has_reconstruction_neighbor(f, msh)) then
          ctx_ok = resolve_reconstruction_context(f, Fn, msh, ctx_up, ctx_down, ctx_x_up, ctx_x_down, &
                                                   ctx_x_face, ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down)
          if (ctx_ok .and. (ctx_up > size(U, 2) .or. ctx_down > size(U, 2))) ctx_ok = .false.
          if (ctx_ok .and. scheme /= SCHEME_QUICK) &
            ctx_far_up = find_far_upwind(ctx_up, ctx_down, ctx_x_up, ctx_up_shift, ctx_dir_hat, size(U, 2), msh)
        end if
      end if

      ! QUICK's LS-fit stencil geometry is also component-independent; compute the fit
      ! for all 3 components in one pass instead of once per component inside the k loop.
      if (ctx_ok .and. scheme == SCHEME_QUICK) then
        if (present(grad_U)) then
          call quick_reconstruction_vec_ctx(U, ctx_up, ctx_down, ctx_x_up, ctx_x_down, &
                    ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down, msh, grad_U, phi_f_quick)
        else
          call quick_reconstruction_vec_ctx(U, ctx_up, ctx_down, ctx_x_up, ctx_x_down, &
                    ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down, msh, phi_f=phi_f_quick)
        end if
      end if

      do k = 1, 3
        if (scheme == SCHEME_KEP) then
          ! KEP: arithmetic mean (KE conservation) plus the midpoint skewness term (2nd order).
          if (nb >= 1) then
            phi_f = 0.5_rp * (U(k, o) + U(k, nb))
            if (present(grad_U)) then
              if (nb <= size(grad_U, 3) .and. o <= size(grad_U, 3)) &
                phi_f = phi_f + kep_midpoint_correction(f, msh, grad_U(:, k, o), grad_U(:, k, nb))
            end if
          else if (nb_per >= 1 .and. nb_per <= size(U, 2)) then
            phi_f = 0.5_rp * (U(k, o) + U(k, nb_per))
          else
            phi_f = Uf(k)  ! wall / inlet / outlet: no off-cell partner
          end if
        else if (nb >= 1 .or. face_has_reconstruction_neighbor(f, msh)) then
          select case (scheme)
          case (SCHEME_QUICK)
            if (.not. ctx_ok) then
              phi_f = face_interpolate_scalar(U(k,:), f, msh)
            else
              phi_f = phi_f_quick(k)
            end if
          case (SCHEME_VANLEER_TVD)
            if (.not. ctx_ok) then
              phi_f = face_interpolate_scalar(U(k,:), f, msh)
            else
              phi_f = vanleer_tvd_reconstruction_ctx(U(k,:), ctx_up, ctx_down, ctx_x_up, &
                        ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down, ctx_far_up, msh)
            end if
          case (SCHEME_KOREN_TVD)
            if (.not. ctx_ok) then
              phi_f = face_interpolate_scalar(U(k,:), f, msh)
            else
              phi_f = koren_tvd_reconstruction_ctx(U(k,:), ctx_up, ctx_down, ctx_x_up, &
                        ctx_up_shift, ctx_dir_hat, ctx_s_face, ctx_s_down, ctx_far_up, msh)
            end if
          case default
            ! Linear (CDS) face interpolation with optional skewness/non-orthogonal
            ! eps-correction (Kim & Choi 2000, Eq. 26) when grad_U is available.
            if (present(grad_U)) then
              phi_f = face_interpolate_scalar(U(k,:), f, msh, grad_phi=grad_U(:, k, :))
              ! Fade the eps-correction out toward refinement interfaces; reuses w_blend from Uf above.
              if (nb >= 1) then
                if (w_blend < 1.0_rp) then
                  phi_f_plain = face_interpolate_scalar(U(k,:), f, msh)
                  phi_f = phi_f_plain + w_blend * (phi_f - phi_f_plain)
                end if
              end if
            else
              phi_f = face_interpolate_scalar(U(k,:), f, msh)
            end if
          end select
        else
          patch_id = msh%faces(f)%patch_id
          if (patch_id >= 1 .and. patch_id <= msh%n_patches) then
            if (msh%patches(patch_id)%bc_id == BC_PERIODIC) then
              ! F5: nb_per was already resolved above for this face; reuse.
              if (nb_per >= 1 .and. nb_per <= size(U, 2)) then
                phi_f = 0.5_rp * (U(k, o) + U(k, nb_per))
              else if (present(U_bnd)) then
                phi_f = U_bnd(k, f_bnd)
              else
                phi_f = U(k, o)
              end if
            else if (present(U_bnd)) then
              phi_f = U_bnd(k, f_bnd)
            else
              phi_f = U(k, o)
            end if
          else if (present(U_bnd)) then
            phi_f = U_bnd(k, f_bnd)
          else
            phi_f = U(k, o)
          end if
        end if

        if (ke_count) ke_dot = ke_dot + phi_f * (U(k, o) - U(k, nb))

        if (o >= 1 .and. o <= msh%n_cells) then
          if (present(inv_vol)) then
            rhs(k, o) = rhs(k, o) - Fn * phi_f * inv_vol(o)
          else
            rhs(k, o) = rhs(k, o) - Fn * phi_f / msh%cells(o)%volume
          end if
        end if

        if (nb >= 1 .and. nb <= msh%n_cells) then
          if (present(inv_vol)) then
            rhs(k, nb) = rhs(k, nb) + Fn * phi_f * inv_vol(nb)
          else
            rhs(k, nb) = rhs(k, nb) + Fn * phi_f / msh%cells(nb)%volume
          end if
        end if
      end do

      if (ke_count) then
        if (msh%faces(f)%is_refine) then
          ke_prod_refine = ke_prod_refine - Fn * ke_dot
        else
          ke_prod_smooth = ke_prod_smooth - Fn * ke_dot
        end if
      else if (ke_budget_on .and. nb < 1 .and. o >= 1 .and. o <= msh%n_cells) then
        ! Open boundaries carry real KE in and out; without them the budget is not a rate at all.
        ke_prod_bnd  = ke_prod_bnd  - Fn * dot_product(Uf, U(:, o))
        ke_bnd_exact = ke_bnd_exact + Fn * (0.5_rp * sum(U(:, o)**2) - dot_product(Uf, U(:, o)))
      end if
    end do
      !$omp end parallel do
    end do

    ! Advective half of the split form: rhs += 0.5*u_c*D_c/V_c cancels the divergence-error
    ! term the divergence form leaves in the KE budget, exactly and for any D_c.
    if (split_form) then
      do idx = 1, msh%n_cells
        if (present(inv_vol)) then
          inv_v = inv_vol(idx)
        else
          inv_v = 1.0_rp / msh%cells(idx)%volume
        end if
        do k = 1, 3
          rhs(k, idx) = rhs(k, idx) + 0.5_rp * U(k, idx) * div_work(idx) * inv_v
        end do
        if (ke_budget_on) &
          ke_prod_split = ke_prod_split + 0.5_rp * sum(U(:, idx)**2) * div_work(idx)
      end do
    end if
  end subroutine advection_rhs_vector

  !> KEP fast path via the flat SoA mesh mirror: fixed owner/neighbour(+periodic) stencil only,
  !> no reconstruction-function calls, so the colored loop maps directly onto !$acc parallel loop.
  !> Skips ke_budget diagnostics; caller only dispatches here when ke_budget_on is off.
  subroutine advection_rhs_vector_kep_soa(U, msh, rhs, U_bnd, face_flux, inv_vol)
    real(rp),    intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    type(Mesh),  intent(in), target :: msh
    real(rp),    intent(inout) :: rhs(:,:)  ! (3, n_cells)
    real(rp),    intent(in)    :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),    intent(in)    :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in)    :: inv_vol(:)  ! (n_cells) cached 1/volume

    integer(ip) :: f, o, nb, k, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd
    real(rp)    :: Fn, Uf(3), phi_f, inv_v, sf(3)
    integer(ip) :: col, idx, n_loops, lo, hi, nc, ng, nf, n_int, nbnd, n_patches, nu
    logical     :: use_color
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:)
    real(rp),    pointer :: sf_soa(:,:), gc_soa(:), cvol_soa(:), pgc_soa(:)

    rhs = 0.0_rp
    nc = msh%n_cells
    ng = msh%n_ghost
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches
    nu = size(U, 2)

    fo_soa     => msh%face_owner_soa
    fn_soa     => msh%face_neighbour_soa
    sf_soa     => msh%face_sf_soa
    gc_soa     => msh%face_gc_soa
    pgc_soa    => msh%periodic_gc_soa
    cvol_soa   => msh%cell_volume_soa
    fpatch_soa => msh%face_patch_id_soa
    pbc_soa    => msh%patch_bc_id_soa
    pneigh_soa => msh%periodic_neighbor_soa
    ppair_soa  => msh%periodic_pair_soa

    if (allocated(div_work)) then
      if (size(div_work) < nc) deallocate(div_work)
    end if
    if (.not. allocated(div_work)) allocate(div_work(nc))
    div_work(1:nc) = 0.0_rp

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
    !$acc data copyin(U(:,1:nu), U_bnd(:,:), face_flux(1:nf), inv_vol(1:nc), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), sf_soa(:,1:nf), gc_soa(1:nf), pgc_soa(1:max(nbnd,1_ip)), &
    !$acc      cvol_soa(1:nc), fpatch_soa(1:nf), pbc_soa(1:max(n_patches,1_ip)), &
    !$acc      pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      fcp(1:nf)) &
    !$acc      copy(rhs(:,1:nc), div_work(1:nc)) if(use_color)
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
      !$acc      pair_face, nb_per, pair_bnd, Fn, Uf, phi_f, sf)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, k, f_bnd, patch_id, bc_id, pair_face, nb_per, pair_bnd, Fn, Uf, phi_f, sf)
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
        f_bnd = 0
        patch_id = 0
        bc_id = -1
        nb_per = -1

        if (nb >= 1) then
          Uf = gc_soa(f) * U(:, o) + (1.0_rp - gc_soa(f)) * U(:, nb)
        else
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          if (patch_id >= 1 .and. patch_id <= n_patches) bc_id = pbc_soa(patch_id)

          if (bc_id == BC_PERIODIC) then
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
            if (nb_per >= 1 .and. nb_per <= nu) then
              Uf = pgc_soa(max(f_bnd,1_ip)) * U(:, o) + (1.0_rp - pgc_soa(max(f_bnd,1_ip))) * U(:, nb_per)
            else
              Uf = U_bnd(:, max(f_bnd,1_ip))
            end if
          else
            Uf = U_bnd(:, max(f_bnd,1_ip))
          end if
        end if

        if (nb >= 1) then
          Fn = face_flux(f)
        else
          if (patch_id >= 1 .and. patch_id <= n_patches) then
            if (bc_id == BC_PERIODIC .or. bc_id == BC_OUTLET) then
              Fn = face_flux(f)
            else
              select case (bc_id)
              case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
                Fn = 0.0_rp
              case default
                Fn = dot_product(Uf, sf)
              end select
            end if
          else
            Fn = dot_product(Uf, sf)
          end if
        end if

        if (o  >= 1 .and. o  <= nc) div_work(o)  = div_work(o)  + Fn
        if (nb >= 1 .and. nb <= nc) div_work(nb) = div_work(nb) - Fn

        do k = 1, 3
          if (nb >= 1) then
            phi_f = 0.5_rp * (U(k, o) + U(k, nb))
          else if (nb_per >= 1 .and. nb_per <= nu) then
            phi_f = 0.5_rp * (U(k, o) + U(k, nb_per))
          else
            phi_f = Uf(k)
          end if

          if (o  >= 1 .and. o  <= nc) rhs(k, o)  = rhs(k, o)  - Fn * phi_f * inv_vol(o)
          if (nb >= 1 .and. nb <= nc) rhs(k, nb) = rhs(k, nb) + Fn * phi_f * inv_vol(nb)
        end do
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

    ! Advective half of the split form (see advection_rhs_vector for the rationale).
    do idx = 1, nc
      inv_v = inv_vol(idx)
      do k = 1, 3
        rhs(k, idx) = rhs(k, idx) + 0.5_rp * U(k, idx) * div_work(idx) * inv_v
      end do
    end do

    nullify(fo_soa, fn_soa, fcp, sf_soa, gc_soa, pgc_soa, cvol_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa)
  end subroutine advection_rhs_vector_kep_soa

  !> GPU/SoA fast path for velocity QUICK/VANLEER_TVD/KOREN_TVD: reuses the scalar SoA
  !> reconstruction device routines per component (see advection_rhs_scalar_recon_soa).
  subroutine advection_rhs_vector_recon_soa(U, msh, scheme, rhs, U_bnd, face_flux, inv_vol, grad_U)
    real(rp),    intent(in)    :: U(:,:)  ! (3, n_cells + n_ghost)
    type(Mesh),  intent(in), target :: msh
    integer(ip), intent(in)    :: scheme
    real(rp),    intent(inout) :: rhs(:,:)  ! (3, n_cells)
    real(rp),    intent(in)    :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),    intent(in)    :: face_flux(:)  ! (n_faces), projected flux
    real(rp),    intent(in)    :: inv_vol(:)  ! (n_cells) cached 1/volume
    real(rp),    intent(in)    :: grad_U(:,:,:)  ! (3, 3, n_cells), grad_U(:,k,c), only used by QUICK

    integer(ip) :: f, o, nb, k, f_bnd, patch_id, bc_id, up, down, far_up
    real(rp)    :: Fn, phi_f, x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3), s_face, s_down
    integer(ip) :: col, idx, n_loops, lo, hi, nc, nf, n_int, nbnd, n_patches, nu
    logical     :: use_color, ok
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), pointer :: pneigh_soa(:), ppair_soa(:), cfo_soa(:), cell_faces_p(:)
    real(rp),    pointer :: cvol_soa(:), fcentre_soa(:,:), ccent_soa(:,:)
    real(rp),    pointer :: ppair_dcf_soa(:), fdcf_soa(:), fnorm_soa(:,:)
    integer(ip) :: kk

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'advection_rhs_vector_recon_soa: msh%soa_ready is false; call mesh%build_soa_mirror first'

    rhs = 0.0_rp
    nc = msh%n_cells
    nf = msh%n_faces
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    n_patches = msh%n_patches
    nu = size(U, 2)

    fo_soa        => msh%face_owner_soa
    fn_soa        => msh%face_neighbour_soa
    cvol_soa      => msh%cell_volume_soa
    fpatch_soa    => msh%face_patch_id_soa
    pbc_soa       => msh%patch_bc_id_soa
    pneigh_soa    => msh%periodic_neighbor_soa
    ppair_soa     => msh%periodic_pair_soa
    ppair_dcf_soa => msh%periodic_pair_dcf_soa
    fdcf_soa      => msh%face_dcf_soa
    fnorm_soa     => msh%face_normal_soa
    fcentre_soa   => msh%face_centre_soa
    ccent_soa     => msh%cell_centroid_soa
    cfo_soa       => msh%cell_face_offset_soa
    cell_faces_p  => msh%cell_faces

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets) .and. msh%soa_ready)
    n_loops = merge(msh%n_colors, 1_ip, use_color)

    if (use_color) then
      fcp => msh%face_color_perm
    else
      call ensure_fcp_identity(fcp_identity_buf, nf)
      fcp => fcp_identity_buf(1:nf)
    end if

    ! U(k,:) and grad_U(:,k,:) are non-contiguous (stride-3) sections; passed straight into an
    ! !$acc routine seq dummy they force a per-thread copy scaling with mesh size, which is what
    ! was overflowing device local memory above ~11k cells. Transpose once here instead so the
    ! per-component slices handed to the device routines below are contiguous.
    ! This host fill is the correctness fallback for the use_color=.false. path below (where the
    ! compute loop's "if(use_color)" makes it run on host); when use_color=.true. it is redundant
    ! with the device-side fill further down but cheap enough (O(nc) scalar copies) not to bother
    ! skipping.
    call ensure_buf_n3(recon_v_Ut_buf, nu)
    ! Explicit loop instead of transpose(U): the intrinsic would allocate a full (nu,3)
    ! temporary before the assignment, on every call.
    do kk = 1, 3
      recon_v_Ut_buf(1:nu, kk) = U(kk, 1:nu)
    end do
    ! gradt is only ever read inside the SCHEME_QUICK branch below, so for other schemes the
    ! buffer just needs a valid (already-allocated) shape, not fresh values.
    call ensure_buf_3n3(recon_v_gradt_buf, nc)
    if (scheme == SCHEME_QUICK) then
      ! grad_U's 3rd extent is n_cells+n_ghost; only owned cells are read back out of this
      ! buffer downstream, so slice to 1:nc (a bare `:` here is a non-conforming assignment
      ! that overruns recon_v_gradt_buf by 3*n_ghost reals).
      do kk = 1, 3
        recon_v_gradt_buf(:, 1:nc, kk) = grad_U(:, kk, 1:nc)
      end do
    end if

#ifdef _OPENACC
    ! recon_v_Ut_buf/recon_v_gradt_buf are `create`, not `copyin`: filled by the device kernel
    ! below instead of uploading the host repack above. U/grad_U stay `copyin`, which is a no-op
    ! refcount bump (no transfer) when this is called from advance_rk3w's outer resident region
    ! (the common case) and a real upload only when this routine is exercised standalone (e.g.
    ! unit tests) without an enclosing device data region.
    !$acc data copyin(U(:,1:nu), grad_U(:,:,1:nc), &
    !$acc      face_flux(1:nf), U_bnd(:,:), inv_vol(1:nc), &
    !$acc      fo_soa(1:nf), fn_soa(1:nf), cvol_soa(1:nc), fpatch_soa(1:nf), &
    !$acc      pbc_soa(1:max(n_patches,1_ip)), pneigh_soa(1:max(nbnd,1_ip)), ppair_soa(1:max(nbnd,1_ip)), &
    !$acc      ppair_dcf_soa(1:max(nbnd,1_ip)), fdcf_soa(1:nf), fnorm_soa(:,1:nf), fcentre_soa(:,1:nf), &
    !$acc      ccent_soa(:,1:nc+msh%n_ghost), cfo_soa(1:nc+1), cell_faces_p(1:size(cell_faces_p)), &
    !$acc      fcp(1:nf)) &
    !$acc      create(recon_v_Ut_buf(1:nu,1:3), recon_v_gradt_buf(:,1:nc,:)) &
    !$acc      copy(rhs(:,1:nc)) if(use_color)

    ! Device-side fill: only executes on-device when use_color (matches the compute loop's own
    ! if(use_color) below); when use_color is false this runs on host and is a harmless redundant
    ! repeat of the host loop above, since in that branch the compute loop reads the host buffers.
    !$acc parallel loop if(use_color)
    do idx = 1, nu
      recon_v_Ut_buf(idx, 1) = U(1, idx)
      recon_v_Ut_buf(idx, 2) = U(2, idx)
      recon_v_Ut_buf(idx, 3) = U(3, idx)
    end do
    !$acc end parallel loop
    if (scheme == SCHEME_QUICK) then
      !$acc parallel loop if(use_color)
      do idx = 1, nc
        recon_v_gradt_buf(1, idx, 1) = grad_U(1, 1, idx)
        recon_v_gradt_buf(2, idx, 1) = grad_U(2, 1, idx)
        recon_v_gradt_buf(3, idx, 1) = grad_U(3, 1, idx)
        recon_v_gradt_buf(1, idx, 2) = grad_U(1, 2, idx)
        recon_v_gradt_buf(2, idx, 2) = grad_U(2, 2, idx)
        recon_v_gradt_buf(3, idx, 2) = grad_U(3, 2, idx)
        recon_v_gradt_buf(1, idx, 3) = grad_U(1, 3, idx)
        recon_v_gradt_buf(2, idx, 3) = grad_U(2, 3, idx)
        recon_v_gradt_buf(3, idx, 3) = grad_U(3, 3, idx)
      end do
      !$acc end parallel loop
    end if
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
      !$acc parallel loop if(use_color) private(f, o, nb, k, f_bnd, patch_id, bc_id, Fn, phi_f, &
      !$acc      up, down, far_up, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down, ok)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, k, f_bnd, patch_id, bc_id, Fn, phi_f, &
      !$omp&         up, down, far_up, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down, ok)
#endif
      do idx = lo, hi
        if (use_color) then
          f = fcp(idx)
        else
          f = idx
        end if
        o  = fo_soa(f)
        nb = fn_soa(f)
        f_bnd = 0
        patch_id = 0
        bc_id = -1

        ! Face-normal flux: identical to the scalar recon SoA fast path.
        if (nb >= 1) then
          Fn = face_flux(f)
        else
          f_bnd = f - n_int
          patch_id = fpatch_soa(f)
          if (patch_id >= 1 .and. patch_id <= n_patches) then
            bc_id = pbc_soa(patch_id)
            if (bc_id == BC_PERIODIC .or. bc_id == BC_OUTLET) then
              Fn = face_flux(f)
            else
              select case (bc_id)
              case (BC_WALL_NOSLIP, BC_WALL_MODEL, BC_WALL_VELOCITY, BC_SLIP, BC_SYMMETRY)
                Fn = 0.0_rp
              case default
                Fn = face_flux(f)
              end select
            end if
          else
            Fn = face_flux(f)
          end if
        end if

        ok = soa_resolve_reconstruction_context_dev(f, Fn, n_int, nbnd, n_patches, nf, nu, &
                  fo_soa, fn_soa, fcentre_soa, ccent_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, &
                  ppair_dcf_soa, fdcf_soa, fnorm_soa, up, down, x_up, x_down, x_face, up_shift, dir_hat, &
                  s_face, s_down)

        if (ok .and. (up > nu .or. down > nu)) ok = .false.

        if (ok .and. scheme /= SCHEME_QUICK) then
          far_up = soa_find_far_upwind_dev(up, down, x_up, up_shift, dir_hat, nu, &
                     n_int, n_patches, nf, nc, cfo_soa, cell_faces_p, fo_soa, fn_soa, fpatch_soa, &
                     pbc_soa, pneigh_soa, ppair_soa, ccent_soa)
        end if

        if (.not. soa_is_finite(Fn)) Fn = 0.0_rp

        do k = 1, 3
          if (ok) then
            select case (scheme)
            case (SCHEME_QUICK)
              phi_f = soa_quick_reconstruction_ctx_dev(recon_v_Ut_buf(1:nu,k), up, down, x_up, x_down, up_shift, dir_hat, &
                        s_face, s_down, nc, cfo_soa, cell_faces_p, fo_soa, fn_soa, ccent_soa, recon_v_gradt_buf(:,1:nc,k), nu)
            case (SCHEME_VANLEER_TVD)
              phi_f = soa_vanleer_tvd_ctx_dev(recon_v_Ut_buf(1:nu,k), up, down, x_up, up_shift, dir_hat, s_face, s_down, &
                        far_up, ccent_soa)
            case (SCHEME_KOREN_TVD)
              phi_f = soa_koren_tvd_ctx_dev(recon_v_Ut_buf(1:nu,k), up, down, x_up, up_shift, dir_hat, s_face, s_down, &
                        far_up, ccent_soa)
            case default
              phi_f = 0.5_rp * (U(k,up) + U(k,down))
            end select

            if (.not. soa_is_finite(phi_f)) then
              if (Fn >= 0.0_rp) then
                phi_f = U(k, o)
              else if (nb >= 1 .and. nb <= nu) then
                phi_f = U(k, nb)
              else
                phi_f = U(k, o)
              end if
            end if
          else
            phi_f = U_bnd(k, max(f_bnd, 1_ip))
          end if

          if (.not. soa_is_finite(phi_f)) phi_f = 0.0_rp

          if (o  >= 1 .and. o  <= nc) rhs(k, o)  = rhs(k, o)  - Fn * phi_f * inv_vol(o)
          if (nb >= 1 .and. nb <= nc) rhs(k, nb) = rhs(k, nb) + Fn * phi_f * inv_vol(nb)
        end do
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

    nullify(fo_soa, fn_soa, fcp, cvol_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, &
            ppair_dcf_soa, fdcf_soa, fnorm_soa, fcentre_soa, ccent_soa, cfo_soa, cell_faces_p)
  end subroutine advection_rhs_vector_recon_soa

  subroutine advection_set_ke_budget(on)
    logical, intent(in) :: on
    ke_budget_on = on
  end subroutine advection_set_ke_budget

  subroutine advection_reset_ke_budget()
    ke_prod_smooth = 0.0_rp
    ke_prod_refine = 0.0_rp
    ke_prod_split  = 0.0_rp
    ke_prod_bnd    = 0.0_rp
    ke_bnd_exact   = 0.0_rp
  end subroutine advection_reset_ke_budget

  subroutine advection_get_ke_budget(smooth, refine, split, bnd, bnd_exact)
    real(rp), intent(out) :: smooth, refine
    real(rp), intent(out), optional :: split, bnd, bnd_exact
    smooth = ke_prod_smooth
    refine = ke_prod_refine
    if (present(split))     split     = ke_prod_split
    if (present(bnd))       bnd       = ke_prod_bnd
    if (present(bnd_exact)) bnd_exact = ke_bnd_exact
  end subroutine advection_get_ke_budget

  logical function face_has_reconstruction_neighbor(f, msh) result(has_neighbor)
    integer(ip), intent(in) :: f
    type(Mesh),  intent(in) :: msh

    integer(ip) :: f_bnd, patch_id, pair_cell

    has_neighbor = .false.
    if (msh%faces(f)%neighbour >= 1) then
      has_neighbor = .true.
      return
    end if

    patch_id = msh%faces(f)%patch_id
    if (patch_id < 1 .or. patch_id > msh%n_patches) return
    if (msh%patches(patch_id)%bc_id /= BC_PERIODIC) return

    f_bnd = f - msh%n_int_faces
    pair_cell = resolve_periodic_nb(msh, f_bnd, msh%faces(f)%owner)
    has_neighbor = pair_cell >= 1
  end function face_has_reconstruction_neighbor

  !> Sample difference with the off-axis part removed, so a 1-D fit along dir_hat stays linear-exact.
  pure real(rp) function axial_delta(phi_far, phi_up, grad_up, dvec, s, dir_hat) result(d)
    !$acc routine seq
    real(rp), intent(in) :: phi_far, phi_up, grad_up(3), dvec(3), s, dir_hat(3)
    d = phi_far - phi_up - dot_product(grad_up, dvec - s * dir_hat)
  end function axial_delta

  !> Midpoint skewness term: the arithmetic mean is phi at 0.5*(x_o+x_nb), not at the face centre.
  pure real(rp) function kep_midpoint_correction(f, msh, grad_o, grad_nb) result(corr)
    integer(ip), intent(in) :: f
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in) :: grad_o(3), grad_nb(3)

    integer(ip) :: o, nb
    real(rp)    :: eps_mid(3)

    ! Unsigned KE source on every skewed face; measured at +5.4e2 against -1.1e1 from tagged faces alone.
    corr = 0.0_rp
    if (KEP_PURE_MEAN) return

    o  = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour
    eps_mid = msh%faces(f)%centre - 0.5_rp * (msh%cells(o)%centroid + msh%cells(nb)%centroid)
    corr = dot_product(0.5_rp * (grad_o + grad_nb), eps_mid)
  end function kep_midpoint_correction

  real(rp) function reconstruct_face_value(phi, f, Fn, msh, scheme, grad_phi) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f, scheme
    real(rp),    intent(in) :: Fn
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in), optional :: grad_phi(:,:)

    select case (scheme)
    case (SCHEME_QUICK)
      phi_f = quick_reconstruction(phi, f, Fn, msh, grad_phi)
    case (SCHEME_VANLEER_TVD)
      phi_f = vanleer_tvd_reconstruction(phi, f, Fn, msh)
    case (SCHEME_KOREN_TVD)
      phi_f = koren_tvd_reconstruction(phi, f, Fn, msh)
    case (SCHEME_KEP)
      phi_f = kep_reconstruction(phi, f, Fn, msh, grad_phi)
    case (SCHEME_CDS)
      phi_f = linear_reconstruction(phi, f, Fn, msh)
    case default
      phi_f = face_interpolate_scalar(phi, f, msh)
    end select
  end function reconstruct_face_value

  !> KEP scalar: arithmetic mean phi_f = 0.5*(phi_up + phi_dn) for scalar variance conservation.
  real(rp) function kep_reconstruction(phi, f, Fn, msh, grad_phi) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    real(rp),    intent(in) :: Fn
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in), optional :: grad_phi(:,:)

    integer(ip) :: up, down, o, nb
    real(rp)    :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp)    :: s_face, s_down
    logical     :: ok

    ok = resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down)
    if (.not. ok .or. up > size(phi) .or. down > size(phi)) then
      phi_f = face_interpolate_scalar(phi, f, msh)
      return
    end if
    phi_f = 0.5_rp * (phi(up) + phi(down))

    ! Interior faces only: on a periodic seam the gradients straddle a discontinuity
    ! and averaging them destroys the energy conservation the arithmetic mean buys.
    if (.not. present(grad_phi)) return
    o  = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour
    if (nb < 1 .or. nb > size(grad_phi, 2) .or. o > size(grad_phi, 2)) return
    phi_f = phi_f + kep_midpoint_correction(f, msh, grad_phi(:, o), grad_phi(:, nb))
  end function kep_reconstruction

  real(rp) function linear_reconstruction(phi, f, Fn, msh) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    real(rp),    intent(in) :: Fn
    type(Mesh),  intent(in) :: msh

    integer(ip) :: up, down
    real(rp)    :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp)    :: s_face, s_down
    logical     :: ok

    ok = resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down)
    if (.not. ok .or. up > size(phi) .or. down > size(phi)) then
      phi_f = face_interpolate_scalar(phi, f, msh)
      return
    end if

    phi_f = phi(up) + (phi(down) - phi(up)) * s_face / max(s_down, SMALL)
  end function linear_reconstruction

  subroutine resolve_periodic_partner(f, msh, pair_cell, pair_face, shift)
    integer(ip), intent(in)  :: f
    type(Mesh),  intent(in)  :: msh
    integer(ip), intent(out) :: pair_cell, pair_face
    real(rp),    intent(out) :: shift(3)

    integer(ip) :: f_bnd, patch_id, pair_bnd, owner
    real(rp)    :: d_per

    pair_cell = 0
    pair_face = 0
    shift = 0.0_rp

    if (msh%faces(f)%neighbour >= 1) return
    patch_id = msh%faces(f)%patch_id
    if (patch_id < 1 .or. patch_id > msh%n_patches) return
    if (msh%patches(patch_id)%bc_id /= BC_PERIODIC) return

    f_bnd = f - msh%n_int_faces
    if (f_bnd < 1 .or. f_bnd > msh%n_bnd_faces) return

    owner = msh%faces(f)%owner
    pair_cell = resolve_periodic_nb(msh, f_bnd, msh%faces(f)%owner)
    if (allocated(msh%periodic_pair)) then
      if (f_bnd > size(msh%periodic_pair)) return
      pair_bnd = msh%periodic_pair(f_bnd)
      if (pair_bnd > 0) then
        pair_face = msh%n_int_faces + pair_bnd
        if (pair_face >= 1 .and. pair_face <= size(msh%faces)) then
          if (pair_cell <= 0) pair_cell = msh%faces(pair_face)%owner
        end if
      end if
    end if

    ! Build periodic image shift for face reconstruction geometry (robust for same/cross-rank pairs).
    if (pair_cell > 0 .and. pair_cell <= size(msh%cells)) then
      if (pair_face >= 1 .and. pair_face <= size(msh%faces)) then
        d_per = msh%faces(f)%d_cf + msh%faces(pair_face)%d_cf
      else
        ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()/
        ! the array element must never be evaluated in the same expression as allocated().
        d_per = 2.0_rp * msh%faces(f)%d_cf
        if (allocated(msh%periodic_pair_dcf)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair_dcf)) then
            if (msh%periodic_pair_dcf(f_bnd) > SMALL) &
              d_per = msh%faces(f)%d_cf + msh%periodic_pair_dcf(f_bnd)
          end if
        end if
      end if
      shift = (msh%cells(owner)%centroid + d_per * msh%faces(f)%normal) - msh%cells(pair_cell)%centroid
    end if
  end subroutine resolve_periodic_partner

  integer(ip) function resolve_periodic_nb(msh, f_bnd, owner) result(nb)
    type(Mesh),  intent(in) :: msh
    integer(ip), intent(in) :: f_bnd, owner

    integer(ip) :: pair_bnd, pair_face

    nb = -1
    if (f_bnd < 1 .or. f_bnd > msh%n_bnd_faces) return

    if (allocated(msh%periodic_neighbor)) then
      if (f_bnd <= size(msh%periodic_neighbor)) then
        nb = msh%periodic_neighbor(f_bnd)
      end if
    end if

    if (nb <= 0 .and. allocated(msh%periodic_pair)) then
      if (f_bnd <= size(msh%periodic_pair)) then
        pair_bnd = msh%periodic_pair(f_bnd)
        if (pair_bnd > 0) then
          pair_face = msh%n_int_faces + pair_bnd
          if (pair_face >= 1 .and. pair_face <= size(msh%faces)) nb = msh%faces(pair_face)%owner
        end if
      end if
    end if

    ! On distributed meshes, periodic partner stored as ghost in faces%neighbour: authoritative fallback.
    if (nb <= 0) then
      pair_face = msh%n_int_faces + f_bnd
      if (pair_face >= 1 .and. pair_face <= size(msh%faces)) nb = msh%faces(pair_face)%neighbour
    end if

    if (nb == owner) nb = -1
  end function resolve_periodic_nb

  logical function resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down) result(ok)
    integer(ip), intent(in)  :: f
    real(rp),    intent(in)  :: Fn
    type(Mesh),  intent(in)  :: msh
    integer(ip), intent(out) :: up, down
    real(rp),    intent(out) :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp),    intent(out) :: s_face, s_down

    integer(ip) :: owner, nb, pair_cell, pair_face
    real(rp)    :: pair_shift(3), dvec(3), dnorm

    ok = .false.
    up = 0
    down = 0
    x_up = 0.0_rp
    x_down = 0.0_rp
    x_face = msh%faces(f)%centre
    up_shift = 0.0_rp
    dir_hat = 0.0_rp
    s_face = 0.0_rp
    s_down = 0.0_rp

    owner = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour

    if (nb >= 1) then
      if (Fn >= 0.0_rp) then
        up = owner
        down = nb
      else
        up = nb
        down = owner
      end if
      x_up = msh%cells(up)%centroid
      x_down = msh%cells(down)%centroid
    else
      call resolve_periodic_partner(f, msh, pair_cell, pair_face, pair_shift)
      if (pair_cell <= 0) return
      if (Fn >= 0.0_rp) then
        up = owner
        down = pair_cell
        x_up = msh%cells(owner)%centroid
        x_down = msh%cells(pair_cell)%centroid + pair_shift
      else
        up = pair_cell
        down = owner
        up_shift = pair_shift
        x_up = msh%cells(pair_cell)%centroid + up_shift
        x_down = msh%cells(owner)%centroid
      end if
    end if

    if (up < 1 .or. down < 1) return

    dvec = x_face - x_up
    dnorm = sqrt(sum(dvec**2))
    if (dnorm <= SMALL) then
      dvec = x_down - x_up
      dnorm = sqrt(sum(dvec**2))
    end if
    if (dnorm <= SMALL) return
    dir_hat = dvec / dnorm
    s_face = dot_product(x_face - x_up, dir_hat)
    s_down = dot_product(x_down - x_up, dir_hat)

    if (s_face <= SMALL .or. s_down <= s_face) then
      dvec = x_down - x_up
      dnorm = sqrt(sum(dvec**2))
      if (dnorm <= SMALL) return
      dir_hat = dvec / dnorm
      s_face = dot_product(x_face - x_up, dir_hat)
      s_down = dot_product(x_down - x_up, dir_hat)
    end if
    if (s_down <= SMALL) return

    s_face = min(max(s_face, 0.0_rp), s_down)
    ok = .true.
  end function resolve_reconstruction_context

  ! QUICK reconstruction

  !> Geometry-aware QUICK: quadratic LS fit on an upwind stencil; public API for tests, hot loop uses quick_reconstruction_ctx instead.
  real(rp) function quick_reconstruction(phi, f, Fn, msh, grad_phi) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    real(rp),    intent(in) :: Fn  ! positive = flow from owner to neighbour
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in), optional :: grad_phi(:,:)

    integer(ip) :: up, down
    real(rp)    :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp)    :: s_face, s_down
    logical     :: ok

    ok = resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down)
    if (.not. ok .or. up > size(phi) .or. down > size(phi)) then
      phi_f = face_interpolate_scalar(phi, f, msh)
      return
    end if

    if (present(grad_phi)) then
      phi_f = quick_reconstruction_ctx(phi, up, down, x_up, x_down, up_shift, dir_hat, s_face, s_down, msh, grad_phi)
    else
      phi_f = quick_reconstruction_ctx(phi, up, down, x_up, x_down, up_shift, dir_hat, s_face, s_down, msh)
    end if
  end function quick_reconstruction

  !> Core of quick_reconstruction with a pre-resolved context; caller must have already confirmed up/down are in range for phi.
  real(rp) function quick_reconstruction_ctx(phi, up, down, x_up, x_down, up_shift, dir_hat, &
                                              s_face, s_down, msh, grad_phi) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down
    real(rp),    intent(in) :: x_up(3), x_down(3), up_shift(3), dir_hat(3)
    real(rp),    intent(in) :: s_face, s_down
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in), optional :: grad_phi(:,:)

    integer(ip) :: fi, fj, nb, nfit
    real(rp)    :: phi_U, phi_D, s, dist2, align, weight
    real(rp)    :: a11, a12, a22, b1, b2, det, a1, a2, lin_num, lin_den
    real(rp)    :: x_nb(3), dvec(3), delta_phi, grad_up(3)

    phi_U = phi(up)
    phi_D = phi(down)

    ! Upwind-cell gradient removes each sample's off-axis part; without it the fit is
    ! only first order wherever the stencil is skewed (refinement faces, snapped cells).
    grad_up = 0.0_rp
    if (present(grad_phi)) then
      if (up <= size(grad_phi, 2)) grad_up = grad_phi(:, up)
    end if

    a11 = 0.0_rp; a12 = 0.0_rp; a22 = 0.0_rp
    b1 = 0.0_rp; b2 = 0.0_rp
    lin_num = 0.0_rp; lin_den = 0.0_rp
    nfit = 0

    dvec = x_down - x_up
    dist2 = max(sum(dvec**2), SMALL)
    weight = 1.0_rp / dist2
    delta_phi = axial_delta(phi_D, phi_U, grad_up, dvec, s_down, dir_hat)
    a11 = a11 + weight * s_down * s_down
    a12 = a12 + weight * s_down * s_down * s_down
    a22 = a22 + weight * s_down * s_down * s_down * s_down
    b1 = b1 + weight * s_down * delta_phi
    b2 = b2 + weight * s_down * s_down * delta_phi
    lin_num = lin_num + weight * s_down * delta_phi
    lin_den = lin_den + weight * s_down * s_down
    nfit = nfit + 1

    do fi = msh%cells(up)%face_offset, msh%cells(up)%face_offset + msh%cells(up)%n_faces - 1
      fj = msh%cell_faces(fi)
      nb = msh%faces(fj)%owner
      if (nb == up) nb = msh%faces(fj)%neighbour
      if (nb < 1 .or. nb == down .or. nb > size(phi)) cycle

      x_nb = msh%cells(nb)%centroid + up_shift
      dvec = x_nb - x_up
      s = dot_product(dvec, dir_hat)
      if (s >= -1.0e-12_rp) cycle

      dist2 = max(sum(dvec**2), SMALL)
      if (s*s < 0.25_rp*dist2) cycle  ! reject before sqrt; s<0 already verified above
      align = min(1.0_rp, -s / sqrt(dist2))
      weight = (0.25_rp + align*align) / dist2
      delta_phi = axial_delta(phi(nb), phi_U, grad_up, dvec, s, dir_hat)

      a11 = a11 + weight * s * s
      a12 = a12 + weight * s * s * s
      a22 = a22 + weight * s * s * s * s
      b1 = b1 + weight * s * delta_phi
      b2 = b2 + weight * s * s * delta_phi
      lin_num = lin_num + weight * s * delta_phi
      lin_den = lin_den + weight * s * s
      nfit = nfit + 1
    end do

    det = a11*a22 - a12*a12
    if (nfit >= 2 .and. abs(det) > SMALL) then
      a1 = ( b1*a22 - b2*a12) / det
      a2 = (-b1*a12 + b2*a11) / det
    else if (lin_den > SMALL) then
      a1 = lin_num / lin_den
      a2 = 0.0_rp
    else
      a1 = (phi_D - phi_U) / max(s_down, SMALL)
      a2 = 0.0_rp
    end if

    phi_f = phi_U + a1*s_face + a2*s_face*s_face
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function quick_reconstruction_ctx

  !> Vector counterpart of quick_reconstruction_ctx: fits all 3 components in one pass instead of 3x-repeating the stencil walk.
  subroutine quick_reconstruction_vec_ctx(U, up, down, x_up, x_down, up_shift, dir_hat, &
                                           s_face, s_down, msh, grad_U, phi_f)
    real(rp),    intent(in)  :: U(:,:)  ! (3, n)
    integer(ip), intent(in)  :: up, down
    real(rp),    intent(in)  :: x_up(3), x_down(3), up_shift(3), dir_hat(3)
    real(rp),    intent(in)  :: s_face, s_down
    type(Mesh),  intent(in)  :: msh
    real(rp),    intent(in), optional :: grad_U(:,:,:)  ! (3, 3, n): grad_U(:,k,c)
    real(rp),    intent(out) :: phi_f(3)

    integer(ip) :: fi, fj, nb, nfit, k
    real(rp)    :: phi_U(3), phi_D(3), s, dist2, align, weight
    real(rp)    :: a11, a12, a22, b1(3), b2(3), det, a1(3), a2(3), lin_num(3), lin_den
    real(rp)    :: x_nb(3), dvec(3), delta_phi(3), grad_up(3,3)

    phi_U = U(:, up)
    phi_D = U(:, down)

    grad_up = 0.0_rp
    if (present(grad_U)) then
      if (up <= size(grad_U, 3)) grad_up = grad_U(:, :, up)
    end if

    a11 = 0.0_rp; a12 = 0.0_rp; a22 = 0.0_rp
    b1 = 0.0_rp; b2 = 0.0_rp
    lin_num = 0.0_rp; lin_den = 0.0_rp
    nfit = 0

    dvec = x_down - x_up
    dist2 = max(sum(dvec**2), SMALL)
    weight = 1.0_rp / dist2
    do k = 1, 3
      delta_phi(k) = axial_delta(phi_D(k), phi_U(k), grad_up(:, k), dvec, s_down, dir_hat)
    end do
    a11 = a11 + weight * s_down * s_down
    a12 = a12 + weight * s_down * s_down * s_down
    a22 = a22 + weight * s_down * s_down * s_down * s_down
    b1 = b1 + weight * s_down * delta_phi
    b2 = b2 + weight * s_down * s_down * delta_phi
    lin_num = lin_num + weight * s_down * delta_phi
    lin_den = lin_den + weight * s_down * s_down
    nfit = nfit + 1

    do fi = msh%cells(up)%face_offset, msh%cells(up)%face_offset + msh%cells(up)%n_faces - 1
      fj = msh%cell_faces(fi)
      nb = msh%faces(fj)%owner
      if (nb == up) nb = msh%faces(fj)%neighbour
      if (nb < 1 .or. nb == down .or. nb > size(U, 2)) cycle

      x_nb = msh%cells(nb)%centroid + up_shift
      dvec = x_nb - x_up
      s = dot_product(dvec, dir_hat)
      if (s >= -1.0e-12_rp) cycle

      dist2 = max(sum(dvec**2), SMALL)
      if (s*s < 0.25_rp*dist2) cycle  ! reject before sqrt; s<0 already verified above
      align = min(1.0_rp, -s / sqrt(dist2))
      weight = (0.25_rp + align*align) / dist2
      do k = 1, 3
        delta_phi(k) = axial_delta(U(k, nb), phi_U(k), grad_up(:, k), dvec, s, dir_hat)
      end do

      a11 = a11 + weight * s * s
      a12 = a12 + weight * s * s * s
      a22 = a22 + weight * s * s * s * s
      b1 = b1 + weight * s * delta_phi
      b2 = b2 + weight * s * s * delta_phi
      lin_num = lin_num + weight * s * delta_phi
      lin_den = lin_den + weight * s * s
      nfit = nfit + 1
    end do

    det = a11*a22 - a12*a12
    do k = 1, 3
      if (nfit >= 2 .and. abs(det) > SMALL) then
        a1(k) = ( b1(k)*a22 - b2(k)*a12) / det
        a2(k) = (-b1(k)*a12 + b2(k)*a11) / det
      else if (lin_den > SMALL) then
        a1(k) = lin_num(k) / lin_den
        a2(k) = 0.0_rp
      else
        a1(k) = (phi_D(k) - phi_U(k)) / max(s_down, SMALL)
        a2(k) = 0.0_rp
      end if
      phi_f(k) = phi_U(k) + a1(k)*s_face + a2(k)*s_face*s_face
      phi_f(k) = clamp_value(phi_f(k), min(phi_U(k), phi_D(k)), max(phi_U(k), phi_D(k)))
    end do
  end subroutine quick_reconstruction_vec_ctx

  !> MUSCL van Leer TVD (2nd order): phi_f = phi_U + 0.5*psi(r)*(phi_D-phi_U); public API, hot loop uses the _ctx variant.
  real(rp) function vanleer_tvd_reconstruction(phi, f, Fn, msh) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    real(rp),    intent(in) :: Fn  ! positive = flow from owner to neighbour
    type(Mesh),  intent(in) :: msh

    integer(ip) :: up, down, far_up
    real(rp)    :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp)    :: s_face, s_down
    logical     :: ok

    ok = resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down)
    if (.not. ok .or. up > size(phi) .or. down > size(phi)) then
      phi_f = face_interpolate_scalar(phi, f, msh)
      return
    end if

    far_up = find_far_upwind(up, down, x_up, up_shift, dir_hat, size(phi), msh)
    phi_f = vanleer_tvd_reconstruction_ctx(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, far_up, msh)
  end function vanleer_tvd_reconstruction

  !> Core of vanleer_tvd_reconstruction, taking an already-resolved reconstruction
  !> context and far-upwind cell instead of re-deriving them from (f, Fn, msh).
  real(rp) function vanleer_tvd_reconstruction_ctx(phi, up, down, x_up, up_shift, dir_hat, &
                                                    s_face, s_down, far_up, msh) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down, far_up
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3)
    real(rp),    intent(in) :: s_face, s_down
    type(Mesh),  intent(in) :: msh

    real(rp)    :: x_far(3)
    real(rp)    :: phi_U, phi_D, phi_UU, r, psi, slope_dn, slope_up, s_up

    phi_U  = phi(up)
    phi_D  = phi(down)

    slope_dn = (phi_D - phi_U) / max(s_down, SMALL)
    if (far_up >= 1 .and. abs(slope_dn) > SMALL) then
      x_far = msh%cells(far_up)%centroid + up_shift
      s_up = -dot_product(x_far - x_up, dir_hat)
      if (s_up <= SMALL) then
        psi = 0.0_rp
      else
      phi_UU = phi(far_up)
        slope_up = (phi_U - phi_UU) / s_up
        r = slope_up / slope_dn
        psi = (r + abs(r)) / (1.0_rp + abs(r))
        psi = max(0.0_rp, min(psi, 2.0_rp))
      end if
    else
      psi = 0.0_rp
      !$omp atomic
      tvd_fallback_count = tvd_fallback_count + 1_ip  ! far-upwind stencil missing
    end if

    phi_f = phi_U + psi * slope_dn * s_face
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function vanleer_tvd_reconstruction_ctx

  !> Find the cell on the far-upwind side of cell `up` given that `down` is
  !> the downwind cell. Returns -1 if not found.
  integer(ip) function find_far_upwind(up, down, x_up, up_shift, dir_hat, nvals, msh) result(far_up)
    integer(ip), intent(in) :: up, down, nvals
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3)
    type(Mesh),  intent(in) :: msh
    integer(ip) :: fi, fj, nb
    integer(ip) :: f_bnd_fu, patch_id_fu
    real(rp)    :: d_up_nb(3), s, dist2, align, perp2
    real(rp)    :: best_align, best_perp2

    far_up = -1
    best_align = -1.0_rp
    best_perp2 = huge(1.0_rp)
    do fi = msh%cells(up)%face_offset, msh%cells(up)%face_offset + msh%cells(up)%n_faces - 1
      fj = msh%cell_faces(fi)
      nb = msh%faces(fj)%owner
      if (nb == up) nb = msh%faces(fj)%neighbour
      ! Resolve periodic partner for TVD far-upwind stencil
      if (nb < 1) then
        f_bnd_fu = fj - msh%n_int_faces
        patch_id_fu = msh%faces(fj)%patch_id
        if (f_bnd_fu >= 1 .and. patch_id_fu >= 1 .and. patch_id_fu <= msh%n_patches) then
          if (msh%patches(patch_id_fu)%bc_id == BC_PERIODIC) then
            nb = resolve_periodic_nb(msh, f_bnd_fu, up)
          end if
        end if
      end if
      if (nb == down .or. nb < 1 .or. nb > nvals) cycle
      d_up_nb = msh%cells(nb)%centroid + up_shift - x_up
      s = dot_product(d_up_nb, dir_hat)
      if (s >= -1.0e-12_rp) cycle
      dist2 = max(sum(d_up_nb**2), SMALL)
      align = -s / sqrt(dist2)
      perp2 = max(0.0_rp, dist2 - s*s)
      if (align > best_align + 1.0e-12_rp .or. &
          (abs(align - best_align) <= 1.0e-12_rp .and. perp2 < best_perp2 - 1.0e-12_rp) .or. &
          (abs(align - best_align) <= 1.0e-12_rp .and. abs(perp2 - best_perp2) <= 1.0e-12_rp .and. nb < far_up)) then
        best_align = align
        best_perp2 = perp2
        far_up = nb
      end if
    end do
  end function find_far_upwind

  ! Koren-TVD reconstruction

  !> MUSCL Koren TVD (3rd order in smooth regions, 1st-order near extrema); public API, hot loop uses the _ctx variant.
  real(rp) function koren_tvd_reconstruction(phi, f, Fn, msh) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    real(rp),    intent(in) :: Fn
    type(Mesh),  intent(in) :: msh

    integer(ip) :: up, down, far_up
    real(rp)    :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp)    :: s_face, s_down
    logical     :: ok

    ok = resolve_reconstruction_context(f, Fn, msh, up, down, x_up, x_down, x_face, up_shift, dir_hat, s_face, s_down)
    if (.not. ok .or. up > size(phi) .or. down > size(phi)) then
      phi_f = face_interpolate_scalar(phi, f, msh)
      return
    end if

    far_up = find_far_upwind(up, down, x_up, up_shift, dir_hat, size(phi), msh)
    phi_f = koren_tvd_reconstruction_ctx(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, far_up, msh)
  end function koren_tvd_reconstruction

  !> Core of koren_tvd_reconstruction, taking an already-resolved reconstruction
  !> context and far-upwind cell instead of re-deriving them from (f, Fn, msh).
  real(rp) function koren_tvd_reconstruction_ctx(phi, up, down, x_up, up_shift, dir_hat, &
                                                  s_face, s_down, far_up, msh) result(phi_f)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down, far_up
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3)
    real(rp),    intent(in) :: s_face, s_down
    type(Mesh),  intent(in) :: msh

    real(rp)    :: x_far(3)
    real(rp)    :: phi_U, phi_D, phi_UU, r, psi, slope_dn, slope_up, s_up

    phi_U  = phi(up)
    phi_D  = phi(down)

    slope_dn = (phi_D - phi_U) / max(s_down, SMALL)
    if (far_up >= 1 .and. abs(slope_dn) > SMALL) then
      x_far = msh%cells(far_up)%centroid + up_shift
      s_up = -dot_product(x_far - x_up, dir_hat)
      if (s_up <= SMALL) then
        psi = 0.0_rp
      else
      phi_UU = phi(far_up)
        slope_up = (phi_U - phi_UU) / s_up
        r = slope_up / slope_dn
        psi = max(0.0_rp, min(2.0_rp*r, (1.0_rp + 2.0_rp*r)/3.0_rp, 2.0_rp))
      end if
    else
      psi = 0.0_rp  ! first-order upwind at extrema or missing stencil
      !$omp atomic
      tvd_fallback_count = tvd_fallback_count + 1_ip  ! far-upwind stencil missing
    end if

    phi_f = phi_U + psi * slope_dn * s_face
    ! Hard TVD bound: face value within [phi_U, phi_D] interval
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function koren_tvd_reconstruction_ctx

  pure real(rp) function clamp_value(x, lo, hi) result(y)
    !$acc routine seq
    real(rp), intent(in) :: x, lo, hi
    real(rp) :: l, h
    l = min(lo, hi)
    h = max(lo, hi)
    y = min(max(x, l), h)
  end function clamp_value

  ! ==========================================================================
  ! Device-callable SoA counterparts of the reconstruction-context helpers
  ! above, for the QUICK/VANLEER_TVD/KOREN_TVD GPU fast path. Diagnostic-only
  ! side effects (tvd_fallback_count) are dropped; everything numerically
  ! relevant is preserved.
  ! ==========================================================================

  !> Device-callable SoA counterpart of resolve_periodic_nb.
  integer(ip) function soa_resolve_periodic_nb_dev(f_bnd, owner, n_int, nf, pneigh_soa, ppair_soa, fo_soa, fn_soa) &
                        result(nb)
    !$acc routine seq
    integer(ip), intent(in) :: f_bnd, owner, n_int, nf
    integer(ip), intent(in) :: pneigh_soa(:), ppair_soa(:), fo_soa(:), fn_soa(:)
    integer(ip) :: pair_bnd, pair_face

    nb = -1
    if (f_bnd < 1 .or. f_bnd > size(pneigh_soa)) return
    nb = pneigh_soa(f_bnd)

    if (nb <= 0 .and. f_bnd <= size(ppair_soa)) then
      pair_bnd = ppair_soa(f_bnd)
      if (pair_bnd > 0) then
        pair_face = n_int + pair_bnd
        if (pair_face >= 1 .and. pair_face <= nf) nb = fo_soa(pair_face)
      end if
    end if

    if (nb <= 0) then
      pair_face = n_int + f_bnd
      if (pair_face >= 1 .and. pair_face <= nf) nb = fn_soa(pair_face)
    end if

    if (nb == owner) nb = -1
  end function soa_resolve_periodic_nb_dev

  !> Device-callable SoA counterpart of resolve_periodic_partner.
  subroutine soa_resolve_periodic_partner_dev(f, n_int, nbnd, nf, n_patches, fo_soa, fn_soa, fpatch_soa, pbc_soa, &
                                               pneigh_soa, ppair_soa, ppair_dcf_soa, fdcf_soa, fnorm_soa, ccent_soa, &
                                               pair_cell, pair_face, shift)
    !$acc routine seq
    integer(ip), intent(in)  :: f, n_int, nbnd, nf, n_patches
    integer(ip), intent(in)  :: fo_soa(:), fn_soa(:), fpatch_soa(:), pbc_soa(:), pneigh_soa(:), ppair_soa(:)
    real(rp),    intent(in)  :: ppair_dcf_soa(:), fdcf_soa(:), fnorm_soa(:,:), ccent_soa(:,:)
    integer(ip), intent(out) :: pair_cell, pair_face
    real(rp),    intent(out) :: shift(3)

    integer(ip) :: f_bnd, patch_id, pair_bnd, owner
    real(rp)    :: d_per

    pair_cell = 0
    pair_face = 0
    shift = 0.0_rp

    if (fn_soa(f) >= 1) return
    patch_id = fpatch_soa(f)
    if (patch_id < 1 .or. patch_id > n_patches) return
    if (pbc_soa(patch_id) /= BC_PERIODIC) return

    f_bnd = f - n_int
    if (f_bnd < 1 .or. f_bnd > nbnd) return

    owner = fo_soa(f)
    pair_cell = soa_resolve_periodic_nb_dev(f_bnd, owner, n_int, nf, pneigh_soa, ppair_soa, fo_soa, fn_soa)
    if (f_bnd <= size(ppair_soa)) then
      pair_bnd = ppair_soa(f_bnd)
      if (pair_bnd > 0) then
        pair_face = n_int + pair_bnd
        if (pair_face >= 1 .and. pair_face <= nf) then
          if (pair_cell <= 0) pair_cell = fo_soa(pair_face)
        end if
      end if
    end if

    if (pair_cell > 0 .and. pair_cell <= size(ccent_soa, 2)) then
      if (pair_face >= 1 .and. pair_face <= nf) then
        d_per = fdcf_soa(f) + fdcf_soa(pair_face)
      else
        if (f_bnd >= 1 .and. f_bnd <= size(ppair_dcf_soa) .and. ppair_dcf_soa(f_bnd) > SMALL) then
          d_per = fdcf_soa(f) + ppair_dcf_soa(f_bnd)
        else
          d_per = 2.0_rp * fdcf_soa(f)
        end if
      end if
      shift = (ccent_soa(:, owner) + d_per * fnorm_soa(:, f)) - ccent_soa(:, pair_cell)
    end if
  end subroutine soa_resolve_periodic_partner_dev

  !> Device-callable SoA counterpart of resolve_reconstruction_context.
  logical function soa_resolve_reconstruction_context_dev(f, Fn, n_int, nbnd, n_patches, nf, nu, &
                        fo_soa, fn_soa, fcentre_soa, ccent_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, &
                        ppair_dcf_soa, fdcf_soa, fnorm_soa, up, down, x_up, x_down, x_face, up_shift, dir_hat, &
                        s_face, s_down) result(ok)
    !$acc routine seq
    integer(ip), intent(in)  :: f, n_int, nbnd, n_patches, nf, nu
    real(rp),    intent(in)  :: Fn
    integer(ip), intent(in)  :: fo_soa(:), fn_soa(:), fpatch_soa(:), pbc_soa(:), pneigh_soa(:), ppair_soa(:)
    real(rp),    intent(in)  :: fcentre_soa(:,:), ccent_soa(:,:), ppair_dcf_soa(:), fdcf_soa(:), fnorm_soa(:,:)
    integer(ip), intent(out) :: up, down
    real(rp),    intent(out) :: x_up(3), x_down(3), x_face(3), up_shift(3), dir_hat(3)
    real(rp),    intent(out) :: s_face, s_down

    integer(ip) :: owner, nb, pair_cell, pair_face
    real(rp)    :: pair_shift(3), dvec(3), dnorm

    ok = .false.
    up = 0
    down = 0
    x_up = 0.0_rp
    x_down = 0.0_rp
    x_face = fcentre_soa(:, f)
    up_shift = 0.0_rp
    dir_hat = 0.0_rp
    s_face = 0.0_rp
    s_down = 0.0_rp

    owner = fo_soa(f)
    nb = fn_soa(f)

    if (nb >= 1) then
      if (Fn >= 0.0_rp) then
        up = owner
        down = nb
      else
        up = nb
        down = owner
      end if
      x_up = ccent_soa(:, up)
      x_down = ccent_soa(:, down)
    else
      call soa_resolve_periodic_partner_dev(f, n_int, nbnd, nf, n_patches, fo_soa, fn_soa, fpatch_soa, pbc_soa, &
                                             pneigh_soa, ppair_soa, ppair_dcf_soa, fdcf_soa, fnorm_soa, ccent_soa, &
                                             pair_cell, pair_face, pair_shift)
      if (pair_cell <= 0) return
      if (Fn >= 0.0_rp) then
        up = owner
        down = pair_cell
        x_up = ccent_soa(:, owner)
        x_down = ccent_soa(:, pair_cell) + pair_shift
      else
        up = pair_cell
        down = owner
        up_shift = pair_shift
        x_up = ccent_soa(:, pair_cell) + up_shift
        x_down = ccent_soa(:, owner)
      end if
    end if

    if (up < 1 .or. down < 1 .or. up > nu .or. down > nu) return

    dvec = x_face - x_up
    dnorm = sqrt(sum(dvec**2))
    if (dnorm <= SMALL) then
      dvec = x_down - x_up
      dnorm = sqrt(sum(dvec**2))
    end if
    if (dnorm <= SMALL) return
    dir_hat = dvec / dnorm
    s_face = dot_product(x_face - x_up, dir_hat)
    s_down = dot_product(x_down - x_up, dir_hat)

    if (s_face <= SMALL .or. s_down <= s_face) then
      dvec = x_down - x_up
      dnorm = sqrt(sum(dvec**2))
      if (dnorm <= SMALL) return
      dir_hat = dvec / dnorm
      s_face = dot_product(x_face - x_up, dir_hat)
      s_down = dot_product(x_down - x_up, dir_hat)
    end if
    if (s_down <= SMALL) return

    s_face = min(max(s_face, 0.0_rp), s_down)
    ok = .true.
  end function soa_resolve_reconstruction_context_dev

  !> Device-callable SoA counterpart of find_far_upwind.
  integer(ip) function soa_find_far_upwind_dev(up, down, x_up, up_shift, dir_hat, nvals, n_int, n_patches, nf, nc, &
                        cfo_soa, cell_faces, fo_soa, fn_soa, fpatch_soa, pbc_soa, pneigh_soa, ppair_soa, ccent_soa) &
                        result(far_up)
    !$acc routine seq
    integer(ip), intent(in) :: up, down, nvals, n_int, n_patches, nf, nc
    integer(ip), intent(in) :: cfo_soa(:), cell_faces(:), fo_soa(:), fn_soa(:), fpatch_soa(:), pbc_soa(:)
    integer(ip), intent(in) :: pneigh_soa(:), ppair_soa(:)
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3), ccent_soa(:,:)

    integer(ip) :: fi, fj, nb
    integer(ip) :: f_bnd_fu, patch_id_fu
    real(rp)    :: d_up_nb(3), s, dist2, align, perp2
    real(rp)    :: best_align, best_perp2

    far_up = -1
    if (up < 1 .or. up > nc) return
    best_align = -1.0_rp
    best_perp2 = huge(1.0_rp)
    do fi = cfo_soa(up), cfo_soa(up+1) - 1
      fj = cell_faces(fi)
      nb = fo_soa(fj)
      if (nb == up) nb = fn_soa(fj)
      if (nb < 1) then
        f_bnd_fu = fj - n_int
        patch_id_fu = fpatch_soa(fj)
        if (f_bnd_fu >= 1 .and. patch_id_fu >= 1 .and. patch_id_fu <= n_patches) then
          if (pbc_soa(patch_id_fu) == BC_PERIODIC) then
            nb = soa_resolve_periodic_nb_dev(f_bnd_fu, up, n_int, nf, pneigh_soa, ppair_soa, fo_soa, fn_soa)
          end if
        end if
      end if
      if (nb == down .or. nb < 1 .or. nb > nvals) cycle
      d_up_nb = ccent_soa(:, nb) + up_shift - x_up
      s = dot_product(d_up_nb, dir_hat)
      if (s >= -1.0e-12_rp) cycle
      dist2 = max(sum(d_up_nb**2), SMALL)
      align = -s / sqrt(dist2)
      perp2 = max(0.0_rp, dist2 - s*s)
      if (align > best_align + 1.0e-12_rp .or. &
          (abs(align - best_align) <= 1.0e-12_rp .and. perp2 < best_perp2 - 1.0e-12_rp) .or. &
          (abs(align - best_align) <= 1.0e-12_rp .and. abs(perp2 - best_perp2) <= 1.0e-12_rp .and. nb < far_up)) then
        best_align = align
        best_perp2 = perp2
        far_up = nb
      end if
    end do
  end function soa_find_far_upwind_dev

  !> Device-callable SoA counterpart of quick_reconstruction_ctx.
  real(rp) function soa_quick_reconstruction_ctx_dev(phi, up, down, x_up, x_down, up_shift, dir_hat, s_face, s_down, &
                        nc, cfo_soa, cell_faces, fo_soa, fn_soa, ccent_soa, grad_phi, nu) result(phi_f)
    !$acc routine seq
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down, nc, nu
    real(rp),    intent(in) :: x_up(3), x_down(3), up_shift(3), dir_hat(3), s_face, s_down
    integer(ip), intent(in) :: cfo_soa(:), cell_faces(:), fo_soa(:), fn_soa(:)
    real(rp),    intent(in) :: ccent_soa(:,:), grad_phi(:,:)

    integer(ip) :: fi, fj, nb, nfit
    real(rp)    :: phi_U, phi_D, s, dist2, align, weight
    real(rp)    :: a11, a12, a22, b1, b2, det, a1, a2, lin_num, lin_den
    real(rp)    :: x_nb(3), dvec(3), delta_phi, grad_up(3)

    phi_U = phi(up)
    phi_D = phi(down)

    grad_up = 0.0_rp
    if (up <= size(grad_phi, 2)) grad_up = grad_phi(:, up)

    a11 = 0.0_rp; a12 = 0.0_rp; a22 = 0.0_rp
    b1 = 0.0_rp; b2 = 0.0_rp
    lin_num = 0.0_rp; lin_den = 0.0_rp
    nfit = 0

    dvec = x_down - x_up
    dist2 = max(sum(dvec**2), SMALL)
    weight = 1.0_rp / dist2
    delta_phi = axial_delta(phi_D, phi_U, grad_up, dvec, s_down, dir_hat)
    a11 = a11 + weight * s_down * s_down
    a12 = a12 + weight * s_down * s_down * s_down
    a22 = a22 + weight * s_down * s_down * s_down * s_down
    b1 = b1 + weight * s_down * delta_phi
    b2 = b2 + weight * s_down * s_down * delta_phi
    lin_num = lin_num + weight * s_down * delta_phi
    lin_den = lin_den + weight * s_down * s_down
    nfit = nfit + 1

    if (up >= 1 .and. up <= nc) then
      do fi = cfo_soa(up), cfo_soa(up+1) - 1
        fj = cell_faces(fi)
        nb = fo_soa(fj)
        if (nb == up) nb = fn_soa(fj)
        if (nb < 1 .or. nb == down .or. nb > nu) cycle

        x_nb = ccent_soa(:, nb) + up_shift
        dvec = x_nb - x_up
        s = dot_product(dvec, dir_hat)
        if (s >= -1.0e-12_rp) cycle

        dist2 = max(sum(dvec**2), SMALL)
        if (s*s < 0.25_rp*dist2) cycle  ! reject before sqrt; s<0 already verified above
        align = min(1.0_rp, -s / sqrt(dist2))
        weight = (0.25_rp + align*align) / dist2
        delta_phi = axial_delta(phi(nb), phi_U, grad_up, dvec, s, dir_hat)

        a11 = a11 + weight * s * s
        a12 = a12 + weight * s * s * s
        a22 = a22 + weight * s * s * s * s
        b1 = b1 + weight * s * delta_phi
        b2 = b2 + weight * s * s * delta_phi
        lin_num = lin_num + weight * s * delta_phi
        lin_den = lin_den + weight * s * s
        nfit = nfit + 1
      end do
    end if

    det = a11*a22 - a12*a12
    if (nfit >= 2 .and. abs(det) > SMALL) then
      a1 = ( b1*a22 - b2*a12) / det
      a2 = (-b1*a12 + b2*a11) / det
    else if (lin_den > SMALL) then
      a1 = lin_num / lin_den
      a2 = 0.0_rp
    else
      a1 = (phi_D - phi_U) / max(s_down, SMALL)
      a2 = 0.0_rp
    end if

    phi_f = phi_U + a1*s_face + a2*s_face*s_face
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function soa_quick_reconstruction_ctx_dev

  !> Device-callable SoA counterpart of vanleer_tvd_reconstruction_ctx (diagnostic fallback counter dropped).
  real(rp) function soa_vanleer_tvd_ctx_dev(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, far_up, ccent_soa) &
                        result(phi_f)
    !$acc routine seq
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down, far_up
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3), s_face, s_down, ccent_soa(:,:)

    real(rp) :: x_far(3)
    real(rp) :: phi_U, phi_D, phi_UU, r, psi, slope_dn, slope_up, s_up

    phi_U = phi(up)
    phi_D = phi(down)

    slope_dn = (phi_D - phi_U) / max(s_down, SMALL)
    if (far_up >= 1 .and. abs(slope_dn) > SMALL) then
      x_far = ccent_soa(:, far_up) + up_shift
      s_up = -dot_product(x_far - x_up, dir_hat)
      if (s_up <= SMALL) then
        psi = 0.0_rp
      else
        phi_UU = phi(far_up)
        slope_up = (phi_U - phi_UU) / s_up
        r = slope_up / slope_dn
        psi = (r + abs(r)) / (1.0_rp + abs(r))
        psi = max(0.0_rp, min(psi, 2.0_rp))
      end if
    else
      psi = 0.0_rp
    end if

    phi_f = phi_U + psi * slope_dn * s_face
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function soa_vanleer_tvd_ctx_dev

  !> Device-callable SoA counterpart of koren_tvd_reconstruction_ctx (diagnostic fallback counter dropped).
  real(rp) function soa_koren_tvd_ctx_dev(phi, up, down, x_up, up_shift, dir_hat, s_face, s_down, far_up, ccent_soa) &
                        result(phi_f)
    !$acc routine seq
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: up, down, far_up
    real(rp),    intent(in) :: x_up(3), up_shift(3), dir_hat(3), s_face, s_down, ccent_soa(:,:)

    real(rp) :: x_far(3)
    real(rp) :: phi_U, phi_D, phi_UU, r, psi, slope_dn, slope_up, s_up

    phi_U = phi(up)
    phi_D = phi(down)

    slope_dn = (phi_D - phi_U) / max(s_down, SMALL)
    if (far_up >= 1 .and. abs(slope_dn) > SMALL) then
      x_far = ccent_soa(:, far_up) + up_shift
      s_up = -dot_product(x_far - x_up, dir_hat)
      if (s_up <= SMALL) then
        psi = 0.0_rp
      else
        phi_UU = phi(far_up)
        slope_up = (phi_U - phi_UU) / s_up
        r = slope_up / slope_dn
        psi = max(0.0_rp, min(2.0_rp*r, (1.0_rp + 2.0_rp*r)/3.0_rp, 2.0_rp))
      end if
    else
      psi = 0.0_rp
    end if

    phi_f = phi_U + psi * slope_dn * s_face
    phi_f = clamp_value(phi_f, min(phi_U, phi_D), max(phi_U, phi_D))
  end function soa_koren_tvd_ctx_dev

  !> Reset the TVD fallback counter (call at the start of each time step).
  subroutine advection_reset_tvd_fallback()
    tvd_fallback_count = 0_ip
  end subroutine advection_reset_tvd_fallback

  !> Return the current TVD fallback count (number of faces using 1st-order upwind this step).
  integer(ip) function advection_get_tvd_fallback_count() result(n)
    n = tvd_fallback_count
  end function advection_get_tvd_fallback_count

  ! ensure_buf_1d/ensure_buf_n3/ensure_buf_3n3/ensure_fcp_identity now live in mod_gradient (used here via `use mod_gradient`).

end module mod_advection
