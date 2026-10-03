! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Germano-Lilly dynamic Smagorinsky (1991/1992), localized to a compact one-ring test filter for unstructured meshes.
module mod_dynamic_smagorinsky
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  use mod_vreman, only: vreman_grid_delta
  implicit none
  private

  public :: dynamic_smagorinsky_sgs
  public :: dynamic_smagorinsky_free

  ! Test filter width / grid filter width; standard choice for a doubled top-hat filter.
  real(rp), parameter :: ALPHA_SQ = 4.0_rp

  ! Same cap rationale as NU_SGS_STAB_FRAC in mod_vreman: bounds nu_sgs by the explicit
  ! diffusive stability limit viscous_rhs_vector relies on.
  real(rp), parameter :: NU_SGS_STAB_FRAC = 0.4_rp

  ! Persistent scratch buffers, grown as needed instead of alloc/dealloc every call.
  real(rp), allocatable, save :: field_buf(:,:)    ! (15, n_cells+n_ghost): U, UiUj, (|S|Sij)
  real(rp), allocatable, save :: field_hat(:,:)    ! (15, n_cells): test-filtered field_buf
  real(rp), allocatable, save :: uhat_buf(:,:)     ! (3, n_cells+n_ghost)
  real(rp), allocatable, save :: lm_mm(:,:)        ! (2, n_cells+n_ghost): LijMij, MijMij
  real(rp), allocatable, save :: lm_mm_hat(:,:)    ! (2, n_cells)
  real(rp), allocatable, save :: delta_buf(:)      ! (n_cells)
  real(rp), allocatable, save :: grad_uhat(:,:,:)   ! (3, 3, n_cells)
  logical,  save :: delta_ready = .false.  ! delta_buf holds the (static) filter width, on host and device
  ! The scratch buffers live on the device (allocate/deallocate mirror automatically), so the
  ! whole procedure -- packing, test filtering, the halo exchanges -- runs there.
  !$acc declare create(field_buf, field_hat, uhat_buf, lm_mm, lm_mm_hat, delta_buf, grad_uhat)

contains

  !> Dynamic Smagorinsky eddy viscosity: nu_t = Cs^2 * Delta^2 * |S|, Cs^2 from Lilly (1992), clipped to >=0.
  subroutine dynamic_smagorinsky_sgs(U, grad_U, msh, nu_sgs, U_bnd, dt_stab)
    real(rp),   intent(in)    :: U(:,:)       ! (3, n_cells+n_ghost)
    real(rp),   intent(in)    :: grad_U(:,:,:)  ! (3, 3, n_cells)
    type(Mesh), intent(inout) :: msh
    real(rp),   intent(inout) :: nu_sgs(:)    ! (n_cells)
    real(rp),   intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),   intent(in), optional :: dt_stab

    integer(ip) :: c, nc, ntot
    real(rp)    :: S(3,3), S_norm, dt_stab_local, nu_cap
    logical     :: has_dt_stab

    nc   = msh%n_cells
    ntot = nc + msh%n_ghost
    call ensure_buffers(nc, ntot)

    has_dt_stab = present(dt_stab)
    dt_stab_local = 1.0_rp
    if (has_dt_stab) dt_stab_local = dt_stab

    ! Filter width is static mesh geometry: computed on the host once, pushed once.
    if (.not. delta_ready) then
      call vreman_grid_delta(msh, delta_buf)
      !$acc update device(delta_buf(1:nc))
      delta_ready = .true.
    end if

    ! Pack U, the six unique UiUj, and the six unique (|S| S_ij) at owned cells.
#ifdef _OPENACC
    !$acc data copyin(U(:,1:nc), grad_U(:,:,1:nc))
    !$acc parallel loop present(field_buf) private(S, S_norm)
#else
    !$omp parallel do schedule(static) private(c, S, S_norm)
#endif
    do c = 1, nc
      field_buf(1:3, c) = U(1:3, c)
      field_buf(4, c) = U(1,c)*U(1,c)
      field_buf(5, c) = U(1,c)*U(2,c)
      field_buf(6, c) = U(1,c)*U(3,c)
      field_buf(7, c) = U(2,c)*U(2,c)
      field_buf(8, c) = U(2,c)*U(3,c)
      field_buf(9, c) = U(3,c)*U(3,c)

      S(1,1) = grad_U(1,1,c);                        S(1,2) = 0.5_rp*(grad_U(1,2,c)+grad_U(2,1,c))
      S(1,3) = 0.5_rp*(grad_U(1,3,c)+grad_U(3,1,c));  S(2,2) = grad_U(2,2,c)
      S(2,3) = 0.5_rp*(grad_U(2,3,c)+grad_U(3,2,c));  S(3,3) = grad_U(3,3,c)
      S_norm = sqrt(2.0_rp * (S(1,1)**2 + S(2,2)**2 + S(3,3)**2 &
                    + 2.0_rp*(S(1,2)**2 + S(1,3)**2 + S(2,3)**2)))

      field_buf(10, c) = S_norm * S(1,1)
      field_buf(11, c) = S_norm * S(1,2)
      field_buf(12, c) = S_norm * S(1,3)
      field_buf(13, c) = S_norm * S(2,2)
      field_buf(14, c) = S_norm * S(2,3)
      field_buf(15, c) = S_norm * S(3,3)
    end do
#ifdef _OPENACC
    !$acc end parallel loop
    !$acc end data
#else
    !$omp end parallel do
#endif

    call msh%exchange_halos_scalars(field_buf(:, 1:ntot))
    call test_filter(msh, field_buf, field_hat, 15_ip)

    ! Grad of the test-filtered velocity needs U_hat at ghost cells too: halo-exchange it.
    !$acc parallel loop present(uhat_buf, field_hat)
    do c = 1, nc
      uhat_buf(1, c) = field_hat(1, c)
      uhat_buf(2, c) = field_hat(2, c)
      uhat_buf(3, c) = field_hat(3, c)
    end do
    call msh%exchange_halos_vec(uhat_buf(:, 1:ntot))
    if (present(U_bnd)) then
      call green_gauss_gradient_vector(uhat_buf(:, 1:ntot), msh, grad_uhat, U_bnd)
    else
      call green_gauss_gradient_vector(uhat_buf(:, 1:ntot), msh, grad_uhat)
    end if

#ifdef _OPENACC
    !$acc parallel loop present(field_hat, grad_uhat, delta_buf, lm_mm)
#else
    !$omp parallel do schedule(static) private(c)
#endif
    do c = 1, nc
      call germano_lilly_terms(field_hat(:, c), grad_uhat(:, :, c), delta_buf(c), &
                                lm_mm(1, c), lm_mm(2, c))
    end do
#ifdef _OPENACC
    !$acc end parallel loop
#else
    !$omp end parallel do
#endif

    ! Local averaging of the numerator/denominator over the compact stencil (Lilly's ratio is noisy pointwise), same rationale as mod_vreman's refinement-interface smoothing.
    call msh%exchange_halos_scalars(lm_mm(:, 1:ntot))
    call test_filter(msh, lm_mm, lm_mm_hat, 2_ip)

#ifdef _OPENACC
    !$acc data present(lm_mm_hat, delta_buf) copyin(grad_U(:,:,1:nc)) copyout(nu_sgs(1:nc))
    !$acc parallel loop private(nu_cap)
#else
    !$omp parallel do schedule(static) private(c, nu_cap)
#endif
    do c = 1, nc
      ! Cs^2 clipped to >=0: a negative ratio corresponds to backscatter, which a purely
      ! dissipative eddy-viscosity closure cannot represent without destabilizing the RK stage.
      nu_sgs(c) = max(lm_mm_hat(1, c) / max(lm_mm_hat(2, c), SMALL), 0.0_rp) &
                * delta_buf(c)**2 * s_norm_of(grad_U(:, :, c))

      if (has_dt_stab) then
        nu_cap = NU_SGS_STAB_FRAC * delta_buf(c)**2 / max(dt_stab_local, SMALL)
        nu_sgs(c) = min(nu_sgs(c), nu_cap)
      end if
    end do
#ifdef _OPENACC
    !$acc end parallel loop
    !$acc end data
#else
    !$omp end parallel do
#endif
  end subroutine dynamic_smagorinsky_sgs

  !> |S| from a raw velocity-gradient tensor (same convention as smagorinsky_sgs).
  pure real(rp) function s_norm_of(grad_U_c) result(S_norm)
    !$acc routine seq
    real(rp), intent(in) :: grad_U_c(3,3)
    real(rp) :: S(3,3)
    integer(ip) :: i, j

    do j = 1, 3
      do i = 1, 3
        S(i,j) = 0.5_rp * (grad_U_c(i,j) + grad_U_c(j,i))
      end do
    end do
    S_norm = sqrt(2.0_rp * sum(S*S))
  end function s_norm_of

  !> Germano identity at one cell: L_ij = (UiUj)^ - Ui^ Uj^, M_ij = 2*[(|S|Sij)^ - alpha^2 |S^| S^_ij],
  !> contracted to the two scalars the Lilly ratio needs (LijMij, MijMij).
  pure subroutine germano_lilly_terms(field_hat_c, grad_uhat_c, delta_c, lij_mij, mij_mij)
    !$acc routine seq
    real(rp), intent(in)  :: field_hat_c(15)  ! packed: Uhat(1:3), (UiUj)^(4:9), (|S|Sij)^(10:15)
    real(rp), intent(in)  :: grad_uhat_c(3,3)
    real(rp), intent(in)  :: delta_c
    real(rp), intent(out) :: lij_mij, mij_mij

    real(rp) :: L(3,3), M(3,3), Shat(3,3), Shat_norm
    integer(ip) :: i, j

    do j = 1, 3
      do i = 1, 3
        Shat(i,j) = 0.5_rp * (grad_uhat_c(i,j) + grad_uhat_c(j,i))
      end do
    end do
    Shat_norm = sqrt(2.0_rp * sum(Shat*Shat))

    L(1,1) = field_hat_c(4) - field_hat_c(1)*field_hat_c(1)
    L(1,2) = field_hat_c(5) - field_hat_c(1)*field_hat_c(2)
    L(1,3) = field_hat_c(6) - field_hat_c(1)*field_hat_c(3)
    L(2,2) = field_hat_c(7) - field_hat_c(2)*field_hat_c(2)
    L(2,3) = field_hat_c(8) - field_hat_c(2)*field_hat_c(3)
    L(3,3) = field_hat_c(9) - field_hat_c(3)*field_hat_c(3)

    M(1,1) = 2.0_rp * (field_hat_c(10) - ALPHA_SQ * Shat_norm * Shat(1,1))
    M(1,2) = 2.0_rp * (field_hat_c(11) - ALPHA_SQ * Shat_norm * Shat(1,2))
    M(1,3) = 2.0_rp * (field_hat_c(12) - ALPHA_SQ * Shat_norm * Shat(1,3))
    M(2,2) = 2.0_rp * (field_hat_c(13) - ALPHA_SQ * Shat_norm * Shat(2,2))
    M(2,3) = 2.0_rp * (field_hat_c(14) - ALPHA_SQ * Shat_norm * Shat(2,3))
    M(3,3) = 2.0_rp * (field_hat_c(15) - ALPHA_SQ * Shat_norm * Shat(3,3))
    M = delta_c**2 * M

    lij_mij = L(1,1)*M(1,1) + L(2,2)*M(2,2) + L(3,3)*M(3,3) &
            + 2.0_rp*(L(1,2)*M(1,2) + L(1,3)*M(1,3) + L(2,3)*M(2,3))
    mij_mij = M(1,1)**2 + M(2,2)**2 + M(3,3)**2 &
            + 2.0_rp*(M(1,2)**2 + M(1,3)**2 + M(2,3)**2)
  end subroutine germano_lilly_terms

  !> Compact top-hat test filter: volume-weighted average of a cell and its immediate
  !> face-neighbours (periodic partners included), applied to n_comp packed fields at once.
  subroutine test_filter(msh, field, field_out, n_comp)
    type(Mesh),  intent(in), target :: msh
    real(rp),    intent(in)  :: field(:,:)      ! (n_comp, n_cells+n_ghost)
    real(rp),    intent(out) :: field_out(:,:)  ! (n_comp, n_cells)
    integer(ip), intent(in)  :: n_comp

    integer(ip) :: c, k, fi, fj, o, nb, f_bnd, ntot, nc, n_int, nbnd
    real(rp)    :: wsum, w
    integer(ip), pointer :: cfo(:), cfaces(:), fo(:), fn(:), pnb(:), ppair(:)
    real(rp),    pointer :: cvol(:)
    integer(ip) :: pair_bnd

    nc = msh%n_cells
    ntot = nc + msh%n_ghost
    n_int = msh%n_int_faces
    nbnd = msh%n_bnd_faces
    cfo => msh%cell_face_offset_soa
    cfaces => msh%cell_faces
    fo => msh%face_owner_soa
    fn => msh%face_neighbour_soa
    pnb => msh%periodic_neighbor_soa
    ppair => msh%periodic_pair_soa
    cvol => msh%cell_volume_soa

    ! field/field_out are the module scratch buffers (device-resident); the mesh SoA arrays are
    ! staged once with the mesh.
#ifdef _OPENACC
    !$acc parallel loop present(field, field_out) private(fi, fj, o, nb, f_bnd, wsum, w, k, pair_bnd)
#else
    !$omp parallel do schedule(static) private(c, k, fi, fj, o, nb, f_bnd, wsum, w, pair_bnd)
#endif
    do c = 1, nc
      wsum = cvol(c)
      do k = 1, n_comp
        field_out(k, c) = wsum * field(k, c)
      end do

      do fi = cfo(c), cfo(c+1) - 1
        fj = cfaces(fi)
        if (fj < 1) cycle

        o  = fo(fj)
        nb = fn(fj)
        nb = merge(o, nb, nb == c)  ! the "other side" of the face relative to c

        if (nb < 1 .and. fj > n_int) then
          f_bnd = fj - n_int
          if (f_bnd >= 1 .and. f_bnd <= nbnd) then
            if (pnb(f_bnd) >= 1) then
              nb = pnb(f_bnd)
            else
              ! Single-rank meshes carry only the same-rank pairing, not periodic_neighbor; without
              ! this fallback the filter silently dropped every periodic partner there.
              pair_bnd = ppair(f_bnd)
              if (pair_bnd >= 1 .and. pair_bnd <= nbnd) nb = fo(n_int + pair_bnd)
            end if
          end if
        end if

        if (nb < 1 .or. nb > ntot .or. nb == c) cycle
        w = cvol(nb)
        wsum = wsum + w
        do k = 1, n_comp
          field_out(k, c) = field_out(k, c) + w * field(k, nb)
        end do
      end do

      do k = 1, n_comp
        field_out(k, c) = field_out(k, c) / max(wsum, SMALL)
      end do
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif
  end subroutine test_filter

  subroutine ensure_buffers(nc, ntot)
    integer(ip), intent(in) :: nc, ntot

    if (allocated(field_buf)) then
      if (size(field_buf, 2) == ntot) return
      deallocate(field_buf, field_hat, uhat_buf, lm_mm, lm_mm_hat, delta_buf, grad_uhat)
    end if
    delta_ready = .false.
    allocate(field_buf(15, ntot), field_hat(15, nc), uhat_buf(3, ntot), &
             lm_mm(2, ntot), lm_mm_hat(2, nc), delta_buf(nc), grad_uhat(3, 3, nc))
  end subroutine ensure_buffers

  !> Free module-level persistent allocations (call on solver teardown).
  subroutine dynamic_smagorinsky_free()
    if (allocated(field_buf))  deallocate(field_buf)
    if (allocated(field_hat))  deallocate(field_hat)
    if (allocated(uhat_buf))   deallocate(uhat_buf)
    if (allocated(lm_mm))      deallocate(lm_mm)
    if (allocated(lm_mm_hat))  deallocate(lm_mm_hat)
    if (allocated(delta_buf))  deallocate(delta_buf)
    if (allocated(grad_uhat))  deallocate(grad_uhat)
    delta_ready = .false.
  end subroutine dynamic_smagorinsky_free

end module mod_dynamic_smagorinsky
