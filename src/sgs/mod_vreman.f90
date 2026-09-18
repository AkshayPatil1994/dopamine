! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Vreman (2004) subgrid-scale eddy-viscosity model.
module mod_vreman
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  implicit none
  private

  public :: vreman_sgs_viscosity
  public :: smagorinsky_sgs
  public :: vreman_free
  public :: vreman_filter_tensor
  public :: vreman_grid_delta

  ! Vreman (2004) Eq.(15): nu_t = C * sqrt(B_beta / alpha_sq)
  real(rp), parameter :: C_v = C_VREMAN

  ! Persistent per-cell filter-width tensor G (Vreman 2004 Eq.5 generalised to polyhedra).
  real(rp), allocatable, save :: g_filt_arr(:,:,:)  ! (3, 3, n_cells + n_ghost)

  ! Persistent pack buffer for exchange_g_halo; grown as needed instead of alloc/dealloc per call.
  real(rp), allocatable, save :: g_halo_buf(:,:)  ! (6, n_cells + n_ghost)

  ! Laplacian passes over G so the filter width ramps instead of stepping across a refinement interface.
  integer(ip), parameter :: G_SMOOTH_PASSES = 2
  ! Neighbour blend weight applied in each smoothing pass.
  real(rp),    parameter :: G_SMOOTH_WEIGHT = 0.5_rp

  ! Cap on nu_sgs as a fraction of the explicit diffusive stability limit nu*dt/h^2 ~ O(1):
  ! viscous_rhs_vector applies (nu+nu_sgs) explicitly in the RK stage, so an uncapped nu_sgs
  ! can violate that limit locally (observed: nu_sgs spiking 30x at the smallest refined
  ! cells, feeding back through the rotational pressure correction into a full blowup).
  real(rp), parameter :: NU_SGS_STAB_FRAC = 0.4_rp

contains

  pure real(rp) function det3x3(A) result(d)
    !$acc routine seq
    real(rp), intent(in) :: A(3,3)
    d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
  end function det3x3

  !> Compute the Vreman SGS viscosity at all interior cells.
  !> @param U velocity field (3, n_cells + n_ghost)
  subroutine vreman_sgs_viscosity(U, msh, nu_sgs, U_bnd, grad_U_ext, dt_stab)
    real(rp),   intent(in)    :: U(:,:)  ! (3, n_cells+n_ghost)
    type(Mesh), intent(inout) :: msh  ! inout: the filter-tensor build halo-exchanges G
    real(rp),   intent(inout) :: nu_sgs(:)  ! (n_cells)
    real(rp),   intent(in), optional :: U_bnd(:,:)  ! (3, n_bnd_faces)
    ! Optionally reuse a pre-computed gradient (e.g. the solver's own grad_U) instead of rebuilding it.
    real(rp),   intent(in), optional, target :: grad_U_ext(:,:,:)  ! (3,3,n_cells[+n_ghost])
    ! RK-stage dt: if present, caps nu_sgs at NU_SGS_STAB_FRAC * V^(2/3) / dt_stab (see
    ! NU_SGS_STAB_FRAC) so an anisotropic-gradient spike at a small refined cell can't
    ! violate the explicit diffusive stability limit viscous_rhs_vector relies on.
    real(rp),   intent(in), optional :: dt_stab

    real(rp), allocatable, target :: grad_U_local(:,:,:)
    real(rp), pointer             :: grad_U(:,:,:)
    integer(ip) :: c, m, n, i, j, nc
    real(rp)    :: alpha(3,3), beta(3,3)
    real(rp)    :: alpha_sq, B_beta, nu_cap, dt_stab_local
    logical     :: has_dt_stab
    real(rp), allocatable :: cvol_local(:)

    if (present(grad_U_ext)) then
      ! Point directly at the caller's gradient; no copy needed.
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

    ! Precompute the filter-width tensor G per cell; rebuild if mesh size changed.
    ! Split into two IFs: .or. does not short-circuit in Fortran, so size(g_filt_arr,...)
    ! must never be evaluated in the same expression as the allocated() check.
    if (.not. allocated(g_filt_arr)) then
      call build_filter_tensor(msh)
    else if (size(g_filt_arr, 3) /= msh%n_cells + msh%n_ghost) then
      call build_filter_tensor(msh)
    end if

    ! GPU offload: flat-array OpenACC kernel with nc hoisted to avoid deep-copying the mesh.
    nc = msh%n_cells
    has_dt_stab = present(dt_stab)
    dt_stab_local = 1.0_rp
    if (has_dt_stab) dt_stab_local = dt_stab
    ! Always-present local volume array (host loop, cheap): avoids referencing an absent
    ! optional or an unbuilt SoA mirror inside the offloaded region below.
    allocate(cvol_local(nc))
    if (has_dt_stab) then
      do c = 1, nc
        cvol_local(c) = msh%cells(c)%volume
      end do
    else
      cvol_local = 0.0_rp
    end if
#ifdef _OPENACC
    !$acc data copyin(grad_U(:,:,1:nc), g_filt_arr(:,:,1:nc), cvol_local(1:nc)) &
    !$acc      copyout(nu_sgs(1:nc))
    !$acc parallel loop private(alpha, beta, alpha_sq, B_beta, nu_cap, i, j, m, n)
#else
    !$omp parallel do private(c, alpha, beta, alpha_sq, B_beta, nu_cap, i, j, m, n) &
    !$omp schedule(static)
#endif
    do c = 1, nc
      ! alpha_ij = d u_j / d x_i
      alpha = grad_U(:, :, c)

      ! Vreman (2004) Eq.(5) in tensor form: beta = alpha^T G alpha. G is symmetric, so is beta;
      ! B_beta below only reads the upper triangle (i<=j), so only that half is computed.
      do j = 1, 3
        do i = 1, j
          beta(i, j) = 0.0_rp
          do n = 1, 3
            do m = 1, 3
              beta(i, j) = beta(i, j) + g_filt_arr(m, n, c) * alpha(m, i) * alpha(n, j)
            end do
          end do
        end do
      end do

      ! alpha_ij * alpha_ij
      alpha_sq = alpha(1,1)**2 + alpha(1,2)**2 + alpha(1,3)**2 &
               + alpha(2,1)**2 + alpha(2,2)**2 + alpha(2,3)**2 &
               + alpha(3,1)**2 + alpha(3,2)**2 + alpha(3,3)**2

      ! B_beta (second invariant of beta)
      B_beta = beta(1,1)*beta(2,2) - beta(1,2)**2 &
             + beta(1,1)*beta(3,3) - beta(1,3)**2 &
             + beta(2,2)*beta(3,3) - beta(2,3)**2

      ! SGS viscosity: nu_t = C * sqrt(B_beta / alpha_sq) (anisotropic delta absorbed in beta)
      if (alpha_sq > SMALL .and. B_beta > 0.0_rp) then
        nu_sgs(c) = C_v * sqrt(B_beta / alpha_sq)
      else
        nu_sgs(c) = 0.0_rp
      end if

      if (has_dt_stab) then
        nu_cap = NU_SGS_STAB_FRAC * cvol_local(c)**(2.0_rp/3.0_rp) / max(dt_stab_local, SMALL)
        nu_sgs(c) = min(nu_sgs(c), nu_cap)
      end if
    end do
#ifdef _OPENACC
    !$acc end parallel loop
    !$acc end data
#else
    !$omp end parallel do
#endif

    deallocate(cvol_local)
    if (allocated(grad_U_local)) deallocate(grad_U_local)
    nullify(grad_U)
  end subroutine vreman_sgs_viscosity

  !> Build the per-cell filter-width tensor G; area weighting makes it invariant to face subdivision.
  subroutine build_filter_tensor(msh)
    type(Mesh), intent(inout) :: msh

    integer(ip) :: c, fi_wm, fj_wm, m, n
    real(rp)    :: g_raw(3,3), dx_wm(3), a_f, d2_iso, det_raw, scale_g, tr3

    if (allocated(g_filt_arr)) then
#ifdef _OPENACC
      !$acc exit data delete(g_filt_arr)
#endif
      deallocate(g_filt_arr)
    end if
    allocate(g_filt_arr(3, 3, msh%n_cells + msh%n_ghost), source=0.0_rp)

    do c = 1, msh%n_cells
      d2_iso = msh%cells(c)%volume**(2.0_rp/3.0_rp)
      g_raw = 0.0_rp
      if (allocated(msh%cell_faces)) then
        do fi_wm = msh%cells(c)%face_offset, &
                   msh%cells(c)%face_offset + msh%cells(c)%n_faces - 1
          fj_wm = msh%cell_faces(fi_wm)
          if (fj_wm < 1 .or. fj_wm > size(msh%faces)) cycle
          dx_wm = msh%faces(fj_wm)%centre - msh%cells(c)%centroid
          ! Area weighting: four coplanar sub-faces must contribute as the one face they replaced.
          a_f = msh%faces(fj_wm)%area
          do n = 1, 3
            do m = 1, 3
              g_raw(m, n) = g_raw(m, n) + a_f * dx_wm(m) * dx_wm(n)
            end do
          end do
        end do
      end if
      ! Normalise so det(G) = V^2; on an axis-aligned hex this recovers G = diag(h_m^2).
      det_raw = det3x3(g_raw)
      ! Scale-free degeneracy guard: area weighting changed the units of g_raw.
      tr3 = ((g_raw(1,1) + g_raw(2,2) + g_raw(3,3)) / 3.0_rp)**3
      if (tr3 > 0.0_rp .and. det_raw > 1.0e-12_rp * tr3) then
        scale_g = (msh%cells(c)%volume**2 / det_raw)**(1.0_rp/3.0_rp)
        g_filt_arr(:, :, c) = scale_g * g_raw
      else
        ! Degenerate cell (coplanar face offsets): isotropic width.
        g_filt_arr(:, :, c) = 0.0_rp
        do m = 1, 3
          g_filt_arr(m, m, c) = d2_iso
        end do
      end if
    end do

    call smooth_filter_tensor(msh)

    ! g_filt_arr is static once built (fixed mesh geometry); stage it once here so every
    ! !$acc kernel's own copyin(g_filt_arr(...)) finds it already present and skips the transfer.
#ifdef _OPENACC
    !$acc enter data copyin(g_filt_arr)
#endif
  end subroutine build_filter_tensor

  !> Blend G toward its area-weighted neighbour mean so nu_t ramps across a refinement step.
  subroutine smooth_filter_tensor(msh)
    type(Mesh), intent(inout) :: msh

    real(rp), allocatable :: g_nb(:,:,:), w_sum(:)
    integer(ip) :: pass, f, o, nb, c, ntot
    real(rp)    :: a_f

    if (G_SMOOTH_PASSES < 1) return
    ntot = msh%n_cells + msh%n_ghost
    allocate(g_nb(3, 3, msh%n_cells), w_sum(msh%n_cells))

    do pass = 1, G_SMOOTH_PASSES
      call exchange_g_halo(msh)
      g_nb  = 0.0_rp
      w_sum = 0.0_rp
      do f = 1, msh%n_faces
        o  = msh%faces(f)%owner
        nb = msh%faces(f)%neighbour
        if (nb < 1 .or. nb > ntot) cycle  ! boundary face: no neighbour to blend against
        a_f = msh%faces(f)%area
        if (o <= msh%n_cells) then
          g_nb(:, :, o) = g_nb(:, :, o) + a_f * g_filt_arr(:, :, nb)
          w_sum(o)      = w_sum(o) + a_f
        end if
        if (nb <= msh%n_cells) then
          g_nb(:, :, nb) = g_nb(:, :, nb) + a_f * g_filt_arr(:, :, o)
          w_sum(nb)      = w_sum(nb) + a_f
        end if
      end do
      do c = 1, msh%n_cells
        if (w_sum(c) > SMALL) &
          g_filt_arr(:, :, c) = (1.0_rp - G_SMOOTH_WEIGHT) * g_filt_arr(:, :, c) &
                              + G_SMOOTH_WEIGHT * g_nb(:, :, c) / w_sum(c)
      end do
    end do

    call exchange_g_halo(msh)
    deallocate(g_nb, w_sum)
  end subroutine smooth_filter_tensor

  !> Halo-exchange the six unique components of G in one fused call so smoothing sees off-rank neighbours.
  subroutine exchange_g_halo(msh)
    type(Mesh), intent(inout) :: msh

    integer(ip) :: m, n, ntot, k
    integer(ip), parameter :: PACK_M(6) = [1, 2, 3, 1, 1, 2]
    integer(ip), parameter :: PACK_N(6) = [1, 2, 3, 2, 3, 3]

    if (msh%n_ghost < 1) return
    ntot = msh%n_cells + msh%n_ghost
    if (.not. allocated(g_halo_buf)) then
      allocate(g_halo_buf(6, ntot))
    else if (size(g_halo_buf, 2) < ntot) then
      deallocate(g_halo_buf)
      allocate(g_halo_buf(6, ntot))
    end if

    do k = 1, 6
      g_halo_buf(k, 1:ntot) = g_filt_arr(PACK_M(k), PACK_N(k), 1:ntot)
    end do
    call msh%exchange_halos_scalars(g_halo_buf(:, 1:ntot))
    do k = 1, 6
      m = PACK_M(k); n = PACK_N(k)
      g_filt_arr(m, n, 1:ntot) = g_halo_buf(k, 1:ntot)
      g_filt_arr(n, m, 1:ntot) = g_halo_buf(k, 1:ntot)
    end do
  end subroutine exchange_g_halo

  !> Constant-coefficient Smagorinsky: nu_t = (Cs*Delta)^2*|S|. Not dynamic; Cs is fixed.
  subroutine smagorinsky_sgs(grad_U, msh, nu_sgs, cs_in, dt_stab)
    real(rp),   intent(in)    :: grad_U(:,:,:)  ! (3, 3, n_cells)
    type(Mesh), intent(inout) :: msh  ! inout: the filter-tensor build halo-exchanges G
    real(rp),   intent(inout) :: nu_sgs(:)
    real(rp),   intent(in), optional :: cs_in  ! override C_s (default 0.1)
    ! RK-stage dt: if present, caps nu_sgs the same way as vreman_sgs_viscosity (see
    ! NU_SGS_STAB_FRAC).
    real(rp),   intent(in), optional :: dt_stab

    integer(ip) :: c, i, j, nc
    real(rp)    :: S(3,3), S_norm, cs_val, delta_f, Sij_Sij, nu_cap, dt_stab_local
    logical     :: has_dt_stab

    cs_val = 0.1_rp
    if (present(cs_in)) cs_val = cs_in
    has_dt_stab = present(dt_stab)
    dt_stab_local = 1.0_rp
    if (has_dt_stab) dt_stab_local = dt_stab

    if (.not. allocated(g_filt_arr)) then
      call build_filter_tensor(msh)
    else if (size(g_filt_arr, 3) /= msh%n_cells + msh%n_ghost) then
      call build_filter_tensor(msh)
    end if

    ! GPU offload: same rationale as vreman_sgs_viscosity above, flat-array kernel only.
    nc = msh%n_cells
#ifdef _OPENACC
    !$acc data copyin(grad_U(:,:,1:nc), g_filt_arr(:,:,1:nc)) &
    !$acc      copyout(nu_sgs(1:nc))
    !$acc parallel loop private(S, S_norm, delta_f, Sij_Sij, nu_cap, i, j)
#else
    !$omp parallel do schedule(static) private(c, S, S_norm, delta_f, Sij_Sij, nu_cap, i, j)
#endif
    do c = 1, nc
      ! Strain rate: S_ij = 0.5 * (dU_j/dx_i + dU_i/dx_j)
      do j = 1, 3
        do i = 1, 3
          S(i, j) = 0.5_rp * (grad_U(i, j, c) + grad_U(j, i, c))
        end do
      end do

      ! Magnitude: |S| = sqrt(2 * S_ij * S_ij)
      Sij_Sij = 0.0_rp
      do j = 1, 3
        do i = 1, 3
          Sij_Sij = Sij_Sij + S(i, j) * S(i, j)
        end do
      end do
      S_norm = sqrt(2.0_rp * Sij_Sij)

      ! Delta from the smoothed filter tensor: det(G)^(1/6) equals V^(1/3) on a uniform hex.
      delta_f = max(max(det3x3(g_filt_arr(:, :, c)), 0.0_rp)**(1.0_rp/6.0_rp), SMALL)

      ! Smagorinsky: nu_t = (C_s * Delta)^2 * |S|
      nu_sgs(c) = (cs_val * delta_f)**2 * S_norm

      if (has_dt_stab) then
        nu_cap = NU_SGS_STAB_FRAC * delta_f**2 / max(dt_stab_local, SMALL)
        nu_sgs(c) = min(nu_sgs(c), nu_cap)
      end if
    end do
#ifdef _OPENACC
    !$acc end parallel loop
    !$acc end data
#else
    !$omp end parallel do
#endif

  end subroutine smagorinsky_sgs

  !> Grid filter width Delta(c) = det(G)^(1/6) at every owned cell; shared with mod_dynamic_smagorinsky.
  subroutine vreman_grid_delta(msh, delta_out)
    type(Mesh), intent(inout) :: msh
    real(rp),   intent(out)   :: delta_out(:)  ! (n_cells)

    integer(ip) :: c

    if (.not. allocated(g_filt_arr)) then
      call build_filter_tensor(msh)
    else if (size(g_filt_arr, 3) /= msh%n_cells + msh%n_ghost) then
      call build_filter_tensor(msh)
    end if

    do c = 1, msh%n_cells
      delta_out(c) = max(max(det3x3(g_filt_arr(:, :, c)), 0.0_rp)**(1.0_rp/6.0_rp), SMALL)
    end do
  end subroutine vreman_grid_delta

  !> Read back the filter tensor of one cell; for diagnostics and tests.
  subroutine vreman_filter_tensor(c, g_out)
    integer(ip), intent(in)  :: c
    real(rp),    intent(out) :: g_out(3,3)

    ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
    ! must never be evaluated in the same expression as the allocated() check.
    g_out = 0.0_rp
    if (allocated(g_filt_arr)) then
      if (c >= 1 .and. c <= size(g_filt_arr, 3)) g_out = g_filt_arr(:, :, c)
    end if
  end subroutine vreman_filter_tensor

  !> Free module-level persistent allocations (call on solver teardown).
  subroutine vreman_free()
    if (allocated(g_filt_arr)) deallocate(g_filt_arr)
    if (allocated(g_halo_buf)) deallocate(g_halo_buf)
  end subroutine vreman_free

end module mod_vreman
