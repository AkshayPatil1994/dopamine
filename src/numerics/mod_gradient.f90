! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Gradient reconstruction: Green-Gauss and weighted least-squares.
module mod_gradient
  use mod_kinds
  use mod_parameters
  use mod_mesh
  implicit none
  private

  public :: green_gauss_gradient_scalar
  public :: green_gauss_gradient_vector
  public :: ls_gradient_scalar
  public :: ls_gradient_vector
  public :: face_interpolate_scalar
  public :: face_interpolate_vector
  public :: face_gradient_normal_scalar
  public :: ensure_buf_1d
  public :: ensure_buf_2d3
  public :: ensure_buf_n3
  public :: ensure_buf_3d33
  public :: ensure_buf_3n3
  public :: ensure_fcp_identity

  ! Persistent scratch buffer grown once instead of allocate/deallocate every RK substage.
  integer(ip), allocatable, target, save :: fcp_identity_buf(:)     ! shared by both Green-Gauss routines

  ! Stand-ins bound to absent optional arguments (phi_bnd, grad_*_prev): an absent optional can't
  ! be named inside an !$acc region, so the kernels always see a real, device-resident array.
  ! They are only ever read behind the matching has_* guard.
  real(rp), target, save :: grad_dummy_s(1)       = 0.0_rp
  real(rp), target, save :: grad_dummy_v(3, 1)    = 0.0_rp
  real(rp), target, save :: grad_dummy_t(3, 3, 1) = 0.0_rp
  !$acc declare create(grad_dummy_s, grad_dummy_v, grad_dummy_t)

contains

  !> Green-Gauss gradient: grad_phi_c = (1/V_c) sum_f phi_f * S_f.
  !> Pass grad_phi_prev to apply skewness correction: phi_f += grad_mean . eps.
  subroutine green_gauss_gradient_scalar(phi, msh, grad_phi, phi_bnd, grad_phi_prev)
    real(rp),    intent(in)    :: phi(:)
    type(Mesh),  intent(in), target :: msh
    real(rp),    intent(inout) :: grad_phi(:,:)  ! (3, n_cells)
    real(rp),    intent(in), optional, target :: phi_bnd(:)  ! (n_bnd_faces)
    real(rp),    intent(in), optional, target :: grad_phi_prev(:,:)  ! (3, n_cells)

    integer(ip) :: f, o, nb, c
    real(rp)    :: phi_f, phi_f_plain, sf(3), gc, eps(3), grad_mean(3), w_blend
    integer(ip) :: col, idx, n_loops, lo, hi, nc, n_prev, nf, n_int, ng
    logical     :: use_color, has_bnd, has_prev
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:)
    real(rp),    pointer :: sf_soa(:,:), gc_soa(:), eps_soa(:,:), cvol_soa(:), skew_soa(:)
    real(rp),    pointer :: pbnd(:), gprev(:,:)

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'green_gauss_gradient_scalar: msh%soa_ready is false; call mesh%build_soa_mirror first'

    nc = msh%n_cells
    nf = msh%n_faces
    n_int = msh%n_int_faces
    ng = size(grad_phi, 2)
    fo_soa   => msh%face_owner_soa
    fn_soa   => msh%face_neighbour_soa
    sf_soa   => msh%face_sf_soa
    gc_soa   => msh%face_gc_soa
    eps_soa  => msh%face_eps_soa
    cvol_soa => msh%cell_volume_soa
    skew_soa => msh%face_skew_soa

    has_bnd = present(phi_bnd)
    pbnd => grad_dummy_s
    if (has_bnd) pbnd => phi_bnd
    has_prev = present(grad_phi_prev)
    gprev => grad_dummy_v
    if (has_prev) gprev => grad_phi_prev
    n_prev = merge(size(gprev, 2), 1_ip, has_prev)

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets) .and. msh%soa_ready)
    n_loops = merge(msh%n_colors, 1_ip, use_color)

    ! fcp always a valid size-nf array (identity permutation when uncolored) so the single
    ! data region below can copy it unconditionally instead of reopening per color.
    if (use_color) then
      fcp => msh%face_color_perm
    else
      call ensure_fcp_identity(fcp_identity_buf, nf)
      fcp => fcp_identity_buf(1:nf)
    end if

    ! The output is zeroed (ghost slots too, filled later by the halo exchange) and scaled on
    ! the device, so the routine works the same whether grad_phi is resident or copied.
    !$acc data copyin(phi, pbnd, gprev, fcp(1:nf)) copyout(grad_phi) if(use_color)
    !$acc parallel loop collapse(2) if(use_color) default(present)
    do c = 1, ng
      do f = 1, 3
        grad_phi(f, c) = 0.0_rp
      end do
    end do
    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) default(present) private(f, o, nb, sf, gc, eps, grad_mean, phi_f, &
      !$acc      phi_f_plain, w_blend)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, sf, gc, eps, grad_mean, phi_f, phi_f_plain, w_blend)
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

        if (nb >= 1) then
          if (has_prev .and. nb <= n_prev) then
            gc        = gc_soa(f)
            eps       = eps_soa(:, f)  ! precomputed at mesh setup: static interior-face geometry
            ! Ghost gradients are halo-exchanged before each sweep, so both sides of a processor
            ! face see the same mean gradient (result independent of the partitioning).
            grad_mean = gc * gprev(:, o) + (1.0_rp - gc) * gprev(:, nb)
            phi_f_plain = gc * phi(o) + (1.0_rp - gc) * phi(nb)
            phi_f       = phi_f_plain + dot_product(grad_mean, eps)
            ! Fade the eps-correction out toward refinement interfaces; see refine_blend_weight.
            if (f <= n_int) then
              w_blend = refine_blend_weight(skew_soa(f))
              if (w_blend < 1.0_rp) phi_f = phi_f_plain + w_blend * (phi_f - phi_f_plain)
            end if
          else
            phi_f = gc_soa(f) * phi(o) + (1.0_rp - gc_soa(f)) * phi(nb)
          end if

          if (o  >= 1 .and. o  <= nc) grad_phi(:, o)  = grad_phi(:, o)  + phi_f * sf
          if (nb >= 1 .and. nb <= nc) grad_phi(:, nb) = grad_phi(:, nb) - phi_f * sf
        else
          if (has_bnd) then
            phi_f = pbnd(f - n_int)
          else
            phi_f = phi(o)
          end if

          if (o >= 1 .and. o <= nc) &
            grad_phi(:, o) = grad_phi(:, o) + phi_f * sf
        end if
      end do
#ifdef _OPENACC
      !$acc end parallel loop
#else
      !$omp end parallel do
#endif
    end do

    !$acc parallel loop if(use_color) default(present)
    do o = 1, nc
      if (cvol_soa(o) > SMALL) &
        grad_phi(:, o) = grad_phi(:, o) * (1.0_rp / cvol_soa(o))
    end do
    !$acc end data
    nullify(fo_soa, fn_soa, fcp, sf_soa, gc_soa, eps_soa, cvol_soa, skew_soa, pbnd, gprev)
  end subroutine green_gauss_gradient_scalar

  !> Green-Gauss gradient of a vector field: grad_U(i,j,c) = dU_j/dx_i.
  !> Pass grad_U_prev to apply skewness correction on each sweep.
  subroutine green_gauss_gradient_vector(U, msh, grad_U, U_bnd, grad_U_prev)
    real(rp),   intent(in)    :: U(:,:)
    type(Mesh), intent(in), target :: msh
    real(rp),   intent(inout) :: grad_U(:,:,:)  ! (3, 3, n_cells)
    real(rp),   intent(in), optional, target :: U_bnd(:,:)  ! (3, n_bnd_faces)
    real(rp),   intent(in), optional, target :: grad_U_prev(:,:,:)  ! (3, 3, n_cells)

    integer(ip) :: f, o, nb, k, c
    real(rp)    :: Uf(3), Uf_plain(3), sf(3), gc, eps(3), grad_mean(3), w_blend
    integer(ip) :: col, idx, n_loops, lo, hi, nc, n_prev, nf, n_int, ng
    logical     :: use_color, has_bnd, has_prev
    integer(ip), pointer :: fo_soa(:), fn_soa(:), fcp(:)
    real(rp),    pointer :: sf_soa(:,:), gc_soa(:), eps_soa(:,:), cvol_soa(:), skew_soa(:)
    real(rp),    pointer :: ubnd(:,:), gprev(:,:,:)

    if (.not. msh%soa_ready .and. msh%n_faces > 0) &
      error stop 'green_gauss_gradient_vector: msh%soa_ready is false; call mesh%build_soa_mirror first'

    nc = msh%n_cells
    nf = msh%n_faces
    n_int = msh%n_int_faces
    ng = size(grad_U, 3)
    fo_soa   => msh%face_owner_soa
    fn_soa   => msh%face_neighbour_soa
    sf_soa   => msh%face_sf_soa
    gc_soa   => msh%face_gc_soa
    eps_soa  => msh%face_eps_soa
    cvol_soa => msh%cell_volume_soa
    skew_soa => msh%face_skew_soa

    has_bnd = present(U_bnd)
    ubnd => grad_dummy_v
    if (has_bnd) ubnd => U_bnd
    has_prev = present(grad_U_prev)
    gprev => grad_dummy_t
    if (has_prev) gprev => grad_U_prev
    n_prev = merge(size(gprev, 3), 1_ip, has_prev)

    use_color = (msh%n_colors > 0 .and. allocated(msh%face_color_perm) &
                 .and. allocated(msh%color_offsets) .and. msh%soa_ready)
    n_loops = merge(msh%n_colors, 1_ip, use_color)

    ! fcp always a valid size-nf array (identity permutation when uncolored) so the single
    ! data region below can copy it unconditionally instead of reopening per color.
    if (use_color) then
      fcp => msh%face_color_perm
    else
      call ensure_fcp_identity(fcp_identity_buf, nf)
      fcp => fcp_identity_buf(1:nf)
    end if

    ! Zeroed and scaled on the device, see green_gauss_gradient_scalar.
    !$acc data copyin(U, ubnd, gprev, fcp(1:nf)) copyout(grad_U) if(use_color)
    !$acc parallel loop collapse(3) if(use_color) default(present)
    do c = 1, ng
      do k = 1, 3
        do f = 1, 3
          grad_U(f, k, c) = 0.0_rp
        end do
      end do
    end do
    do col = 1, n_loops
      if (use_color) then
        lo = msh%color_offsets(col)
        hi = msh%color_offsets(col+1) - 1_ip
      else
        lo = 1_ip
        hi = nf
      end if
#ifdef _OPENACC
      !$acc parallel loop if(use_color) default(present) private(f, o, nb, k, Uf, Uf_plain, sf, gc, eps, &
      !$acc      grad_mean, w_blend)
#else
      !$omp parallel do if(use_color) default(shared) &
      !$omp& private(f, idx, o, nb, k, Uf, Uf_plain, sf, gc, eps, grad_mean, w_blend)
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

        if (nb >= 1) then
          if (has_prev .and. nb <= n_prev) then
            gc   = gc_soa(f)
            eps  = eps_soa(:, f)  ! precomputed at mesh setup: static interior-face geometry
            do k = 1, 3
              grad_mean = gc * gprev(:, k, o) + (1.0_rp - gc) * gprev(:, k, nb)
              Uf(k) = gc * U(k, o) + (1.0_rp - gc) * U(k, nb) + dot_product(grad_mean, eps)
            end do
            ! Fade the eps-correction out toward refinement interfaces; see refine_blend_weight.
            if (f <= n_int) then
              w_blend = refine_blend_weight(skew_soa(f))
              if (w_blend < 1.0_rp) then
                Uf_plain = gc * U(:, o) + (1.0_rp - gc) * U(:, nb)
                Uf = Uf_plain + w_blend * (Uf - Uf_plain)
              end if
            end if
          else
            Uf = gc_soa(f) * U(:, o) + (1.0_rp - gc_soa(f)) * U(:, nb)
          end if
        else
          if (has_bnd) then
            Uf = ubnd(:, f - n_int)
          else
            Uf = U(:, o)
          end if
        end if

        do k = 1, 3
          if (o  >= 1 .and. o  <= nc) grad_U(:, k, o)  = grad_U(:, k, o)  + Uf(k) * sf
          if (nb >= 1 .and. nb <= nc) grad_U(:, k, nb) = grad_U(:, k, nb) - Uf(k) * sf
        end do
      end do
#ifdef _OPENACC
      !$acc end parallel loop
#else
      !$omp end parallel do
#endif
    end do

    !$acc parallel loop if(use_color) default(present)
    do o = 1, nc
      if (cvol_soa(o) > SMALL) &
        grad_U(:, :, o) = grad_U(:, :, o) * (1.0_rp / cvol_soa(o))
    end do
    !$acc end data
    nullify(fo_soa, fn_soa, fcp, sf_soa, gc_soa, eps_soa, cvol_soa, skew_soa, ubnd, gprev)
  end subroutine green_gauss_gradient_vector

  !> Weighted least-squares gradient using precomputed stencil.
  subroutine ls_gradient_scalar(phi, msh, grad_phi, phi_bnd)
    real(rp),   intent(in)    :: phi(:)
    type(Mesh), intent(in), target :: msh
    real(rp),   intent(inout) :: grad_phi(:,:)  ! (3, n_cells)
    real(rp),   intent(in), optional, target :: phi_bnd(:)  ! (n_bnd_faces)

    integer(ip) :: c, k, nb, f, bnd_idx
    integer(ip) :: nvals, nc, n_int_faces, n_bnd, n_stencil
    real(rp)    :: dphi, phi_nb
#ifdef _OPENACC
    real(rp)    :: gx, gy, gz  ! GPU-path-only per-cell reduction accumulators, see below
    integer(ip) :: ng
#endif
    logical     :: has_bnd
    integer(ip), pointer :: stencil_offset(:), stencil(:), stencil_face(:)
    real(rp),    pointer :: ls_w(:,:), pbnd(:)

    if (.not. allocated(msh%ls_stencil_face)) &
      error stop 'ls_gradient_scalar: ls_stencil_face not allocated; call mesh%build_ls_stencil first'

    nvals = size(phi)
    nc = msh%n_cells
    n_int_faces = msh%n_int_faces
    stencil_offset => msh%ls_stencil_offset
    stencil        => msh%ls_stencil
    stencil_face   => msh%ls_stencil_face
    ls_w           => msh%ls_weights
    n_stencil = size(stencil)

    ! An absent optional can't be named inside an !$acc region: bind a resident stand-in.
    has_bnd = present(phi_bnd)
    pbnd => grad_dummy_s
    if (has_bnd) pbnd => phi_bnd
    n_bnd = size(pbnd)

#ifdef _OPENACC
    ! GPU path only: spreads the stencil walk across a gang's vector threads with a manual reduction (nsys: 97.5% of kernel time otherwise, one thread crawling the whole list).
    ! Owned cells are fully overwritten below and ghost slots zeroed (filled later by the halo
    ! exchange), so copyout -- not copy -- avoids uploading a value the kernel never uses.
    ng = size(grad_phi, 2)
    !$acc data copyin(phi(1:nvals), pbnd, stencil_offset(1:nc+1), &
    !$acc              stencil(1:n_stencil), stencil_face(1:n_stencil), ls_w(:,1:n_stencil)) &
    !$acc      copyout(grad_phi)
    !$acc parallel loop default(present)
    do c = nc + 1, ng
      grad_phi(1, c) = 0.0_rp
      grad_phi(2, c) = 0.0_rp
      grad_phi(3, c) = 0.0_rp
    end do
    !$acc parallel loop gang private(gx, gy, gz)
    do c = 1, nc
      gx = 0.0_rp
      gy = 0.0_rp
      gz = 0.0_rp
      !$acc loop vector reduction(+:gx,gy,gz) private(k, nb, f, bnd_idx, dphi, phi_nb)
      do k = stencil_offset(c), stencil_offset(c+1) - 1
        nb = stencil(k)
        if (nb < 1) then
          f = stencil_face(k)
          if (f > 0) then
            bnd_idx = f - n_int_faces
          else
            bnd_idx = -1
          end if
          if (has_bnd .and. bnd_idx >= 1 .and. bnd_idx <= n_bnd) then
            phi_nb = pbnd(bnd_idx)
          else
            phi_nb = phi(c)
          end if
          dphi = phi_nb - phi(c)
        else if (nb > nvals) then
          dphi = 0.0_rp  ! out-of-range stencil entry: contributes nothing (replaces the CPU path's 'cycle')
        else
          phi_nb = phi(nb)
          dphi = phi_nb - phi(c)
        end if
        gx = gx + ls_w(1, k) * dphi
        gy = gy + ls_w(2, k) * dphi
        gz = gz + ls_w(3, k) * dphi
      end do
      grad_phi(1, c) = gx
      grad_phi(2, c) = gy
      grad_phi(3, c) = gz
    end do
    !$acc end parallel loop
    !$acc end data
#else
    grad_phi = 0.0_rp
    !$omp parallel do private(c, k, nb, f, bnd_idx, dphi, phi_nb) schedule(static)
    do c = 1, nc
      do k = stencil_offset(c), stencil_offset(c+1) - 1
        nb   = stencil(k)
        if (nb < 1) then
          f = stencil_face(k)
          if (f > 0) then
            bnd_idx = f - n_int_faces
          else
            bnd_idx = -1
          end if
          if (has_bnd .and. bnd_idx >= 1 .and. bnd_idx <= n_bnd) then
            phi_nb = pbnd(bnd_idx)
          else
            phi_nb = phi(c)
          end if
        else if (nb > nvals) then
          cycle
        else
          phi_nb = phi(nb)
        end if
        dphi = phi_nb - phi(c)
        grad_phi(:, c) = grad_phi(:, c) + ls_w(:, k) * dphi
      end do
    end do
    !$omp end parallel do
#endif

    nullify(stencil_offset, stencil, stencil_face, ls_w, pbnd)
  end subroutine ls_gradient_scalar

  !> Weighted least-squares gradient for each component of a vector field.
  subroutine ls_gradient_vector(U, msh, grad_U, U_bnd)
    real(rp),   intent(in)    :: U(:,:)
    type(Mesh), intent(in), target :: msh
    real(rp),   intent(inout) :: grad_U(:,:,:)  ! (3, 3, n_cells)
    real(rp),   intent(in), optional, target :: U_bnd(:,:)  ! (3, n_bnd_faces)

    integer(ip) :: c, k, nb, j, f, bnd_idx
    integer(ip) :: nvals, nc, n_int_faces, n_bnd, n_stencil
    real(rp)    :: dU(3), U_nb(3)
#ifdef _OPENACC
    real(rp)    :: g11, g21, g31, g12, g22, g32, g13, g23, g33  ! GPU-path-only reduction accumulators
    integer(ip) :: ng
#endif
    logical     :: has_bnd
    integer(ip), pointer :: stencil_offset(:), stencil(:), stencil_face(:)
    real(rp),    pointer :: ls_w(:,:), ubnd(:,:)

    if (.not. allocated(msh%ls_stencil_face)) &
      error stop 'ls_gradient_vector: ls_stencil_face not allocated; call mesh%build_ls_stencil first'

    nvals = size(U, 2)
    nc = msh%n_cells
    n_int_faces = msh%n_int_faces
    stencil_offset => msh%ls_stencil_offset
    stencil        => msh%ls_stencil
    stencil_face   => msh%ls_stencil_face
    ls_w           => msh%ls_weights
    n_stencil = size(stencil)

    ! An absent optional can't be named inside an !$acc region: bind a resident stand-in.
    has_bnd = present(U_bnd)
    ubnd => grad_dummy_v
    if (has_bnd) ubnd => U_bnd
    n_bnd = size(ubnd, 2)

#ifdef _OPENACC
    ! GPU path only, same rationale as ls_gradient_scalar above (97.5% of measured GPU kernel time); CPU loop below is untouched.
    ! Unrolled into 9 scalar accumulators, not a g(3,3) array reduction: the array form hit an nvfortran 26.5 shared-memory codegen bug.
    ! Owned cells are fully overwritten below and ghost slots zeroed (filled later by the halo
    ! exchange), so copyout -- not copy -- avoids uploading a value the kernel never uses.
    ng = size(grad_U, 3)
    !$acc data copyin(U(:,1:nvals), ubnd, stencil_offset(1:nc+1), &
    !$acc              stencil(1:n_stencil), stencil_face(1:n_stencil), ls_w(:,1:n_stencil)) &
    !$acc      copyout(grad_U)
    !$acc parallel loop collapse(2) default(present)
    do c = nc + 1, ng
      do j = 1, 3
        grad_U(1, j, c) = 0.0_rp
        grad_U(2, j, c) = 0.0_rp
        grad_U(3, j, c) = 0.0_rp
      end do
    end do
    !$acc parallel loop gang private(g11,g21,g31,g12,g22,g32,g13,g23,g33)
    do c = 1, nc
      g11 = 0.0_rp; g21 = 0.0_rp; g31 = 0.0_rp
      g12 = 0.0_rp; g22 = 0.0_rp; g32 = 0.0_rp
      g13 = 0.0_rp; g23 = 0.0_rp; g33 = 0.0_rp
      !$acc loop vector reduction(+:g11,g21,g31,g12,g22,g32,g13,g23,g33) &
      !$acc      private(k, nb, f, bnd_idx, dU, U_nb)
      do k = stencil_offset(c), stencil_offset(c+1) - 1
        nb = stencil(k)
        if (nb < 1) then
          f = stencil_face(k)
          if (f > 0) then
            bnd_idx = f - n_int_faces
          else
            bnd_idx = -1
          end if
          if (has_bnd .and. bnd_idx >= 1 .and. bnd_idx <= n_bnd) then
            U_nb = ubnd(:, bnd_idx)
          else
            U_nb = U(:, c)
          end if
          dU = U_nb - U(:, c)
        else if (nb > nvals) then
          dU = 0.0_rp  ! out-of-range stencil entry: contributes nothing (replaces the CPU path's 'cycle')
        else
          U_nb = U(:, nb)
          dU = U_nb - U(:, c)
        end if
        g11 = g11 + ls_w(1, k) * dU(1)
        g21 = g21 + ls_w(2, k) * dU(1)
        g31 = g31 + ls_w(3, k) * dU(1)
        g12 = g12 + ls_w(1, k) * dU(2)
        g22 = g22 + ls_w(2, k) * dU(2)
        g32 = g32 + ls_w(3, k) * dU(2)
        g13 = g13 + ls_w(1, k) * dU(3)
        g23 = g23 + ls_w(2, k) * dU(3)
        g33 = g33 + ls_w(3, k) * dU(3)
      end do
      grad_U(1, 1, c) = g11
      grad_U(2, 1, c) = g21
      grad_U(3, 1, c) = g31
      grad_U(1, 2, c) = g12
      grad_U(2, 2, c) = g22
      grad_U(3, 2, c) = g32
      grad_U(1, 3, c) = g13
      grad_U(2, 3, c) = g23
      grad_U(3, 3, c) = g33
    end do
    !$acc end parallel loop
    !$acc end data
#else
    grad_U = 0.0_rp
    !$omp parallel do private(c, k, nb, j, f, bnd_idx, dU, U_nb) schedule(static)
    do c = 1, nc
      do k = stencil_offset(c), stencil_offset(c+1) - 1
        nb = stencil(k)
        if (nb < 1) then
          f = stencil_face(k)
          if (f > 0) then
            bnd_idx = f - n_int_faces
          else
            bnd_idx = -1
          end if
          if (has_bnd .and. bnd_idx >= 1 .and. bnd_idx <= n_bnd) then
            U_nb = ubnd(:, bnd_idx)
          else
            U_nb = U(:, c)
          end if
        else if (nb > nvals) then
          cycle
        else
          U_nb = U(:, nb)
        end if
        dU = U_nb - U(:, c)
        do j = 1, 3
          grad_U(:, j, c) = grad_U(:, j, c) + ls_w(:, k) * dU(j)
        end do
      end do
    end do
    !$omp end parallel do
#endif

    nullify(stencil_offset, stencil, stencil_face, ls_w, ubnd)
  end subroutine ls_gradient_vector

  !> Face interpolation of a scalar with optional gradient correction (O(h^2) on skewed grids).
  pure real(rp) function face_interpolate_scalar(phi, f, msh, grad_phi) result(phi_f)
    real(rp),   intent(in) :: phi(:)
    integer(ip),intent(in) :: f
    type(Mesh), intent(in) :: msh
    real(rp),   intent(in), optional :: grad_phi(:,:)  ! (3, n_cells)

    integer(ip) :: o, nb, f_bnd
    real(rp)    :: gc, phi_linear
    real(rp)    :: d_on(3), epsilon(3), grad_mean(3)

    f_bnd = 0
    o  = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour

    ! Periodic boundary face: resolve partner via periodic_neighbor or faces%neighbour.
    if (nb < 1 .and. f > msh%n_int_faces) then
      f_bnd = f - msh%n_int_faces
      if (allocated(msh%periodic_neighbor)) then
        if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) then
          if (msh%periodic_neighbor(f_bnd) >= 1) nb = msh%periodic_neighbor(f_bnd)
        end if
      end if
      if (nb < 1) nb = msh%faces(f)%neighbour
    end if

    if (nb >= 1 .and. nb <= size(phi)) then
      ! Linear interpolation at face center. faces(f)%gc is hard-coded to 1 on boundary/periodic
      ! faces (see mesh_compute_geometry), so periodic faces must use the true periodic-pair weight.
      if (f > msh%n_int_faces) then
        gc = msh%periodic_gc_soa(f_bnd)
      else
        gc = msh%faces(f)%gc
      end if
      phi_linear = gc * phi(o) + (1.0_rp - gc) * phi(nb)

      ! O(h^2) skewness correction: phi_f += grad_mean . (face_centre - gc-interp point).
      if (present(grad_phi)) then
        if (f <= msh%n_int_faces) then
          ! Genuine interior face: nb above is msh%faces(f)%neighbour, matching the cache.
          epsilon = msh%faces(f)%eps
        else
          ! Periodic face: nb was resolved via periodic_neighbor above, so the cache
          ! (built from msh%faces(f)%neighbour only) does not apply; recompute directly.
          d_on = msh%cells(nb)%centroid - msh%cells(o)%centroid
          epsilon = msh%faces(f)%centre - (msh%cells(o)%centroid + (1.0_rp - gc) * d_on)
        end if
        grad_mean = gc * grad_phi(:, o) + (1.0_rp - gc) * grad_phi(:, nb)
        phi_f = phi_linear + dot_product(grad_mean, epsilon)
      else
        phi_f = phi_linear
      end if
    else
      phi_f = phi(o)
    end if
  end function face_interpolate_scalar

  !> Face interpolation of a vector with optional gradient correction (O(h^2) on skewed grids).
  pure function face_interpolate_vector(U, f, msh, grad_U) result(Uf)
    real(rp),   intent(in) :: U(:,:)
    integer(ip),intent(in) :: f
    type(Mesh), intent(in) :: msh
    real(rp),   intent(in), optional :: grad_U(:,:,:)  ! (3, 3, n_cells)

    real(rp)    :: Uf(3)
    integer(ip) :: o, nb, j, f_bnd
    real(rp)    :: gc
    real(rp)    :: d_on(3), epsilon(3), grad_mean(3)

    f_bnd = 0

    o  = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour

    ! Periodic boundary face: resolve partner via periodic_neighbor or faces%neighbour.
    if (nb < 1 .and. f > msh%n_int_faces) then
      f_bnd = f - msh%n_int_faces
      if (allocated(msh%periodic_neighbor)) then
        if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) then
          if (msh%periodic_neighbor(f_bnd) >= 1) nb = msh%periodic_neighbor(f_bnd)
        end if
      end if
      if (nb < 1) nb = msh%faces(f)%neighbour
    end if

    if (nb >= 1 .and. nb <= size(U,2)) then
      ! faces(f)%gc is hard-coded to 1 on boundary/periodic faces (see mesh_compute_geometry),
      ! so periodic faces must use the true periodic-pair weight.
      if (f > msh%n_int_faces) then
        gc = msh%periodic_gc_soa(f_bnd)
      else
        gc = msh%faces(f)%gc
      end if
      Uf = gc * U(:,o) + (1.0_rp - gc) * U(:,nb)

      ! O(h^2) skewness correction (Kim & Choi Eq. 26): Uf_j += grad_mean_j . epsilon.
      if (present(grad_U)) then
        if (f <= msh%n_int_faces) then
          ! Genuine interior face: nb above is msh%faces(f)%neighbour, matching the cache.
          epsilon = msh%faces(f)%eps
        else
          ! Periodic face: nb was resolved via periodic_neighbor above, so the cache
          ! (built from msh%faces(f)%neighbour only) does not apply; recompute directly.
          d_on = msh%cells(nb)%centroid - msh%cells(o)%centroid
          epsilon = msh%faces(f)%centre - (msh%cells(o)%centroid + (1.0_rp - gc) * d_on)
        end if
        do j = 1, 3
          grad_mean = gc * grad_U(:, j, o) + (1.0_rp - gc) * grad_U(:, j, nb)
          Uf(j) = Uf(j) + dot_product(grad_mean, epsilon)
        end do
      end if
    else
      Uf = U(:, o)
    end if
  end function face_interpolate_vector

  !> Face-normal gradient of a scalar: d(phi)/dn at face f.
  !> Uses compact stencil: (phi_nb - phi_o) / d_cf
  pure real(rp) function face_gradient_normal_scalar(phi, f, msh, &
                                                       grad_phi_o, grad_phi_nb)
    real(rp),    intent(in) :: phi(:)
    integer(ip), intent(in) :: f
    type(Mesh),  intent(in) :: msh
    real(rp),    intent(in) :: grad_phi_o(3), grad_phi_nb(3)  ! cell gradients

    integer(ip) :: o, nb
    real(rp)    :: d_on(3), phi_corr, grad_mean(3), d_on_len, d_dot_n

    o  = msh%faces(f)%owner
    nb = msh%faces(f)%neighbour

    if (nb >= 1) then
      d_on = msh%cells(nb)%centroid - msh%cells(o)%centroid
      d_on_len = msh%faces(f)%d_cf  ! == sqrt(sum(d_on**2)), already cached at mesh setup
      ! Floor denominator to 10% of |d_on| to cap amplification at 10x on highly
      ! non-orthogonal faces (theta > ~84 deg); consistent over-relaxed k-vector.
      d_dot_n = max(dot_product(d_on, msh%faces(f)%normal), 0.1_rp * d_on_len) + SMALL
      face_gradient_normal_scalar = (phi(nb) - phi(o)) / d_dot_n
      grad_mean = 0.5_rp * (grad_phi_o + grad_phi_nb)
      phi_corr  = dot_product(grad_mean, msh%faces(f)%normal - d_on / d_dot_n)
      face_gradient_normal_scalar = face_gradient_normal_scalar + phi_corr
    else
      face_gradient_normal_scalar = dot_product(grad_phi_o, msh%faces(f)%normal)
    end if
  end function face_gradient_normal_scalar

  !> Grow buf to at least n (never shrinks); mesh topology is static so this allocates once.
  subroutine ensure_buf_1d(buf, n)
    real(rp), allocatable, intent(inout) :: buf(:)
    integer(ip), intent(in) :: n
    if (allocated(buf)) then
      if (size(buf) >= n) return
      deallocate(buf)
    end if
    allocate(buf(n))
  end subroutine ensure_buf_1d

  !> Grow buf to shape (3, >=n) (never shrinks).
  subroutine ensure_buf_2d3(buf, n)
    real(rp), allocatable, intent(inout) :: buf(:,:)
    integer(ip), intent(in) :: n
    if (allocated(buf)) then
      if (size(buf, 2) >= n) return
      deallocate(buf)
    end if
    allocate(buf(3, n))
  end subroutine ensure_buf_2d3

  !> Grow buf to shape (3, 3, >=n) (never shrinks).
  subroutine ensure_buf_3d33(buf, n)
    real(rp), allocatable, intent(inout) :: buf(:,:,:)
    integer(ip), intent(in) :: n
    if (allocated(buf)) then
      if (size(buf, 3) >= n) return
      deallocate(buf)
    end if
    allocate(buf(3, 3, n))
  end subroutine ensure_buf_3d33

  !> Grow buf to shape (>=n, 3) (never shrinks).
  subroutine ensure_buf_n3(buf, n)
    real(rp), allocatable, intent(inout) :: buf(:,:)
    integer(ip), intent(in) :: n
    if (allocated(buf)) then
      if (size(buf, 1) >= n) return
      deallocate(buf)
    end if
    allocate(buf(n, 3))
  end subroutine ensure_buf_n3

  !> Grow buf to shape (3, >=n, 3) (never shrinks).
  subroutine ensure_buf_3n3(buf, n)
    real(rp), allocatable, intent(inout) :: buf(:,:,:)
    integer(ip), intent(in) :: n
    if (allocated(buf)) then
      if (size(buf, 2) >= n) return
      deallocate(buf)
    end if
    allocate(buf(3, n, 3))
  end subroutine ensure_buf_3n3

  !> Grow buf to an identity permutation of length >=n (never shrinks).
  subroutine ensure_fcp_identity(buf, n)
    integer(ip), allocatable, intent(inout) :: buf(:)
    integer(ip), intent(in) :: n
    integer(ip) :: i
    if (allocated(buf)) then
      if (size(buf) >= n) return
      deallocate(buf)
    end if
    allocate(buf(n))
    do i = 1, n
      buf(i) = i
    end do
  end subroutine ensure_fcp_identity

end module mod_gradient
