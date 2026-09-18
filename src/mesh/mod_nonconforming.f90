! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Detects T-junction hanging-node boundary-face pairs and converts them to internal faces.
module mod_nonconforming
  use mod_kinds
  use mod_parameters
  use mod_mesh
  implicit none
  private
  public :: enforce_face_conformity

  ! Tolerances - all relative to characteristic length sqrt(area_A) or area_A.
  real(rp), parameter :: NC_NRML_TOL  = 0.98_rp  ! antiparallel dot-product threshold (>0.98)
  real(rp), parameter :: NC_PLANE_TOL = 0.05_rp   ! coplanarity tol: fraction of sqrt(area_A)
  real(rp), parameter :: NC_TANG_TOL  = 0.75_rp   ! tangential distance tol: fraction of sqrt(area_A)
  real(rp), parameter :: NC_AREA_TOL  = 0.02_rp   ! relative area-sum balance tolerance (2%)

contains

  !> Convert T-junction boundary-face groups to internal faces.
  subroutine enforce_face_conformity(msh)
    type(Mesh), intent(inout) :: msh

    integer(ip) :: n_bnd, a, b, j, ia, ib, f_old, ptr
    integer(ip) :: n_added_int, n_removed_bnd, n_t_junctions, new_n_faces
    integer(ip), allocatable :: bnd_idx(:)       ! face index of each boundary face (1..n_bnd)
    integer(ip), allocatable :: n_fine_of(:)     ! how many fine faces belong to coarse a
    real(rp),    allocatable :: area_sum_of(:)   ! accumulated fine area for coarse a
    integer(ip), allocatable :: fine_off(:)      ! CSR offsets for fine_list (size n_bnd+1)
    integer(ip), allocatable :: fine_list(:)     ! concatenated fine face indices
    integer(ip), allocatable :: assigned_to(:)   ! 0=unassigned, -1=coarse, >0=coarse index
    logical                  :: valid
    real(rp)     :: n_A(3), n_B(3), c_A(3), c_B(3), area_A, area_B, char_len
    real(rp)     :: d_vec(3), d_plane, tang(3), d_tang
    type(FaceData), allocatable :: new_faces(:)
    integer(ip), allocatable :: pcount(:)
    integer(ip) :: p, cnt_p

    n_bnd = msh%n_bnd_faces
    if (n_bnd < 2) return

    ! Map boundary face b (1-based local index) to global face index in msh%faces
    allocate(bnd_idx(n_bnd))
    do b = 1, n_bnd
      bnd_idx(b) = msh%n_int_faces + b
    end do

    allocate(n_fine_of(n_bnd), source=0_ip)
    allocate(area_sum_of(n_bnd), source=0.0_rp)

    ! --- Pass 1: count fine-face candidates for each boundary face ---
    do a = 1, n_bnd
      ia       = bnd_idx(a)
      n_A      = msh%faces(ia)%normal
      c_A      = msh%faces(ia)%centre
      area_A   = msh%faces(ia)%area
      char_len = sqrt(max(area_A, SMALL))

      do b = 1, n_bnd
        if (b == a) cycle

        ib     = bnd_idx(b)
        area_B = msh%faces(ib)%area

        ! (4) B must be strictly smaller than A
        if (area_B >= area_A * (1.0_rp - NC_AREA_TOL)) cycle

        ! (1) Antiparallel normals
        n_B = msh%faces(ib)%normal
        if (dot_product(n_A, n_B) > -NC_NRML_TOL) cycle

        ! (2) Coplanarity
        c_B    = msh%faces(ib)%centre
        d_vec  = c_B - c_A
        d_plane = abs(dot_product(n_A, d_vec))
        if (d_plane > char_len * NC_PLANE_TOL) cycle

        ! (3) Tangential proximity (c_B projects near c_A on the face plane)
        tang  = d_vec - dot_product(d_vec, n_A) * n_A
        d_tang = sqrt(tang(1)**2 + tang(2)**2 + tang(3)**2)
        if (d_tang > char_len * NC_TANG_TOL) cycle

        n_fine_of(a)   = n_fine_of(a) + 1
        area_sum_of(a) = area_sum_of(a) + area_B
      end do
    end do

    ! Build CSR offsets for fine_list
    allocate(fine_off(n_bnd + 1))
    fine_off(1) = 1
    do a = 1, n_bnd
      fine_off(a+1) = fine_off(a) + n_fine_of(a)
    end do

    if (fine_off(n_bnd+1) == 1) then
      ! No candidate T-junctions found; nothing to do.
      deallocate(bnd_idx, n_fine_of, area_sum_of, fine_off)
      return
    end if

    allocate(fine_list(fine_off(n_bnd+1) - 1), source=0_ip)

    ! --- Pass 2: fill fine_list (same geometry checks as pass 1) ---
    n_fine_of = 0  ! reuse as write pointer
    do a = 1, n_bnd
      ia       = bnd_idx(a)
      n_A      = msh%faces(ia)%normal
      c_A      = msh%faces(ia)%centre
      area_A   = msh%faces(ia)%area
      char_len = sqrt(max(area_A, SMALL))

      do b = 1, n_bnd
        if (b == a) cycle
        ib     = bnd_idx(b)
        area_B = msh%faces(ib)%area
        if (area_B >= area_A * (1.0_rp - NC_AREA_TOL)) cycle
        n_B = msh%faces(ib)%normal
        if (dot_product(n_A, n_B) > -NC_NRML_TOL) cycle
        c_B    = msh%faces(ib)%centre
        d_vec  = c_B - c_A
        d_plane = abs(dot_product(n_A, d_vec))
        if (d_plane > char_len * NC_PLANE_TOL) cycle
        tang  = d_vec - dot_product(d_vec, n_A) * n_A
        d_tang = sqrt(tang(1)**2 + tang(2)**2 + tang(3)**2)
        if (d_tang > char_len * NC_TANG_TOL) cycle

        n_fine_of(a) = n_fine_of(a) + 1
        fine_list(fine_off(a) + n_fine_of(a) - 1) = b
      end do
    end do

    ! --- Validate groups and build assignment map ---
    allocate(assigned_to(n_bnd), source=0_ip)
    n_added_int   = 0
    n_removed_bnd = 0
    n_t_junctions = 0

    do a = 1, n_bnd
      if (n_fine_of(a) == 0) cycle

      ia     = bnd_idx(a)
      area_A = msh%faces(ia)%area

      ! Reject if area sum does not balance
      if (abs(area_sum_of(a) - area_A) / max(area_A, SMALL) > NC_AREA_TOL) cycle

      ! Reject if any fine face is already assigned to a different coarse face
      valid = .true.
      do j = fine_off(a), fine_off(a) + n_fine_of(a) - 1
        b = fine_list(j)
        if (assigned_to(b) /= 0) then; valid = .false.; exit; end if
      end do
      if (.not. valid) cycle

      ! Accept this T-junction group
      do j = fine_off(a), fine_off(a) + n_fine_of(a) - 1
        b = fine_list(j)
        assigned_to(b) = a   ! b is a fine face belonging to coarse face a
      end do
      assigned_to(a)  = -1   ! a is a coarse face
      n_t_junctions   = n_t_junctions + 1
      n_added_int     = n_added_int   + n_fine_of(a)
      n_removed_bnd   = n_removed_bnd + 1 + n_fine_of(a)  ! remove 1 coarse + N fine
    end do

    if (n_t_junctions == 0) then
      deallocate(bnd_idx, n_fine_of, area_sum_of, fine_off, fine_list, assigned_to)
      return
    end if

    write(*, '(2X,A,I0,A,I0,A)') &
      'enforce_face_conformity: ', n_t_junctions, &
      ' T-junction(s); replacing ', n_removed_bnd, &
      ' boundary faces with internal faces.'

    ! ---- Rebuild msh%faces ----
    ! New layout: [original internal | new internal | kept boundary]
    new_n_faces = msh%n_int_faces + n_added_int + (n_bnd - n_removed_bnd)
    allocate(new_faces(new_n_faces))

    ! Copy original internal faces unchanged
    do f_old = 1, msh%n_int_faces
      new_faces(f_old) = msh%faces(f_old)
    end do

    ! New internal faces: geometry from fine face, owner = coarse cell.
    ptr = msh%n_int_faces
    do a = 1, n_bnd
      if (assigned_to(a) /= -1) cycle   ! not a coarse face
      ia = bnd_idx(a)
      do j = fine_off(a), fine_off(a) + n_fine_of(a) - 1
        b   = fine_list(j)
        ib  = bnd_idx(b)
        ptr = ptr + 1
        ! Copy fine face geometry (centre, area, normal, d_cf, gc)
        new_faces(ptr)           = msh%faces(ib)
        ! Override topology: coarse cell owns, fine cell is neighbour
        new_faces(ptr)%owner     = msh%faces(ia)%owner    ! coarse cell
        new_faces(ptr)%neighbour = msh%faces(ib)%owner    ! fine cell
        new_faces(ptr)%patch_id  = -1                     ! internal (no patch)
        ! sf oriented owner->neighbour = n_A direction.
        new_faces(ptr)%normal = msh%faces(ia)%normal
        new_faces(ptr)%sf     = msh%faces(ib)%area * msh%faces(ia)%normal
      end do
    end do

    ! Copy kept boundary faces (those not involved in any T-junction)
    do b = 1, n_bnd
      if (assigned_to(b) /= 0) cycle    ! skip coarse (=-1) and fine (=a>0) faces
      ib  = bnd_idx(b)
      ptr = ptr + 1
      new_faces(ptr) = msh%faces(ib)
    end do

    if (ptr /= new_n_faces) &
      error stop 'enforce_face_conformity: face count mismatch - check tolerances'

    ! Update mesh counters
    msh%n_int_faces = msh%n_int_faces + n_added_int
    msh%n_bnd_faces = n_bnd - n_removed_bnd
    msh%n_faces     = new_n_faces

    ! Replace face array
    deallocate(msh%faces)
    allocate(msh%faces(new_n_faces))
    msh%faces = new_faces
    deallocate(new_faces)

    ! Face-vertex topology arrays are now stale; deallocate so callers know to skip them.
    if (allocated(msh%face_offsets))  deallocate(msh%face_offsets)
    if (allocated(msh%face_vertices)) deallocate(msh%face_vertices)

    ! Recompute patch start/count from face patch_id labels.
    if (allocated(msh%patches) .and. msh%n_patches > 0) then
      allocate(pcount(msh%n_patches), source=0_ip)
      do b = 1, msh%n_bnd_faces
        p = msh%faces(msh%n_int_faces + b)%patch_id
        if (p >= 1 .and. p <= msh%n_patches) pcount(p) = pcount(p) + 1
      end do
      cnt_p = 0
      do p = 1, msh%n_patches
        msh%patches(p)%start_face = cnt_p + 1
        msh%patches(p)%n_faces    = pcount(p)
        cnt_p = cnt_p + pcount(p)
      end do
      deallocate(pcount)
    end if

    deallocate(bnd_idx, n_fine_of, area_sum_of, fine_off, fine_list, assigned_to)
  end subroutine enforce_face_conformity

end module mod_nonconforming
