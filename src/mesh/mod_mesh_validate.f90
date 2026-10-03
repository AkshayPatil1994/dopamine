!> mod_mesh_validate
!> Cross-rank consistency self-tests for the local Mesh. Run once at startup
module mod_mesh_validate
  use mpi_f08
  use mod_kinds, only: ip, lp, rp
  use mod_parameters, only: SMALL, BC_PERIODIC
  use mod_mesh, only: Mesh, REFINE_SKEW_TOL
  use mod_warn_log, only: warn_log_message
  implicit none
  private

  public :: validate_processor_faces
  public :: validate_periodic_ghost_centroids
  public :: report_mesh_quality
  public :: report_mesh_integrity

contains

  !> For every processor face (interior face whose neighbour is a ghost cell),
  !> exchange (gc, sf, area, d_cf) with the partner rank and assert:
  subroutine validate_processor_faces(msh, tol, abort_on_fail, n_fail_out)
    type(Mesh), intent(in) :: msh
    real(rp), intent(in), optional :: tol
    logical,  intent(in), optional :: abort_on_fail
    integer(ip), intent(out), optional :: n_fail_out

    real(rp)        :: rtol
    logical         :: do_abort
    integer(ip)     :: nh, h, f, nb, i, ierr
    integer(ip)     :: nc, ng, n_proc_faces_local, n_proc_faces_global
    integer(ip)     :: n_fail_local, n_fail_global
    integer(ip), allocatable :: ghost_halo(:)  ! ghost local-cell -> halo idx (0 if none)
    integer(ip), allocatable :: face_count(:)  ! per-halo processor-face count
    integer(ip), allocatable :: face_off(:)  ! per-halo packing offset
    integer(ip), allocatable :: face_idx(:)  ! flattened face indices, grouped by halo
    integer(ip), allocatable :: order(:)  ! sort permutation per halo
    integer(lp), allocatable :: keys(:,:)  ! (2, nf_h) sorted keys (lo, hi)
    real(rp),    allocatable :: sbuf(:), rbuf(:)
    integer,     parameter   :: NPAY = 6  ! gc, sf(3), area, d_cf
    integer(ip)              :: tag

    if (msh%nranks <= 1) return
    nh = 0
    if (allocated(msh%halos)) nh = size(msh%halos)
    if (nh == 0) return

    rtol     = 1.0e-10_rp
    do_abort = .true.
    if (present(tol)) rtol = tol
    if (present(abort_on_fail)) do_abort = abort_on_fail
    if (present(n_fail_out)) n_fail_out = 0

    nc = msh%n_cells
    ng = msh%n_ghost
    if (ng == 0) return

    ! Build ghost -> halo lookup
    allocate(ghost_halo(nc + ng), source=0_ip)
    do h = 1, nh
      do i = 1, msh%halos(h)%n_recv
        ghost_halo(msh%halos(h)%recv_cells(i)) = h
      end do
    end do

    ! Count processor faces per halo
    allocate(face_count(nh), source=0_ip)
    do f = 1, msh%n_int_faces
      nb = msh%faces(f)%neighbour
      if (nb > nc) then
        h = ghost_halo(nb)
        if (h > 0) face_count(h) = face_count(h) + 1
      end if
    end do

    allocate(face_off(nh + 1))
    face_off(1) = 0
    do h = 1, nh
      face_off(h + 1) = face_off(h) + face_count(h)
    end do
    n_proc_faces_local = face_off(nh + 1)

    allocate(face_idx(max(1, n_proc_faces_local)))
    face_count = 0
    do f = 1, msh%n_int_faces
      nb = msh%faces(f)%neighbour
      if (nb > nc) then
        h = ghost_halo(nb)
        if (h > 0) then
          face_count(h) = face_count(h) + 1
          face_idx(face_off(h) + face_count(h)) = f
        end if
      end if
    end do

    n_fail_local = 0
    tag = 19427

    ! Per-halo: sort by (lo,hi) gid pair, then exchange payload
    do h = 1, nh
      block
        integer(ip) :: nf_h, k, kk
        integer(lp) :: gid_o, gid_n, lo, hi
        real(rp)    :: gc_l, gc_r, area_l, area_r, dcf_l, dcf_r
        real(rp)    :: sf_l(3), sf_r(3)
        type(MPI_Request) :: req(2)

        nf_h = face_count(h)
        if (nf_h == 0) cycle

        allocate(keys(2, nf_h))
        allocate(order(nf_h))
        allocate(sbuf(NPAY * nf_h))
        allocate(rbuf(NPAY * nf_h))

        do k = 1, nf_h
          f     = face_idx(face_off(h) + k)
          gid_o = msh%global_id(msh%faces(f)%owner)
          gid_n = msh%global_id(msh%faces(f)%neighbour)
          lo    = min(gid_o, gid_n)
          hi    = max(gid_o, gid_n)
          keys(1, k) = lo
          keys(2, k) = hi
          order(k)   = k
        end do
        call sort_pairs(keys, order, nf_h)

        do k = 1, nf_h
          f = face_idx(face_off(h) + order(k))
          kk = NPAY * (k - 1)
          sbuf(kk + 1) = msh%faces(f)%gc
          sbuf(kk + 2) = msh%faces(f)%sf(1)
          sbuf(kk + 3) = msh%faces(f)%sf(2)
          sbuf(kk + 4) = msh%faces(f)%sf(3)
          sbuf(kk + 5) = msh%faces(f)%area
          sbuf(kk + 6) = msh%faces(f)%d_cf
        end do

        call MPI_Irecv(rbuf, NPAY * nf_h, MPI_DOUBLE_PRECISION, &
                       msh%halos(h)%rank, tag, msh%comm, req(1), ierr)
        call MPI_Isend(sbuf, NPAY * nf_h, MPI_DOUBLE_PRECISION, &
                       msh%halos(h)%rank, tag, msh%comm, req(2), ierr)
        call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)

        do k = 1, nf_h
          kk     = NPAY * (k - 1)
          gc_l   = sbuf(kk + 1); gc_r   = rbuf(kk + 1)
          sf_l   = sbuf(kk + 2 : kk + 4)
          sf_r   = rbuf(kk + 2 : kk + 4)
          area_l = sbuf(kk + 5); area_r = rbuf(kk + 5)
          dcf_l  = sbuf(kk + 6); dcf_r  = rbuf(kk + 6)

          if (abs((gc_l + gc_r) - 1.0_rp) > rtol .or. &
              maxval(abs(sf_l + sf_r))    > rtol .or. &
              abs(area_l - area_r)        > rtol .or. &
              abs(dcf_l  - dcf_r)         > rtol) then
            n_fail_local = n_fail_local + 1
            if (n_fail_local <= 5) then
              f = face_idx(face_off(h) + order(k))
              write(*, '(A,I0,A,I0,A,I0,A,I0)') &
                '[mesh-validate] rank ', msh%rank, &
                ' face ', f, ' partner-rank ', msh%halos(h)%rank, &
                ' face_idx_in_halo ', k
              write(*, '(A,ES14.6,A,ES14.6,A,ES14.6)') &
                '   gc_local=', gc_l, ' gc_remote=', gc_r, ' sum-1=', (gc_l+gc_r)-1.0_rp
              write(*, '(A,3ES14.6)') '   sf_local + sf_remote = ', sf_l + sf_r
              write(*, '(A,ES14.6,A,ES14.6)') '   area_l=', area_l, ' area_r=', area_r
              write(*, '(A,ES14.6,A,ES14.6)') '   dcf_l =', dcf_l,  ' dcf_r =', dcf_r
            end if
          end if
        end do

        deallocate(keys, order, sbuf, rbuf)
      end block
    end do

    ! Global reduction
    call MPI_Allreduce(n_fail_local,        n_fail_global,        1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(n_proc_faces_local,  n_proc_faces_global,  1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)

    if (msh%rank == 0) then
      if (n_fail_global == 0) then
        write(*, '(A,I0,A)') ' [mesh-validate] processor-face geometry consistent (', &
                              n_proc_faces_global, ' face-pairs checked)'
      else
        write(*, '(A,I0,A,I0,A)') ' [mesh-validate] FAILED: ', n_fail_global, &
                                  ' inconsistent processor faces of ', n_proc_faces_global
      end if
    end if

    if (n_fail_global > 0 .and. do_abort) then
      call MPI_Abort(msh%comm, 71, ierr)
    end if
    if (present(n_fail_out)) n_fail_out = n_fail_global

    deallocate(ghost_halo, face_count, face_off, face_idx)
  end subroutine validate_processor_faces

  !> In-place insertion sort of the permutation `order` so that
  !> keys(:, order(:)) is non-decreasing in (key1, key2). Insertion sort is
  subroutine sort_pairs(keys, order, n)
    integer(lp), intent(in)    :: keys(:,:)
    integer(ip), intent(inout) :: order(:)
    integer(ip), intent(in)    :: n
    integer(ip) :: i, j, tmp
    integer(lp) :: k1, k2
    do i = 2, n
      tmp = order(i)
      k1  = keys(1, tmp)
      k2  = keys(2, tmp)
      j   = i - 1
      do while (j >= 1)
        if (keys(1, order(j)) < k1) exit
        if (keys(1, order(j)) == k1 .and. keys(2, order(j)) <= k2) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = tmp
    end do
  end subroutine sort_pairs

  !> Periodic ghost centroid storage assertion.
  !> For every periodic boundary face we recover the periodic translation
  subroutine validate_periodic_ghost_centroids(msh, tol, abort_on_fail, n_fail_out)
    type(Mesh), intent(in) :: msh
    real(rp), intent(in), optional :: tol
    logical,  intent(in), optional :: abort_on_fail
    integer(ip), intent(out), optional :: n_fail_out

    real(rp)    :: rtol, mag_ref, mag_diff
    real(rp)    :: ref_shift(3), shift(3)
    logical     :: do_abort, have_ref
    integer(ip) :: p, f_loc, f_bnd, f_glob, owner, pair_bnd, pair_face, pair_cell
    integer(ip) :: n_fail_local, n_fail_global, n_checked_local, n_checked_global, ierr

    rtol     = 1.0e-9_rp
    do_abort = .true.
    if (present(tol)) rtol = tol
    if (present(abort_on_fail)) do_abort = abort_on_fail
    if (present(n_fail_out)) n_fail_out = 0

    n_fail_local    = 0
    n_checked_local = 0

    if (msh%n_patches <= 0) goto 900

    do p = 1, msh%n_patches
      if (msh%patches(p)%bc_id /= BC_PERIODIC) cycle
      have_ref = .false.
      ref_shift = 0.0_rp
      mag_ref   = 0.0_rp

      do f_loc = 1, msh%patches(p)%n_faces
        f_bnd  = msh%patches(p)%start_face + f_loc - 1
        f_glob = msh%n_int_faces + f_bnd
        if (f_glob < 1 .or. f_glob > msh%n_faces) cycle
        owner = msh%faces(f_glob)%owner
        if (owner < 1 .or. owner > size(msh%cells)) cycle

        ! Resolve partner cell + partner face d_cf.
        pair_bnd = 0; pair_face = 0; pair_cell = 0
        if (allocated(msh%periodic_pair)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) pair_bnd = msh%periodic_pair(f_bnd)
        end if
        if (pair_bnd > 0) then
          pair_face = msh%n_int_faces + pair_bnd
          if (pair_face >= 1 .and. pair_face <= msh%n_faces) pair_cell = msh%faces(pair_face)%owner
        end if
        if (pair_cell <= 0 .and. allocated(msh%periodic_neighbor)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) pair_cell = msh%periodic_neighbor(f_bnd)
        end if
        ! For preprocessed distributed meshes the decomposer stores the periodic partner
        ! as a ghost cell directly in faces%neighbour; use it as the authoritative fallback.
        if (pair_cell <= 0) then
          if (f_glob >= 1 .and. f_glob <= msh%n_faces) pair_cell = msh%faces(f_glob)%neighbour
        end if
        if (pair_cell <= 0 .or. pair_cell > size(msh%cells)) cycle

        ! Face-centre translation is strictly constant within a periodic patch for any
        ! valid mesh (graded or not). 
        if (pair_face < 1 .or. pair_face > msh%n_faces) cycle

        shift = msh%faces(pair_face)%centre - msh%faces(f_glob)%centre
        n_checked_local = n_checked_local + 1

        if (.not. have_ref) then
          ref_shift = shift
          mag_ref   = sqrt(sum(ref_shift**2))
          have_ref  = .true.
          cycle
        end if

        mag_diff = sqrt(sum((shift - ref_shift)**2))
        if (mag_diff > rtol * max(mag_ref, 1.0_rp)) then
          n_fail_local = n_fail_local + 1
          if (n_fail_local <= 5) then
            write(*, '(A,I0,A,I0,A,I0,A)') &
              '[mesh-validate] rank ', msh%rank, &
              ' periodic patch ', p, ' face ', f_glob, &
              ' has inconsistent ghost-centroid shift vs patch reference:'
            write(*, '(A,3ES14.6)') '   ref_shift = ', ref_shift
            write(*, '(A,3ES14.6)') '   shift     = ', shift
            write(*, '(A,ES14.6)')  '   |diff|    = ', mag_diff
          end if
        end if
      end do
    end do

900 continue
    if (msh%nranks > 1) then
      call MPI_Allreduce(n_fail_local,    n_fail_global,    1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
      call MPI_Allreduce(n_checked_local, n_checked_global, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    else
      n_fail_global    = n_fail_local
      n_checked_global = n_checked_local
    end if

    if (msh%rank == 0) then
      if (n_fail_global == 0) then
        write(*, '(A,I0,A)') ' [mesh-validate] periodic ghost centroids consistent (', &
                              n_checked_global, ' periodic faces checked)'
      else
        write(*, '(A,I0,A,I0,A)') ' [mesh-validate] FAILED: ', n_fail_global, &
                                  ' inconsistent periodic ghost centroids of ', n_checked_global
      end if
    end if

    if (n_fail_global > 0 .and. do_abort) then
      call MPI_Abort(msh%comm, 73, ierr)
    end if
    if (present(n_fail_out)) n_fail_out = n_fail_global
  end subroutine validate_periodic_ghost_centroids

  subroutine report_mesh_quality(msh, nonortho_max_deg_out, n_refine_out)
    use mpi_f08, only: MPI_DOUBLE_PRECISION, MPI_MAX, MPI_Allreduce, MPI_INTEGER8, MPI_SUM
    type(Mesh), intent(in) :: msh
    real(rp), optional, intent(out) :: nonortho_max_deg_out
    integer(lp), optional, intent(out) :: n_refine_out

    integer(ip) :: f, o, nb
    real(rp)    :: d_on(3), cos_th, sin2_th, skew_mag, d_mid(3)
    real(rp)    :: nonortho_max_local, skew_max_local, nonortho_max, skew_max, h_ref
    real(rp)    :: skew_mid_local, skew_mid_max
    integer(lp) :: n_refine_local, n_refine, n_int_local, n_int_total
    real(rp)    :: refine_frac
    integer     :: ierr
    logical     :: count_face

    nonortho_max_local = 0.0_rp
    skew_max_local     = 0.0_rp
    n_refine_local     = 0_lp
    skew_mid_local     = 0.0_rp
    n_int_local        = 0_lp

    do f = 1, msh%n_int_faces
      o  = msh%faces(f)%owner
      nb = msh%faces(f)%neighbour
      if (o < 1 .or. o > size(msh%cells)) cycle
      if (nb < 1 .or. nb > size(msh%cells)) cycle
      ! A processor face exists on both ranks; count it once, on the rank owning the lower global id.
      count_face = .true.
      if (nb > msh%n_cells .or. o > msh%n_cells) then
        if (allocated(msh%global_id)) count_face = (msh%global_id(o) < msh%global_id(nb))
      end if
      if (count_face) n_int_local = n_int_local + 1_lp

      d_on  = msh%cells(nb)%centroid - msh%cells(o)%centroid
      h_ref = sqrt(sum(d_on**2))
      if (h_ref < 1.0e-30_rp) cycle

      cos_th = dot_product(msh%faces(f)%sf, d_on) / (msh%faces(f)%area * h_ref)
      cos_th = max(-1.0_rp, min(1.0_rp, cos_th))
      sin2_th = 1.0_rp - cos_th**2
      nonortho_max_local = max(nonortho_max_local, acos(abs(cos_th)))

      ! gc is the owner weight, so the gc-interpolation point is x_o + (1-gc)*d_on.
      d_mid  = msh%faces(f)%centre - (msh%cells(o)%centroid + (1.0_rp - msh%faces(f)%gc) * d_on)
      skew_mag = sqrt(sum(d_mid**2)) / max(msh%cells(o)%volume**(1.0_rp/3.0_rp), 1.0e-30_rp)
      skew_max_local = max(skew_max_local, skew_mag)

      if (msh%faces(f)%is_refine .and. count_face) n_refine_local = n_refine_local + 1_lp
      skew_mid_local = max(skew_mid_local, msh%faces(f)%skew)
    end do

    call MPI_Allreduce(nonortho_max_local, nonortho_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)
    call MPI_Allreduce(skew_max_local,     skew_max,     1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)
    call MPI_Allreduce(n_refine_local, n_refine,    1, MPI_INTEGER8, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(n_int_local,    n_int_total, 1, MPI_INTEGER8, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(skew_mid_local, skew_mid_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)

    refine_frac = 0.0_rp
    if (n_int_total > 0_lp) refine_frac = 100.0_rp * real(n_refine, rp) / real(n_int_total, rp)

    if (msh%rank == 0) then
      write(*, '(A,F6.1,A)') ' [mesh-quality] max non-orthogonality = ', nonortho_max * 180.0_rp / acos(-1.0_rp), ' deg'
      write(*, '(A,F8.4)')   ' [mesh-quality] max face skewness      = ', skew_max
      write(*, '(A,I0,A,F6.2,A)') ' [mesh-quality] refinement-interface faces = ', n_refine, &
                                  ' (', refine_frac, '% of internal)'
      ! A zero count with visible refinement means the metric or the tag path is wrong, not the mesh.
      if (n_refine == 0_lp) &
        write(*, '(A,F8.4,A,F8.4)') ' [mesh-quality] none tagged: max midpoint skew = ', skew_mid_max, &
                                    ' vs threshold ', REFINE_SKEW_TOL
      if (nonortho_max * 180.0_rp / acos(-1.0_rp) > 70.0_rp) &
        call warn_log_message('[mesh-quality]', 'non-orthogonality > 70 deg; mesh quality may be insufficient')
      if (nonortho_max * 180.0_rp / acos(-1.0_rp) > 50.0_rp) &
        write(*, '(A)') ' [mesh-quality] INFO: non-orthogonality > 50 deg; recommend n_nonortho_corr >= 2'
      ! Green-Gauss is only advisable when no side of a cell is split; see test_gradient_refinement.
      if (skew_max > 0.5_rp .and. n_refine == 0_lp) &
        write(*, '(A)') ' [mesh-quality] INFO: skewness > 0.5; enable GREEN_GAUSS with gg_iterations >= 2'
    end if
    if (present(nonortho_max_deg_out)) nonortho_max_deg_out = nonortho_max * 180.0_rp / acos(-1.0_rp)
    if (present(n_refine_out)) n_refine_out = n_refine
  end subroutine report_mesh_quality

  !> Check the geometric preconditions the FV operators assume: face closure, positive Sf.d_on, no slivers.
  subroutine report_mesh_integrity(msh)
    use mpi_f08, only: MPI_DOUBLE_PRECISION, MPI_MAX, MPI_SUM, MPI_INTEGER8, MPI_Allreduce
    type(Mesh), intent(in) :: msh

    integer(ip) :: f, c, o, nb
    real(rp)    :: closure(3), d_on(3), s_dot_d, h3, ar
    real(rp)    :: closure_max_local, closure_max
    integer(lp) :: n_neg_local, n_neg, n_sliver_local, n_sliver
    integer     :: ierr

    closure_max_local = 0.0_rp
    n_neg_local       = 0_lp
    n_sliver_local    = 0_lp

    ! Discrete Gauss closure: sum of outward Sf over a cell must vanish relative to its area scale.
    do c = 1, msh%n_cells
      closure = 0.0_rp
      ar      = 0.0_rp
      do f = msh%cells(c)%face_offset, msh%cells(c)%face_offset + msh%cells(c)%n_faces - 1
        if (f < 1 .or. f > size(msh%cell_faces)) cycle
        block
          integer(ip) :: fg
          fg = msh%cell_faces(f)
          if (fg < 1 .or. fg > msh%n_faces) cycle
          ar = ar + msh%faces(fg)%area
          if (msh%faces(fg)%owner == c) then
            closure = closure + msh%faces(fg)%sf
          else
            closure = closure - msh%faces(fg)%sf
          end if
        end block
      end do
      if (ar > SMALL) closure_max_local = max(closure_max_local, sqrt(sum(closure**2)) / ar)

      ! Sliver detection: cell volume far below the cube of its own area scale.
      h3 = (ar / 6.0_rp)**1.5_rp
      if (h3 > SMALL) then
        if (msh%cells(c)%volume / h3 < 0.01_rp) n_sliver_local = n_sliver_local + 1_lp
      end if
    end do

    ! Sf.d_on <= 0 makes the Poisson coefficient non-positive and silently trips the d_cf fallback.
    do f = 1, msh%n_int_faces
      o  = msh%faces(f)%owner
      nb = msh%faces(f)%neighbour
      if (o < 1 .or. nb < 1 .or. o > size(msh%cells) .or. nb > size(msh%cells)) cycle
      d_on    = msh%cells(nb)%centroid - msh%cells(o)%centroid
      s_dot_d = dot_product(msh%faces(f)%sf, d_on)
      if (s_dot_d <= SMALL) n_neg_local = n_neg_local + 1_lp
    end do

    call MPI_Allreduce(closure_max_local, closure_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)
    call MPI_Allreduce(n_neg_local,       n_neg,       1, MPI_INTEGER8, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(n_sliver_local,    n_sliver,    1, MPI_INTEGER8, MPI_SUM, msh%comm, ierr)

    if (msh%rank == 0) then
      write(*, '(A,ES10.3)') ' [mesh-integrity] max face-closure residual  = ', closure_max
      write(*, '(A,I0)')     ' [mesh-integrity] faces with Sf.d_on <= 0    = ', n_neg
      write(*, '(A,I0)')     ' [mesh-integrity] sliver cells (V/h^3<0.01)  = ', n_sliver
      if (closure_max > 1.0e-10_rp) &
        call warn_log_message('[mesh-integrity]', 'face closure residual > 1e-10; FV operators are inconsistent')
      if (n_neg > 0_lp) &
        call warn_log_message('[mesh-integrity]', 'Sf.d_on <= 0 faces; Poisson falls back to a different operator')
      if (n_sliver > 0_lp) &
        call warn_log_message('[mesh-integrity]', 'sliver cells present; 1/V amplifies every error there')
    end if
  end subroutine report_mesh_integrity

end module mod_mesh_validate
