! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Mesh partitioning for MPI runs.
module mod_partitioner
  use mod_kinds
  use mod_mesh
  use mod_parameters, only: SMALL
  use mpi_f08
  implicit none
  private

  public :: partition_mesh, build_hilbert_partition
  public :: pack_global_mesh, sort_and_offset_partition, build_local_partition_mesh
  public :: compute_hilbert_key

contains

  subroutine partition_mesh(mesh_global, mesh_local, comm)
    type(Mesh),     intent(inout) :: mesh_global
    type(Mesh),     intent(inout) :: mesh_local
    type(MPI_Comm), intent(in)    :: comm

    integer(ip) :: nranks, ierr

    call MPI_Comm_size(comm, nranks, ierr)

    if (nranks == 1) then
      mesh_local = mesh_global
      mesh_local%n_cells_global = int(mesh_global%n_cells, lp)  ! For single-rank: global = local
      return
    end if

    call partition_mesh_hilbert(mesh_global, mesh_local, comm)

  end subroutine partition_mesh


  !> Build a Hilbert space-filling curve partition; cell_centroid must be broadcast to every rank first.
  subroutine build_hilbert_partition(nc, nranks, rank, cell_centroid, part, cell_volume, weight_by_volume)
    integer(ip), intent(in)  :: nc, nranks, rank
    real(rp), intent(in)     :: cell_centroid(3, nc)
    real(rp), intent(in), optional :: cell_volume(nc)
    integer(ip), intent(out) :: part(nc)
    logical, intent(in), optional :: weight_by_volume  ! default: balance cell count, not volume

    integer(ip) :: i, p, min_cnt, max_cnt
    real(rp)    :: avg_cnt, xmin, xmax, ymin, ymax, zmin, zmax, scale
    real(rp)    :: total_volume, target_volume, cumul_volume
    integer(ip), allocatable :: idx(:), part_counts(:)
    integer(lp), allocatable :: hilbert_keys(:)
    logical :: use_load_balance

    ! Guard against zero-size array (nc == 0)
    if (nc == 0) return

    ! Equal cell count is the default: solver cost scales with cells, not volume, so
    ! volume weighting badly imbalances graded meshes (3.5x at 1024 ranks on a tanh mesh).
    use_load_balance = .false.
    if (present(weight_by_volume)) use_load_balance = weight_by_volume .and. present(cell_volume)

    ! Find bounding box
    xmin = minval(cell_centroid(1, :))
    xmax = maxval(cell_centroid(1, :))
    ymin = minval(cell_centroid(2, :))
    ymax = maxval(cell_centroid(2, :))
    zmin = minval(cell_centroid(3, :))
    zmax = maxval(cell_centroid(3, :))

    ! Expand slightly to avoid boundary cells mapping to 2^order
    scale = max(xmax - xmin, ymax - ymin, zmax - zmin) * 1.001_rp + 1.0e-14_rp

    ! Compute Hilbert keys for all cells
    allocate(hilbert_keys(nc), idx(nc))
    do i = 1, nc
      hilbert_keys(i) = compute_hilbert_key(cell_centroid(1:3, i), &
                                            [xmin, ymin, zmin], scale)
    end do

    ! Argsort by Hilbert key
    call argsort_int64(hilbert_keys, idx)

    ! Assign sorted cells to ranks
    if (use_load_balance) then
      ! Volume-weighted load balancing: assign cells so each rank gets roughly equal total volume
      total_volume = sum(cell_volume)
      target_volume = total_volume / real(nranks, rp)
      cumul_volume = 0.0_rp
      p = 0
      do i = 1, nc
        cumul_volume = cumul_volume + cell_volume(idx(i))
        ! Assign to current rank until cumulative volume exceeds target
        if (p < nranks - 1 .and. cumul_volume >= target_volume * real(p + 1, rp)) then
          p = p + 1
        end if
        part(idx(i)) = p
      end do
    else
      ! Equal cell count per rank (naive partitioning)
      do i = 1, nc
        p = idx(i)
        part(p) = int(int(i - 1, lp) * int(nranks, lp) / int(nc, lp), ip)
      end do
    end if

    ! Report statistics (rank 0 only)
    if (rank == 0) then
      allocate(part_counts(nranks), source=0_ip)
      do i = 1, nc
        if (part(i) >= 0 .and. part(i) < nranks) then
          part_counts(part(i) + 1) = part_counts(part(i) + 1) + 1
        end if
      end do
      min_cnt = minval(part_counts)
      max_cnt = maxval(part_counts)
      avg_cnt = real(nc, rp) / real(max(1_ip, nranks), rp)
      write(*, '(2X,A,I0,A,I0,A,F8.3)') 'Cells/rank  min=', min_cnt, '  max=', max_cnt, &
                                         '  imbalance=', real(max_cnt, rp) / max(avg_cnt, 1.0_rp)
      if (use_load_balance) then
        write(*, '(2X,A)') 'Load balance: volume-weighted'
      else
        write(*, '(2X,A)') 'Load balance: equal cell count'
      end if
      deallocate(part_counts)
    end if

    deallocate(idx, hilbert_keys)

  end subroutine build_hilbert_partition

  !> Compute 3D Hilbert curve index for a point.
  !> Normalized coordinates in [xmin, xmin+scale) mapped to integer space.
  function compute_hilbert_key(xyz, origin, scale) result(key)
    real(rp), intent(in)  :: xyz(3), origin(3), scale
    integer(lp)           :: key

    integer(ip), parameter :: ORDER = 20
    integer(ip) :: ix, iy, iz, i
    integer(ip) :: x, y, z, rx, ry, rz
    integer(lp) :: h

    ! Normalize to [0, 2^ORDER)
    ix = max(0_ip, min(2_ip**ORDER - 1_ip, &
             int((xyz(1) - origin(1)) / scale * real(2_ip**ORDER, rp), ip)))
    iy = max(0_ip, min(2_ip**ORDER - 1_ip, &
             int((xyz(2) - origin(2)) / scale * real(2_ip**ORDER, rp), ip)))
    iz = max(0_ip, min(2_ip**ORDER - 1_ip, &
             int((xyz(3) - origin(3)) / scale * real(2_ip**ORDER, rp), ip)))

    x = ix
    y = iy
    z = iz
    h = 0_lp

    ! 3D Hilbert curve: process bits from MSB to LSB with rotation
    do i = ORDER - 1, 0, -1
      rx = ibits(x, i, 1)
      ry = ibits(y, i, 1)
      rz = ibits(z, i, 1)

      ! Compute Hilbert index bit using rotation
      h = h * 8_lp + int(hilbert_combine_bits(rx, ry, rz), lp)
    end do

    key = h

  end function compute_hilbert_key

  !> Combine 3 bits into Hilbert index with rotation.
  pure function hilbert_combine_bits(rx, ry, rz) result(h_idx)
    integer(ip), intent(in) :: rx, ry, rz
    integer(ip) :: h_idx

    ! Map (rx, ry, rz) to Hilbert octant index [0..7]
    h_idx = ior(ior(ishft(rx, 2), ishft(ry, 1)), rz)

  end function hilbert_combine_bits

  !> Argsort array of 64-bit integers using quicksort, return indices in ascending order.
  subroutine argsort_int64(arr, idx)
    integer(lp), intent(in)    :: arr(:)
    integer(ip), intent(out)   :: idx(:)

    integer(ip) :: i, n

    n = size(arr)
    do i = 1, n
      idx(i) = i
    end do

    ! Use radix sort for Hilbert keys (5-10x faster than quicksort for uniformly distributed keys)
    call radix_sort_int64(arr, idx, n)

  end subroutine argsort_int64

  !> Radix sort for 64-bit integers using counting sort by bytes (256 radices).
  !> O(8*n) complexity, much faster than quicksort for uniformly distributed keys like Hilbert curves.
  subroutine radix_sort_int64(arr, idx, n)
    integer(ip), intent(in)    :: n
    integer(lp), intent(in)    :: arr(:)
    integer(ip), intent(inout) :: idx(:)

    integer(ip), allocatable :: idx_tmp(:), count(:), offset(:)
    integer(ip) :: i, byte_idx, key_byte
    integer(ip) :: pos

    if (n <= 1) return

    allocate(idx_tmp(n), count(0:255), offset(0:255))

    ! Process each byte of the 64-bit key (8 passes, 256 radices)
    do byte_idx = 0, 7
      count = 0
      ! Count occurrences of each byte value
      do i = 1, n
        key_byte = int(ibits(arr(idx(i)), byte_idx*8, 8), ip)
        count(key_byte) = count(key_byte) + 1
      end do
      ! Compute cumulative offsets
      offset(0) = 1
      do i = 1, 255
        offset(i) = offset(i-1) + count(i-1)
      end do
      ! Place elements into sorted order
      do i = 1, n
        key_byte = int(ibits(arr(idx(i)), byte_idx*8, 8), ip)
        pos = offset(key_byte)
        idx_tmp(pos) = idx(i)
        offset(key_byte) = offset(key_byte) + 1
      end do
      ! Swap arrays for next pass
      idx(1:n) = idx_tmp(1:n)
    end do

    deallocate(idx_tmp, count, offset)

  end subroutine radix_sort_int64

  !> Recursive quicksort on index array, sorting by array values.
  recursive subroutine quicksort_indices(arr, idx, left, right)
    integer(lp), intent(in)    :: arr(:)
    integer(ip), intent(inout) :: idx(:)
    integer(ip), intent(in)    :: left, right

    integer(ip) :: i, j, pivot_idx
    integer(lp) :: pivot_val
    integer(ip) :: temp_idx

    if (left >= right) return

    ! Choose middle element as pivot
    pivot_idx = (left + right) / 2
    pivot_val = arr(idx(pivot_idx))

    ! Partition
    i = left
    j = right
    do while (i <= j)
      do while (arr(idx(i)) < pivot_val)
        i = i + 1
      end do
      do while (arr(idx(j)) > pivot_val)
        j = j - 1
      end do
      if (i <= j) then
        temp_idx = idx(i)
        idx(i) = idx(j)
        idx(j) = temp_idx
        i = i + 1
        j = j - 1
      end if
    end do

    ! Recursive calls
    if (left < j) call quicksort_indices(arr, idx, left, j)
    if (i < right) call quicksort_indices(arr, idx, i, right)

  end subroutine quicksort_indices

  subroutine pack_global_mesh(msh, nc, nf, np, cell_centroid, cell_volume, &
                              face_owner, face_nei, face_patch, face_area, face_centre, &
                              face_normal, face_sf, face_dcf, face_gc, patch_start, patch_nfaces, &
                              patch_bcid, patch_name, patch_type)
    type(Mesh), intent(in) :: msh
    integer(ip), intent(in) :: nc, nf, np
    real(rp), intent(out) :: cell_centroid(3, nc), cell_volume(nc)
    integer(ip), intent(out) :: face_owner(nf), face_nei(nf), face_patch(nf)
    real(rp), intent(out) :: face_area(nf), face_centre(3, nf), face_normal(3, nf), face_sf(3, nf), face_dcf(nf), face_gc(nf)
    integer(ip), intent(out) :: patch_start(np), patch_nfaces(np), patch_bcid(np)
    character(len=128), intent(out) :: patch_name(np)
    character(len=64), intent(out)  :: patch_type(np)
    integer(ip) :: i

    do i = 1, nc
      cell_centroid(:, i) = msh%cells(i)%centroid
      cell_volume(i)      = msh%cells(i)%volume
    end do

    do i = 1, nf
      face_owner(i)     = msh%faces(i)%owner
      face_nei(i)       = msh%faces(i)%neighbour
      face_patch(i)     = msh%faces(i)%patch_id
      face_area(i)      = msh%faces(i)%area
      face_centre(:, i) = msh%faces(i)%centre
      face_normal(:, i) = msh%faces(i)%normal
      face_sf(:, i)     = msh%faces(i)%sf
      face_dcf(i)       = msh%faces(i)%d_cf
      face_gc(i)        = msh%faces(i)%gc
    end do

    do i = 1, np
      patch_start(i)  = msh%patches(i)%start_face
      patch_nfaces(i) = msh%patches(i)%n_faces
      patch_bcid(i)   = msh%patches(i)%bc_id
      patch_name(i)   = msh%patches(i)%name
      patch_type(i)   = msh%patches(i)%bc_type
    end do
  end subroutine pack_global_mesh

  subroutine build_local_partition_mesh(msh, comm, target_rank, nranks, nc, nf, np, part, &
                                        cell_centroid, cell_volume, &
                                        face_owner, face_nei, face_patch, &
                                        face_area, face_centre, face_normal, face_dcf, &
                                        patch_start, patch_nfaces, patch_bcid, patch_name, patch_type, &
                                        periodic_pair_global, &
                                        topo_flag, n_vertices_global, vertices_global, cell_offsets_global, &
                                        cell_vertices_global, cell_types_global, &
                                        face_topo_flag, face_offsets_global, face_vertices_global, &
                                        cell_order, cell_start, cell_end, face_chunk_size, face_sf, face_gc, staged_geom_unit, &
                                        int_face_list, n_int_faces_global)
    type(Mesh), intent(inout) :: msh
    type(MPI_Comm), intent(in) :: comm
    integer(ip), intent(in) :: target_rank, nranks, nc, nf, np
    integer(ip), intent(in) :: part(nc)
    integer(ip), intent(in), optional :: cell_order(:), cell_start, cell_end
    integer(ip), intent(in), optional :: face_chunk_size
    real(rp), intent(in) :: cell_centroid(3, nc), cell_volume(nc)
    integer(ip), intent(in) :: face_owner(nf), face_nei(nf), face_patch(nf)
    real(rp), intent(in) :: face_area(nf), face_dcf(nf)
    real(rp), intent(in), optional :: face_centre(3, nf), face_normal(3, nf)
    real(rp), intent(in), optional :: face_sf(3, nf)
    real(rp), intent(in), optional :: face_gc(nf)
    integer(ip), intent(in), optional :: staged_geom_unit
    integer(ip), intent(in), optional :: int_face_list(:)
    integer(ip), intent(in), optional :: n_int_faces_global
    integer(ip), intent(in) :: patch_start(np), patch_nfaces(np), patch_bcid(np)
    character(len=128), intent(in) :: patch_name(np)
    character(len=64), intent(in)  :: patch_type(np)
    integer(ip), intent(in) :: periodic_pair_global(:)
    integer(ip), intent(in) :: topo_flag, n_vertices_global
    real(rp),    intent(in) :: vertices_global(3, n_vertices_global)
    integer(ip), intent(in) :: cell_offsets_global(nc + 1), cell_vertices_global(:), cell_types_global(nc)
    integer(ip), intent(in) :: face_topo_flag
    integer(ip), intent(in) :: face_offsets_global(:), face_vertices_global(:)

    integer(ip) :: i, f, p, o, nb, n_int_global, f_bnd
    integer(ip) :: g_bnd, g_pair, pair_face, nb_global
    integer(ip) :: n_owned, n_ghost, n_int_local, n_bnd_local, n_faces_local
    integer(ip) :: lo, ln, cursor, bcursor
    integer(ip) :: gid, n0, n1, nconn, pos, vlocal, v, g0, g1, nfv
    integer(ip) :: f0, f1, face_chunk
    integer(ip) :: staged_cache_size, staged_cache_start, staged_cache_end
    real(rp) :: sf_local(3), gc_local, d_on(3), d_oc(3), d_len
    real(rp) :: centre_local(3), normal_local(3)
    logical :: flip
    logical :: have_sf_gc
    integer(ip), allocatable :: owned_ids(:), ghost_ids(:), patch_local_count(:), bnd_order(:), g2lbnd(:)
    integer(ip), allocatable :: ghost_candidates(:), map_gid(:), map_lid(:)
    integer(ip) :: n_ghost_candidates
    integer(ip), allocatable :: vlist(:), local_face_global(:)
    integer(ip) :: n_vlist, vpos
    real(rp), allocatable :: staged_geom_block(:,:)

    if (present(face_chunk_size)) then
      face_chunk = face_chunk_size
    else
      face_chunk = nf
    end if
    face_chunk = max(1_ip, face_chunk)

    staged_cache_size = 0
    staged_cache_start = 1
    staged_cache_end = 0
    if (present(staged_geom_unit)) then
      staged_cache_size = max(256_ip, min(face_chunk, nf))
      ! 10 reals per face: centre[3] + normal[3] + sf[3] + gc[1].
      ! See fetch_face_geom; the streaming decomposer writes the same layout.
      allocate(staged_geom_block(10, staged_cache_size))
    end if

    ! If cell_order provided (from sorted partition), use contiguous range; else count
    if (present(cell_order) .and. present(cell_start) .and. present(cell_end)) then
      n_owned = cell_end - cell_start + 1
    else
      n_owned = count(part == target_rank)
    end if
    ! Global constant: take it from the caller when available, it costs an O(nf) scan per rank.
    if (present(n_int_faces_global)) then
      n_int_global = n_int_faces_global
    else
      n_int_global = count(face_nei(1:nf) >= 1)
    end if
    allocate(owned_ids(n_owned), patch_local_count(np))
    patch_local_count = 0

    ! Build owned_ids from sorted cell_order or from full iteration
    if (present(cell_order) .and. present(cell_start) .and. present(cell_end)) then
      do i = cell_start, cell_end
        owned_ids(i - cell_start + 1) = cell_order(i)
      end do
    else
      cursor = 0
      do i = 1, nc
        if (part(i) == target_rank) then
          cursor = cursor + 1
          owned_ids(cursor) = i
        end if
      end do
    end if

    n_ghost_candidates = 0
    call init_dynamic_int_buffer(ghost_candidates, max(1024_ip, n_owned / 2_ip))
    if (present(int_face_list)) then
      do i = 1, size(int_face_list)
        f = int_face_list(i)
        o = face_owner(f)
        nb = face_nei(f)
        if (part(o) == target_rank .and. part(nb) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, nb)
        if (part(nb) == target_rank .and. part(o) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, o)
      end do
    else
      do f0 = 1, nf, face_chunk
        f1 = min(nf, f0 + face_chunk - 1)
        do f = f0, f1
          o = face_owner(f)
          nb = face_nei(f)
          if (o < 1) cycle
          if (nb >= 1) then
            if (part(o) == target_rank .and. part(nb) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, nb)
            if (part(nb) == target_rank .and. part(o) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, o)
          end if
        end do
      end do
    end if

    do g_bnd = 1, size(periodic_pair_global)
      g_pair = periodic_pair_global(g_bnd)
      if (g_pair <= 0) cycle
      f = n_int_global + g_bnd
      pair_face = n_int_global + g_pair
      if (f < 1 .or. f > nf .or. pair_face < 1 .or. pair_face > nf) cycle
      o = face_owner(f)
      nb_global = face_owner(pair_face)
      if (o < 1 .or. nb_global < 1) cycle
      if (part(o) == target_rank .and. part(nb_global) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, nb_global)
      if (part(nb_global) == target_rank .and. part(o) /= target_rank) call append_dynamic_int(ghost_candidates, n_ghost_candidates, o)
    end do

    if (n_ghost_candidates > 0) then
      call sort_int_values(ghost_candidates, n_ghost_candidates)
      n_ghost = unique_sorted_int(ghost_candidates, n_ghost_candidates)
    else
      n_ghost = 0
    end if

    allocate(ghost_ids(max(1_ip, n_ghost)))
    if (n_ghost > 0) ghost_ids(1:n_ghost) = ghost_candidates(1:n_ghost)
    if (allocated(ghost_candidates)) deallocate(ghost_candidates)

    call build_gid_to_lid_lookup(owned_ids, n_owned, ghost_ids, n_ghost, map_gid, map_lid)

    if (present(int_face_list)) then
      n_int_local = size(int_face_list)
    else
      n_int_local = 0
      do f0 = 1, nf, face_chunk
        f1 = min(nf, f0 + face_chunk - 1)
        do f = f0, f1
          o = face_owner(f)
          nb = face_nei(f)
          if (nb < 1) cycle
          if (part(o) == target_rank .or. part(nb) == target_rank) n_int_local = n_int_local + 1
        end do
      end do
    end if

    do p = 1, np
      do i = 0, patch_nfaces(p) - 1
        f_bnd = patch_start(p) + i
        f = n_int_global + f_bnd
        if (f < 1 .or. f > nf) cycle
        o = face_owner(f)
        nb = face_nei(f)
        if (nb >= 1) cycle
        if (o >= 1 .and. part(o) == target_rank) patch_local_count(p) = patch_local_count(p) + 1
      end do
    end do

    n_bnd_local = sum(patch_local_count)
    n_faces_local = n_int_local + n_bnd_local

    msh%comm = comm
    msh%rank = target_rank
    msh%nranks = nranks
    msh%n_cells = n_owned
    msh%n_ghost = n_ghost
    msh%n_faces = n_faces_local
    msh%n_int_faces = n_int_local
    msh%n_bnd_faces = n_bnd_local
    msh%n_patches = np
    msh%n_vertices = 0
    msh%n_cells_global = int(nc, lp)

    allocate(msh%cells(n_owned + n_ghost))
    allocate(msh%faces(n_faces_local))
    allocate(msh%patches(np))
    allocate(msh%global_id(n_owned + n_ghost))
    allocate(local_face_global(max(1_ip, n_faces_local)), source=0_ip)

    do i = 1, n_owned
      msh%cells(i)%centroid = cell_centroid(:, owned_ids(i))
      msh%cells(i)%volume   = cell_volume(owned_ids(i))
      msh%global_id(i)      = int(owned_ids(i), lp)
    end do
    do i = 1, n_ghost
      msh%cells(n_owned + i)%centroid = cell_centroid(:, ghost_ids(i))
      msh%cells(n_owned + i)%volume   = cell_volume(ghost_ids(i))
      msh%global_id(n_owned + i)      = int(ghost_ids(i), lp)
    end do

    if (topo_flag == 1 .and. n_vertices_global > 0) then
      ! Gather only the vertices this rank touches, then sort-unique. Ascending global
      ! order reproduces the old dense-scan numbering exactly, at O(n_owned) not O(n_vertices_global).
      n_vlist = 0
      do i = 1, n_owned
        gid = owned_ids(i)
        n_vlist = n_vlist + cell_offsets_global(gid + 1) - cell_offsets_global(gid)
      end do
      allocate(vlist(max(1_ip, n_vlist)))
      vpos = 0
      do i = 1, n_owned
        gid = owned_ids(i)
        n0 = cell_offsets_global(gid)
        n1 = cell_offsets_global(gid + 1) - 1
        do v = n0, n1
          vpos = vpos + 1
          vlist(vpos) = cell_vertices_global(v)
        end do
      end do
      call sort_int_values(vlist, vpos)
      vlocal = unique_sorted_int(vlist, vpos)

      msh%n_vertices = vlocal
      allocate(msh%vertices(3, msh%n_vertices))
      do v = 1, vlocal
        msh%vertices(:, v) = vertices_global(:, vlist(v))
      end do

      allocate(msh%cell_types(n_owned), msh%cell_offsets(n_owned + 1))
      msh%cell_offsets(1) = 1
      do i = 1, n_owned
        gid = owned_ids(i)
        nconn = cell_offsets_global(gid + 1) - cell_offsets_global(gid)
        msh%cell_types(i) = cell_types_global(gid)
        msh%cell_offsets(i + 1) = msh%cell_offsets(i) + nconn
      end do

      allocate(msh%cell_vertices(msh%cell_offsets(n_owned + 1) - 1))
      pos = 1
      do i = 1, n_owned
        gid = owned_ids(i)
        n0 = cell_offsets_global(gid)
        n1 = cell_offsets_global(gid + 1) - 1
        nconn = n1 - n0 + 1
        do v = n0, n1
          msh%cell_vertices(pos + v - n0) = bsearch_sorted_int(vlist, vlocal, cell_vertices_global(v))
        end do
        pos = pos + nconn
      end do

    end if

    cursor = 0
    if (present(int_face_list)) then
      do i = 1, size(int_face_list)
        f = int_face_list(i)
        o = face_owner(f)
        nb = face_nei(f)
        if (part(o) == target_rank) then
          lo = g2l_lookup_gid(o, map_gid, map_lid)
          ln = g2l_lookup_gid(nb, map_gid, map_lid)
          flip = .false.
        else
          lo = g2l_lookup_gid(nb, map_gid, map_lid)
          ln = g2l_lookup_gid(o, map_gid, map_lid)
          flip = .true.
        end if
        if (lo <= 0 .or. ln <= 0) cycle

        cursor = cursor + 1
        local_face_global(cursor) = f
        call fetch_face_geom(f, nf, face_centre, face_normal, face_sf, face_gc, &
                 staged_geom_unit, staged_geom_block, staged_cache_start, staged_cache_end, &
                 centre_local, normal_local, sf_out=sf_local, gc_out=gc_local, has_sf_gc=have_sf_gc)
        msh%faces(cursor)%owner = lo
        msh%faces(cursor)%neighbour = ln
        msh%faces(cursor)%patch_id = -1
        msh%faces(cursor)%area = face_area(f)
        msh%faces(cursor)%centre = centre_local
        if (.not. have_sf_gc) then
          ! Fallback: orthogonal reconstruction. This loses non-orthogonal
          ! correction on skewed meshes - see -04-28.md F1.
          sf_local = face_area(f) * normal_local
        end if
        if (flip) then
          msh%faces(cursor)%normal = -normal_local
          msh%faces(cursor)%sf     = -sf_local
        else
          msh%faces(cursor)%normal = normal_local
          msh%faces(cursor)%sf     = sf_local
        end if
        msh%faces(cursor)%d_cf = face_dcf(f)
        if (.not. have_sf_gc) then
          d_on = cell_centroid(:, nb) - cell_centroid(:, o)
          d_len = sqrt(sum(d_on**2)) + SMALL
          d_oc = centre_local - cell_centroid(:, o)
          ! Projection form, matching mesh_compute_geometry; |d_oc|/d_len differs on skewed faces.
          gc_local = max(0.0_rp, min(1.0_rp, 1.0_rp - dot_product(d_oc, d_on) / (d_len * d_len)))
        end if
        ! gc above is computed/stored from the GLOBAL owner perspective
        ! (weight of the global owner cell). When the local face owner has
        if (flip) gc_local = 1.0_rp - gc_local
        msh%faces(cursor)%gc = gc_local
      end do
    else
      do f0 = 1, nf, face_chunk
        f1 = min(nf, f0 + face_chunk - 1)
        do f = f0, f1
          o = face_owner(f)
          nb = face_nei(f)
          if (o < 1 .or. nb < 1) cycle
          if (part(o) /= target_rank .and. part(nb) /= target_rank) cycle

          if (part(o) == target_rank) then
            lo = g2l_lookup_gid(o, map_gid, map_lid)
            ln = g2l_lookup_gid(nb, map_gid, map_lid)
            flip = .false.
          else
            lo = g2l_lookup_gid(nb, map_gid, map_lid)
            ln = g2l_lookup_gid(o, map_gid, map_lid)
            flip = .true.
          end if
          if (lo <= 0 .or. ln <= 0) cycle

          cursor = cursor + 1
          local_face_global(cursor) = f
          call fetch_face_geom(f, nf, face_centre, face_normal, face_sf, face_gc, &
                   staged_geom_unit, staged_geom_block, staged_cache_start, staged_cache_end, &
                   centre_local, normal_local, sf_out=sf_local, gc_out=gc_local, has_sf_gc=have_sf_gc)
          msh%faces(cursor)%owner = lo
          msh%faces(cursor)%neighbour = ln
          msh%faces(cursor)%patch_id = -1
          msh%faces(cursor)%area = face_area(f)
          msh%faces(cursor)%centre = centre_local
          if (.not. have_sf_gc) then
            ! Fallback: orthogonal reconstruction. This loses non-orthogonal
            ! correction on skewed meshes - see -04-28.md F1.
            sf_local = face_area(f) * normal_local
          end if
          if (flip) then
            msh%faces(cursor)%normal = -normal_local
            msh%faces(cursor)%sf     = -sf_local
          else
            msh%faces(cursor)%normal = normal_local
            msh%faces(cursor)%sf     = sf_local
          end if
          msh%faces(cursor)%d_cf = face_dcf(f)
          if (.not. have_sf_gc) then
            d_on = cell_centroid(:, nb) - cell_centroid(:, o)
            d_len = sqrt(sum(d_on**2)) + SMALL
            d_oc = centre_local - cell_centroid(:, o)
            ! Projection form, matching mesh_compute_geometry; |d_oc|/d_len differs on skewed faces.
            gc_local = max(0.0_rp, min(1.0_rp, 1.0_rp - dot_product(d_oc, d_on) / (d_len * d_len)))
          end if
          ! gc above is computed/stored from the GLOBAL owner perspective
          ! (weight of the global owner cell). When the local face owner has
          if (flip) gc_local = 1.0_rp - gc_local
          msh%faces(cursor)%gc = gc_local
        end do
      end do
    end if

    ! Boundary faces are appended patch-by-patch so patch ranges are contiguous.
    allocate(bnd_order(max(1_ip, n_bnd_local)))
    bcursor = 0
    do p = 1, np
      do i = 0, patch_nfaces(p) - 1
        f_bnd = patch_start(p) + i
        f = n_int_global + f_bnd
        if (f < 1 .or. f > nf) cycle
        o = face_owner(f)
        nb = face_nei(f)
        if (nb >= 1) cycle
        if (o < 1 .or. part(o) /= target_rank) cycle
        bcursor = bcursor + 1
        bnd_order(bcursor) = f
      end do
    end do

    do i = 1, n_bnd_local
      f = bnd_order(i)
      cursor = cursor + 1
      local_face_global(cursor) = f
      call fetch_face_geom(f, nf, face_centre, face_normal, face_sf, face_gc, &
               staged_geom_unit, staged_geom_block, staged_cache_start, staged_cache_end, &
               centre_local, normal_local, sf_out=sf_local, gc_out=gc_local, has_sf_gc=have_sf_gc)
      msh%faces(cursor)%owner = g2l_lookup_gid(face_owner(f), map_gid, map_lid)
      msh%faces(cursor)%neighbour = -1
      msh%faces(cursor)%patch_id = face_patch(f)
      msh%faces(cursor)%area = face_area(f)
      msh%faces(cursor)%centre = centre_local
      msh%faces(cursor)%normal = normal_local
      if (have_sf_gc) then
        msh%faces(cursor)%sf = sf_local
        msh%faces(cursor)%gc = gc_local
      else
        msh%faces(cursor)%sf = face_area(f) * normal_local
        msh%faces(cursor)%gc = 1.0_rp
      end if
      msh%faces(cursor)%d_cf = face_dcf(f)
    end do

    if (topo_flag == 1 .and. face_topo_flag == 1 .and. msh%n_vertices > 0 .and. allocated(vlist)) then
      allocate(msh%face_offsets(n_faces_local + 1))
      msh%face_offsets(1) = 1
      do i = 1, n_faces_local
        f = local_face_global(i)
        g0 = face_offsets_global(f)
        g1 = face_offsets_global(f + 1) - 1
        nfv = g1 - g0 + 1
        msh%face_offsets(i + 1) = msh%face_offsets(i) + nfv
      end do

      allocate(msh%face_vertices(msh%face_offsets(n_faces_local + 1) - 1))
      pos = 1
      do i = 1, n_faces_local
        f = local_face_global(i)
        g0 = face_offsets_global(f)
        g1 = face_offsets_global(f + 1) - 1
        nfv = g1 - g0 + 1

        ! Global face vertex winding is owner->neighbour. During partitioning,
        ! local face owner may be flipped so owner is always local. Keep face
        if (msh%faces(i)%owner == g2l_lookup_gid(face_owner(f), map_gid, map_lid)) then
          do v = g0, g1
            msh%face_vertices(pos + v - g0) = bsearch_sorted_int(vlist, vlocal, face_vertices_global(v))
          end do
        else
          do v = g0, g1
            msh%face_vertices(pos + g1 - v) = bsearch_sorted_int(vlist, vlocal, face_vertices_global(v))
          end do
        end if

        pos = pos + nfv
      end do
    end if

    if (allocated(vlist)) deallocate(vlist)

    bcursor = 1
    do p = 1, np
      msh%patches(p)%name = patch_name(p)
      msh%patches(p)%bc_type = patch_type(p)
      msh%patches(p)%bc_id = patch_bcid(p)
      msh%patches(p)%start_face = bcursor
      msh%patches(p)%n_faces = patch_local_count(p)
      bcursor = bcursor + patch_local_count(p)
    end do

    allocate(msh%periodic_pair(msh%n_bnd_faces), source=0_ip)
    allocate(msh%periodic_neighbor(msh%n_bnd_faces), source=0_ip)
    allocate(msh%periodic_pair_dcf(msh%n_bnd_faces), source=0.0_rp)
    allocate(g2lbnd(max(1_ip, size(periodic_pair_global))), source=0_ip)

    do i = 1, n_bnd_local
      g_bnd = bnd_order(i) - n_int_global
      if (g_bnd >= 1 .and. g_bnd <= size(g2lbnd)) g2lbnd(g_bnd) = i
    end do

    do i = 1, n_bnd_local
      g_bnd = bnd_order(i) - n_int_global
      if (g_bnd < 1 .or. g_bnd > size(periodic_pair_global)) cycle

      g_pair = periodic_pair_global(g_bnd)
      if (g_pair <= 0) cycle

      if (g_pair >= 1 .and. g_pair <= size(g2lbnd)) msh%periodic_pair(i) = g2lbnd(g_pair)

      pair_face = n_int_global + g_pair
      if (pair_face < 1 .or. pair_face > nf) cycle
      msh%periodic_pair_dcf(i) = face_dcf(pair_face)
      nb_global = face_owner(pair_face)
      if (nb_global < 1 .or. nb_global > nc) cycle
      msh%periodic_neighbor(i) = g2l_lookup_gid(nb_global, map_gid, map_lid)
    end do

    call build_halos_from_partition(msh, target_rank, nranks, nc, part, face_owner, face_nei, map_gid, map_lid, &
                                    periodic_pair_global, n_int_global, int_face_list=int_face_list)

    call remap_periodic_neighbors_to_ghosts(msh, n_bnd_local, bnd_order, &
                                            periodic_pair_global, face_owner, &
                                            n_int_global, nf, nc)

    deallocate(owned_ids, ghost_ids, g2lbnd, patch_local_count, bnd_order, local_face_global)
    if (allocated(map_gid)) deallocate(map_gid)
    if (allocated(map_lid)) deallocate(map_lid)
    if (allocated(staged_geom_block)) deallocate(staged_geom_block)
  end subroutine build_local_partition_mesh

  !> Remap unresolved periodic neighbors after halo build and before deallocating bnd_order/face_owner.
  subroutine remap_periodic_neighbors_to_ghosts(msh, n_bnd_local, bnd_order, &
                                                periodic_pair_global, face_owner, &
                                                n_int_global, nf, nc)
    type(Mesh),  intent(inout) :: msh
    integer(ip), intent(in)    :: n_bnd_local, n_int_global, nf, nc
    integer(ip), intent(in)    :: bnd_order(:)
    integer(ip), intent(in)    :: periodic_pair_global(:)
    integer(ip), intent(in)    :: face_owner(:)

    integer(ip) :: i, j, g_bnd, g_pair, pair_face, nb_global
    integer(lp) :: nb_global_lp
    integer(ip) :: n_remapped

    if (.not. allocated(msh%periodic_neighbor)) return
    if (msh%n_ghost <= 0) return

    n_remapped = 0
    do i = 1, n_bnd_local
      if (msh%periodic_neighbor(i) > 0) cycle  ! already resolved

      g_bnd = bnd_order(i) - n_int_global
      if (g_bnd < 1 .or. g_bnd > size(periodic_pair_global)) cycle
      g_pair = periodic_pair_global(g_bnd)
      if (g_pair <= 0) cycle

      pair_face = n_int_global + g_pair
      if (pair_face < 1 .or. pair_face > nf) cycle
      nb_global = face_owner(pair_face)
      if (nb_global < 1 .or. nb_global > nc) cycle

      nb_global_lp = int(nb_global, lp)
      do j = msh%n_cells + 1, msh%n_cells + msh%n_ghost
        if (msh%global_id(j) == nb_global_lp) then
          msh%periodic_neighbor(i) = j
          n_remapped = n_remapped + 1
          exit
        end if
      end do
    end do

    if (n_remapped > 0) then
      write(*,'(A,I0,A,I0,A)') ' [remap_periodic] rank ', msh%rank, &
        ': remapped ', n_remapped, ' cross-rank periodic neighbors'
    end if
  end subroutine remap_periodic_neighbors_to_ghosts

  !> Fetch per-face centre, normal, sf, and gc - either from in-memory arrays
  !> or, when the streaming decomposer staged geometry on disk, by reading a
  subroutine fetch_face_geom(f, nf, face_centre, face_normal, face_sf_in, face_gc_in, &
                             staged_geom_unit, staged_geom_block, &
                             cache_start, cache_end, centre, normal, sf_out, gc_out, has_sf_gc)
    integer(ip), intent(in) :: f, nf
    real(rp), intent(in), optional :: face_centre(:,:), face_normal(:,:)
    real(rp), intent(in), optional :: face_sf_in(:,:), face_gc_in(:)
    integer(ip), intent(in), optional :: staged_geom_unit
    real(rp), intent(inout), optional :: staged_geom_block(:,:)
    integer(ip), intent(inout), optional :: cache_start, cache_end
    real(rp), intent(out) :: centre(3), normal(3)
    real(rp), intent(out), optional :: sf_out(3), gc_out
    logical,  intent(out), optional :: has_sf_gc

    integer(lp) :: pos_bytes
    integer(ip) :: rbytes, nload, idx
    integer(ip) :: nrec  ! number of reals per staged record

    if (present(has_sf_gc)) has_sf_gc = .false.

    if (present(face_centre) .and. present(face_normal)) then
      centre = face_centre(:, f)
      normal = face_normal(:, f)
      if (present(sf_out) .and. present(face_sf_in)) sf_out = face_sf_in(:, f)
      if (present(gc_out) .and. present(face_gc_in)) gc_out = face_gc_in(f)
      if (present(has_sf_gc)) has_sf_gc = present(face_sf_in) .and. present(face_gc_in)
      return
    end if

    if (.not. present(staged_geom_unit)) then
      error stop 'fetch_face_geom: neither in-memory nor staged geometry available'
    end if

    if (.not. present(staged_geom_block) .or. .not. present(cache_start) .or. .not. present(cache_end)) then
      error stop 'fetch_face_geom: staged geometry cache not provided'
    end if

    nrec = size(staged_geom_block, 1)  ! 6 (legacy) or 10 (with sf+gc)

    if (f < cache_start .or. f > cache_end) then
      cache_start = f
      nload = min(size(staged_geom_block, 2), nf - cache_start + 1)
      cache_end = cache_start + nload - 1
      rbytes = storage_size(0.0_rp) / 8
      pos_bytes = int(cache_start - 1_ip, lp) * int(nrec * rbytes, lp) + 1_lp
      read(staged_geom_unit, pos=pos_bytes) staged_geom_block(:, 1:nload)
    end if

    idx = f - cache_start + 1
    centre = staged_geom_block(1:3, idx)
    normal = staged_geom_block(4:6, idx)
    if (nrec >= 10) then
      if (present(sf_out)) sf_out = staged_geom_block(7:9,  idx)
      if (present(gc_out)) gc_out = staged_geom_block(10,   idx)
      if (present(has_sf_gc)) has_sf_gc = .true.
    end if
  end subroutine fetch_face_geom

  subroutine build_halos_from_partition(msh, rank, nranks, nc, part, face_owner, face_nei, map_gid, map_lid, &
                                        periodic_pair_global, n_int_global, int_face_list)
    type(Mesh), intent(inout) :: msh
    integer(ip), intent(in) :: rank, nranks, nc
    integer(ip), intent(in) :: part(nc), face_owner(:), face_nei(:)
    integer(ip), intent(in) :: map_gid(:), map_lid(:)
    integer(ip), intent(in) :: periodic_pair_global(:), n_int_global
    integer(ip), intent(in), optional :: int_face_list(:)

    integer(ip) :: r, f, o, nb, nsend, nrecv, nh, ihf
    integer(ip) :: f_bnd, pair_bnd, pair_face
    integer(ip), allocatable :: send_candidates(:), recv_candidates(:)
    integer(ip) :: nsend_candidates, nrecv_candidates
    type(HaloData), allocatable :: tmp(:)

    allocate(tmp(max(1_ip, nranks-1)))
    nh = 0

    do r = 0, nranks - 1
      if (r == rank) cycle

      nsend_candidates = 0
      nrecv_candidates = 0
      call init_dynamic_int_buffer(send_candidates, 1024_ip)
      call init_dynamic_int_buffer(recv_candidates, 1024_ip)

      if (present(int_face_list)) then
        do ihf = 1, size(int_face_list)
          f = int_face_list(ihf)
          o = face_owner(f)
          nb = face_nei(f)
          if (part(o) == rank .and. part(nb) == r) then
            call append_dynamic_int(send_candidates, nsend_candidates, o)
            call append_dynamic_int(recv_candidates, nrecv_candidates, nb)
          else if (part(nb) == rank .and. part(o) == r) then
            call append_dynamic_int(send_candidates, nsend_candidates, nb)
            call append_dynamic_int(recv_candidates, nrecv_candidates, o)
          end if
        end do
      else
        do f = 1, size(face_owner)
          o = face_owner(f)
          nb = face_nei(f)
          if (o < 1 .or. nb < 1) cycle

          if (part(o) == rank .and. part(nb) == r) then
            call append_dynamic_int(send_candidates, nsend_candidates, o)
            call append_dynamic_int(recv_candidates, nrecv_candidates, nb)
          else if (part(nb) == rank .and. part(o) == r) then
            call append_dynamic_int(send_candidates, nsend_candidates, nb)
            call append_dynamic_int(recv_candidates, nrecv_candidates, o)
          end if
        end do
      end if

      do f_bnd = 1, size(periodic_pair_global)
        pair_bnd = periodic_pair_global(f_bnd)
        if (pair_bnd <= 0) cycle
        f = n_int_global + f_bnd
        pair_face = n_int_global + pair_bnd
        if (f < 1 .or. f > size(face_owner) .or. pair_face < 1 .or. pair_face > size(face_owner)) cycle
        o = face_owner(f)
        nb = face_owner(pair_face)
        if (o < 1 .or. nb < 1) cycle

        if (part(o) == rank .and. part(nb) == r) then
          call append_dynamic_int(send_candidates, nsend_candidates, o)
          call append_dynamic_int(recv_candidates, nrecv_candidates, nb)
        else if (part(nb) == rank .and. part(o) == r) then
          call append_dynamic_int(send_candidates, nsend_candidates, nb)
          call append_dynamic_int(recv_candidates, nrecv_candidates, o)
        end if
      end do

      if (nsend_candidates > 0) then
        call sort_int_values(send_candidates, nsend_candidates)
        nsend = unique_sorted_int(send_candidates, nsend_candidates)
      else
        nsend = 0
      end if
      if (nrecv_candidates > 0) then
        call sort_int_values(recv_candidates, nrecv_candidates)
        nrecv = unique_sorted_int(recv_candidates, nrecv_candidates)
      else
        nrecv = 0
      end if

      if (nsend > 0 .or. nrecv > 0) then
        nh = nh + 1
        tmp(nh)%rank = r
        tmp(nh)%n_send = nsend
        tmp(nh)%n_recv = nrecv
        allocate(tmp(nh)%send_cells(max(1_ip, nsend)))
        allocate(tmp(nh)%recv_cells(max(1_ip, nrecv)))
        call fill_local_indices_from_gids(send_candidates, nsend, map_gid, map_lid, tmp(nh)%send_cells)
        call fill_local_indices_from_gids(recv_candidates, nrecv, map_gid, map_lid, tmp(nh)%recv_cells)
      end if

      if (allocated(send_candidates)) deallocate(send_candidates)
      if (allocated(recv_candidates)) deallocate(recv_candidates)
    end do

    if (allocated(msh%halos)) deallocate(msh%halos)
    if (nh > 0) then
      allocate(msh%halos(nh))
      msh%halos = tmp(1:nh)
    else
      allocate(msh%halos(0))
    end if

    deallocate(tmp)
  end subroutine build_halos_from_partition

  subroutine init_dynamic_int_buffer(buf, initial_capacity)
    integer(ip), allocatable, intent(inout) :: buf(:)
    integer(ip), intent(in) :: initial_capacity

    integer(ip) :: cap

    cap = max(1_ip, initial_capacity)
    if (allocated(buf)) deallocate(buf)
    allocate(buf(cap))
  end subroutine init_dynamic_int_buffer

  subroutine append_dynamic_int(buf, n, value)
    integer(ip), allocatable, intent(inout) :: buf(:)
    integer(ip), intent(inout) :: n
    integer(ip), intent(in) :: value

    integer(ip), allocatable :: tmp(:)
    integer(ip) :: new_cap

    if (.not. allocated(buf)) then
      allocate(buf(1024))
    end if
    if (n >= size(buf)) then
      new_cap = max(2_ip * size(buf), size(buf) + 1024_ip)
      allocate(tmp(new_cap))
      if (n > 0) tmp(1:n) = buf(1:n)
      call move_alloc(tmp, buf)
    end if
    n = n + 1
    buf(n) = value
  end subroutine append_dynamic_int

  recursive subroutine quicksort_int(arr, left, right)
    integer(ip), intent(inout) :: arr(:)
    integer(ip), intent(in) :: left, right

    integer(ip) :: i, j, pivot, tmp

    if (left >= right) return

    i = left
    j = right
    pivot = arr((left + right) / 2)

    do while (i <= j)
      do while (arr(i) < pivot)
        i = i + 1
      end do
      do while (arr(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = arr(i)
        arr(i) = arr(j)
        arr(j) = tmp
        i = i + 1
        j = j - 1
      end if
    end do

    if (left < j) call quicksort_int(arr, left, j)
    if (i < right) call quicksort_int(arr, i, right)
  end subroutine quicksort_int

  subroutine sort_int_values(arr, n)
    integer(ip), intent(inout) :: arr(:)
    integer(ip), intent(in) :: n

    if (n <= 1) return
    call quicksort_int(arr, 1, n)
  end subroutine sort_int_values

  integer(ip) function unique_sorted_int(arr, n) result(n_unique)
    integer(ip), intent(inout) :: arr(:)
    integer(ip), intent(in) :: n

    integer(ip) :: i

    if (n <= 0) then
      n_unique = 0
      return
    end if

    n_unique = 1
    do i = 2, n
      if (arr(i) /= arr(n_unique)) then
        n_unique = n_unique + 1
        arr(n_unique) = arr(i)
      end if
    end do
  end function unique_sorted_int

  !> Position of gid in a sorted unique array; 0 if absent.
  pure integer(ip) function bsearch_sorted_int(arr, n, gid) result(pos)
    integer(ip), intent(in) :: arr(:), n, gid
    integer(ip) :: lo, hi, mid

    pos = 0
    lo = 1
    hi = n
    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (arr(mid) == gid) then
        pos = mid
        return
      else if (arr(mid) < gid) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function bsearch_sorted_int

  recursive subroutine quicksort_gid_lid(gid, lid, left, right)
    integer(ip), intent(inout) :: gid(:), lid(:)
    integer(ip), intent(in) :: left, right

    integer(ip) :: i, j, pivot, tmp

    if (left >= right) return

    i = left
    j = right
    pivot = gid((left + right) / 2)

    do while (i <= j)
      do while (gid(i) < pivot)
        i = i + 1
      end do
      do while (gid(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = gid(i); gid(i) = gid(j); gid(j) = tmp
        tmp = lid(i); lid(i) = lid(j); lid(j) = tmp
        i = i + 1
        j = j - 1
      end if
    end do

    if (left < j) call quicksort_gid_lid(gid, lid, left, j)
    if (i < right) call quicksort_gid_lid(gid, lid, i, right)
  end subroutine quicksort_gid_lid

  subroutine build_gid_to_lid_lookup(owned_ids, n_owned, ghost_ids, n_ghost, map_gid, map_lid)
    integer(ip), intent(in) :: owned_ids(:), ghost_ids(:)
    integer(ip), intent(in) :: n_owned, n_ghost
    integer(ip), allocatable, intent(out) :: map_gid(:), map_lid(:)

    integer(ip) :: i, n_total

    n_total = n_owned + n_ghost
    allocate(map_gid(max(1_ip, n_total)), map_lid(max(1_ip, n_total)))
    map_gid = huge(1_ip)
    map_lid = 0

    do i = 1, n_owned
      map_gid(i) = owned_ids(i)
      map_lid(i) = i
    end do
    do i = 1, n_ghost
      map_gid(n_owned + i) = ghost_ids(i)
      map_lid(n_owned + i) = n_owned + i
    end do

    if (n_total > 1) call quicksort_gid_lid(map_gid, map_lid, 1, n_total)
  end subroutine build_gid_to_lid_lookup

  integer(ip) function g2l_lookup_gid(gid, map_gid, map_lid) result(lid)
    integer(ip), intent(in) :: gid
    integer(ip), intent(in) :: map_gid(:), map_lid(:)

    integer(ip) :: lo, hi, mid

    lid = 0
    lo = 1
    hi = size(map_gid)
    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (map_gid(mid) == gid) then
        lid = map_lid(mid)
        return
      else if (map_gid(mid) < gid) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function g2l_lookup_gid

  subroutine fill_local_indices_from_gids(gids, n, map_gid, map_lid, out_idx)
    integer(ip), intent(in) :: gids(:)
    integer(ip), intent(in) :: n
    integer(ip), intent(in) :: map_gid(:), map_lid(:)
    integer(ip), intent(out) :: out_idx(:)

    integer(ip) :: i

    if (n <= 0) return
    do i = 1, n
      out_idx(i) = g2l_lookup_gid(gids(i), map_gid, map_lid)
    end do
  end subroutine fill_local_indices_from_gids

  !> Quicksort cell indices by partition value for efficient rank-contiguous processing.
  !> After sorting, all cells belonging to rank r are contiguous in cell_order array.
  recursive subroutine quicksort_cells_by_partition(cell_order, part, left, right)
    integer(ip), intent(inout) :: cell_order(:)
    integer(ip), intent(in)    :: part(:)
    integer(ip), intent(in)    :: left, right
    integer(ip) :: pivot, i, j, tmp

    if (left >= right) return

    ! Partition: pivot on middle cell's partition value
    i = left
    j = right
    pivot = part(cell_order((left + right) / 2))

    do while (i <= j)
      do while (part(cell_order(i)) < pivot)
        i = i + 1
      end do
      do while (part(cell_order(j)) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = cell_order(i)
        cell_order(i) = cell_order(j)
        cell_order(j) = tmp
        i = i + 1
        j = j - 1
      end if
    end do

    if (left < j) call quicksort_cells_by_partition(cell_order, part, left, j)
    if (i < right) call quicksort_cells_by_partition(cell_order, part, i, right)
  end subroutine quicksort_cells_by_partition

  !> Compute rank offsets from sorted partition array.
  !> rank_offsets(r) = first index in cell_order where part(cell_order(i)) == r
  subroutine compute_rank_offsets(cell_order, part, nc, nranks, rank_offsets)
    integer(ip), intent(in)  :: cell_order(:), part(:), nc, nranks
    integer(ip), intent(out) :: rank_offsets(0:nranks)
    integer(ip) :: i, r

    rank_offsets = nc + 1  ! Initialize beyond array bounds

    do i = 1, nc
      r = part(cell_order(i))
      if (rank_offsets(r) > nc) then
        rank_offsets(r) = i
      end if
    end do

    ! Last rank's offset should point to nc+1
    do r = 0, nranks - 1
      if (rank_offsets(r) == nc + 1) then
        rank_offsets(r) = rank_offsets(r + 1)  ! Share next rank's start if empty
      end if
    end do
    rank_offsets(nranks) = nc + 1
  end subroutine compute_rank_offsets

  !> Thin wrapper: sort cells by partition and compute rank offsets.
  !> Calls quicksort_cells_by_partition + compute_rank_offsets in sequence.
  subroutine sort_and_offset_partition(nc, nranks, part, cell_order, rank_offsets)
    integer(ip), intent(in)  :: nc, nranks
    integer(ip), intent(in)  :: part(:)
    integer(ip), intent(out) :: cell_order(:), rank_offsets(0:nranks)
    integer(ip) :: i

    ! Initialize cell_order: cell_order(i) = i
    do i = 1, nc
      cell_order(i) = i
    end do

    ! Sort by partition
    call quicksort_cells_by_partition(cell_order, part, 1, nc)

    ! Compute offsets
    call compute_rank_offsets(cell_order, part, nc, nranks, rank_offsets)
  end subroutine sort_and_offset_partition

  !> Hilbert space-filling curve partitioning with memory-efficient local mesh distribution.
  !> Only rank 0 builds/holds global mesh. Non-zero ranks receive their local mesh via MPI.
  subroutine partition_mesh_hilbert(mesh_global, mesh_local, comm)
    type(Mesh), intent(inout) :: mesh_global, mesh_local
    type(MPI_Comm), intent(in) :: comm

    integer(ip) :: rank, nranks, ierr, r, i
    integer(ip) :: nc, nf, np, ng, nintf, nbndf, nv
    integer(ip) :: topo_flag, total_conn, face_topo_flag, total_face_conn
    integer(ip), allocatable :: part(:), cell_order(:), rank_offsets(:)
    type(Mesh) :: tmp_mesh

    ! Broadcast-packed global mesh arrays (only on rank 0)
    integer(ip), allocatable :: face_owner(:), face_nei(:), face_patch(:)
    integer(ip), allocatable :: patch_start(:), patch_nfaces(:), patch_bcid(:)
    integer(ip), allocatable :: periodic_pair_global(:)
    character(len=128), allocatable :: patch_name(:)
    character(len=64),  allocatable :: patch_type(:)
    real(rp), allocatable :: cell_centroid(:,:), cell_volume(:)
    real(rp), allocatable :: face_area(:), face_centre(:,:), face_normal(:,:), face_sf(:,:), face_dcf(:), face_gc(:)
    real(rp), allocatable :: vertices(:,:)
    integer(ip), allocatable :: cell_offsets(:), cell_vertices(:), cell_types(:)
    integer(ip), allocatable :: face_offsets(:), face_vertices(:)

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nranks, ierr)

    ! --- Broadcast scalar dimensions to all ranks ---
    if (rank == 0) then
      nc    = mesh_global%n_cells
      ng    = mesh_global%n_ghost
      nf    = mesh_global%n_faces
      nintf = mesh_global%n_int_faces
      nbndf = mesh_global%n_bnd_faces
      np    = mesh_global%n_patches
      nv    = mesh_global%n_vertices
      topo_flag = 0
      total_conn = 0
      face_topo_flag = 0
      total_face_conn = 0
      if (allocated(mesh_global%vertices) .and. allocated(mesh_global%cell_offsets) .and. &
          allocated(mesh_global%cell_vertices) .and. allocated(mesh_global%cell_types)) then
        topo_flag = 1
        total_conn = mesh_global%cell_offsets(nc + 1) - 1
      end if
      if (allocated(mesh_global%face_offsets) .and. allocated(mesh_global%face_vertices)) then
        face_topo_flag = 1
        total_face_conn = mesh_global%face_offsets(nf + 1) - 1
      end if
    end if

    call MPI_Bcast(nc,    1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(ng,    1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(nf,    1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(nintf, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(nbndf, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(np,    1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(nv,    1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(topo_flag, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(total_conn, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(face_topo_flag, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Bcast(total_face_conn, 1, MPI_INTEGER, 0, comm, ierr)

    allocate(periodic_pair_global(nbndf), source=0_ip)
    if (rank == 0) then
      if (allocated(mesh_global%periodic_pair)) then
        periodic_pair_global = mesh_global%periodic_pair
      end if
    end if
    call MPI_Bcast(periodic_pair_global, nbndf, MPI_INTEGER, 0, comm, ierr)

    if (topo_flag == 1) then
      allocate(vertices(3, nv), cell_offsets(nc + 1), cell_vertices(total_conn), cell_types(nc))
      if (rank == 0) then
        vertices = mesh_global%vertices
        cell_offsets = mesh_global%cell_offsets
        cell_vertices = mesh_global%cell_vertices
        cell_types = mesh_global%cell_types
      end if
      call MPI_Bcast(vertices, 3*nv, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(cell_offsets, nc + 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(cell_vertices, total_conn, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(cell_types, nc, MPI_INTEGER, 0, comm, ierr)
    end if

    if (face_topo_flag == 1) then
      allocate(face_offsets(nf + 1), face_vertices(total_face_conn))
      if (rank == 0) then
        face_offsets = mesh_global%face_offsets
        face_vertices = mesh_global%face_vertices
      end if
      call MPI_Bcast(face_offsets, nf + 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(face_vertices, total_face_conn, MPI_INTEGER, 0, comm, ierr)
    end if

    ! --- Broadcast cell centroids only (needed for Hilbert SFC) ---
    allocate(cell_centroid(3, nc))
    if (rank == 0) then
      do i = 1, nc
        cell_centroid(:, i) = mesh_global%cells(i)%centroid
      end do
    end if

    call MPI_Bcast(cell_centroid, 3*nc, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    if (rank == 0) write(*, '(A)') ' [HILBERT] Computing partition from centroids...'

    ! --- All ranks compute Hilbert partition ---
    allocate(part(nc))
    call build_hilbert_partition(nc, nranks, rank, cell_centroid, part)

    ! --- b: Rank 0 sorts cells by partition for efficient rank-contiguous processing ---
    if (rank == 0) then
      allocate(cell_order(nc), rank_offsets(0:nranks))
      do i = 1, nc
        cell_order(i) = i
      end do
      if (rank == 0) write(*, '(A)') ' [HILBERT] Sorting cells by partition (quicksort)...'
      call quicksort_cells_by_partition(cell_order, part, 1, nc)
      if (rank == 0) write(*, '(A)') ' [HILBERT] Computing rank offsets...'
      call compute_rank_offsets(cell_order, part, nc, nranks, rank_offsets)
    end if

    ! --- Prepare geometry on rank 0 (non-zero ranks skip and go directly to recv_local_mesh) ---

    ! Non-zero ranks exit early to avoid memory pressure from collective MPI operations
    ! They receive their complete mesh via point-to-point recv_local_mesh
    if (rank /= 0) then
      ! Deallocate temporary data before exiting to free memory under pressure
      deallocate(cell_centroid, periodic_pair_global)
      if (allocated(vertices)) deallocate(vertices)
      if (allocated(cell_offsets)) deallocate(cell_offsets)
      if (allocated(cell_vertices)) deallocate(cell_vertices)
      if (allocated(cell_types)) deallocate(cell_types)
      if (allocated(face_offsets)) deallocate(face_offsets)
      if (allocated(face_vertices)) deallocate(face_vertices)

      call recv_local_mesh(mesh_local, 0, comm)
      deallocate(part)
      if (allocated(cell_order)) deallocate(cell_order)
      if (allocated(rank_offsets)) deallocate(rank_offsets)

      ! Synchronize with rank 0 before returning
      call MPI_Barrier(comm, ierr)
      return
    end if

    ! --- Rank 0 only - Allocate and pack geometry ---
    if (.not. allocated(face_owner)) then
      allocate(face_owner(nf), face_nei(nf), face_patch(nf))
      allocate(face_area(nf), face_centre(3, nf), face_normal(3, nf), face_sf(3, nf), face_dcf(nf), face_gc(nf))
      allocate(patch_start(np), patch_nfaces(np), patch_bcid(np), patch_name(np), patch_type(np))
      allocate(cell_volume(nc))
      call pack_global_mesh(mesh_global, nc, nf, np, cell_centroid, cell_volume, &
                            face_owner, face_nei, face_patch, face_area, face_centre, &
                            face_normal, face_sf, face_dcf, face_gc, patch_start, patch_nfaces, &
                            patch_bcid, patch_name, patch_type)
    end if

    ! --- Rank 0 destroys global mesh and builds/sends local meshes ---
    if (rank == 0) then
      call mesh_global%destroy()
    end if

    do r = 0, nranks - 1
      ! Skip rank if it has no cells (empty offset range)
      if (rank_offsets(r) >= rank_offsets(r + 1)) cycle

      ! Clear tmp_mesh from previous iteration before building new one
      if (r > 0) call tmp_mesh%destroy()

      tmp_mesh%comm = comm
      tmp_mesh%rank = r
      tmp_mesh%nranks = nranks

      call build_local_partition_mesh(tmp_mesh, comm, r, nranks, &
              nc, nf, np, part, &
              cell_centroid=cell_centroid, cell_volume=cell_volume, &
              face_owner=face_owner, face_nei=face_nei, face_patch=face_patch, &
              face_area=face_area, face_dcf=face_dcf, &
              patch_start=patch_start, patch_nfaces=patch_nfaces, patch_bcid=patch_bcid, &
              patch_name=patch_name, patch_type=patch_type, &
              periodic_pair_global=periodic_pair_global, &
              topo_flag=topo_flag, n_vertices_global=nv, vertices_global=vertices, &
              cell_offsets_global=cell_offsets, cell_vertices_global=cell_vertices, &
              cell_types_global=cell_types, &
              face_topo_flag=face_topo_flag, face_offsets_global=face_offsets, &
              face_vertices_global=face_vertices, &
              cell_order=cell_order, cell_start=rank_offsets(r), &
              cell_end=rank_offsets(r + 1) - 1, &
              face_centre=face_centre, face_normal=face_normal, &
              face_sf=face_sf, face_gc=face_gc)

      if (r == 0) then
        mesh_local = tmp_mesh
      else
        call send_local_mesh(tmp_mesh, r, comm)
      end if
    end do

    ! Cleanup
    deallocate(part, cell_centroid, cell_volume)
    deallocate(face_owner, face_nei, face_patch)
    deallocate(face_area, face_centre, face_normal, face_sf, face_dcf, face_gc)
    deallocate(patch_start, patch_nfaces, patch_bcid, patch_name, patch_type)
    deallocate(periodic_pair_global, cell_order, rank_offsets)
    if (allocated(vertices)) deallocate(vertices)
    if (allocated(cell_offsets)) deallocate(cell_offsets)
    if (allocated(cell_vertices)) deallocate(cell_vertices)
    if (allocated(cell_types)) deallocate(cell_types)
    if (allocated(face_offsets)) deallocate(face_offsets)
    if (allocated(face_vertices)) deallocate(face_vertices)

    call MPI_Barrier(comm, ierr)

    if (rank == 0) write(*, '(A)') ' [HILBERT] Partitioning complete.'

  end subroutine partition_mesh_hilbert

  !> Send a local mesh from rank 0 to a destination rank.
  !> Serializes mesh data and sends via MPI.
  subroutine send_local_mesh(msh, dest_rank, comm)
    type(Mesh), intent(in)       :: msh
    integer(ip), intent(in)      :: dest_rank
    type(MPI_Comm), intent(in)   :: comm

    integer(ip) :: nc, ng, nf, nif, nbf, np, nv, nh, i, ierr, tag, rank
    integer(ip) :: header(13)
    integer(ip) :: total_cell_verts, total_face_verts

    call MPI_Comm_rank(comm, rank, ierr)

    nc  = msh%n_cells
    ng  = msh%n_ghost
    nf  = msh%n_faces
    nif = msh%n_int_faces
    nbf = msh%n_bnd_faces
    np  = msh%n_patches
    nv  = msh%n_vertices
    nh  = 0
    if (allocated(msh%halos)) nh = size(msh%halos)

    ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
    ! must never be evaluated in the same expression as the allocated() check.
    total_cell_verts = 0
    total_face_verts = 0
    if (allocated(msh%cell_offsets)) then
      if (size(msh%cell_offsets) > nc) total_cell_verts = msh%cell_offsets(nc + 1) - 1
    end if
    if (allocated(msh%face_offsets)) then
      if (size(msh%face_offsets) > nf) total_face_verts = msh%face_offsets(nf + 1) - 1
    end if

    ! Pack header
    header(1)  = nc
    header(2)  = ng
    header(3)  = nf
    header(4)  = nif
    header(5)  = nbf
    header(6)  = np
    header(7)  = nv
    header(8)  = int(msh%n_cells_global, ip)
    header(9)  = 0
    if (allocated(msh%vertices)) header(9) = 1
    header(10) = 0
    if (allocated(msh%face_offsets)) header(10) = 1
    header(11) = nh
    header(12) = total_cell_verts
    header(13) = total_face_verts

    tag = 9000
    call MPI_Send(header, 13, MPI_INTEGER, dest_rank, tag, comm, ierr)

    ! Send cell data (flat arrays)
    if (nc + ng > 0) then
      tag = 9001
      call send_cell_data(msh, nc, ng, dest_rank, comm, tag)
    end if

    ! Send face data (flat arrays)
    if (nf > 0) then
      tag = 9010
      call send_face_data(msh, nf, dest_rank, comm, tag)
    end if

    ! Send patch data
    if (np > 0) then
      tag = 9020
      call send_patch_data(msh, np, dest_rank, comm, tag)
    end if

    ! Send periodic pairing (local same-rank pairs)
    ! Must always be allocated if nbf > 0
    if (nbf > 0) then
      if (.not. allocated(msh%periodic_pair)) then
        error stop 'send_local_mesh: periodic_pair not allocated but nbf > 0'
      end if
      tag = 9030
      call MPI_Send(msh%periodic_pair, nbf, MPI_INTEGER, dest_rank, tag, comm, ierr)
    end if

    ! Send periodic neighbors (cross-rank periodic mappings)
    ! Must always be allocated if nbf > 0
    if (nbf > 0) then
      if (.not. allocated(msh%periodic_neighbor)) then
        error stop 'send_local_mesh: periodic_neighbor not allocated but nbf > 0'
      end if
      tag = 9032
      call MPI_Send(msh%periodic_neighbor, nbf, MPI_INTEGER, dest_rank, tag, comm, ierr)
    end if

    ! Send periodic partner distances used for decomposition-invariant metrics.
    if (nbf > 0) then
      tag = 9033
      if (allocated(msh%periodic_pair_dcf)) then
        call MPI_Send(msh%periodic_pair_dcf, nbf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
      else
        call MPI_Send([ (msh%faces(msh%n_int_faces + i)%d_cf, i=1,nbf) ], nbf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
      end if
    end if

    ! Send global IDs
    if (nc + ng > 0) then
      tag = 9031
      if (allocated(msh%global_id)) then
        call MPI_Send(msh%global_id, nc + ng, MPI_INT64_T, dest_rank, tag, comm, ierr)
      else
        call MPI_Send([1_lp], 1, MPI_INT64_T, dest_rank, tag, comm, ierr)
      end if
    end if

    ! Send optional topology
    if (header(9) == 1) then
      tag = 9040
      call MPI_Send(msh%vertices, 3 * nv, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
      tag = 9041
      call MPI_Send(msh%cell_offsets, nc + 1, MPI_INTEGER, dest_rank, tag, comm, ierr)
      tag = 9042
      call MPI_Send(msh%cell_vertices, total_cell_verts, MPI_INTEGER, dest_rank, tag, comm, ierr)
      tag = 9043
      call MPI_Send(msh%cell_types, nc, MPI_INTEGER, dest_rank, tag, comm, ierr)
    end if

    if (header(10) == 1) then
      tag = 9044
      call MPI_Send(msh%face_offsets, nf + 1, MPI_INTEGER, dest_rank, tag, comm, ierr)
      tag = 9045
      call MPI_Send(msh%face_vertices, total_face_verts, MPI_INTEGER, dest_rank, tag, comm, ierr)
    end if

    ! Send cell_faces if allocated
    if (allocated(msh%cell_faces)) then
      tag = 9046
      call MPI_Send(size(msh%cell_faces), 1, MPI_INTEGER, dest_rank, tag, comm, ierr)
      tag = 9047
      call MPI_Send(msh%cell_faces, size(msh%cell_faces), MPI_INTEGER, dest_rank, tag, comm, ierr)
    else
      tag = 9046
      call MPI_Send([0_ip], 1, MPI_INTEGER, dest_rank, tag, comm, ierr)
    end if

    ! Send halos
    if (nh > 0) then
      do i = 1, nh
        tag = 9050 + i
        call MPI_Send([msh%halos(i)%rank, msh%halos(i)%n_send, msh%halos(i)%n_recv], &
                      3, MPI_INTEGER, dest_rank, tag, comm, ierr)
        if (msh%halos(i)%n_send > 0) then
          tag = 9100 + i
          call MPI_Send(msh%halos(i)%send_cells, msh%halos(i)%n_send, MPI_INTEGER, dest_rank, tag, comm, ierr)
        end if
        if (msh%halos(i)%n_recv > 0) then
          tag = 9150 + i
          call MPI_Send(msh%halos(i)%recv_cells, msh%halos(i)%n_recv, MPI_INTEGER, dest_rank, tag, comm, ierr)
        end if
      end do
    end if


  end subroutine send_local_mesh

  !> Helper: send cell data arrays
  subroutine send_cell_data(msh, nc, ng, dest_rank, comm, tag_start)
    type(Mesh), intent(in)      :: msh
    integer(ip), intent(in)     :: nc, ng, dest_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    real(rp), allocatable :: centroid(:,:), volume(:)
    integer(ip), allocatable :: nfaces(:), offsets(:)

    allocate(centroid(3, nc+ng), volume(nc+ng), nfaces(nc+ng), offsets(nc+ng))

    do i = 1, nc + ng
      centroid(:, i) = msh%cells(i)%centroid
      volume(i) = msh%cells(i)%volume
      nfaces(i) = msh%cells(i)%n_faces
      offsets(i) = msh%cells(i)%face_offset
    end do

    tag = tag_start
    call MPI_Send(centroid, 3*(nc+ng), MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(volume, nc+ng, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(nfaces, nc+ng, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(offsets, nc+ng, MPI_INTEGER, dest_rank, tag, comm, ierr)

    deallocate(centroid, volume, nfaces, offsets)

  end subroutine send_cell_data

  !> Helper: send face data arrays
  subroutine send_face_data(msh, nf, dest_rank, comm, tag_start)
    type(Mesh), intent(in)      :: msh
    integer(ip), intent(in)     :: nf, dest_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    integer(ip), allocatable :: owner(:), nei(:), patch(:)
    real(rp), allocatable :: area(:), centre(:,:), normal(:,:), sf(:,:), dcf(:), gc(:)

    allocate(owner(nf), nei(nf), patch(nf), area(nf), centre(3,nf), normal(3,nf), sf(3,nf), dcf(nf), gc(nf))

    do i = 1, nf
      owner(i) = msh%faces(i)%owner
      nei(i) = msh%faces(i)%neighbour
      patch(i) = msh%faces(i)%patch_id
      area(i) = msh%faces(i)%area
      centre(:, i) = msh%faces(i)%centre
      normal(:, i) = msh%faces(i)%normal
      sf(:, i) = msh%faces(i)%sf
      dcf(i) = msh%faces(i)%d_cf
      gc(i) = msh%faces(i)%gc
    end do

    tag = tag_start
    call MPI_Send(owner, nf, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(nei, nf, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(patch, nf, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(area, nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(centre, 3*nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(normal, 3*nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(sf, 3*nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(dcf, nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(gc, nf, MPI_DOUBLE_PRECISION, dest_rank, tag, comm, ierr)

    deallocate(owner, nei, patch, area, centre, normal, sf, dcf, gc)

  end subroutine send_face_data

  !> Helper: send patch data arrays
  subroutine send_patch_data(msh, np, dest_rank, comm, tag_start)
    type(Mesh), intent(in)      :: msh
    integer(ip), intent(in)     :: np, dest_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    integer(ip), allocatable :: bcid(:), start(:), nfaces(:)
    character(len=128), allocatable :: names(:)
    character(len=64), allocatable  :: types(:)

    allocate(bcid(np), start(np), nfaces(np), names(np), types(np))

    do i = 1, np
      bcid(i) = msh%patches(i)%bc_id
      start(i) = msh%patches(i)%start_face
      nfaces(i) = msh%patches(i)%n_faces
      names(i) = msh%patches(i)%name
      types(i) = msh%patches(i)%bc_type
    end do

    tag = tag_start
    call MPI_Send(bcid, np, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(start, np, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(nfaces, np, MPI_INTEGER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(names, 128*np, MPI_CHARACTER, dest_rank, tag, comm, ierr)
    tag = tag + 1
    call MPI_Send(types, 64*np, MPI_CHARACTER, dest_rank, tag, comm, ierr)

    deallocate(bcid, start, nfaces, names, types)

  end subroutine send_patch_data

  !> Receive a local mesh from rank 0 on a non-zero rank.
  !> Deserializes mesh data received via MPI.
  subroutine recv_local_mesh(msh, src_rank, comm)
    type(Mesh), intent(inout)    :: msh
    integer(ip), intent(in)      :: src_rank
    type(MPI_Comm), intent(in)   :: comm

    integer(ip) :: nc, ng, nf, nif, nbf, np, nv, nh, i, ierr, tag, rank, nranks
    integer(ip) :: header(13)
    integer(ip) :: total_cell_verts, total_face_verts, ncf
    integer(ip), allocatable :: halo_header(:)
    type(MPI_Status) :: status

    ! Get current rank and size for mesh initialization
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nranks, ierr)

    ! Receive header
    tag = 9000
    call MPI_Recv(header, 13, MPI_INTEGER, src_rank, tag, comm, status, ierr)

    nc  = header(1)
    ng  = header(2)
    nf  = header(3)
    nif = header(4)
    nbf = header(5)
    np  = header(6)
    nv  = header(7)
    nh  = header(11)
    total_cell_verts = header(12)
    total_face_verts = header(13)

    ! Initialize mesh metadata (needed for BC setup and other operations)
    msh%comm = comm
    msh%rank = rank
    msh%nranks = nranks

    ! Allocate mesh arrays
    msh%n_cells = nc
    msh%n_ghost = ng
    msh%n_faces = nf
    msh%n_int_faces = nif
    msh%n_bnd_faces = nbf
    msh%n_patches = np
    msh%n_vertices = nv
    msh%n_cells_global = int(header(8), lp)

    if (nc + ng > 0) then
      allocate(msh%cells(nc + ng))
      tag = 9001
      call recv_cell_data(msh, nc, ng, src_rank, comm, tag)
    end if

    if (nf > 0) then
      allocate(msh%faces(nf))
      tag = 9010
      call recv_face_data(msh, nf, src_rank, comm, tag)
    end if

    if (np > 0) then
      allocate(msh%patches(np))
      tag = 9020
      call recv_patch_data(msh, np, src_rank, comm, tag)
    end if

    ! Receive periodic pairing (local same-rank pairs)
    if (nbf > 0) then
      allocate(msh%periodic_pair(nbf), source=0_ip)
      allocate(msh%periodic_neighbor(nbf), source=0_ip)
      allocate(msh%periodic_pair_dcf(nbf), source=0.0_rp)
      tag = 9030
      call MPI_Recv(msh%periodic_pair, nbf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    end if

    ! Receive periodic neighbors (cross-rank periodic mappings)
    if (nbf > 0) then
      tag = 9032
      call MPI_Recv(msh%periodic_neighbor, nbf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    end if

    if (nbf > 0) then
      tag = 9033
      call MPI_Recv(msh%periodic_pair_dcf, nbf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    end if

    ! Receive global IDs
    if (nc + ng > 0) then
      allocate(msh%global_id(nc + ng))
      tag = 9031
      call MPI_Recv(msh%global_id, nc + ng, MPI_INT64_T, src_rank, tag, comm, status, ierr)
    end if

    ! Receive optional topology
    if (header(9) == 1 .and. nv > 0) then
      allocate(msh%vertices(3, nv))
      allocate(msh%cell_offsets(nc + 1))
      allocate(msh%cell_vertices(total_cell_verts))
      allocate(msh%cell_types(nc))

      tag = 9040
      call MPI_Recv(msh%vertices, 3*nv, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
      tag = 9041
      call MPI_Recv(msh%cell_offsets, nc+1, MPI_INTEGER, src_rank, tag, comm, status, ierr)
      tag = 9042
      call MPI_Recv(msh%cell_vertices, total_cell_verts, MPI_INTEGER, src_rank, tag, comm, status, ierr)
      tag = 9043
      call MPI_Recv(msh%cell_types, nc, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    end if

    if (header(10) == 1 .and. nf > 0) then
      allocate(msh%face_offsets(nf + 1))
      allocate(msh%face_vertices(total_face_verts))

      tag = 9044
      call MPI_Recv(msh%face_offsets, nf+1, MPI_INTEGER, src_rank, tag, comm, status, ierr)
      tag = 9045
      call MPI_Recv(msh%face_vertices, total_face_verts, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    end if

    ! Receive cell_faces if sent
    tag = 9046
    call MPI_Recv(ncf, 1, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    if (ncf > 0) then
      allocate(msh%cell_faces(ncf))
      tag = 9047
      call MPI_Recv(msh%cell_faces, ncf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    end if

    ! Receive halos
    if (nh > 0) then
      allocate(msh%halos(nh))
      allocate(halo_header(3))
      do i = 1, nh
        tag = 9050 + i
        call MPI_Recv(halo_header, 3, MPI_INTEGER, src_rank, tag, comm, status, ierr)
        msh%halos(i)%rank = halo_header(1)
        msh%halos(i)%n_send = halo_header(2)
        msh%halos(i)%n_recv = halo_header(3)

        if (msh%halos(i)%n_send > 0) then
          allocate(msh%halos(i)%send_cells(msh%halos(i)%n_send))
          tag = 9100 + i
          call MPI_Recv(msh%halos(i)%send_cells, msh%halos(i)%n_send, MPI_INTEGER, src_rank, tag, comm, status, ierr)
        end if

        if (msh%halos(i)%n_recv > 0) then
          allocate(msh%halos(i)%recv_cells(msh%halos(i)%n_recv))
          tag = 9150 + i
          call MPI_Recv(msh%halos(i)%recv_cells, msh%halos(i)%n_recv, MPI_INTEGER, src_rank, tag, comm, status, ierr)
        end if
      end do
      deallocate(halo_header)
    end if


  end subroutine recv_local_mesh

  !> Helper: receive cell data arrays
  subroutine recv_cell_data(msh, nc, ng, src_rank, comm, tag_start)
    type(Mesh), intent(inout)   :: msh
    integer(ip), intent(in)     :: nc, ng, src_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    real(rp), allocatable :: centroid(:,:), volume(:)
    integer(ip), allocatable :: nfaces(:), offsets(:)
    type(MPI_Status) :: status

    allocate(centroid(3, nc+ng), volume(nc+ng), nfaces(nc+ng), offsets(nc+ng))

    tag = tag_start
    call MPI_Recv(centroid, 3*(nc+ng), MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(volume, nc+ng, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(nfaces, nc+ng, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(offsets, nc+ng, MPI_INTEGER, src_rank, tag, comm, status, ierr)

    do i = 1, nc + ng
      msh%cells(i)%centroid = centroid(:, i)
      msh%cells(i)%volume = volume(i)
      msh%cells(i)%n_faces = nfaces(i)
      msh%cells(i)%face_offset = offsets(i)
    end do

    deallocate(centroid, volume, nfaces, offsets)

  end subroutine recv_cell_data

  !> Helper: receive face data arrays
  subroutine recv_face_data(msh, nf, src_rank, comm, tag_start)
    type(Mesh), intent(inout)   :: msh
    integer(ip), intent(in)     :: nf, src_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    integer(ip), allocatable :: owner(:), nei(:), patch(:)
    real(rp), allocatable :: area(:), centre(:,:), normal(:,:), sf(:,:), dcf(:), gc(:)
    type(MPI_Status) :: status

    allocate(owner(nf), nei(nf), patch(nf), area(nf), centre(3,nf), normal(3,nf), sf(3,nf), dcf(nf), gc(nf))

    tag = tag_start
    call MPI_Recv(owner, nf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(nei, nf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(patch, nf, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(area, nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(centre, 3*nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(normal, 3*nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(sf, 3*nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(dcf, nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(gc, nf, MPI_DOUBLE_PRECISION, src_rank, tag, comm, status, ierr)

    do i = 1, nf
      msh%faces(i)%owner = owner(i)
      msh%faces(i)%neighbour = nei(i)
      msh%faces(i)%patch_id = patch(i)
      msh%faces(i)%area = area(i)
      msh%faces(i)%centre = centre(:, i)
      msh%faces(i)%normal = normal(:, i)
      msh%faces(i)%sf = sf(:, i)
      msh%faces(i)%d_cf = dcf(i)
      msh%faces(i)%gc = gc(i)
    end do

    deallocate(owner, nei, patch, area, centre, normal, sf, dcf, gc)

  end subroutine recv_face_data

  !> Helper: receive patch data arrays
  subroutine recv_patch_data(msh, np, src_rank, comm, tag_start)
    type(Mesh), intent(inout)   :: msh
    integer(ip), intent(in)     :: np, src_rank
    type(MPI_Comm), intent(in)  :: comm
    integer(ip), intent(in)     :: tag_start

    integer(ip) :: i, tag, ierr
    integer(ip), allocatable :: bcid(:), start(:), nfaces(:)
    character(len=128), allocatable :: names(:)
    character(len=64), allocatable  :: types(:)
    type(MPI_Status) :: status

    allocate(bcid(np), start(np), nfaces(np), names(np), types(np))

    tag = tag_start
    call MPI_Recv(bcid, np, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(start, np, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(nfaces, np, MPI_INTEGER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(names, 128*np, MPI_CHARACTER, src_rank, tag, comm, status, ierr)
    tag = tag + 1
    call MPI_Recv(types, 64*np, MPI_CHARACTER, src_rank, tag, comm, status, ierr)

    do i = 1, np
      msh%patches(i)%bc_id = bcid(i)
      msh%patches(i)%start_face = start(i)
      msh%patches(i)%n_faces = nfaces(i)
      msh%patches(i)%name = names(i)
      msh%patches(i)%bc_type = types(i)
    end do

    deallocate(bcid, start, nfaces, names, types)

  end subroutine recv_patch_data

end module mod_partitioner
