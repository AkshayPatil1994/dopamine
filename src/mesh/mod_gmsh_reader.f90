! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Gmsh v4 .msh reader (ASCII): tet4, hex8, prism6, pyr5, tri3, quad4.
module mod_gmsh_reader
  use mod_kinds
  use mod_parameters
  use mod_mesh
  implicit none
  private

  public :: gmsh_read

  ! Gmsh element type codes and their node counts
  integer(ip), parameter :: GMSH_LINE2  = 1
  integer(ip), parameter :: GMSH_TRI3   = 2
  integer(ip), parameter :: GMSH_QUAD4  = 3
  integer(ip), parameter :: GMSH_TET4   = 4
  integer(ip), parameter :: GMSH_HEX8   = 5
  integer(ip), parameter :: GMSH_PRISM6 = 6
  integer(ip), parameter :: GMSH_PYR5   = 7
  integer(ip), parameter :: GMSH_POINT  = 15
  integer(ip), parameter :: GMSH_TET10  = 11
  integer(ip), parameter :: GMSH_HEX20  = 17

  type :: PhysGroup
    integer(ip)   :: dim
    integer(ip)   :: tag
    character(128):: name
  end type PhysGroup

  type :: GmshElement
    integer(ip)              :: etype
    integer(ip)              :: phys_tag
    integer(ip), allocatable :: nodes(:)
  end type GmshElement

  type :: EntityPhysMap
    integer(ip) :: dim
    integer(ip) :: entity_tag
    integer(ip) :: phys_tag
  end type EntityPhysMap

contains

  subroutine gmsh_read(filename, msh)
    character(len=*), intent(in)    :: filename
    type(Mesh),       intent(inout) :: msh

    integer(ip) :: iunit, ios
    character(256) :: line

    type(PhysGroup),  allocatable :: phys_groups(:)
    type(EntityPhysMap), allocatable :: entity_phys(:)
    real(rp),         allocatable :: node_coords(:,:)  ! (3, n_nodes)
    type(GmshElement),allocatable :: elements(:)
    integer(lp),      allocatable :: node_tags(:)

    integer(ip) :: n_nodes, n_elements, n_phys, n_entity_phys
    logical :: got_phys, got_nodes, got_elements

    n_nodes = 0
    n_elements = 0
    n_phys = 0
    n_entity_phys = 0
    got_phys = .false.
    got_nodes = .false.
    got_elements = .false.

    allocate(phys_groups(0))
    allocate(entity_phys(0))
    allocate(node_coords(3,0), node_tags(0), elements(0))

    open(newunit=iunit, file=trim(filename), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) error stop 'gmsh_read: cannot open file: '//trim(filename)

    ! Scan sections
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      select case (trim(line))
      case ('$PhysicalNames')
        if (allocated(phys_groups)) deallocate(phys_groups)
        call read_physical_names(iunit, phys_groups, n_phys)
        got_phys = .true.
      case ('$Entities')
        if (allocated(entity_phys)) deallocate(entity_phys)
        call read_entities_v4(iunit, entity_phys, n_entity_phys)
      case ('$Nodes')
        if (allocated(node_coords)) deallocate(node_coords)
        if (allocated(node_tags)) deallocate(node_tags)
        call read_nodes_v4(iunit, node_coords, node_tags, n_nodes)
        got_nodes = .true.
      case ('$Elements')
        if (allocated(elements)) deallocate(elements)
        call read_elements_v4(iunit, entity_phys, n_entity_phys, elements, n_elements)
        got_elements = .true.
      end select
    end do
    close(iunit)

    if (.not. got_nodes) error stop 'gmsh_read: missing $Nodes section'
    if (.not. got_elements) error stop 'gmsh_read: missing $Elements section'
    if (.not. got_phys) then
      if (allocated(phys_groups)) deallocate(phys_groups)
      allocate(phys_groups(0))
      n_phys = 0
    end if

    ! Build msh from raw data
    call build_mesh_from_gmsh(msh, node_coords, node_tags, elements, &
                               phys_groups, n_nodes, n_elements, n_phys)

    deallocate(node_coords, node_tags, elements, phys_groups, entity_phys)
  end subroutine gmsh_read

  subroutine read_physical_names(iunit, groups, ngroups)
    integer(ip),                        intent(in)  :: iunit
    type(PhysGroup), allocatable,       intent(out) :: groups(:)
    integer(ip),                        intent(out) :: ngroups
    integer(ip) :: i
    character(256) :: line

    read(iunit, *) ngroups
    allocate(groups(ngroups))
    do i = 1, ngroups
      read(iunit, *) groups(i)%dim, groups(i)%tag, groups(i)%name
      ! Strip surrounding quotes from name
      groups(i)%name = trim(adjustl(groups(i)%name))
      if (groups(i)%name(1:1) == '"') &
        groups(i)%name = groups(i)%name(2:len_trim(groups(i)%name)-1)
    end do
    read(iunit, '(A)') line  ! $EndPhysicalNames
  end subroutine read_physical_names

  subroutine read_nodes_v4(iunit, coords, tags, n_nodes)
    integer(ip),          intent(in)  :: iunit
    real(rp), allocatable,intent(out) :: coords(:,:)
    integer(lp),allocatable,intent(out):: tags(:)
    integer(ip),          intent(out) :: n_nodes

    integer(ip)  :: n_entity_blocks, min_tag, max_tag
    integer(ip)  :: ib, n_in_block, entity_dim, entity_tag, parametric
    integer(ip)  :: i, offset
    character(256) :: line

    read(iunit, *) n_entity_blocks, n_nodes, min_tag, max_tag
    allocate(coords(3, n_nodes), tags(n_nodes))

    offset = 0
    do ib = 1, n_entity_blocks
      read(iunit, *) entity_dim, entity_tag, parametric, n_in_block
      do i = 1, n_in_block
        read(iunit, *) tags(offset + i)
      end do
      do i = 1, n_in_block
        read(iunit, *) coords(1, offset+i), coords(2, offset+i), coords(3, offset+i)
      end do
      offset = offset + n_in_block
    end do
    read(iunit, '(A)') line  ! $EndNodes
  end subroutine read_nodes_v4

  subroutine read_entities_v4(iunit, entity_phys, n_entity_phys)
    integer(ip), intent(in) :: iunit
    type(EntityPhysMap), allocatable, intent(out) :: entity_phys(:)
    integer(ip), intent(out) :: n_entity_phys

    integer(ip) :: n_points, n_curves, n_surfaces, n_volumes
    integer(ip) :: i
    character(2048) :: line

    n_entity_phys = 0
    allocate(entity_phys(0))

    read(iunit, *) n_points, n_curves, n_surfaces, n_volumes

    do i = 1, n_points
      read(iunit, '(A)') line
      call parse_entity_line_for_phys(line, 0_ip, entity_phys, n_entity_phys)
    end do
    do i = 1, n_curves
      read(iunit, '(A)') line
      call parse_entity_line_for_phys(line, 1_ip, entity_phys, n_entity_phys)
    end do
    do i = 1, n_surfaces
      read(iunit, '(A)') line
      call parse_entity_line_for_phys(line, 2_ip, entity_phys, n_entity_phys)
    end do
    do i = 1, n_volumes
      read(iunit, '(A)') line
      call parse_entity_line_for_phys(line, 3_ip, entity_phys, n_entity_phys)
    end do

    read(iunit, '(A)') line  ! $EndEntities
  end subroutine read_entities_v4

  subroutine read_elements_v4(iunit, entity_phys, n_entity_phys, elements, n_elem)
    integer(ip),               intent(in)  :: iunit
    type(EntityPhysMap),       intent(in)  :: entity_phys(:)
    integer(ip),               intent(in)  :: n_entity_phys
    type(GmshElement),allocatable,intent(out) :: elements(:)
    integer(ip),               intent(out) :: n_elem

    integer(ip) :: n_entity_blocks, min_tag, max_tag
    integer(ip) :: ib, n_in_block, entity_dim, entity_tag, etype
    integer(ip) :: i, offset, nn
    integer(lp) :: etag
    character(256) :: line

    read(iunit, *) n_entity_blocks, n_elem, min_tag, max_tag
    allocate(elements(n_elem))

    offset = 0
    do ib = 1, n_entity_blocks
      read(iunit, *) entity_dim, entity_tag, etype, n_in_block
      nn = gmsh_nodes_per_element(etype)
      ! A guessed node count would desynchronise the list-directed read of every later element.
      if (nn < 1) then
        write(*, '(A,I0,A)') 'ERROR: unsupported Gmsh element type ', etype, &
          ' (supported: first-order line, tri, quad, tet, hex, prism, pyramid, point)'
        error stop 'gmsh_read: unsupported element type; export first-order elements'
      end if
      do i = 1, n_in_block
        elements(offset+i)%etype    = etype
        ! Untagged entities keep phys_tag = -1: falling back to the entity tag could silently match an
        ! unrelated physical group that happens to carry the same number.
        elements(offset+i)%phys_tag = lookup_entity_phys_tag(entity_phys, n_entity_phys, entity_dim, entity_tag)
        allocate(elements(offset+i)%nodes(nn))
        read(iunit, *) etag, elements(offset+i)%nodes(1:nn)
      end do
      offset = offset + n_in_block
    end do
    read(iunit, '(A)') line  ! $EndElements
  end subroutine read_elements_v4

  !> Assemble Mesh from parsed Gmsh data.
  subroutine build_mesh_from_gmsh(msh, coords, node_tags, elements, &
                                   groups, n_nodes, n_elem, n_phys)
    type(Mesh),        intent(inout) :: msh
    real(rp),          intent(in)    :: coords(:,:)
    integer(lp),       intent(in)    :: node_tags(:)
    type(GmshElement), intent(in)    :: elements(:)
    type(PhysGroup),   intent(in)    :: groups(:)
    integer(ip),       intent(in)    :: n_nodes, n_elem, n_phys

    ! Build a tag->local index map for nodes
    integer(ip), allocatable :: tag_to_local(:)
    integer(ip) :: max_tag, i, e, c
    integer(ip) :: n_vol, n_bnd
    integer(ip), allocatable :: vol_indices(:), bnd_indices(:)  ! heap: n_elem can exceed the OpenMP stack
    integer(ip), allocatable :: conn(:)
    integer(ip) :: total_conn, vtk_code, nconn

    ! Map node global tags to 1..n_nodes local indices
    max_tag = int(maxval(node_tags))
    allocate(tag_to_local(max_tag))
    tag_to_local = 0
    do i = 1, n_nodes
      tag_to_local(node_tags(i)) = i
    end do

    ! Separate volume and boundary elements
    allocate(vol_indices(max(n_elem, 1)), bnd_indices(max(n_elem, 1)))
    n_vol = 0; n_bnd = 0
    do e = 1, n_elem
      if (is_volume_element(elements(e)%etype)) then
        n_vol = n_vol + 1
        vol_indices(n_vol) = e
      else if (is_bnd_element(elements(e)%etype)) then
        n_bnd = n_bnd + 1
        bnd_indices(n_bnd) = e
      end if
    end do

    msh%n_cells = n_vol
    msh%n_ghost  = 0
    allocate(msh%cells(n_vol))
    allocate(msh%global_id(n_vol))

    ! Optional full topology for geometry-aware output (serial path).
    msh%n_vertices = n_nodes
    allocate(msh%vertices(3, n_nodes))
    msh%vertices = coords
    allocate(msh%cell_types(n_vol))
    allocate(msh%cell_offsets(n_vol + 1))

    total_conn = 0
    msh%cell_offsets(1) = 1
    do c = 1, n_vol
      e = vol_indices(c)
      call extract_volume_connectivity(elements(e), tag_to_local, vtk_code, conn)
      nconn = size(conn)
      msh%cell_types(c) = vtk_code
      total_conn = total_conn + nconn
      msh%cell_offsets(c + 1) = msh%cell_offsets(c) + nconn
      deallocate(conn)
    end do
    allocate(msh%cell_vertices(total_conn))

    total_conn = 1
    do c = 1, n_vol
      e = vol_indices(c)
      call extract_volume_connectivity(elements(e), tag_to_local, vtk_code, conn)
      nconn = size(conn)
      msh%cell_vertices(total_conn:total_conn+nconn-1) = conn
      total_conn = total_conn + nconn
      deallocate(conn)
    end do

    ! Compute cell centroids and volumes
    do c = 1, n_vol
      e = vol_indices(c)
      msh%global_id(c) = int(c, lp)
      call compute_cell_geometry(elements(e), coords, tag_to_local, &
                                  msh%cells(c)%centroid, msh%cells(c)%volume)
    end do

    ! Build faces from cell connectivity (shared faces detected via hash)
    call build_faces_from_cells(msh, elements, vol_indices, bnd_indices, &
                                 n_vol, n_bnd, coords, tag_to_local, &
                                 groups, n_phys)

    deallocate(tag_to_local, vol_indices, bnd_indices)
  end subroutine build_mesh_from_gmsh

  !> Extract low-order volume connectivity and matching VTK cell type code.
  subroutine extract_volume_connectivity(elem, t2l, vtk_code, conn)
    type(GmshElement), intent(in) :: elem
    integer(ip),       intent(in) :: t2l(:)
    integer(ip),       intent(out):: vtk_code
    integer(ip), allocatable, intent(out) :: conn(:)

    select case (elem%etype)
    case (GMSH_TET4, GMSH_TET10)
      vtk_code = 10_ip
      allocate(conn(4))
      conn = [t2l(elem%nodes(1)), t2l(elem%nodes(2)), t2l(elem%nodes(3)), t2l(elem%nodes(4))]
    case (GMSH_HEX8, GMSH_HEX20)
      vtk_code = 12_ip
      allocate(conn(8))
      conn = [t2l(elem%nodes(1)), t2l(elem%nodes(2)), t2l(elem%nodes(3)), t2l(elem%nodes(4)), &
              t2l(elem%nodes(5)), t2l(elem%nodes(6)), t2l(elem%nodes(7)), t2l(elem%nodes(8))]
    case (GMSH_PRISM6)
      vtk_code = 13_ip
      allocate(conn(6))
      conn = [t2l(elem%nodes(1)), t2l(elem%nodes(2)), t2l(elem%nodes(3)), &
              t2l(elem%nodes(4)), t2l(elem%nodes(5)), t2l(elem%nodes(6))]
    case (GMSH_PYR5)
      vtk_code = 14_ip
      allocate(conn(5))
      conn = [t2l(elem%nodes(1)), t2l(elem%nodes(2)), t2l(elem%nodes(3)), t2l(elem%nodes(4)), t2l(elem%nodes(5))]
    case default
      vtk_code = 10_ip
      allocate(conn(4))
      conn = [1_ip, 1_ip, 1_ip, 1_ip]
    end select
  end subroutine extract_volume_connectivity

  !> Build all faces from cell element connectivity.
  !> Internal faces are shared between two cells; boundary faces come from
  subroutine build_faces_from_cells(msh, elements, vol_idx, bnd_idx, &
                                     n_vol, n_bnd, coords, t2l, groups, n_phys)
    type(Mesh),        intent(inout) :: msh
    type(GmshElement), intent(in)    :: elements(:)
    integer(ip),       intent(in)    :: vol_idx(:), bnd_idx(:)
    integer(ip),       intent(in)    :: n_vol, n_bnd
    real(rp),          intent(in)    :: coords(:,:)
    integer(ip),       intent(in)    :: t2l(:)
    type(PhysGroup),   intent(in)    :: groups(:)
    integer(ip),       intent(in)    :: n_phys

    ! Hash-based face matching for production-scale meshes.
    integer(ip), parameter :: MAX_FACES_PER_CELL = 6
    integer(ip), parameter :: MAX_NODES_PER_FACE = 4

    integer(ip), allocatable :: cell_face_nodes(:,:,:)  ! (MAX_NODES_PER_FACE, MAX_FACES_PER_CELL, n_vol)
    integer(ip), allocatable :: cell_nfaces(:)
    integer(ip), allocatable :: face_owner(:), face_nb(:), face_patch(:)
    real(rp),    allocatable :: face_centre(:,:), face_normal(:,:), face_area(:)

    integer(ip) :: i, j, p, c, e, f, g, fi
    integer(ip) :: fnodes(MAX_NODES_PER_FACE)
    integer(ip) :: nf
    integer(ip) :: max_faces, n_int_f, n_bnd_f, n_total_f

    integer(ip), allocatable :: h_nodes(:,:), h_owner(:), h_state(:)
    integer(ip), allocatable :: patch_count(:), patch_map(:), patch_cursor(:), patch_offset(:), compact_count(:)
    integer(ip), allocatable :: bnd_owner(:), bnd_nb(:), bnd_patch(:)
    real(rp),    allocatable :: bnd_centre(:,:), bnd_normal(:,:), bnd_area(:)
    integer(ip) :: sig(MAX_NODES_PER_FACE), slot, probe, owner
    integer(ip) :: table_size
    integer(ip) :: npatches, pb, pidx

    ! Count faces per cell (first pass) - needed for max_faces estimate
    allocate(cell_nfaces(n_vol))
    allocate(cell_face_nodes(MAX_NODES_PER_FACE, MAX_FACES_PER_CELL, 1))  ! Rolling buffer (rank 3, size 1 in 3rd dim)
    do c = 1, n_vol
      e = vol_idx(c)
      call get_element_faces(elements(e), t2l, cell_face_nodes(:,:,1), cell_nfaces(c))
    end do

    ! Upper bound on faces
    max_faces = sum(cell_nfaces) / 2 + n_bnd + 100
    allocate(face_owner(max_faces), face_nb(max_faces))
    allocate(face_patch(max_faces), source=-1)
    allocate(face_centre(3, max_faces), face_normal(3, max_faces))
    allocate(face_area(max_faces))

    ! Hash table size: estimate unique faces with 0.75 load factor, cap at 500M entries.
    ! For typical meshes: ~n_vol * avg_faces_per_cell / 2 unique faces.
    table_size = min(500000000_ip, max(65536_ip, nint(real(max_faces, rp) * 2.0_rp / 0.75_rp)))
    ! Round up to next power of 2 for efficient hashing
    do while (iand(table_size, table_size - 1) /= 0)
      table_size = table_size + 1
    end do
    table_size = min(536870912_ip, table_size)  ! Cap at 2^29

    allocate(h_nodes(MAX_NODES_PER_FACE, table_size), source=0)
    allocate(h_owner(table_size), source=0)
    allocate(h_state(table_size), source=0)

    n_int_f = 0

    ! Match faces between cells using open-addressing hash table.
    do c = 1, n_vol
      ! Refill rolling buffer for current cell
      e = vol_idx(c)
      cell_face_nodes(:,:,1) = 0
      call get_element_faces(elements(e), t2l, cell_face_nodes(:,:,1), cell_nfaces(c))

      do fi = 1, cell_nfaces(c)
        fnodes = cell_face_nodes(:, fi, 1)
        nf = count(fnodes > 0)
        call canonical_face_nodes(fnodes(1:nf), nf, sig)
        slot = face_hash(sig, table_size)
        do probe = 1, table_size
          if (h_state(slot) == 0) then
            h_nodes(:, slot) = sig
            h_owner(slot) = c
            h_state(slot) = 1
            exit
          else if (all(h_nodes(:, slot) == sig)) then
            if (h_state(slot) == 1) then
              owner = h_owner(slot)
              n_int_f = n_int_f + 1
              face_owner(n_int_f) = owner
              face_nb(n_int_f)    = c
              call compute_face_geometry(fnodes(1:nf), coords, &
                                          face_centre(:,n_int_f), &
                                          face_normal(:,n_int_f), &
                                          face_area(n_int_f))
              call orient_normal(face_normal(:,n_int_f), &
                                  msh%cells(owner)%centroid, &
                                  msh%cells(c)%centroid)
              h_state(slot) = 2
            end if
            exit
          else
            if (slot == table_size) then
              slot = 1
            else
              slot = slot + 1
            end if
          end if
        end do
      end do
    end do

    msh%n_int_faces = n_int_f

    ! Boundary faces from 2-D elements
    n_bnd_f = 0
    do i = 1, n_bnd
      e   = bnd_idx(i)
      nf  = count(elements(e)%nodes > 0)
      n_bnd_f = n_bnd_f + 1
      f       = n_int_f + n_bnd_f

      ! Find owner cell for this boundary face via the same hash table
      face_nb(f) = -1
      do g = 1, nf
        fnodes(g) = t2l(elements(e)%nodes(g))
      end do
      call canonical_face_nodes(fnodes(1:nf), nf, sig)
      slot = face_hash(sig, table_size)
      face_owner(f) = -1
      do probe = 1, table_size
        if (h_state(slot) == 0) exit
        if (all(h_nodes(:, slot) == sig)) then
          face_owner(f) = h_owner(slot)
          exit
        end if
        if (slot == table_size) then
          slot = 1
        else
          slot = slot + 1
        end if
      end do
      if (face_owner(f) < 1) error stop 'build_faces_from_cells: boundary owner not found'

      call compute_face_geometry(fnodes(1:nf), coords, &
                                  face_centre(:,f), face_normal(:,f), face_area(f))
      ! Boundary normal points outward (away from owner centroid)
      call orient_outward_normal(face_normal(:,f), &
                                  msh%cells(face_owner(f))%centroid, &
                                  face_centre(:,f))
      ! Find patch
      do j = 1, n_phys
        if (groups(j)%dim == 2 .and. groups(j)%tag == elements(e)%phys_tag) then
          face_patch(f) = j
          exit
        end if
      end do
      if (n_phys > 0 .and. face_patch(f) < 1) then
        write(*, '(A,I0,A)') 'ERROR: boundary element with physical tag ', elements(e)%phys_tag, &
          ' belongs to no 2-D physical group'
        error stop 'build_faces_from_cells: boundary face without a physical surface group'
      end if
    end do

    msh%n_bnd_faces = n_bnd_f
    n_total_f        = n_int_f + n_bnd_f
    msh%n_faces     = n_total_f

    ! Reorder boundary faces by patch so each patch has a contiguous face range.
    allocate(patch_count(n_phys), source=0)
    do f = n_int_f + 1, n_total_f
      if (face_patch(f) >= 1 .and. face_patch(f) <= n_phys) then
        patch_count(face_patch(f)) = patch_count(face_patch(f)) + 1
      end if
    end do

    npatches = count(patch_count > 0)
    allocate(patch_map(n_phys), source=0)
    allocate(compact_count(max(1, npatches)), source=0)
    pidx = 0
    do p = 1, n_phys
      if (patch_count(p) > 0) then
        pidx = pidx + 1
        patch_map(p) = pidx
        compact_count(pidx) = patch_count(p)
      end if
    end do

    if (n_bnd_f > 0 .and. npatches > 0) then
      allocate(patch_offset(npatches), patch_cursor(npatches), source=0)
      patch_offset(1) = 1
      do p = 2, npatches
        patch_offset(p) = patch_offset(p-1) + compact_count(p-1)
      end do
      patch_cursor = patch_offset

      allocate(bnd_owner(n_bnd_f), bnd_nb(n_bnd_f), bnd_patch(n_bnd_f))
      allocate(bnd_centre(3, n_bnd_f), bnd_normal(3, n_bnd_f), bnd_area(n_bnd_f))

      do f = n_int_f + 1, n_total_f
        pb = patch_map(face_patch(f))
        if (pb < 1) cycle
        i = patch_cursor(pb)
        bnd_owner(i)   = face_owner(f)
        bnd_nb(i)      = face_nb(f)
        bnd_patch(i)   = pb
        bnd_centre(:,i)= face_centre(:,f)
        bnd_normal(:,i)= face_normal(:,f)
        bnd_area(i)    = face_area(f)
        patch_cursor(pb) = patch_cursor(pb) + 1
      end do

      do i = 1, n_bnd_f
        f = n_int_f + i
        face_owner(f) = bnd_owner(i)
        face_nb(f)    = bnd_nb(i)
        face_patch(f) = bnd_patch(i)
        face_centre(:,f) = bnd_centre(:,i)
        face_normal(:,f) = bnd_normal(:,i)
        face_area(f)     = bnd_area(i)
      end do

      deallocate(patch_offset, patch_cursor)
      deallocate(bnd_owner, bnd_nb, bnd_patch)
      deallocate(bnd_centre, bnd_normal, bnd_area)
    else if (n_bnd_f > 0) then
      ! Keep boundary faces in-place if no patch mapping was found.
      face_patch(n_int_f+1:n_total_f) = 1
    end if

    ! Allocate and fill msh faces
    allocate(msh%faces(n_total_f))
    do f = 1, n_total_f
      msh%faces(f)%owner      = face_owner(f)
      msh%faces(f)%neighbour  = face_nb(f)
      msh%faces(f)%area       = face_area(f)
      msh%faces(f)%centre     = face_centre(:, f)
      msh%faces(f)%normal     = face_normal(:, f)
      msh%faces(f)%sf         = face_area(f) * face_normal(:, f)
      if (f > n_int_f) msh%faces(f)%patch_id = face_patch(f)
    end do

    ! Build boundary patches
    call build_patches(msh, face_patch, n_int_f+1, n_total_f, groups, n_phys)

    deallocate(cell_nfaces, cell_face_nodes, face_owner, face_nb, face_patch)
    deallocate(face_centre, face_normal, face_area)
    deallocate(h_nodes, h_owner, h_state)
    deallocate(patch_count, patch_map, compact_count)
  end subroutine build_faces_from_cells

  subroutine build_patches(msh, face_patch, first_bnd, last_bnd, groups, n_phys)
    type(Mesh),      intent(inout) :: msh
    integer(ip),     intent(in)    :: face_patch(:)
    integer(ip),     intent(in)    :: first_bnd, last_bnd
    type(PhysGroup), intent(in)    :: groups(:)
    integer(ip),     intent(in)    :: n_phys

    integer(ip) :: f, p, npatches, start_f
    integer(ip), allocatable :: patch_count(:)

    allocate(patch_count(n_phys), source=0)
    do f = first_bnd, last_bnd
      p = face_patch(f)
      if (p >= 1 .and. p <= n_phys) patch_count(p) = patch_count(p) + 1
    end do

    npatches = count(patch_count > 0)
    msh%n_patches = npatches
    allocate(msh%patches(npatches))

    ! Assign contiguous face ranges per patch (boundary faces are reordered by patch).
    p = 0
    start_f = 1
    do f = 1, n_phys
      if (patch_count(f) == 0) cycle
      p = p + 1
      msh%patches(p)%name    = trim(groups(f)%name)
      msh%patches(p)%n_faces = patch_count(f)
      msh%patches(p)%start_face = start_f
      start_f = start_f + patch_count(f)
    end do
  end subroutine build_patches

  subroutine parse_entity_line_for_phys(line, dim, entity_phys, n_entity_phys)
    character(len=*), intent(in) :: line
    integer(ip),      intent(in) :: dim
    type(EntityPhysMap), allocatable, intent(inout) :: entity_phys(:)
    integer(ip), intent(inout) :: n_entity_phys

    character(64) :: tok(512)
    integer(ip) :: ntok, entity_tag, nphys, phys_tag

    call split_tokens(line, tok, ntok)
    if (ntok < 8) return

    read(tok(1), *, err=900) entity_tag
    read(tok(8), *, err=900) nphys
    if (nphys <= 0) return
    if (ntok < 8 + nphys) return

    read(tok(9), *, err=900) phys_tag
    if (nphys > 1) write(*, '(A,I0,A,I0,A,I0,A)') ' WARNING: Gmsh entity ', entity_tag, ' belongs to ', nphys, &
      ' physical groups; only the first (', phys_tag, ') is used'
    call append_entity_phys_map(entity_phys, n_entity_phys, dim, entity_tag, phys_tag)
    return

900 continue
    return
  end subroutine parse_entity_line_for_phys

  subroutine split_tokens(line, tokens, ntok)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: tokens(:)
    integer(ip), intent(out) :: ntok

    integer(ip) :: i, n, s, e

    ntok = 0
    n = len_trim(line)
    i = 1
    do while (i <= n)
      do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
        i = i + 1
      end do
      if (i > n) exit
      s = i
      do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= achar(9))
        i = i + 1
      end do
      e = i - 1
      if (ntok < size(tokens)) then
        ntok = ntok + 1
        tokens(ntok) = line(s:e)
      end if
    end do
  end subroutine split_tokens

  subroutine append_entity_phys_map(entity_phys, n_entity_phys, dim, entity_tag, phys_tag)
    type(EntityPhysMap), allocatable, intent(inout) :: entity_phys(:)
    integer(ip), intent(inout) :: n_entity_phys
    integer(ip), intent(in) :: dim, entity_tag, phys_tag

    type(EntityPhysMap), allocatable :: tmp(:)

    allocate(tmp(n_entity_phys + 1))
    if (n_entity_phys > 0) tmp(1:n_entity_phys) = entity_phys
    tmp(n_entity_phys + 1)%dim = dim
    tmp(n_entity_phys + 1)%entity_tag = entity_tag
    tmp(n_entity_phys + 1)%phys_tag = phys_tag

    call move_alloc(tmp, entity_phys)
    n_entity_phys = n_entity_phys + 1
  end subroutine append_entity_phys_map

  integer(ip) function lookup_entity_phys_tag(entity_phys, n_entity_phys, dim, entity_tag) result(phys_tag)
    type(EntityPhysMap), intent(in) :: entity_phys(:)
    integer(ip), intent(in) :: n_entity_phys, dim, entity_tag
    integer(ip) :: i

    phys_tag = -1
    do i = 1, n_entity_phys
      if (entity_phys(i)%dim == dim .and. entity_phys(i)%entity_tag == entity_tag) then
        phys_tag = entity_phys(i)%phys_tag
        return
      end if
    end do
  end function lookup_entity_phys_tag

  !> Compute cell centroid and volume via divergence-theorem tet decomposition (FPS Sect.A.2).
  subroutine compute_cell_geometry(elem, coords, t2l, centroid, volume)
    type(GmshElement), intent(in)  :: elem
    real(rp),          intent(in)  :: coords(:,:)
    integer(ip),       intent(in)  :: t2l(:)
    real(rp),          intent(out) :: centroid(3), volume

    select case (elem%etype)
    case (GMSH_TET4, GMSH_TET10)
      call tet_geometry(coords(:, t2l(elem%nodes(1))), &
                        coords(:, t2l(elem%nodes(2))), &
                        coords(:, t2l(elem%nodes(3))), &
                        coords(:, t2l(elem%nodes(4))), centroid, volume)
    case (GMSH_HEX8, GMSH_HEX20)
      call hex_geometry(elem%nodes(1:8), coords, t2l, centroid, volume)
    case (GMSH_PRISM6)
      call prism_geometry(elem%nodes(1:6), coords, t2l, centroid, volume)
    case (GMSH_PYR5)
      call pyr_geometry(elem%nodes(1:5), coords, t2l, centroid, volume)
    case default
      ! Fallback: vertex arithmetic mean (element type not recognised)
      block
        integer(ip) :: i, n, nn
        real(rp)    :: xc(3)
        nn = size(elem%nodes)
        xc = 0.0_rp
        do i = 1, nn
          n = t2l(elem%nodes(i))
          xc = xc + coords(:, n)
        end do
        centroid = xc / real(nn, rp)
        volume   = 1.0_rp
      end block
    end select
  end subroutine compute_cell_geometry

  !> Exact tet centroid and volume.
  subroutine tet_geometry(p1, p2, p3, p4, centroid, volume)
    real(rp), intent(in)  :: p1(3), p2(3), p3(3), p4(3)
    real(rp), intent(out) :: centroid(3), volume
    real(rp) :: a(3), b(3), c(3)
    a = p2 - p1; b = p3 - p1; c = p4 - p1
    volume   = abs(a(1)*(b(2)*c(3)-b(3)*c(2)) - a(2)*(b(1)*c(3)-b(3)*c(1)) + &
                   a(3)*(b(1)*c(2)-b(2)*c(1))) / 6.0_rp
    centroid = (p1 + p2 + p3 + p4) / 4.0_rp
  end subroutine tet_geometry

  !> Hex centroid and volume via 5-tet decomposition (Freudenthal partition).
  subroutine hex_geometry(nodes, coords, t2l, centroid, volume)
    integer(ip), intent(in)  :: nodes(8)
    real(rp),    intent(in)  :: coords(:,:)
    integer(ip), intent(in)  :: t2l(:)
    real(rp),    intent(out) :: centroid(3), volume
    real(rp)    :: p(3,8), tc(3), tv
    integer(ip) :: i
    integer(ip), parameter :: conn(4,5) = reshape( &
      [1,2,4,5, 2,3,4,7, 4,5,7,8, 2,5,6,7, 2,4,5,7], [4,5])
    do i = 1, 8; p(:,i) = coords(:, t2l(nodes(i))); end do
    volume = 0.0_rp; centroid = 0.0_rp
    do i = 1, 5
      call tet_geometry(p(:,conn(1,i)), p(:,conn(2,i)), &
                        p(:,conn(3,i)), p(:,conn(4,i)), tc, tv)
      volume = volume + tv; centroid = centroid + tv * tc
    end do
    if (volume > SMALL) centroid = centroid / volume
  end subroutine hex_geometry

  !> Prism centroid and volume via 3-tet decomposition.
  subroutine prism_geometry(nodes, coords, t2l, centroid, volume)
    integer(ip), intent(in)  :: nodes(6)
    real(rp),    intent(in)  :: coords(:,:)
    integer(ip), intent(in)  :: t2l(:)
    real(rp),    intent(out) :: centroid(3), volume
    real(rp)    :: p(3,6), tc(3), tv
    integer(ip) :: i
    integer(ip), parameter :: conn(4,3) = reshape( &
      [1,2,3,4, 2,4,5,3, 3,4,5,6], [4,3])
    do i = 1, 6; p(:,i) = coords(:, t2l(nodes(i))); end do
    volume = 0.0_rp; centroid = 0.0_rp
    do i = 1, 3
      call tet_geometry(p(:,conn(1,i)), p(:,conn(2,i)), &
                        p(:,conn(3,i)), p(:,conn(4,i)), tc, tv)
      volume = volume + tv; centroid = centroid + tv * tc
    end do
    if (volume > SMALL) centroid = centroid / volume
  end subroutine prism_geometry

  !> Pyramid centroid and volume via 2-tet decomposition.
  subroutine pyr_geometry(nodes, coords, t2l, centroid, volume)
    integer(ip), intent(in)  :: nodes(5)
    real(rp),    intent(in)  :: coords(:,:)
    integer(ip), intent(in)  :: t2l(:)
    real(rp),    intent(out) :: centroid(3), volume
    real(rp)    :: p(3,5), tc(3), tv
    integer(ip) :: i
    integer(ip), parameter :: conn(4,2) = reshape( &
      [1,2,4,5, 2,3,4,5], [4,2])
    do i = 1, 5; p(:,i) = coords(:, t2l(nodes(i))); end do
    volume = 0.0_rp; centroid = 0.0_rp
    do i = 1, 2
      call tet_geometry(p(:,conn(1,i)), p(:,conn(2,i)), &
                        p(:,conn(3,i)), p(:,conn(4,i)), tc, tv)
      volume = volume + tv; centroid = centroid + tv * tc
    end do
    if (volume > SMALL) centroid = centroid / volume
  end subroutine pyr_geometry

  subroutine compute_face_geometry(fnodes, coords, centre, normal, area)
    integer(ip), intent(in)  :: fnodes(:)
    real(rp),    intent(in)  :: coords(:,:)
    real(rp),    intent(out) :: centre(3), normal(3), area
    integer(ip) :: nf, i
    real(rp)    :: p0(3), pi(3), pip1(3), cross(3), n_acc(3)

    nf = size(fnodes)
    ! Centroid
    centre = 0.0_rp
    do i = 1, nf
      centre = centre + coords(:, fnodes(i))
    end do
    centre = centre / real(nf, rp)

    ! Area vector via fan triangulation from centroid
    n_acc = 0.0_rp
    p0 = centre
    do i = 1, nf
      pi   = coords(:, fnodes(i))
      pip1 = coords(:, fnodes(mod(i, nf) + 1))
      call cross_product(pi - p0, pip1 - p0, cross)
      n_acc = n_acc + 0.5_rp * cross
    end do
    area   = sqrt(sum(n_acc**2))
    normal = n_acc / (area + SMALL)
  end subroutine compute_face_geometry

  subroutine orient_normal(normal, owner_c, nb_c)
    real(rp), intent(inout) :: normal(3)
    real(rp), intent(in)    :: owner_c(3), nb_c(3)
    real(rp) :: d(3)
    d = nb_c - owner_c
    if (dot_product(normal, d) < 0.0_rp) normal = -normal
  end subroutine orient_normal

  subroutine orient_outward_normal(normal, owner_c, face_c)
    real(rp), intent(inout) :: normal(3)
    real(rp), intent(in)    :: owner_c(3), face_c(3)
    real(rp) :: d(3)
    d = face_c - owner_c
    if (dot_product(normal, d) < 0.0_rp) normal = -normal
  end subroutine orient_outward_normal

  subroutine get_element_faces(elem, t2l, face_nodes, nfaces)
    type(GmshElement), intent(in)  :: elem
    integer(ip),       intent(in)  :: t2l(:)
    integer(ip),       intent(out) :: face_nodes(4,6)  ! (max_nodes_per_face, max_faces)
    integer(ip),       intent(out) :: nfaces

    integer(ip) :: n(8)

    face_nodes = 0
    nfaces = 0

    select case (elem%etype)
    case (GMSH_TET4, GMSH_TET10)
      nfaces = 4
      n(1:4) = t2l(elem%nodes(1:4))
      face_nodes(:,1) = [n(1),n(2),n(3),0]
      face_nodes(:,2) = [n(1),n(2),n(4),0]
      face_nodes(:,3) = [n(2),n(3),n(4),0]
      face_nodes(:,4) = [n(1),n(3),n(4),0]
    case (GMSH_HEX8, GMSH_HEX20)
      nfaces = 6
      n(1:8) = t2l(elem%nodes(1:8))
      face_nodes(:,1) = [n(1),n(2),n(3),n(4)]
      face_nodes(:,2) = [n(5),n(6),n(7),n(8)]
      face_nodes(:,3) = [n(1),n(2),n(6),n(5)]
      face_nodes(:,4) = [n(3),n(4),n(8),n(7)]
      face_nodes(:,5) = [n(1),n(4),n(8),n(5)]
      face_nodes(:,6) = [n(2),n(3),n(7),n(6)]
    case (GMSH_PRISM6)
      nfaces = 5
      n(1:6) = t2l(elem%nodes(1:6))
      face_nodes(:,1) = [n(1),n(2),n(3),0]
      face_nodes(:,2) = [n(4),n(5),n(6),0]
      face_nodes(:,3) = [n(1),n(2),n(5),n(4)]
      face_nodes(:,4) = [n(2),n(3),n(6),n(5)]
      face_nodes(:,5) = [n(1),n(3),n(6),n(4)]
    case (GMSH_PYR5)
      nfaces = 5
      n(1:5) = t2l(elem%nodes(1:5))
      face_nodes(:,1) = [n(1),n(2),n(3),n(4)]
      face_nodes(:,2) = [n(1),n(2),n(5),0]
      face_nodes(:,3) = [n(2),n(3),n(5),0]
      face_nodes(:,4) = [n(3),n(4),n(5),0]
      face_nodes(:,5) = [n(1),n(4),n(5),0]
    end select
  end subroutine get_element_faces

  pure logical function faces_match(a, b)
    integer(ip), intent(in) :: a(:), b(:)
    integer(ip) :: i, j, nmatch
    nmatch = 0
    do i = 1, size(a)
      do j = 1, size(b)
        if (a(i) == b(j)) nmatch = nmatch + 1
      end do
    end do
    faces_match = (nmatch == size(a))
  end function faces_match

  pure subroutine canonical_face_nodes(nodes_in, n, nodes_out)
    integer(ip), intent(in)  :: nodes_in(:)
    integer(ip), intent(in)  :: n
    integer(ip), intent(out) :: nodes_out(4)
    integer(ip) :: i, j, tmp

    nodes_out = 0
    nodes_out(1:n) = nodes_in(1:n)
    do i = 2, n
      tmp = nodes_out(i)
      j = i - 1
      do while (j >= 1)
        if (nodes_out(j) <= tmp) exit
        nodes_out(j+1) = nodes_out(j)
        j = j - 1
      end do
      nodes_out(j+1) = tmp
    end do
  end subroutine canonical_face_nodes

  pure integer(ip) function face_hash(nodes, table_size)
    integer(ip), intent(in) :: nodes(4)
    integer(ip), intent(in) :: table_size
    integer(lp) :: h

    h = int(nodes(1), lp) * 73856093_lp + &
        int(nodes(2), lp) * 19349663_lp + &
        int(nodes(3), lp) * 83492791_lp + &
        int(nodes(4), lp) * 2654435761_lp
    if (h < 0_lp) h = -h
    face_hash = int(mod(h, int(table_size, lp)), ip) + 1
  end function face_hash

  pure subroutine cross_product(a, b, c)
    real(rp), intent(in)  :: a(3), b(3)
    real(rp), intent(out) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine cross_product

  pure integer function gmsh_nodes_per_element(etype)
    integer(ip), intent(in) :: etype
    select case (etype)
    case (GMSH_LINE2);  gmsh_nodes_per_element = 2
    case (GMSH_TRI3);   gmsh_nodes_per_element = 3
    case (GMSH_QUAD4);  gmsh_nodes_per_element = 4
    case (GMSH_TET4);   gmsh_nodes_per_element = 4
    case (GMSH_HEX8);   gmsh_nodes_per_element = 8
    case (GMSH_PRISM6); gmsh_nodes_per_element = 6
    case (GMSH_PYR5);   gmsh_nodes_per_element = 5
    case (GMSH_TET10);  gmsh_nodes_per_element = 10
    case (GMSH_HEX20);  gmsh_nodes_per_element = 20
    case (GMSH_POINT);  gmsh_nodes_per_element = 1
    case default;       gmsh_nodes_per_element = -1  ! unsupported: caller must stop
    end select
  end function gmsh_nodes_per_element

  pure logical function is_volume_element(etype)
    integer(ip), intent(in) :: etype
    is_volume_element = (etype == GMSH_TET4 .or. etype == GMSH_HEX8 .or. &
                         etype == GMSH_PRISM6 .or. etype == GMSH_PYR5 .or. &
                         etype == GMSH_TET10 .or. etype == GMSH_HEX20)
  end function is_volume_element

  pure logical function is_bnd_element(etype)
    integer(ip), intent(in) :: etype
    is_bnd_element = (etype == GMSH_TRI3 .or. etype == GMSH_QUAD4)
  end function is_bnd_element

end module mod_gmsh_reader
