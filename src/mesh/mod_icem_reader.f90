! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> ICEM/Fluent .msh polyhedral mesh reader (ASCII and binary).
module mod_icem_reader
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use iso_fortran_env, only: int8, int32, int64, real32, real64
  implicit none
  private

  public :: icem_read

  ! Fluent face-zone bc-type numeric codes (subset relevant to classification/mapping).
  integer(ip), parameter :: FL_BC_INTERIOR        = 2
  integer(ip), parameter :: FL_BC_WALL            = 3
  integer(ip), parameter :: FL_BC_PRESSURE_INLET  = 4
  integer(ip), parameter :: FL_BC_PRESSURE_OUTLET = 5
  integer(ip), parameter :: FL_BC_SYMMETRY        = 7
  integer(ip), parameter :: FL_BC_PRESSURE_FARFIELD = 9
  integer(ip), parameter :: FL_BC_VELOCITY_INLET  = 10
  integer(ip), parameter :: FL_BC_PERIODIC        = 12
  integer(ip), parameter :: FL_BC_MASS_FLOW_INLET = 20
  integer(ip), parameter :: FL_BC_INTERFACE       = 24
  integer(ip), parameter :: FL_BC_OUTFLOW         = 36
  integer(ip), parameter :: FL_BC_AXIS            = 37

  ! Fluent cell-type tags -> VTK cell codes.
  integer(ip), parameter :: FL_CELL_TET  = 2, FL_CELL_HEX  = 4
  integer(ip), parameter :: FL_CELL_PYR  = 5, FL_CELL_WEDGE = 6
  integer(ip), parameter :: FL_CELL_POLY = 7

  !> Header/data span for a single (10)/(12)/(13)-style zone section.
  type :: ZoneRec
    integer(ip) :: zone_id = 0, first_idx = 0, last_idx = 0, f4 = -1, f5 = -1
    integer(ip) :: data_a = 0, data_b = -1  ! byte span of the data list; data_b < data_a means no data
    logical     :: is_binary = .false.      ! section index was the binary variant (30XX), not ASCII (XX)
  end type ZoneRec

contains

  subroutine icem_read(filename, msh)
    character(len=*), intent(in)    :: filename
    type(Mesh),        intent(inout) :: msh

    character(len=:), allocatable :: buf
    integer(ip) :: dimension_decl
    integer(ip) :: n_zone_names, n_node_recs, n_cell_recs, n_face_recs
    integer(ip), allocatable :: zn_id(:)
    character(len=64), allocatable :: zn_type(:), zn_name(:)
    type(ZoneRec), allocatable :: node_recs(:), cell_recs(:), face_recs(:)
    real(rp),    allocatable :: points(:,:)
    integer(ip), allocatable :: cell_type_tag(:)
    integer(ip), allocatable :: face_pts(:), face_offset(:), face_owner(:), face_neigh(:)
    integer(ip) :: n_points, n_cells, n_faces_total, n_int_faces
    integer(ip) :: z, f, c

    call read_whole_file(filename, buf)

    ! Pass 1: count zone-name records and node/cell/face zone headers.
    call scan_sections(buf, dimension_decl, n_zone_names, n_node_recs, n_cell_recs, n_face_recs, &
                        zn_id, zn_type, zn_name, node_recs, cell_recs, face_recs, .false.)
    allocate(zn_id(max(n_zone_names,1)), zn_type(max(n_zone_names,1)), zn_name(max(n_zone_names,1)))
    allocate(node_recs(max(n_node_recs,1)), cell_recs(max(n_cell_recs,1)), face_recs(max(n_face_recs,1)))

    ! Pass 2: fill the zone tables.
    call scan_sections(buf, dimension_decl, n_zone_names, n_node_recs, n_cell_recs, n_face_recs, &
                        zn_id, zn_type, zn_name, node_recs, cell_recs, face_recs, .true.)

    if (dimension_decl /= 3) error stop 'icem_read: only 3-D Fluent meshes are supported'

    ! -- Nodes --
    n_points = 0
    do z = 1, n_node_recs
      if (node_recs(z)%zone_id == 0) cycle
      n_points = max(n_points, node_recs(z)%last_idx)
    end do
    if (n_points <= 0) error stop 'icem_read: no node data found in mesh file'
    allocate(points(3, n_points))
    points = 0.0_rp
    do z = 1, n_node_recs
      if (node_recs(z)%zone_id == 0) cycle
      call fill_node_zone(buf, node_recs(z), points)
    end do

    ! -- Cells --
    n_cells = 0
    do z = 1, n_cell_recs
      if (cell_recs(z)%zone_id == 0) cycle
      n_cells = max(n_cells, cell_recs(z)%last_idx)
    end do
    if (n_cells <= 0) error stop 'icem_read: no cell data found in mesh file'
    allocate(cell_type_tag(n_cells), source=0_ip)
    do z = 1, n_cell_recs
      if (cell_recs(z)%zone_id == 0) cycle
      call fill_cell_zone(buf, cell_recs(z), cell_type_tag)
    end do
    if (any(cell_type_tag == 0)) error stop 'icem_read: not all cells were assigned a type'

    ! -- Faces: classify zones as interior vs boundary, size arrays --
    n_int_faces = 0
    n_faces_total = 0
    do z = 1, n_face_recs
      if (face_recs(z)%zone_id == 0) cycle
      n_faces_total = n_faces_total + (face_recs(z)%last_idx - face_recs(z)%first_idx + 1)
      if (face_recs(z)%f4 == FL_BC_INTERIOR) &
        n_int_faces = n_int_faces + (face_recs(z)%last_idx - face_recs(z)%first_idx + 1)
    end do
    if (n_faces_total <= 0) error stop 'icem_read: no face data found in mesh file'

    allocate(face_owner(n_faces_total), face_neigh(n_faces_total))
    allocate(face_offset(n_faces_total + 1))
    msh%n_faces     = n_faces_total
    msh%n_int_faces = n_int_faces
    msh%n_bnd_faces = n_faces_total - n_int_faces
    msh%n_cells     = n_cells
    msh%n_ghost     = 0
    allocate(msh%cells(n_cells))
    allocate(msh%faces(n_faces_total))
    allocate(msh%global_id(n_cells))
    do c = 1, n_cells
      msh%cells(c)%centroid = 0.0_rp
      msh%cells(c)%volume   = 0.0_rp
      msh%global_id(c)      = int(c, lp)
    end do

    ! Boundary patches: one per non-interior face zone, in file order.
    msh%n_patches = count(face_recs(1:n_face_recs)%zone_id /= 0 .and. &
                           face_recs(1:n_face_recs)%f4 /= FL_BC_INTERIOR)
    allocate(msh%patches(msh%n_patches))

    call fill_faces(buf, face_recs, n_face_recs, n_int_faces, n_faces_total, &
                     zn_id, zn_type, zn_name, n_zone_names, &
                     face_owner, face_neigh, face_pts, face_offset, msh)

    do f = 1, n_faces_total
      msh%faces(f)%owner     = face_owner(f)
      msh%faces(f)%neighbour = face_neigh(f)
    end do

    call build_cell_topology(msh, points, face_pts, face_offset, face_owner, face_neigh, &
                              n_points, n_faces_total, n_cells, cell_type_tag)
    allocate(msh%face_offsets(n_faces_total + 1))
    allocate(msh%face_vertices(size(face_pts)))
    msh%face_offsets  = face_offset
    msh%face_vertices = face_pts

    do f = 1, n_faces_total
      call compute_polygon_face_geometry(f, face_pts, face_offset, points, &
                                          msh%faces(f)%centre, msh%faces(f)%normal, msh%faces(f)%area)
      msh%faces(f)%sf = msh%faces(f)%area * msh%faces(f)%normal
    end do
    call compute_polyhedron_cell_moments(msh, n_cells, n_faces_total, face_pts, face_offset, points)

    call validate_icem_boundary_connectivity(msh)
  end subroutine icem_read

  !> Read a face zone's node/owner/neighbour records, split into interior-first + per-patch boundary layout.
  subroutine fill_faces(buf, face_recs, n_face_recs, n_int_faces, n_faces_total, &
                         zn_id, zn_type, zn_name, n_zone_names, &
                         face_owner, face_neigh, face_pts, face_offset, msh)
    character(len=*), intent(in) :: buf
    type(ZoneRec),     intent(in) :: face_recs(:)
    integer(ip),       intent(in) :: n_face_recs, n_int_faces, n_faces_total, n_zone_names
    integer(ip),       intent(in) :: zn_id(:)
    character(len=*),  intent(in) :: zn_type(:), zn_name(:)
    integer(ip), allocatable, intent(out) :: face_owner(:), face_neigh(:)
    integer(ip), allocatable, intent(out) :: face_pts(:), face_offset(:)
    type(Mesh),        intent(inout) :: msh

    integer(ip), allocatable :: rec_npts(:)
    integer(ip) :: z, i, nf, slot_int, slot_bnd, p, tot_pts

    allocate(face_owner(n_faces_total), face_neigh(n_faces_total))
    allocate(rec_npts(n_faces_total))

    ! First sub-pass: extract node-count-per-face for every zone, to build the CSR offset table
    ! in final (interior-first, then per-patch boundary) slot order.
    slot_int = 0
    slot_bnd = n_int_faces
    do z = 1, n_face_recs
      if (face_recs(z)%zone_id == 0) cycle
      nf = face_recs(z)%last_idx - face_recs(z)%first_idx + 1
      if (face_recs(z)%f4 == FL_BC_INTERIOR) then
        call read_face_zone_records(buf, face_recs(z), nf, rec_npts, slot_int, face_owner, face_neigh)
        slot_int = slot_int + nf
      end if
    end do

    p = 0
    do z = 1, n_face_recs
      if (face_recs(z)%zone_id == 0) cycle
      if (face_recs(z)%f4 == FL_BC_INTERIOR) cycle
      nf = face_recs(z)%last_idx - face_recs(z)%first_idx + 1
      p = p + 1
      call read_face_zone_records(buf, face_recs(z), nf, rec_npts, slot_bnd, face_owner, face_neigh)
      msh%patches(p)%name       = zone_display_name(face_recs(z)%zone_id, zn_id, zn_name, n_zone_names)
      msh%patches(p)%bc_type    = zone_bc_keyword(face_recs(z)%zone_id, face_recs(z)%f4, zn_id, zn_type, n_zone_names)
      msh%patches(p)%bc_id      = fluent_bc_to_id(face_recs(z)%f4, msh%patches(p)%bc_type)
      msh%patches(p)%start_face = (slot_bnd - n_int_faces) + 1
      msh%patches(p)%n_faces    = nf
      do i = slot_bnd + 1, slot_bnd + nf
        msh%faces(i)%patch_id = p
      end do
      slot_bnd = slot_bnd + nf
    end do

    allocate(face_offset(n_faces_total + 1))
    face_offset(1) = 1
    do i = 1, n_faces_total
      face_offset(i + 1) = face_offset(i) + rec_npts(i)
    end do
    tot_pts = face_offset(n_faces_total + 1) - 1
    allocate(face_pts(tot_pts))

    ! Second sub-pass: fill node lists at the now-known offsets.
    slot_int = 0
    slot_bnd = n_int_faces
    do z = 1, n_face_recs
      if (face_recs(z)%zone_id == 0) cycle
      nf = face_recs(z)%last_idx - face_recs(z)%first_idx + 1
      if (face_recs(z)%f4 == FL_BC_INTERIOR) then
        call read_face_zone_points(buf, face_recs(z), nf, face_offset, slot_int, face_pts)
        slot_int = slot_int + nf
      end if
    end do
    do z = 1, n_face_recs
      if (face_recs(z)%zone_id == 0) cycle
      if (face_recs(z)%f4 == FL_BC_INTERIOR) cycle
      nf = face_recs(z)%last_idx - face_recs(z)%first_idx + 1
      call read_face_zone_points(buf, face_recs(z), nf, face_offset, slot_bnd, face_pts)
      slot_bnd = slot_bnd + nf
    end do
  end subroutine fill_faces

  !> First sub-pass over a face zone's records: node-count-per-face + owner/neighbour, into slots [slot0+1..].
  subroutine read_face_zone_records(buf, zr, nf, rec_npts, slot0, face_owner, face_neigh)
    character(len=*), intent(in)    :: buf
    type(ZoneRec),     intent(in)   :: zr
    integer(ip),       intent(in)   :: nf, slot0
    integer(ip),       intent(inout):: rec_npts(:), face_owner(:), face_neigh(:)

    integer(ip) :: is_mixed, fixed_n
    integer(ip) :: i, pos, nn, c0, c1

    is_mixed = merge(1, 0, zr%f5 == 0 .or. zr%f5 == 5)
    fixed_n  = 0
    select case (zr%f5)
    case (2); fixed_n = 2
    case (3); fixed_n = 3
    case (4); fixed_n = 4
    end select

    if (zr%data_b < zr%data_a) error stop 'icem_read: face zone declares no connectivity data'

    if (zr%is_binary) then
      pos = zr%data_a
      do i = 1, nf
        if (is_mixed == 1) then
          nn = read_bin_i32(buf, pos)
        else
          nn = fixed_n
        end if
        rec_npts(slot0 + i) = nn
        pos = pos + nn * 4
        c0 = read_bin_i32(buf, pos)
        c1 = read_bin_i32(buf, pos)
        face_owner(slot0 + i) = c0
        face_neigh(slot0 + i) = merge(-1_ip, c1, c1 == 0)
      end do
    else
      pos = zr%data_a
      do i = 1, nf
        if (is_mixed == 1) then
          call next_int_ascii(buf, pos, zr%data_b, nn)
        else
          nn = fixed_n
        end if
        rec_npts(slot0 + i) = nn
        call skip_n_tokens_ascii(buf, pos, zr%data_b, nn)
        call next_int_ascii(buf, pos, zr%data_b, c0)
        call next_int_ascii(buf, pos, zr%data_b, c1)
        face_owner(slot0 + i) = c0
        face_neigh(slot0 + i) = merge(-1_ip, c1, c1 == 0)
      end do
    end if
  end subroutine read_face_zone_records

  !> Second sub-pass: fill node-id lists for a face zone at known CSR offsets.
  subroutine read_face_zone_points(buf, zr, nf, face_offset, slot0, face_pts)
    character(len=*), intent(in)    :: buf
    type(ZoneRec),     intent(in)   :: zr
    integer(ip),       intent(in)   :: nf, slot0
    integer(ip),       intent(in)   :: face_offset(:)
    integer(ip),       intent(inout):: face_pts(:)

    integer(ip) :: is_mixed, fixed_n
    integer(ip) :: i, j, pos, nn, c0, c1, w

    is_mixed = merge(1, 0, zr%f5 == 0 .or. zr%f5 == 5)
    fixed_n  = 0
    select case (zr%f5)
    case (2); fixed_n = 2
    case (3); fixed_n = 3
    case (4); fixed_n = 4
    end select

    if (zr%is_binary) then
      pos = zr%data_a
      do i = 1, nf
        if (is_mixed == 1) then
          nn = read_bin_i32(buf, pos)
        else
          nn = fixed_n
        end if
        w = face_offset(slot0 + i)
        do j = 1, nn
          face_pts(w + j - 1) = read_bin_i32(buf, pos)
        end do
        c0 = read_bin_i32(buf, pos)
        c1 = read_bin_i32(buf, pos)
      end do
    else
      pos = zr%data_a
      do i = 1, nf
        if (is_mixed == 1) then
          call next_int_ascii(buf, pos, zr%data_b, nn)
        else
          nn = fixed_n
        end if
        w = face_offset(slot0 + i)
        do j = 1, nn
          call next_int_ascii(buf, pos, zr%data_b, face_pts(w + j - 1))
        end do
        call next_int_ascii(buf, pos, zr%data_b, c0)
        call next_int_ascii(buf, pos, zr%data_b, c1)
      end do
    end if
  end subroutine read_face_zone_points

  subroutine fill_node_zone(buf, zr, points)
    character(len=*), intent(in)    :: buf
    type(ZoneRec),     intent(in)   :: zr
    real(rp),          intent(inout):: points(:,:)

    integer(ip) :: nd, i, d, pos
    real(rp)    :: v

    nd = zr%f5
    if (nd /= 2 .and. nd /= 3) nd = 3
    if (zr%data_b < zr%data_a) error stop 'icem_read: node zone declares no coordinate data'

    if (zr%is_binary) then
      pos = zr%data_a
      do i = zr%first_idx, zr%last_idx
        do d = 1, nd
          points(d, i) = read_bin_f64(buf, pos)
        end do
        if (nd == 2) points(3, i) = 0.0_rp
      end do
    else
      pos = zr%data_a
      do i = zr%first_idx, zr%last_idx
        do d = 1, nd
          call next_real_ascii(buf, pos, zr%data_b, v)
          points(d, i) = v
        end do
        if (nd == 2) points(3, i) = 0.0_rp
      end do
    end if
  end subroutine fill_node_zone

  subroutine fill_cell_zone(buf, zr, cell_type_tag)
    character(len=*), intent(in)    :: buf
    type(ZoneRec),     intent(in)   :: zr
    integer(ip),       intent(inout):: cell_type_tag(:)

    integer(ip) :: i, pos, tag

    if (zr%f5 /= 0) then
      ! Uniform zone: every cell shares the element-type declared in the header.
      do i = zr%first_idx, zr%last_idx
        cell_type_tag(i) = zr%f5
      end do
      return
    end if

    if (zr%data_b < zr%data_a) error stop 'icem_read: mixed cell zone declares no type data'
    if (zr%is_binary) then
      pos = zr%data_a
      do i = zr%first_idx, zr%last_idx
        tag = read_bin_i32(buf, pos)
        cell_type_tag(i) = tag
      end do
    else
      pos = zr%data_a
      do i = zr%first_idx, zr%last_idx
        call next_int_ascii(buf, pos, zr%data_b, tag)
        cell_type_tag(i) = tag
      end do
    end if
  end subroutine fill_cell_zone

  !> Build msh%vertices/cell_offsets/cell_vertices/cell_types from face-cell incidence (mirrors
  !> mod_foam_reader's build_cell_topology_metadata, but preserves native per-cell VTK types).
  subroutine build_cell_topology(msh, points, face_pts, face_offset, owner, neighbour, &
                                  n_points, n_faces, n_cells, cell_type_tag)
    type(Mesh),  intent(inout) :: msh
    real(rp),    intent(in)    :: points(:,:)
    integer(ip), intent(in)    :: face_pts(:), face_offset(:)
    integer(ip), intent(in)    :: owner(:), neighbour(:)
    integer(ip), intent(in)    :: n_points, n_faces, n_cells
    integer(ip), intent(in)    :: cell_type_tag(:)

    integer(ip) :: f, c, idx, pos, total_incidence, write_pos, nverts
    integer(ip), allocatable :: cell_face_count(:), cell_face_offset(:), cell_face_ids(:)
    integer(ip), allocatable :: fill_count(:), cell_vertex_count(:), vertex_marker(:)

    msh%n_vertices = n_points
    allocate(msh%vertices(3, n_points))
    msh%vertices = points

    allocate(cell_face_count(n_cells), source=0_ip)
    do f = 1, n_faces
      cell_face_count(owner(f)) = cell_face_count(owner(f)) + 1
      if (neighbour(f) >= 1) cell_face_count(neighbour(f)) = cell_face_count(neighbour(f)) + 1
    end do

    allocate(cell_face_offset(n_cells + 1))
    cell_face_offset(1) = 1
    do c = 1, n_cells
      cell_face_offset(c + 1) = cell_face_offset(c) + cell_face_count(c)
    end do
    total_incidence = cell_face_offset(n_cells + 1) - 1
    allocate(cell_face_ids(total_incidence), fill_count(n_cells), source=0_ip)

    do f = 1, n_faces
      c = owner(f)
      pos = cell_face_offset(c) + fill_count(c)
      cell_face_ids(pos) = f
      fill_count(c) = fill_count(c) + 1
      if (neighbour(f) >= 1) then
        c = neighbour(f)
        pos = cell_face_offset(c) + fill_count(c)
        cell_face_ids(pos) = f
        fill_count(c) = fill_count(c) + 1
      end if
    end do
    deallocate(cell_face_count, fill_count)

    allocate(cell_vertex_count(n_cells), vertex_marker(n_points), source=0_ip)
    do c = 1, n_cells
      nverts = 0
      do idx = cell_face_offset(c), cell_face_offset(c + 1) - 1
        f = cell_face_ids(idx)
        do pos = face_offset(f), face_offset(f + 1) - 1
          if (vertex_marker(face_pts(pos)) /= c) then
            vertex_marker(face_pts(pos)) = c
            nverts = nverts + 1
          end if
        end do
      end do
      cell_vertex_count(c) = nverts
    end do

    allocate(msh%cell_offsets(n_cells + 1))
    msh%cell_offsets(1) = 1
    do c = 1, n_cells
      msh%cell_offsets(c + 1) = msh%cell_offsets(c) + cell_vertex_count(c)
    end do
    allocate(msh%cell_vertices(msh%cell_offsets(n_cells + 1) - 1))
    allocate(msh%cell_types(n_cells))
    vertex_marker = 0

    do c = 1, n_cells
      write_pos = msh%cell_offsets(c)
      nverts = 0
      do idx = cell_face_offset(c), cell_face_offset(c + 1) - 1
        f = cell_face_ids(idx)
        do pos = face_offset(f), face_offset(f + 1) - 1
          if (vertex_marker(face_pts(pos)) /= c) then
            vertex_marker(face_pts(pos)) = c
            msh%cell_vertices(write_pos + nverts) = face_pts(pos)
            nverts = nverts + 1
          end if
        end do
      end do
      msh%cell_types(c) = fluent_celltype_to_vtk(cell_type_tag(c))
    end do

    deallocate(cell_face_offset, cell_face_ids, cell_vertex_count, vertex_marker)
  end subroutine build_cell_topology

  pure integer(ip) function fluent_celltype_to_vtk(tag) result(vtk)
    integer(ip), intent(in) :: tag
    select case (tag)
    case (FL_CELL_TET);   vtk = 10
    case (FL_CELL_HEX);   vtk = 12
    case (FL_CELL_WEDGE); vtk = 13
    case (FL_CELL_PYR);   vtk = 14
    case default;         vtk = 42  ! polyhedral (7), or any unrecognized tag
    end select
  end function fluent_celltype_to_vtk

  pure integer(ip) function fluent_bc_to_id(bc_code, bc_keyword) result(bc_id)
    integer(ip),      intent(in) :: bc_code
    character(len=*), intent(in) :: bc_keyword

    ! Prefer the named BC keyword (from a (45)/(39) zone-name section) when available.
    select case (trim(bc_keyword))
    case ('wall');                             bc_id = BC_WALL_MODEL;  return
    case ('velocity-inlet', 'mass-flow-inlet', 'pressure-inlet'); bc_id = BC_INLET; return
    case ('pressure-outlet', 'outflow');       bc_id = BC_OUTLET;      return
    case ('symmetry');                         bc_id = BC_SYMMETRY;    return
    case ('axis');                             bc_id = BC_SYMMETRY;    return
    case ('periodic');                         bc_id = BC_PERIODIC;    return
    end select

    select case (bc_code)
    case (FL_BC_WALL);                          bc_id = BC_WALL_MODEL
    case (FL_BC_VELOCITY_INLET, FL_BC_MASS_FLOW_INLET, FL_BC_PRESSURE_INLET, FL_BC_PRESSURE_FARFIELD)
      bc_id = BC_INLET
    case (FL_BC_PRESSURE_OUTLET, FL_BC_OUTFLOW); bc_id = BC_OUTLET
    case (FL_BC_SYMMETRY, FL_BC_AXIS);           bc_id = BC_SYMMETRY
    case (FL_BC_PERIODIC);                       bc_id = BC_PERIODIC
    case default;                                bc_id = BC_ZEROGRADIENT
    end select
  end function fluent_bc_to_id

  pure function zone_display_name(zone_id, zn_id, zn_name, n_zone_names) result(nm)
    integer(ip),      intent(in) :: zone_id, n_zone_names
    integer(ip),      intent(in) :: zn_id(:)
    character(len=*), intent(in) :: zn_name(:)
    character(len=128) :: nm

    integer(ip) :: k
    nm = ''
    do k = 1, n_zone_names
      if (zn_id(k) == zone_id) then
        nm = trim(zn_name(k))
        return
      end if
    end do
    write(nm, '(A,I0)') 'zone_', zone_id
  end function zone_display_name

  pure function zone_bc_keyword(zone_id, bc_code, zn_id, zn_type, n_zone_names) result(kw)
    integer(ip),      intent(in) :: zone_id, bc_code, n_zone_names
    integer(ip),      intent(in) :: zn_id(:)
    character(len=*), intent(in) :: zn_type(:)
    character(len=64) :: kw

    integer(ip) :: k
    kw = ''
    do k = 1, n_zone_names
      if (zn_id(k) == zone_id) then
        kw = trim(zn_type(k))
        return
      end if
    end do
    select case (bc_code)
    case (FL_BC_WALL);            kw = 'wall'
    case (FL_BC_VELOCITY_INLET);  kw = 'velocity-inlet'
    case (FL_BC_MASS_FLOW_INLET); kw = 'mass-flow-inlet'
    case (FL_BC_PRESSURE_INLET);  kw = 'pressure-inlet'
    case (FL_BC_PRESSURE_OUTLET); kw = 'pressure-outlet'
    case (FL_BC_PRESSURE_FARFIELD); kw = 'pressure-far-field'
    case (FL_BC_SYMMETRY);        kw = 'symmetry'
    case (FL_BC_AXIS);            kw = 'axis'
    case (FL_BC_PERIODIC);        kw = 'periodic'
    case (FL_BC_OUTFLOW);         kw = 'outflow'
    case default;                 kw = 'wall'
    end select
  end function zone_bc_keyword

  subroutine validate_icem_boundary_connectivity(msh)
    type(Mesh), intent(in) :: msh
    integer(ip) :: p, f_bnd, f_glob, n_assigned

    n_assigned = 0
    do f_bnd = 1, msh%n_bnd_faces
      f_glob = msh%n_int_faces + f_bnd
      if (msh%faces(f_glob)%patch_id < 1 .or. msh%faces(f_glob)%patch_id > msh%n_patches) cycle
      n_assigned = n_assigned + 1
    end do
    if (n_assigned /= msh%n_bnd_faces) &
      error stop 'icem_read: boundary connectivity invalid (unmapped boundary faces)'

    do p = 1, msh%n_patches
      if (msh%patches(p)%start_face < 1 .or. msh%patches(p)%n_faces < 0) &
        error stop 'icem_read: boundary connectivity invalid (bad patch start_face/n_faces)'
      if (msh%patches(p)%start_face + msh%patches(p)%n_faces - 1 > msh%n_bnd_faces) &
        error stop 'icem_read: boundary connectivity invalid (patch range exceeds boundary faces)'
    end do

    write(*, '(A,I0,A,I0)') ' ICEM/Fluent mesh: patches=', msh%n_patches, ', boundary_faces=', msh%n_bnd_faces
    do p = 1, msh%n_patches
      write(*, '(A,I0,A,A,A,A,A,I0)') '   patch ', p, ': name=', trim(msh%patches(p)%name), &
          ', type=', trim(msh%patches(p)%bc_type), ', bc_id=', msh%patches(p)%bc_id
    end do
  end subroutine validate_icem_boundary_connectivity

  ! ==========================================================================
  ! Low-level Fluent .msh section scanner
  ! ==========================================================================

  subroutine read_whole_file(filename, buf)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable, intent(out) :: buf
    integer(ip) :: iunit, ios
    integer(int64) :: nbytes

    open(newunit=iunit, file=trim(filename), status='old', access='stream', form='unformatted', &
         action='read', iostat=ios)
    if (ios /= 0) error stop 'icem_read: cannot open '//trim(filename)
    inquire(unit=iunit, size=nbytes)
    if (nbytes <= 0_int64) error stop 'icem_read: empty file '//trim(filename)
    allocate(character(len=int(nbytes,ip)) :: buf)
    read(iunit, pos=1, iostat=ios) buf
    if (ios /= 0) error stop 'icem_read: failed reading '//trim(filename)
    close(iunit)
  end subroutine read_whole_file

  !> Walk all top-level "(<index> ...)" sections; count (do_collect=F) or fill (do_collect=T) tables.
  subroutine scan_sections(buf, dimension_decl, n_zone_names, n_node_recs, n_cell_recs, n_face_recs, &
                            zn_id, zn_type, zn_name, node_recs, cell_recs, face_recs, do_collect)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(out)  :: dimension_decl
    integer(ip),       intent(inout):: n_zone_names, n_node_recs, n_cell_recs, n_face_recs
    integer(ip), allocatable, intent(inout) :: zn_id(:)
    character(len=*),  allocatable, intent(inout) :: zn_type(:), zn_name(:)
    type(ZoneRec), allocatable, intent(inout) :: node_recs(:), cell_recs(:), face_recs(:)
    logical,           intent(in)   :: do_collect

    integer(ip) :: n, pos, sec_open, sec_close, idx, base_idx, cpos
    integer(ip) :: cnt_zn, cnt_node, cnt_cell, cnt_face

    n = len(buf)
    pos = 1
    dimension_decl = 0
    cnt_zn = 0; cnt_node = 0; cnt_cell = 0; cnt_face = 0

    do
      call skip_ws(buf, pos)
      if (pos > n) exit
      if (buf(pos:pos) /= '(') error stop 'icem_read: expected top-level "(" in mesh file'
      sec_open = pos
      cpos = sec_open + 1
      call skip_ws(buf, cpos)
      call parse_dec_int(buf, cpos, idx)

      ! Zone sections (10/12/13, or a binary variant such as 3010/2013 -- observed prefixes vary
      ! by field/precision) may contain raw binary payload bytes that coincidentally equal '(' or
      ! ')', so their extent must NOT be found by generic paren-matching; handle_zone_section
      ! locates it from the declared record count/size instead. Binary variants are recognised by
      ! base_idx = mod(idx, 1000) falling in the known zone-section set, for any idx >= 1000.
      base_idx = idx
      if (idx >= 1000) base_idx = mod(idx, 1000)

      select case (base_idx)
      case (2)
        call skip_ws(buf, cpos)
        call parse_dec_int(buf, cpos, dimension_decl)
        sec_close = matching_paren(buf, sec_open)
      case (10, 12, 13)
        call handle_zone_section(buf, idx, cpos, &
                                  cnt_node, cnt_cell, cnt_face, node_recs, cell_recs, face_recs, do_collect, sec_close)
      case (39, 45)
        sec_close = matching_paren(buf, sec_open)
        cnt_zn = cnt_zn + 1
        if (do_collect) call parse_zone_name_section(buf, cpos, sec_close, zn_id(cnt_zn), zn_type(cnt_zn), zn_name(cnt_zn))
      case default
        ! comment(0), header(1), and any other/unhandled section: skip.
        sec_close = matching_paren(buf, sec_open)
      end select

      pos = sec_close + 1
    end do

    n_zone_names = cnt_zn
    n_node_recs  = cnt_node
    n_cell_recs  = cnt_cell
    n_face_recs  = cnt_face
  end subroutine scan_sections

  !> Parse a (10)/(12)/(13)-style "(header)(data)" or "(header)" section (ASCII or binary variant)
  !> into a ZoneRec, and return the position of its closing ")" via `sec_close`.
  subroutine handle_zone_section(buf, idx, cpos_in, cnt_node, cnt_cell, cnt_face, &
                                  node_recs, cell_recs, face_recs, do_collect, sec_close)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(in)   :: idx, cpos_in
    integer(ip),       intent(inout):: cnt_node, cnt_cell, cnt_face
    type(ZoneRec), allocatable, intent(inout) :: node_recs(:), cell_recs(:), face_recs(:)
    logical,            intent(in)  :: do_collect
    integer(ip),        intent(out) :: sec_close

    integer(ip) :: cpos, hdr_open, hdr_close, data_open, data_close, base_idx
    integer(ip), allocatable :: fields(:)
    type(ZoneRec) :: zr
    logical :: is_bin

    is_bin = (idx >= 1000)
    base_idx = idx
    if (is_bin) base_idx = mod(idx, 1000)

    cpos = cpos_in
    call skip_ws(buf, cpos)
    if (buf(cpos:cpos) /= '(') error stop 'icem_read: expected zone header "("'
    hdr_open = cpos
    hdr_close = matching_paren(buf, hdr_open)  ! header tuple is plain ASCII hex fields: safe to paren-match
    call parse_hex_fields(buf, hdr_open + 1, hdr_close - 1, fields)

    zr%zone_id   = fields(1)
    zr%first_idx = merge(fields(2), 0, size(fields) >= 2)
    zr%last_idx  = merge(fields(3), 0, size(fields) >= 3)
    zr%f4        = merge(fields(4), -1, size(fields) >= 4)
    zr%f5        = merge(fields(5), -1, size(fields) >= 5)
    zr%data_a = 0
    zr%data_b = -1
    zr%is_binary = is_bin

    cpos = hdr_close + 1
    call skip_ws(buf, cpos)
    if (buf(cpos:cpos) == '(') then
      data_open = cpos
      if (is_bin) then
        data_close = binary_data_close(buf, base_idx, data_open, zr)
      else
        data_close = matching_paren(buf, data_open)  ! ASCII payload: digits/whitespace only, safe
      end if
      if (data_close - data_open > 1) then
        zr%data_a = data_open + 1
        zr%data_b = data_close - 1
      end if
      cpos = data_close + 1
    end if

    if (is_bin) then
      ! Binary sections carry a literal "End of Binary Section <idx>" comment between the data
      ! payload and the closing ")" (not whitespace). The payload itself was already consumed by
      ! binary_data_close above (which validates its own end against a real ")"), so everything
      ! from here to the next ")" is plain ASCII trailer text -- safe to scan for literally.
      cpos = find_char_fwd(buf, ')', cpos)
      if (cpos < 0) error stop 'icem_read: could not find closing ")" after binary section trailer'
    else
      call skip_ws(buf, cpos)
      if (buf(cpos:cpos) /= ')') error stop 'icem_read: expected zone section closing ")"'
    end if
    sec_close = cpos

    select case (base_idx)
    case (10)
      cnt_node = cnt_node + 1
      if (do_collect) node_recs(cnt_node) = zr
    case (12)
      cnt_cell = cnt_cell + 1
      if (do_collect) cell_recs(cnt_cell) = zr
    case (13)
      cnt_face = cnt_face + 1
      if (do_collect) face_recs(cnt_face) = zr
    end select
  end subroutine handle_zone_section

  !> Compute the byte position of a binary zone's data-list closing ")" from the declared record
  !> count/size (uniform zones) or by walking variable-length records (mixed/polygonal faces).
  integer(ip) function binary_data_close(buf, base_idx, data_open, zr) result(data_close)
    character(len=*), intent(in) :: buf
    integer(ip),       intent(in) :: base_idx, data_open
    type(ZoneRec),     intent(in) :: zr
    integer(ip) :: nrec, nbytes, is_mixed, fixed_n, pos, i, nn, nd

    nrec = zr%last_idx - zr%first_idx + 1
    select case (base_idx)
    case (10)
      nd = zr%f5
      if (nd /= 2 .and. nd /= 3) nd = 3
      nbytes = nrec * nd * 8
    case (12)
      nbytes = nrec * 4  ! mixed cell-type tag list (one int32 per cell); only reached when f5==0
    case (13)
      is_mixed = merge(1, 0, zr%f5 == 0 .or. zr%f5 == 5)
      fixed_n = 0
      select case (zr%f5)
      case (2); fixed_n = 2
      case (3); fixed_n = 3
      case (4); fixed_n = 4
      end select
      if (is_mixed == 0) then
        nbytes = nrec * (fixed_n + 2) * 4
      else
        pos = data_open + 1
        do i = 1, nrec
          nn = read_bin_i32(buf, pos)
          pos = pos + nn * 4 + 8
        end do
        nbytes = pos - (data_open + 1)
      end if
    case default
      error stop 'icem_read: unsupported binary zone section index'
    end select

    ! The writer inserts a separator byte (observed: a single newline) between the last payload
    ! byte and the closing ")", so locate it by skipping whitespace rather than assuming adjacency;
    ! validate against the computed length instead of trusting a bare byte offset.
    data_close = data_open + 1 + nbytes
    call skip_ws(buf, data_close)
    if (data_close > len(buf) .or. buf(data_close:data_close) /= ')') &
      error stop 'icem_read: binary zone data length mismatch (no closing ")" at computed offset)'
  end function binary_data_close

  !> Parse "(zone-id zone-type zone-name)()" from a (39)/(45) section body.
  subroutine parse_zone_name_section(buf, cpos_in, sec_close, id_out, type_out, name_out)
    character(len=*), intent(in)  :: buf
    integer(ip),       intent(in) :: cpos_in, sec_close
    integer(ip),       intent(out):: id_out
    character(len=*),  intent(out):: type_out, name_out

    integer(ip) :: cpos, hdr_open, hdr_close, tok_a, tok_b

    cpos = cpos_in
    call skip_ws(buf, cpos)
    if (cpos > sec_close .or. buf(cpos:cpos) /= '(') error stop 'icem_read: expected zone-name header "("'
    hdr_open = cpos
    hdr_close = matching_paren(buf, hdr_open)

    cpos = hdr_open + 1
    call skip_ws(buf, cpos)
    call parse_dec_int(buf, cpos, id_out)
    call skip_ws(buf, cpos)
    call next_token_span(buf, cpos, hdr_close - 1, tok_a, tok_b)
    type_out = buf(tok_a:tok_b)
    call skip_ws(buf, cpos)
    call next_token_span(buf, cpos, hdr_close - 1, tok_a, tok_b)
    name_out = buf(tok_a:tok_b)
  end subroutine parse_zone_name_section

  !> Find the position of the ')' matching the '(' at `open_pos`, honouring quoted strings.
  integer(ip) function matching_paren(buf, open_pos) result(close_pos)
    character(len=*), intent(in) :: buf
    integer(ip),       intent(in) :: open_pos
    integer(ip) :: pos, depth, n
    logical :: in_quote

    n = len(buf)
    depth = 0
    in_quote = .false.
    pos = open_pos
    do while (pos <= n)
      if (buf(pos:pos) == '"') then
        in_quote = .not. in_quote
      else if (.not. in_quote) then
        if (buf(pos:pos) == '(') depth = depth + 1
        if (buf(pos:pos) == ')') then
          depth = depth - 1
          if (depth == 0) then
            close_pos = pos
            return
          end if
        end if
      end if
      pos = pos + 1
    end do
    error stop 'icem_read: unbalanced parentheses in mesh file'
  end function matching_paren

  !> First occurrence of `ch` at or after `start_pos` (plain forward scan, no depth-tracking).
  integer(ip) function find_char_fwd(buf, ch, start_pos) result(pos)
    character(len=*), intent(in) :: buf
    character(1),      intent(in) :: ch
    integer(ip),        intent(in) :: start_pos
    integer(ip) :: n
    n = len(buf)
    pos = start_pos
    do while (pos <= n)
      if (buf(pos:pos) == ch) return
      pos = pos + 1
    end do
    pos = -1
  end function find_char_fwd

  subroutine skip_ws(buf, pos)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip) :: n
    n = len(buf)
    do while (pos <= n)
      select case (buf(pos:pos))
      case (' ', char(9), char(10), char(13))
        pos = pos + 1
      case default
        return
      end select
    end do
  end subroutine skip_ws

  subroutine parse_dec_int(buf, pos, val)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(out)  :: val
    integer(ip) :: a
    call skip_ws(buf, pos)
    a = pos
    if (pos <= len(buf) .and. buf(pos:pos) == '-') pos = pos + 1
    do while (pos <= len(buf))
      if (buf(pos:pos) < '0' .or. buf(pos:pos) > '9') exit
      pos = pos + 1
    end do
    if (pos == a) error stop 'icem_read: expected decimal integer'
    read(buf(a:pos-1), *) val
  end subroutine parse_dec_int

  !> Parse whitespace-separated hexadecimal integers between [a,b] (Fluent zone header convention).
  subroutine parse_hex_fields(buf, a, b, fields)
    character(len=*), intent(in) :: buf
    integer(ip),       intent(in) :: a, b
    integer(ip), allocatable, intent(out) :: fields(:)
    integer(ip) :: pos, cnt, ta, tb, tmp(8)

    pos = a
    cnt = 0
    do
      call skip_ws_bounded(buf, pos, b)
      if (pos > b) exit
      ta = pos
      do while (pos <= b .and. is_hex_char(buf(pos:pos)))
        pos = pos + 1
      end do
      tb = pos - 1
      if (tb < ta) error stop 'icem_read: malformed zone header field'
      cnt = cnt + 1
      if (cnt > size(tmp)) error stop 'icem_read: too many zone header fields'
      tmp(cnt) = hex_to_int(buf(ta:tb))
    end do
    allocate(fields(cnt))
    fields = tmp(1:cnt)
  end subroutine parse_hex_fields

  logical function is_hex_char(c) result(r)
    character(1), intent(in) :: c
    r = (c >= '0' .and. c <= '9') .or. (c >= 'a' .and. c <= 'f') .or. (c >= 'A' .and. c <= 'F')
  end function is_hex_char

  integer(ip) function hex_to_int(tok) result(v)
    character(len=*), intent(in) :: tok
    integer(ip) :: i
    character(1) :: c
    v = 0
    do i = 1, len_trim(tok)
      c = tok(i:i)
      v = v * 16
      if (c >= '0' .and. c <= '9') then
        v = v + (iachar(c) - iachar('0'))
      else if (c >= 'a' .and. c <= 'f') then
        v = v + (iachar(c) - iachar('a') + 10)
      else
        v = v + (iachar(c) - iachar('A') + 10)
      end if
    end do
  end function hex_to_int

  subroutine skip_ws_bounded(buf, pos, b)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(in)   :: b
    do while (pos <= b)
      select case (buf(pos:pos))
      case (' ', char(9), char(10), char(13))
        pos = pos + 1
      case default
        return
      end select
    end do
  end subroutine skip_ws_bounded

  subroutine next_token_span(buf, pos, b, ta, tb)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(in)   :: b
    integer(ip),       intent(out)  :: ta, tb
    call skip_ws_bounded(buf, pos, b)
    ta = pos
    do while (pos <= b .and. buf(pos:pos) /= ' ' .and. buf(pos:pos) /= char(9) &
              .and. buf(pos:pos) /= char(10) .and. buf(pos:pos) /= char(13))
      pos = pos + 1
    end do
    tb = pos - 1
    if (tb < ta) tb = ta - 1  ! empty token (e.g. trailing "()")
  end subroutine next_token_span

  subroutine next_int_ascii(buf, pos, b, val)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(in)   :: b
    integer(ip),       intent(out)  :: val
    integer(ip) :: ta, tb
    call next_token_span(buf, pos, b, ta, tb)
    if (tb < ta) error stop 'icem_read: expected integer token, found none'
    read(buf(ta:tb), *) val
  end subroutine next_int_ascii

  subroutine next_real_ascii(buf, pos, b, val)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(in)   :: b
    real(rp),          intent(out)  :: val
    integer(ip) :: ta, tb
    call next_token_span(buf, pos, b, ta, tb)
    if (tb < ta) error stop 'icem_read: expected real token, found none'
    read(buf(ta:tb), *) val
  end subroutine next_real_ascii

  subroutine skip_n_tokens_ascii(buf, pos, b, n)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(ip),       intent(in)   :: b, n
    integer(ip) :: i, ta, tb
    do i = 1, n
      call next_token_span(buf, pos, b, ta, tb)
    end do
  end subroutine skip_n_tokens_ascii


  !> Little-endian 32-bit integer at `pos`; advances `pos` by 4 bytes.
  integer(ip) function read_bin_i32(buf, pos) result(val)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(int8)  :: raw(4)
    integer(int32) :: v32
    integer(ip) :: i
    do i = 1, 4
      raw(i) = int(iachar(buf(pos+i-1:pos+i-1)), int8)
    end do
    v32 = transfer(raw, v32)
    val = int(v32, ip)
    pos = pos + 4
  end function read_bin_i32

  !> Little-endian 64-bit double at `pos`; advances `pos` by 8 bytes.
  real(rp) function read_bin_f64(buf, pos) result(val)
    character(len=*), intent(in)    :: buf
    integer(ip),       intent(inout):: pos
    integer(int8)  :: raw(8)
    real(real64)   :: v64
    integer(ip) :: i
    do i = 1, 8
      raw(i) = int(iachar(buf(pos+i-1:pos+i-1)), int8)
    end do
    v64 = transfer(raw, v64)
    val = real(v64, rp)
    pos = pos + 8
  end function read_bin_f64

end module mod_icem_reader
