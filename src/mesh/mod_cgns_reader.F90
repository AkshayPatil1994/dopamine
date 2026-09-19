! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> CGNS polyhedral mesh reader (NGON_n/NFACE_n and standard element sections).
module mod_cgns_reader
#include "cgnstypes_f03.h"
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use CGNS
  implicit none
  private

  public :: cgns_read

  !> Pad value for triangle keys; large so it sorts after every real vertex id.
  integer(ip), parameter :: NO_VERT = huge(0_ip)

  !> Face-signature record: sorted vertex-set key plus owning cell (or boundary element id).

  type :: FaceSig
    integer(ip) :: key(4) = 0   ! sorted vertex ids, NO_VERT-padded for tri faces
    integer(ip) :: cell   = 0
  end type FaceSig

  !> A boundary (TRI_3/QUAD_4) element section: its name and global element-id range.
  type :: BndSection
    character(len=64) :: name = ''
    integer(ip) :: estart = 0, eend = 0
  end type BndSection

contains

  subroutine cgns_read(filename, msh)
    character(len=*), intent(in)    :: filename
    type(Mesh),        intent(inout) :: msh

    integer :: fn, ier, nbases, cell_dim, phys_dim
    integer :: nzones, iz, zonetype
    integer(cgsize_t) :: isize(1,3)
    character(len=32) :: basename, zonename
    integer(ip) :: n_points, n_cells
    integer(ip), allocatable :: zone_node_offset(:), zone_cell_offset(:)
    real(rp), allocatable :: points(:,:)
    integer(ip), allocatable :: cell_vertex_offsets(:), cell_vertex_ids(:), cell_types(:)
    integer(ip), allocatable :: face_offsets(:), face_vertices(:), face_owner(:), face_neigh(:)
    integer(ip) :: n_faces_total, n_int_faces
    integer(cgsize_t) :: ngon_estart
    integer(ip), allocatable :: orig_to_mesh_face(:)
    type(BndSection), allocatable :: bnd_secs(:)

    call cg_open_f(trim(filename), CG_MODE_READ, fn, ier)
    call cgns_check(ier, 'cg_open_f')

    call cg_nbases_f(fn, nbases, ier)
    call cgns_check(ier, 'cg_nbases_f')
    if (nbases < 1) error stop 'cgns_read: no CGNS base found in file'

    call cg_base_read_f(fn, 1, basename, cell_dim, phys_dim, ier)
    call cgns_check(ier, 'cg_base_read_f')
    if (cell_dim /= 3) error stop 'cgns_read: only 3-D CGNS bases are supported'

    call cg_nzones_f(fn, 1, nzones, ier)
    call cgns_check(ier, 'cg_nzones_f')
    if (nzones < 1) error stop 'cgns_read: no zones found in CGNS base'

    allocate(zone_node_offset(nzones + 1), zone_cell_offset(nzones + 1))
    zone_node_offset(1) = 0
    zone_cell_offset(1) = 0
    do iz = 1, nzones
      call cg_zone_type_f(fn, 1, iz, zonetype, ier)
      call cgns_check(ier, 'cg_zone_type_f')
      if (zonetype /= CGNS_ENUMV(Unstructured)) &
        error stop 'cgns_read: only Unstructured zones are supported'
      call cg_zone_read_f(fn, 1, iz, zonename, isize, ier)
      call cgns_check(ier, 'cg_zone_read_f')
      zone_node_offset(iz + 1) = zone_node_offset(iz) + int(isize(1,1), ip)
      zone_cell_offset(iz + 1) = zone_cell_offset(iz) + int(isize(1,2), ip)
    end do
    n_points = zone_node_offset(nzones + 1)
    n_cells  = zone_cell_offset(nzones + 1)
    if (n_points <= 0) error stop 'cgns_read: zero total vertices across zones'
    if (n_cells  <= 0) error stop 'cgns_read: zero total cells across zones'

    allocate(points(3, n_points))
    allocate(cell_types(n_cells), source=0_ip)

    do iz = 1, nzones
      call read_zone_coords(fn, iz, zone_node_offset(iz), points)
    end do

    call read_zone_elements(fn, nzones, zone_node_offset, zone_cell_offset, n_points, n_cells, &
                             cell_types, cell_vertex_offsets, cell_vertex_ids, &
                             face_offsets, face_vertices, face_owner, face_neigh, &
                             n_faces_total, n_int_faces, points, ngon_estart, orig_to_mesh_face, bnd_secs)

    msh%n_cells     = n_cells
    msh%n_ghost     = 0
    msh%n_faces     = n_faces_total
    msh%n_int_faces = n_int_faces
    msh%n_bnd_faces = n_faces_total - n_int_faces
    msh%n_vertices  = n_points
    allocate(msh%cells(n_cells))
    allocate(msh%global_id(n_cells))
    block
      integer(ip) :: c
      do c = 1, n_cells
        msh%cells(c)%centroid = 0.0_rp
        msh%cells(c)%volume   = 0.0_rp
        msh%global_id(c)      = int(c, lp)
      end do
    end block

    ! Determine each boundary face's patch from ZoneBC_t/BC_t, then regroup the boundary segment
    ! of the face arrays contiguously by patch (the Mesh convention every reader must satisfy) --
    ! this must happen before geometry is computed, since it reorders face_offsets/face_vertices.
    call read_boundary_conditions(fn, nzones, ngon_estart, orig_to_mesh_face, n_int_faces, n_faces_total, &
                                   bnd_secs, face_offsets, face_vertices, face_owner, face_neigh, msh)

    allocate(msh%faces(n_faces_total))
    block
      integer(ip) :: f, p, lf
      do f = 1, n_faces_total
        msh%faces(f)%owner     = face_owner(f)
        msh%faces(f)%neighbour = face_neigh(f)
      end do
      do p = 1, msh%n_patches
        do lf = msh%patches(p)%start_face, msh%patches(p)%start_face + msh%patches(p)%n_faces - 1
          msh%faces(n_int_faces + lf)%patch_id = p
        end do
      end do
    end block

    allocate(msh%vertices(3, n_points))
    msh%vertices = points
    call move_alloc(cell_vertex_offsets, msh%cell_offsets)
    call move_alloc(cell_vertex_ids, msh%cell_vertices)
    call move_alloc(cell_types, msh%cell_types)
    call move_alloc(face_offsets, msh%face_offsets)
    call move_alloc(face_vertices, msh%face_vertices)

    block
      integer(ip) :: f
      do f = 1, n_faces_total
        call compute_polygon_face_geometry(f, msh%face_vertices, msh%face_offsets, points, &
                                            msh%faces(f)%centre, msh%faces(f)%normal, msh%faces(f)%area)
        msh%faces(f)%sf = msh%faces(f)%area * msh%faces(f)%normal
      end do
    end block
    call compute_polyhedron_cell_moments(msh, n_cells, n_faces_total, msh%face_vertices, &
                                          msh%face_offsets, points)

    call cg_close_f(fn, ier)
    call cgns_check(ier, 'cg_close_f')

    call validate_cgns_boundary_connectivity(msh)
  end subroutine cgns_read

  subroutine read_zone_coords(fn, iz, node_offset, points)
    integer,     intent(in)    :: fn, iz
    integer(ip), intent(in)    :: node_offset
    real(rp),    intent(inout) :: points(:,:)

    integer(cgsize_t) :: isize(1,3), rmin(3), rmax(3)
    character(len=32) :: zonename
    integer :: ier
    integer(ip) :: n, i
    real(kind=8), allocatable :: xc(:), yc(:), zc(:)

    call cg_zone_read_f(fn, 1, iz, zonename, isize, ier)
    call cgns_check(ier, 'cg_zone_read_f (coords)')
    n = int(isize(1,1), ip)
    allocate(xc(n), yc(n), zc(n))
    rmin = 1_cgsize_t
    rmax(1) = int(n, cgsize_t)

    call cg_coord_read_f(fn, 1, iz, 'CoordinateX', CGNS_ENUMV(RealDouble), rmin, rmax, xc, ier)
    call cgns_check(ier, 'cg_coord_read_f (X)')
    call cg_coord_read_f(fn, 1, iz, 'CoordinateY', CGNS_ENUMV(RealDouble), rmin, rmax, yc, ier)
    call cgns_check(ier, 'cg_coord_read_f (Y)')
    call cg_coord_read_f(fn, 1, iz, 'CoordinateZ', CGNS_ENUMV(RealDouble), rmin, rmax, zc, ier)
    call cgns_check(ier, 'cg_coord_read_f (Z)')

    do i = 1, n
      points(1, node_offset + i) = real(xc(i), rp)
      points(2, node_offset + i) = real(yc(i), rp)
      points(3, node_offset + i) = real(zc(i), rp)
    end do
    deallocate(xc, yc, zc)
  end subroutine read_zone_coords

  !> Read every element section across all zones, splitting into the polyhedral NGON_n/NFACE_n
  !> path (used directly as face/cell CSR connectivity) and the standard-element path (tet/hex/
  !> prism/pyra, whose faces are derived by sorting vertex-set signatures and merging duplicates).
  subroutine read_zone_elements(fn, nzones, zone_node_offset, zone_cell_offset, n_points, n_cells, &
                                 cell_types, cell_vertex_offsets, cell_vertex_ids, &
                                 face_offsets, face_vertices, face_owner, face_neigh, &
                                 n_faces_total, n_int_faces, points, ngon_estart, orig_to_mesh_face, bnd_secs)
    integer,     intent(in) :: fn, nzones
    integer(ip), intent(in) :: zone_node_offset(:), zone_cell_offset(:), n_points, n_cells
    integer(ip), intent(inout) :: cell_types(:)
    integer(ip), allocatable, intent(out) :: cell_vertex_offsets(:), cell_vertex_ids(:)
    integer(ip), allocatable, intent(out) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(out) :: face_owner(:), face_neigh(:)
    integer(ip), intent(out) :: n_faces_total, n_int_faces
    real(rp),    intent(in) :: points(:,:)
    integer(cgsize_t), intent(out) :: ngon_estart
    integer(ip), allocatable, intent(out) :: orig_to_mesh_face(:)
    type(BndSection), allocatable, intent(out) :: bnd_secs(:)

    integer :: iz, ier, nsections, is, etype, nbndry, parent_flag
    integer(cgsize_t) :: estart, eend
    type(FaceSig), allocatable :: standard_faces(:), bnd_sigs(:)
    integer(ip) :: n_std_faces, n_bnd_sigs, n_bsec, emin, emax
    logical :: seen_ngon, seen_nface

    ! Polyhedral CSR connectivity accumulated directly from a single NGON_n/NFACE_n section pair
    ! (the standard ICEM/Fluent-Meshing convention: one face section, one cell section per zone).
    integer(ip), allocatable :: ngon_offsets(:), ngon_verts(:)
    integer(ip), allocatable :: nface_offsets(:), nface_faceids(:)
    logical :: have_poly

    have_poly = .false.
    seen_ngon = .false.
    seen_nface = .false.
    n_std_faces = 0
    ngon_estart = 1_cgsize_t
    allocate(standard_faces(0), bnd_sigs(0), bnd_secs(0))
    n_bnd_sigs = 0

    if (nzones > 1) &
      error stop 'cgns_read: multi-zone CGNS files are not yet supported for the polyhedral BC path'

    do iz = 1, nzones
      call cg_nsections_f(fn, 1, iz, nsections, ier)
      call cgns_check(ier, 'cg_nsections_f')
      do is = 1, nsections
        call section_element_type(fn, iz, is, etype, estart, eend, nbndry, parent_flag)
        select case (etype)
        case (CGNS_ENUMV(NGON_n))
          if (seen_ngon) error stop 'cgns_read: more than one NGON_n section is not yet supported'
          seen_ngon = .true.
          have_poly = .true.
          ngon_estart = estart
          call read_ngon_section(fn, iz, is, zone_node_offset(iz), estart, eend, ngon_offsets, ngon_verts)
        case (CGNS_ENUMV(NFACE_n))
          if (seen_nface) error stop 'cgns_read: more than one NFACE_n section is not yet supported'
          seen_nface = .true.
          have_poly = .true.
          call read_nface_section(fn, iz, is, zone_cell_offset(iz), estart, eend, &
                                   nface_offsets, nface_faceids, cell_types)
        case (CGNS_ENUMV(TETRA_4), CGNS_ENUMV(HEXA_8), CGNS_ENUMV(PENTA_6), CGNS_ENUMV(PYRA_5))
          call read_standard_section(fn, iz, is, zone_node_offset(iz), zone_cell_offset(iz), &
                                      etype, estart, eend, cell_types, standard_faces, n_std_faces)
        case (CGNS_ENUMV(TRI_3), CGNS_ENUMV(QUAD_4))
          call read_boundary_section(fn, iz, is, zone_node_offset(iz), etype, estart, eend, &
                                      bnd_sigs, n_bnd_sigs, bnd_secs)
        case default
          ! Other section types (e.g. BAR_2 edges) carry no boundary faces and are ignored.
        end select
      end do
    end do

    if (have_poly .and. n_std_faces > 0) &
      error stop 'cgns_read: mixing NGON_n/NFACE_n and standard element sections is not supported'

    if (have_poly) then
      if (.not. allocated(nface_offsets)) &
        error stop 'cgns_read: NGON_n faces present without a matching NFACE_n cell section'
      call assemble_polyhedral_mesh(n_cells, ngon_offsets, ngon_verts, nface_offsets, nface_faceids, &
                                     points, cell_vertex_offsets, cell_vertex_ids, &
                                     face_offsets, face_vertices, face_owner, face_neigh, &
                                     n_faces_total, n_int_faces, orig_to_mesh_face)
    else
      call assemble_standard_mesh(n_cells, standard_faces(1:n_std_faces), points, &
                                   cell_vertex_offsets, cell_vertex_ids, cell_types, &
                                   face_offsets, face_vertices, face_owner, face_neigh, &
                                   n_faces_total, n_int_faces)
      n_bsec = size(bnd_secs)
      if (n_bsec > 0) then
        emin = minval(bnd_secs(:)%estart)
        emax = maxval(bnd_secs(:)%eend)
        ngon_estart = int(emin, cgsize_t)
        call map_boundary_elements(bnd_sigs(1:n_bnd_sigs), face_offsets, face_vertices, n_int_faces, &
                                    n_faces_total, emin, emax, orig_to_mesh_face)
      else
        allocate(orig_to_mesh_face(0))
      end if
    end if
  end subroutine read_zone_elements

  subroutine section_element_type(fn, iz, is, etype, estart, eend, nbndry, parent_flag)
    integer,            intent(in)  :: fn, iz, is
    integer,            intent(out) :: etype, nbndry, parent_flag
    integer(cgsize_t),  intent(out) :: estart, eend
    character(len=32) :: sname
    integer :: ier
    call cg_section_read_f(fn, 1, iz, is, sname, etype, estart, eend, nbndry, parent_flag, ier)
    call cgns_check(ier, 'cg_section_read_f')
  end subroutine section_element_type

  !> NGON_n section: face-vertex CSR connectivity, directly usable as Mesh face_vertices/offsets.
  subroutine read_ngon_section(fn, iz, is, node_offset, estart, eend, ngon_offsets, ngon_verts)
    integer,           intent(in) :: fn, iz, is
    integer(ip),       intent(in) :: node_offset
    integer(cgsize_t), intent(in) :: estart, eend
    integer(ip), allocatable, intent(inout) :: ngon_offsets(:), ngon_verts(:)

    integer :: ier
    integer(cgsize_t) :: data_size
    integer(cgsize_t), allocatable :: conn(:), off(:)
    integer(cgsize_t) :: parent_data(1)
    integer(ip) :: nf, i

    nf = int(eend - estart + 1, ip)
    call cg_ElementDataSize_f(fn, 1, iz, is, data_size, ier)
    call cgns_check(ier, 'cg_ElementDataSize_f (NGON_n)')

    allocate(conn(data_size))
    allocate(off(nf + 1))
    call cg_poly_elements_read_f(fn, 1, iz, is, conn, off, parent_data, ier)
    call cgns_check(ier, 'cg_poly_elements_read_f (NGON_n)')

    allocate(ngon_offsets(nf + 1))
    ngon_offsets(1) = 1
    do i = 1, nf
      ngon_offsets(i + 1) = ngon_offsets(i) + int(off(i + 1) - off(i), ip)
    end do
    allocate(ngon_verts(ngon_offsets(nf + 1) - 1))
    do i = 1, size(conn)
      ngon_verts(i) = int(conn(i), ip) + node_offset
    end do

    deallocate(conn, off)
  end subroutine read_ngon_section

  !> NFACE_n section: cell-face CSR connectivity (signed face ids: sign gives orientation).
  subroutine read_nface_section(fn, iz, is, cell_offset, estart, eend, &
                                 nface_offsets, nface_faceids, cell_types)
    integer,           intent(in) :: fn, iz, is
    integer(ip),       intent(in) :: cell_offset
    integer(cgsize_t), intent(in) :: estart, eend
    integer(ip), allocatable, intent(inout) :: nface_offsets(:), nface_faceids(:)
    integer(ip), intent(inout) :: cell_types(:)

    integer :: ier
    integer(cgsize_t) :: data_size
    integer(cgsize_t), allocatable :: conn(:), off(:)
    integer(cgsize_t) :: parent_data(1)
    integer(ip) :: nc, i, c

    nc = int(eend - estart + 1, ip)
    call cg_ElementDataSize_f(fn, 1, iz, is, data_size, ier)
    call cgns_check(ier, 'cg_ElementDataSize_f (NFACE_n)')

    allocate(conn(data_size))
    allocate(off(nc + 1))
    call cg_poly_elements_read_f(fn, 1, iz, is, conn, off, parent_data, ier)
    call cgns_check(ier, 'cg_poly_elements_read_f (NFACE_n)')

    allocate(nface_offsets(nc + 1))
    nface_offsets(1) = 1
    do i = 1, nc
      nface_offsets(i + 1) = nface_offsets(i) + int(off(i + 1) - off(i), ip)
    end do
    allocate(nface_faceids(nface_offsets(nc + 1) - 1))
    do i = 1, size(conn)
      nface_faceids(i) = int(conn(i), ip)
    end do

    do c = 1, nc
      cell_types(cell_offset + c) = 42_ip  ! polyhedron
    end do

    deallocate(conn, off)
  end subroutine read_nface_section

  !> Assign owner/neighbour per NGON_n face from the NFACE_n cell->face incidence (signed face id:
  !> positive = this cell owns the face with its stored orientation, negative = neighbour side).
  subroutine assemble_polyhedral_mesh(n_cells, ngon_offsets, ngon_verts, nface_offsets, nface_faceids, &
                                       points, cell_vertex_offsets, cell_vertex_ids, &
                                       face_offsets, face_vertices, face_owner, face_neigh, &
                                       n_faces_total, n_int_faces, orig_to_mesh_face)
    integer(ip), intent(in) :: n_cells
    integer(ip), intent(in) :: ngon_offsets(:), ngon_verts(:)
    integer(ip), intent(in) :: nface_offsets(:), nface_faceids(:)
    real(rp),    intent(in) :: points(:,:)
    integer(ip), allocatable, intent(out) :: cell_vertex_offsets(:), cell_vertex_ids(:)
    integer(ip), allocatable, intent(out) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(out) :: face_owner(:), face_neigh(:)
    integer(ip), intent(out) :: n_faces_total, n_int_faces
    integer(ip), allocatable, intent(out) :: orig_to_mesh_face(:)

    integer(ip) :: n_faces, c, i, fid, sgn, afid
    integer(ip), allocatable :: owner0(:), neigh0(:)

    n_faces = size(ngon_offsets) - 1
    allocate(owner0(n_faces), neigh0(n_faces), source=-1_ip)

    do c = 1, n_cells
      do i = nface_offsets(c), nface_offsets(c + 1) - 1
        fid = nface_faceids(i)
        sgn = sign(1_ip, fid)
        afid = abs(fid)
        if (sgn > 0) then
          if (owner0(afid) /= -1) error stop 'cgns_read: NFACE_n face referenced as owner twice'
          owner0(afid) = c
        else
          if (neigh0(afid) /= -1) error stop 'cgns_read: NFACE_n face referenced as neighbour twice'
          neigh0(afid) = c
        end if
      end do
    end do
    if (any(owner0 == -1)) error stop 'cgns_read: NGON_n face with no owning cell in NFACE_n'

    ! Reorder: internal faces first, then boundary faces (Mesh convention), preserving CSR shape.
    call reorder_faces_internal_first(n_faces, owner0, neigh0, ngon_offsets, ngon_verts, &
                                       face_offsets, face_vertices, face_owner, face_neigh, &
                                       n_faces_total, n_int_faces, orig_to_mesh_face)

    call build_cell_vertices_from_faces(n_cells, face_offsets, face_vertices, face_owner, face_neigh, &
                                         cell_vertex_offsets, cell_vertex_ids)
  end subroutine assemble_polyhedral_mesh

  !> Reorder a face set (arbitrary order in) into internal-first, boundary-after CSR order.
  !> `orig_to_mesh_face(f)` gives the new Mesh face index for original (NGON_n-local) face `f`,
  !> so boundary-condition point lists (which reference the original CGNS element numbering) can
  !> be resolved to Mesh face indices after this reordering.
  subroutine reorder_faces_internal_first(n_faces, owner0, neigh0, offsets0, verts0, &
                                           face_offsets, face_vertices, face_owner, face_neigh, &
                                           n_faces_total, n_int_faces, orig_to_mesh_face)
    integer(ip), intent(in) :: n_faces
    integer(ip), intent(in) :: owner0(:), neigh0(:), offsets0(:), verts0(:)
    integer(ip), allocatable, intent(out) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(out) :: face_owner(:), face_neigh(:)
    integer(ip), intent(out) :: n_faces_total, n_int_faces
    integer(ip), allocatable, intent(out) :: orig_to_mesh_face(:)

    integer(ip) :: f, wi, wpos, nn

    n_int_faces = count(neigh0 /= -1)
    n_faces_total = n_faces
    allocate(face_owner(n_faces_total), face_neigh(n_faces_total))
    allocate(face_offsets(n_faces_total + 1))
    allocate(orig_to_mesh_face(n_faces))

    face_offsets(1) = 1
    wi = 0
    do f = 1, n_faces
      if (neigh0(f) /= -1) then
        wi = wi + 1
        face_owner(wi) = owner0(f)
        face_neigh(wi) = neigh0(f)
        face_offsets(wi + 1) = face_offsets(wi) + (offsets0(f + 1) - offsets0(f))
        orig_to_mesh_face(f) = wi
      end if
    end do
    do f = 1, n_faces
      if (neigh0(f) == -1) then
        wi = wi + 1
        face_owner(wi) = owner0(f)
        face_neigh(wi) = -1
        face_offsets(wi + 1) = face_offsets(wi) + (offsets0(f + 1) - offsets0(f))
        orig_to_mesh_face(f) = wi
      end if
    end do

    allocate(face_vertices(face_offsets(n_faces_total + 1) - 1))
    wi = 0
    do f = 1, n_faces
      if (neigh0(f) /= -1) then
        wi = wi + 1
        wpos = face_offsets(wi)
        nn = offsets0(f + 1) - offsets0(f)
        face_vertices(wpos:wpos + nn - 1) = verts0(offsets0(f):offsets0(f) + nn - 1)
      end if
    end do
    do f = 1, n_faces
      if (neigh0(f) == -1) then
        wi = wi + 1
        wpos = face_offsets(wi)
        nn = offsets0(f + 1) - offsets0(f)
        face_vertices(wpos:wpos + nn - 1) = verts0(offsets0(f):offsets0(f) + nn - 1)
      end if
    end do
  end subroutine reorder_faces_internal_first

  !> Build cell_offsets/cell_vertices (dedup'd vertex list per cell) from face incidence, mirroring
  !> mod_foam_reader.f90's build_cell_topology_metadata.
  subroutine build_cell_vertices_from_faces(n_cells, face_offsets, face_vertices, face_owner, face_neigh, &
                                             cell_vertex_offsets, cell_vertex_ids)
    integer(ip), intent(in) :: n_cells
    integer(ip), intent(in) :: face_offsets(:), face_vertices(:), face_owner(:), face_neigh(:)
    integer(ip), allocatable, intent(out) :: cell_vertex_offsets(:), cell_vertex_ids(:)

    integer(ip) :: n_faces, f, c, pos, nverts, write_pos
    integer(ip), allocatable :: cell_face_count(:), cell_face_offset(:), cell_face_ids(:), fill_count(:)
    integer(ip), allocatable :: cell_vertex_count(:), vertex_marker(:)
    integer(ip) :: n_points

    n_faces = size(face_owner)
    n_points = maxval(face_vertices)

    allocate(cell_face_count(n_cells), source=0_ip)
    do f = 1, n_faces
      cell_face_count(face_owner(f)) = cell_face_count(face_owner(f)) + 1
      if (face_neigh(f) >= 1) cell_face_count(face_neigh(f)) = cell_face_count(face_neigh(f)) + 1
    end do
    allocate(cell_face_offset(n_cells + 1))
    cell_face_offset(1) = 1
    do c = 1, n_cells
      cell_face_offset(c + 1) = cell_face_offset(c) + cell_face_count(c)
    end do
    allocate(cell_face_ids(cell_face_offset(n_cells + 1) - 1), fill_count(n_cells), source=0_ip)
    do f = 1, n_faces
      c = face_owner(f)
      cell_face_ids(cell_face_offset(c) + fill_count(c)) = f
      fill_count(c) = fill_count(c) + 1
      if (face_neigh(f) >= 1) then
        c = face_neigh(f)
        cell_face_ids(cell_face_offset(c) + fill_count(c)) = f
        fill_count(c) = fill_count(c) + 1
      end if
    end do

    allocate(cell_vertex_count(n_cells), vertex_marker(n_points), source=0_ip)
    do c = 1, n_cells
      nverts = 0
      do pos = cell_face_offset(c), cell_face_offset(c + 1) - 1
        f = cell_face_ids(pos)
        block
          integer(ip) :: k
          do k = face_offsets(f), face_offsets(f + 1) - 1
            if (vertex_marker(face_vertices(k)) /= c) then
              vertex_marker(face_vertices(k)) = c
              nverts = nverts + 1
            end if
          end do
        end block
      end do
      cell_vertex_count(c) = nverts
    end do

    allocate(cell_vertex_offsets(n_cells + 1))
    cell_vertex_offsets(1) = 1
    do c = 1, n_cells
      cell_vertex_offsets(c + 1) = cell_vertex_offsets(c) + cell_vertex_count(c)
    end do
    allocate(cell_vertex_ids(cell_vertex_offsets(n_cells + 1) - 1))
    vertex_marker = 0
    do c = 1, n_cells
      write_pos = cell_vertex_offsets(c)
      nverts = 0
      do pos = cell_face_offset(c), cell_face_offset(c + 1) - 1
        f = cell_face_ids(pos)
        block
          integer(ip) :: k
          do k = face_offsets(f), face_offsets(f + 1) - 1
            if (vertex_marker(face_vertices(k)) /= c) then
              vertex_marker(face_vertices(k)) = c
              cell_vertex_ids(write_pos + nverts) = face_vertices(k)
              nverts = nverts + 1
            end if
          end do
        end block
      end do
    end do
  end subroutine build_cell_vertices_from_faces

  !> Standard (tet/hex/prism/pyra) element section: read flat connectivity, append each cell's
  !> face vertex-set signatures (see FaceSig) for later sort-merge into internal/boundary faces.
  subroutine read_standard_section(fn, iz, is, node_offset, cell_offset, etype, estart, eend, &
                                    cell_types, standard_faces, n_std_faces)
    integer,           intent(in) :: fn, iz, is, etype
    integer(ip),       intent(in) :: node_offset, cell_offset
    integer(cgsize_t), intent(in) :: estart, eend
    integer(ip), intent(inout) :: cell_types(:)
    type(FaceSig), allocatable, intent(inout) :: standard_faces(:)
    integer(ip), intent(inout) :: n_std_faces

    integer :: ier, npe
    integer(cgsize_t), allocatable :: elems(:)
    integer(cgsize_t) :: parent_data(1)
    integer(ip) :: ne, e, c, gv(8), vtk_type, nfaces_this
    type(FaceSig), allocatable :: tmp(:)

    call cg_npe_f(etype, npe, ier)
    call cgns_check(ier, 'cg_npe_f')
    ne = int(eend - estart + 1, ip)
    allocate(elems(int(ne, cgsize_t) * npe))
    call cg_elements_read_f(fn, 1, iz, is, elems, parent_data, ier)
    call cgns_check(ier, 'cg_elements_read_f')

    select case (etype)
    case (CGNS_ENUMV(TETRA_4)); vtk_type = 10; nfaces_this = 4
    case (CGNS_ENUMV(HEXA_8));  vtk_type = 12; nfaces_this = 6
    case (CGNS_ENUMV(PENTA_6)); vtk_type = 13; nfaces_this = 5
    case (CGNS_ENUMV(PYRA_5));  vtk_type = 14; nfaces_this = 5
    case default
      error stop 'cgns_read: unsupported standard element type'
    end select

    allocate(tmp(n_std_faces + ne * nfaces_this))
    tmp(1:n_std_faces) = standard_faces
    call move_alloc(tmp, standard_faces)

    do e = 1, ne
      c = cell_offset + int(estart, ip) - 1 + (e - 1) + 1  ! 1-based global cell index
      gv(1:npe) = int(elems((e-1)*npe+1 : e*npe), ip) + node_offset
      cell_types(c) = vtk_type
      call append_cell_faces(etype, gv, c, standard_faces, n_std_faces)
    end do
    deallocate(elems)
  end subroutine read_standard_section

  subroutine append_cell_faces(etype, gv, c, standard_faces, n_std_faces)
    integer,       intent(in) :: etype
    integer(ip),   intent(in) :: gv(8), c
    type(FaceSig), intent(inout) :: standard_faces(:)
    integer(ip),   intent(inout) :: n_std_faces

    select case (etype)
    case (CGNS_ENUMV(TETRA_4))
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(3),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(4),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(1),gv(3),gv(4),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(2),gv(3),gv(4),NO_VERT], c)
    case (CGNS_ENUMV(HEXA_8))
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(3),gv(4)], c)
      call push_face(standard_faces, n_std_faces, [gv(5),gv(6),gv(7),gv(8)], c)
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(6),gv(5)], c)
      call push_face(standard_faces, n_std_faces, [gv(2),gv(3),gv(7),gv(6)], c)
      call push_face(standard_faces, n_std_faces, [gv(3),gv(4),gv(8),gv(7)], c)
      call push_face(standard_faces, n_std_faces, [gv(4),gv(1),gv(5),gv(8)], c)
    case (CGNS_ENUMV(PENTA_6))
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(3),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(4),gv(5),gv(6),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(5),gv(4)], c)
      call push_face(standard_faces, n_std_faces, [gv(2),gv(3),gv(6),gv(5)], c)
      call push_face(standard_faces, n_std_faces, [gv(3),gv(1),gv(4),gv(6)], c)
    case (CGNS_ENUMV(PYRA_5))
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(3),gv(4)], c)
      call push_face(standard_faces, n_std_faces, [gv(1),gv(2),gv(5),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(2),gv(3),gv(5),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(3),gv(4),gv(5),NO_VERT], c)
      call push_face(standard_faces, n_std_faces, [gv(4),gv(1),gv(5),NO_VERT], c)
    end select
  end subroutine append_cell_faces

  subroutine push_face(standard_faces, n_std_faces, verts, c)
    type(FaceSig), intent(inout) :: standard_faces(:)
    integer(ip),   intent(inout) :: n_std_faces
    integer(ip),   intent(in)    :: verts(4), c
    integer(ip) :: k(4)
    k = verts
    call sort4(k)
    n_std_faces = n_std_faces + 1
    standard_faces(n_std_faces)%key  = k
    standard_faces(n_std_faces)%cell = c
  end subroutine push_face

  !> Read a TRI_3/QUAD_4 boundary section into keyed signatures (cell field = global element id).
  subroutine read_boundary_section(fn, iz, is, node_offset, etype, estart, eend, bnd_sigs, n_bnd_sigs, bnd_secs)
    integer,           intent(in) :: fn, iz, is, etype
    integer(ip),       intent(in) :: node_offset
    integer(cgsize_t), intent(in) :: estart, eend
    type(FaceSig), allocatable, intent(inout) :: bnd_sigs(:)
    integer(ip),       intent(inout) :: n_bnd_sigs
    type(BndSection), allocatable, intent(inout) :: bnd_secs(:)

    integer :: ier, npe, nbndry, pf, et
    integer(cgsize_t) :: es, ee, parent_data(1)
    integer(cgsize_t), allocatable :: elems(:)
    integer(ip) :: ne, e, k(4)
    type(FaceSig), allocatable :: tmp(:)
    type(BndSection), allocatable :: tmps(:)
    character(len=64) :: nm

    call cg_section_read_f(fn, 1, iz, is, nm, et, es, ee, nbndry, pf, ier)
    call cgns_check(ier, 'cg_section_read_f (boundary)')
    call cg_npe_f(etype, npe, ier)
    call cgns_check(ier, 'cg_npe_f')
    ne = int(eend - estart + 1, ip)
    allocate(elems(int(ne, cgsize_t) * npe))
    call cg_elements_read_f(fn, 1, iz, is, elems, parent_data, ier)
    call cgns_check(ier, 'cg_elements_read_f (boundary)')

    allocate(tmp(n_bnd_sigs + ne))
    tmp(1:n_bnd_sigs) = bnd_sigs(1:n_bnd_sigs)
    call move_alloc(tmp, bnd_sigs)
    do e = 1, ne
      k = NO_VERT
      k(1:npe) = int(elems((e-1)*npe+1 : e*npe), ip) + node_offset
      call sort4(k)
      n_bnd_sigs = n_bnd_sigs + 1
      bnd_sigs(n_bnd_sigs)%key  = k
      bnd_sigs(n_bnd_sigs)%cell = int(estart, ip) + e - 1
    end do

    allocate(tmps(size(bnd_secs) + 1))
    tmps(1:size(bnd_secs)) = bnd_secs
    tmps(size(bnd_secs) + 1) = BndSection(trim(nm), int(estart, ip), int(eend, ip))
    call move_alloc(tmps, bnd_secs)
  end subroutine read_boundary_section

  !> Map each boundary element id to the Mesh boundary face with the same vertex set (0 if none).
  subroutine map_boundary_elements(bsig, face_offsets, face_vertices, n_int_faces, n_faces_total, &
                                    emin, emax, orig_to_mesh_face)
    type(FaceSig), intent(inout) :: bsig(:)
    integer(ip),   intent(in)    :: face_offsets(:), face_vertices(:), n_int_faces, n_faces_total, emin, emax
    integer(ip), allocatable, intent(out) :: orig_to_mesh_face(:)

    integer(ip) :: f, nn, lo, hi, mid, k(4)
    type(FaceSig) :: probe

    call sort_facesigs(bsig)
    allocate(orig_to_mesh_face(emax - emin + 1), source=0_ip)
    do f = n_int_faces + 1, n_faces_total
      nn = face_offsets(f+1) - face_offsets(f)
      k = NO_VERT
      k(1:nn) = face_vertices(face_offsets(f):face_offsets(f+1)-1)
      call sort4(k)
      probe%key = k
      lo = 1; hi = size(bsig)
      do while (lo <= hi)
        mid = (lo + hi) / 2
        if (all(bsig(mid)%key == k)) then
          orig_to_mesh_face(bsig(mid)%cell - emin + 1) = f
          exit
        else if (facesig_lt(bsig(mid), probe)) then
          lo = mid + 1
        else
          hi = mid - 1
        end if
      end do
    end do
  end subroutine map_boundary_elements

  subroutine sort4(k)
    integer(ip), intent(inout) :: k(4)
    integer(ip) :: i, j, t
    do i = 1, 3
      do j = i + 1, 4
        if (k(j) < k(i)) then
          t = k(i); k(i) = k(j); k(j) = t
        end if
      end do
    end do
  end subroutine sort4

  logical function facesig_lt(a, b) result(r)
    type(FaceSig), intent(in) :: a, b
    integer(ip) :: i
    r = .false.
    do i = 1, 4
      if (a%key(i) /= b%key(i)) then
        r = a%key(i) < b%key(i)
        return
      end if
    end do
  end function facesig_lt

  !> Sort-merge standard-element face signatures into internal (2 cells sharing a key) and
  !> boundary (1 cell) faces, then build CSR face_vertices/offsets and orient each face so its
  !> normal points from owner toward neighbour (internal) or outward from owner (boundary).
  subroutine assemble_standard_mesh(n_cells, sig, points, &
                                     cell_vertex_offsets, cell_vertex_ids, cell_types, &
                                     face_offsets, face_vertices, face_owner, face_neigh, &
                                     n_faces_total, n_int_faces)
    integer(ip),   intent(in) :: n_cells
    type(FaceSig), intent(inout) :: sig(:)
    real(rp),      intent(in) :: points(:,:)
    integer(ip), allocatable, intent(out) :: cell_vertex_offsets(:), cell_vertex_ids(:)
    integer(ip), intent(in) :: cell_types(:)
    integer(ip), allocatable, intent(out) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(out) :: face_owner(:), face_neigh(:)
    integer(ip), intent(out) :: n_faces_total, n_int_faces

    integer(ip) :: n, i, j, ngrp, nv
    integer(ip), allocatable :: grp_owner(:), grp_neigh(:), grp_nv(:), grp_verts(:,:)

    n = size(sig)
    call sort_facesigs(sig)

    allocate(grp_owner(n), grp_neigh(n), grp_nv(n), grp_verts(4, n))
    ngrp = 0
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (any(sig(j+1)%key /= sig(i)%key)) exit
        j = j + 1
      end do
      if (j - i + 1 > 2) error stop 'cgns_read: a face is shared by more than two cells'
      ngrp = ngrp + 1
      nv = merge(3, 4, sig(i)%key(4) == NO_VERT)
      grp_nv(ngrp) = nv
      grp_verts(1:nv, ngrp) = sig(i)%key(1:nv)
      grp_owner(ngrp) = sig(i)%cell
      if (j > i) then
        grp_neigh(ngrp) = sig(j)%cell
      else
        grp_neigh(ngrp) = -1
      end if
      i = j + 1
    end do

    n_int_faces = count(grp_neigh(1:ngrp) /= -1)
    n_faces_total = ngrp
    allocate(face_owner(ngrp), face_neigh(ngrp), face_offsets(ngrp + 1))
    face_offsets(1) = 1

    block
      integer(ip) :: wi, g
      wi = 0
      do g = 1, ngrp
        if (grp_neigh(g) /= -1) then
          wi = wi + 1
          face_owner(wi) = grp_owner(g); face_neigh(wi) = grp_neigh(g)
          face_offsets(wi + 1) = face_offsets(wi) + grp_nv(g)
        end if
      end do
      do g = 1, ngrp
        if (grp_neigh(g) == -1) then
          wi = wi + 1
          face_owner(wi) = grp_owner(g); face_neigh(wi) = -1
          face_offsets(wi + 1) = face_offsets(wi) + grp_nv(g)
        end if
      end do

      allocate(face_vertices(face_offsets(ngrp + 1) - 1))
      wi = 0
      do g = 1, ngrp
        if (grp_neigh(g) /= -1) then
          wi = wi + 1
          face_vertices(face_offsets(wi):face_offsets(wi+1)-1) = grp_verts(1:grp_nv(g), g)
        end if
      end do
      do g = 1, ngrp
        if (grp_neigh(g) == -1) then
          wi = wi + 1
          face_vertices(face_offsets(wi):face_offsets(wi+1)-1) = grp_verts(1:grp_nv(g), g)
        end if
      end do
    end block

    call orient_faces_by_cell_centroid(n_cells, face_offsets, face_vertices, face_owner, face_neigh, points)
    call build_cell_vertices_from_faces(n_cells, face_offsets, face_vertices, face_owner, face_neigh, &
                                         cell_vertex_offsets, cell_vertex_ids)
  end subroutine assemble_standard_mesh

  !> Bottom-up merge sort of FaceSig by key, O(n log n).
  subroutine sort_facesigs(sig)
    type(FaceSig), intent(inout) :: sig(:)
    type(FaceSig), allocatable :: buf(:)
    integer(ip) :: n, width, lo, mid, hi, i, j, k
    logical :: take_left
    n = size(sig)
    if (n < 2) return
    allocate(buf(n))
    width = 1
    do while (width < n)
      lo = 1
      do while (lo <= n)
        mid = min(lo + width - 1, n)
        hi = min(lo + 2*width - 1, n)
        i = lo; j = mid + 1
        do k = lo, hi
          if (i > mid) then
            take_left = .false.
          else if (j > hi) then
            take_left = .true.
          else
            take_left = .not. facesig_lt(sig(j), sig(i))
          end if
          if (take_left) then
            buf(k) = sig(i); i = i + 1
          else
            buf(k) = sig(j); j = j + 1
          end if
        end do
        lo = lo + 2*width
      end do
      sig = buf
      width = 2*width
    end do
  end subroutine sort_facesigs

  !> Self-correcting orientation: doesn't rely on a memorized CGNS/VTK face winding convention.
  !> Computes each face's geometric centroid/normal from an approximate cell centroid (vertex
  !> average) and flips vertex order (boundary) or swaps owner/neighbour (internal) as needed so
  !> internal normals point owner->neighbour and boundary normals point outward from the owner.
  subroutine orient_faces_by_cell_centroid(n_cells, face_offsets, face_vertices, face_owner, face_neigh, points)
    integer(ip), intent(in)    :: n_cells
    integer(ip), intent(in)    :: face_offsets(:)
    integer(ip), intent(inout) :: face_vertices(:), face_owner(:), face_neigh(:)
    real(rp),    intent(in)    :: points(:,:)

    real(rp), allocatable :: approx_centroid(:,:)
    integer(ip), allocatable :: vcount(:)
    integer(ip) :: c, f, tmp
    real(rp) :: centre(3), normal(3), area, d

    allocate(approx_centroid(3, n_cells), source=0.0_rp)
    allocate(vcount(n_cells), source=0_ip)
    ! Approximate cell centroid via mean of its faces' vertex-average points (good enough to
    ! disambiguate sign; the real centroid/volume are computed later by compute_polyhedron_cell_moments).
    do f = 1, size(face_owner)
      block
        real(rp) :: fv(3)
        integer(ip) :: nn, i
        fv = 0.0_rp
        nn = face_offsets(f+1) - face_offsets(f)
        do i = face_offsets(f), face_offsets(f+1)-1
          fv = fv + points(:, face_vertices(i))
        end do
        fv = fv / real(nn, rp)
        approx_centroid(:, face_owner(f)) = approx_centroid(:, face_owner(f)) + fv
        vcount(face_owner(f)) = vcount(face_owner(f)) + 1
        if (face_neigh(f) >= 1) then
          approx_centroid(:, face_neigh(f)) = approx_centroid(:, face_neigh(f)) + fv
          vcount(face_neigh(f)) = vcount(face_neigh(f)) + 1
        end if
      end block
    end do
    do c = 1, n_cells
      if (vcount(c) > 0) approx_centroid(:, c) = approx_centroid(:, c) / real(vcount(c), rp)
    end do

    do f = 1, size(face_owner)
      call compute_polygon_face_geometry(f, face_vertices, face_offsets, points, centre, normal, area)
      if (face_neigh(f) >= 1) then
        d = dot_product(approx_centroid(:, face_neigh(f)) - approx_centroid(:, face_owner(f)), normal)
        if (d < 0.0_rp) then
          tmp = face_owner(f); face_owner(f) = face_neigh(f); face_neigh(f) = tmp
        end if
      else
        d = dot_product(centre - approx_centroid(:, face_owner(f)), normal)
        if (d < 0.0_rp) then
          associate(a => face_offsets(f), b => face_offsets(f+1)-1)
            face_vertices(a:b) = face_vertices(b:a:-1)
          end associate
        end if
      end if
    end do
  end subroutine orient_faces_by_cell_centroid

  !> Boundary patches from ZoneBC_t/BC_t point lists of NGON_n face ids (FaceCenter GridLocation).
  !> Resolves each BC's face ids to Mesh face indices via `orig_to_mesh_face`, then regroups the
  !> boundary segment of the face CSR arrays contiguously by patch (the layout every Mesh reader
  !> must produce) before geometry is computed by the caller. `patch_id_of_face` gives the final
  !> patch id per (regrouped) face index, for the caller to copy into msh%faces(:)%patch_id once
  !> msh%faces has been allocated.
  subroutine read_boundary_conditions(fn, nzones, ngon_estart, orig_to_mesh_face, n_int_faces, n_faces_total, &
                                       bnd_secs, face_offsets, face_vertices, face_owner, face_neigh, msh)
    integer,           intent(in)    :: fn, nzones
    integer(cgsize_t), intent(in)    :: ngon_estart
    integer(ip),       intent(in)    :: orig_to_mesh_face(:)
    integer(ip),       intent(in)    :: n_int_faces, n_faces_total
    type(BndSection),  intent(in)    :: bnd_secs(:)
    integer(ip), allocatable, intent(inout) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(inout) :: face_owner(:), face_neigh(:)
    type(Mesh),        intent(inout) :: msh

    integer :: iz, ier, nbocos, ibc, bocotype, ptset_type, normalindex(3), gridloc
    integer(cgsize_t) :: npnts, normallistsize
    integer :: normaldatatype, ndataset
    character(len=32) :: boconame
    integer(cgsize_t), allocatable :: pnts(:)
    integer(ip) :: total_bocos, patch_id, n_bnd
    integer(ip), allocatable :: patch_of_bnd(:)  ! patch id per CURRENT (pre-regroup) boundary-local face

    n_bnd = n_faces_total - n_int_faces
    allocate(patch_of_bnd(n_bnd), source=0_ip)

    if (size(bnd_secs) > 0) then
      call read_bc_from_sections(fn, nzones, ngon_estart, orig_to_mesh_face, n_int_faces, n_bnd, bnd_secs, &
                                  face_offsets, face_vertices, face_owner, face_neigh, msh)
      return
    end if

    total_bocos = 0
    do iz = 1, nzones
      call cg_nbocos_f(fn, 1, iz, nbocos, ier)
      call cgns_check(ier, 'cg_nbocos_f')
      total_bocos = total_bocos + nbocos
    end do
    if (total_bocos == 0) &
      error stop 'cgns_read: no ZoneBC_t/BC_t boundary conditions found (Family_t-based BCs are not yet supported)'

    msh%n_patches = total_bocos
    allocate(msh%patches(total_bocos))
    patch_id = 0

    do iz = 1, nzones
      call cg_nbocos_f(fn, 1, iz, nbocos, ier)
      call cgns_check(ier, 'cg_nbocos_f')
      do ibc = 1, nbocos
        call cg_boco_info_f(fn, 1, iz, ibc, boconame, bocotype, ptset_type, npnts, &
                             normalindex, normallistsize, normaldatatype, ndataset, ier)
        call cgns_check(ier, 'cg_boco_info_f')

        call cg_boco_gridlocation_read_f(fn, 1, iz, ibc, gridloc, ier)
        call cgns_check(ier, 'cg_boco_gridlocation_read_f')
        if (gridloc /= CGNS_ENUMV(FaceCenter)) &
          error stop 'cgns_read: only FaceCenter-located BCs are supported (got a different GridLocation)'

        patch_id = patch_id + 1
        allocate(pnts(npnts))
        call cg_boco_read_f(fn, 1, iz, ibc, pnts, CG_NULL, ier)
        call cgns_check(ier, 'cg_boco_read_f')

        msh%patches(patch_id)%name    = trim(boconame)
        msh%patches(patch_id)%bc_id   = cgns_bctype_to_id(bocotype)
        msh%patches(patch_id)%bc_type = cgns_bctype_name(bocotype)

        call mark_patch_faces(patch_id, ptset_type, pnts, ngon_estart, orig_to_mesh_face, n_int_faces, patch_of_bnd)
        deallocate(pnts)
      end do
    end do

    if (any(patch_of_bnd == 0)) &
      error stop 'cgns_read: some boundary faces are not covered by any ZoneBC_t/BC_t patch'

    call regroup_boundary_by_patch(msh, n_int_faces, n_bnd, patch_of_bnd, &
                                    face_offsets, face_vertices, face_owner, face_neigh)
  end subroutine read_boundary_conditions

  !> Standard-element path: one patch per TRI_3/QUAD_4 section, BC type from the matching BC_t (by vertex containment or element range).
  subroutine read_bc_from_sections(fn, nzones, emin, elem_to_face, n_int_faces, n_bnd, bnd_secs, &
                                    face_offsets, face_vertices, face_owner, face_neigh, msh)
    integer,           intent(in)    :: fn, nzones
    integer(cgsize_t), intent(in)    :: emin
    integer(ip),       intent(in)    :: elem_to_face(:), n_int_faces, n_bnd
    type(BndSection),  intent(in)    :: bnd_secs(:)
    integer(ip), allocatable, intent(inout) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(inout) :: face_owner(:), face_neigh(:)
    type(Mesh),        intent(inout) :: msh

    integer :: iz, ier, nbocos, ibc, bocotype, ptset_type, normalindex(3), gridloc, normaldatatype, ndataset
    integer(cgsize_t) :: npnts, normallistsize
    character(len=32) :: boconame
    integer(cgsize_t), allocatable :: pnts(:)
    integer(ip), allocatable :: patch_of_bnd(:)
    logical, allocatable :: in_sec(:,:)
    integer(ip) :: ns, s, e, f, i, k, sec, nverts, n_hit, n_pts
    integer(ip) :: lo, hi, best_size

    ns = size(bnd_secs)
    allocate(patch_of_bnd(n_bnd), source=0_ip)
    allocate(in_sec(msh%n_vertices, ns), source=.false.)

    do s = 1, ns
      do e = bnd_secs(s)%estart, bnd_secs(s)%eend
        f = elem_to_face(e - int(emin, ip) + 1)
        if (f < n_int_faces + 1) cycle
        if (patch_of_bnd(f - n_int_faces) /= 0) cycle
        patch_of_bnd(f - n_int_faces) = s
        do i = face_offsets(f), face_offsets(f+1) - 1
          in_sec(face_vertices(i), s) = .true.
        end do
      end do
    end do
    if (any(patch_of_bnd == 0)) &
      error stop 'cgns_read: some boundary faces are not covered by any TRI_3/QUAD_4 section'

    msh%n_patches = ns
    allocate(msh%patches(ns))
    do s = 1, ns
      msh%patches(s)%name    = trim(bnd_secs(s)%name)
      msh%patches(s)%bc_id   = BC_WALL_MODEL
      msh%patches(s)%bc_type = 'wall'
    end do

    do iz = 1, nzones
      call cg_nbocos_f(fn, 1, iz, nbocos, ier)
      call cgns_check(ier, 'cg_nbocos_f')
      do ibc = 1, nbocos
        call cg_boco_info_f(fn, 1, iz, ibc, boconame, bocotype, ptset_type, npnts, &
                             normalindex, normallistsize, normaldatatype, ndataset, ier)
        call cgns_check(ier, 'cg_boco_info_f')
        call cg_boco_gridlocation_read_f(fn, 1, iz, ibc, gridloc, ier)
        call cgns_check(ier, 'cg_boco_gridlocation_read_f')
        allocate(pnts(npnts))
        call cg_boco_read_f(fn, 1, iz, ibc, pnts, CG_NULL, ier)
        call cgns_check(ier, 'cg_boco_read_f')

        sec = 0
        n_pts = int(npnts, ip)
        if (ptset_type == CGNS_ENUMV(PointRange)) then
          lo = int(pnts(1), ip); hi = int(pnts(2), ip)
        else
          lo = int(minval(pnts), ip); hi = int(maxval(pnts), ip)
        end if
        if (gridloc == CGNS_ENUMV(FaceCenter)) then
          do s = 1, ns
            if (lo >= bnd_secs(s)%estart .and. hi <= bnd_secs(s)%eend) sec = s
          end do
        else if (gridloc == CGNS_ENUMV(Vertex) .and. ptset_type == CGNS_ENUMV(PointList)) then
          best_size = huge(0_ip)
          do s = 1, ns
            n_hit = 0
            do k = 1, n_pts
              if (in_sec(int(pnts(k), ip), s)) n_hit = n_hit + 1
            end do
            nverts = count(in_sec(:, s))
            if (n_hit == n_pts .and. nverts < best_size) then
              best_size = nverts
              sec = s
            end if
          end do
        end if
        if (sec > 0) then
          msh%patches(sec)%bc_id   = cgns_bctype_to_id(bocotype)
          msh%patches(sec)%bc_type = cgns_bctype_name(bocotype)
        end if
        deallocate(pnts)
      end do
    end do

    call regroup_boundary_by_patch(msh, n_int_faces, n_bnd, patch_of_bnd, &
                                    face_offsets, face_vertices, face_owner, face_neigh)
  end subroutine read_bc_from_sections

  !> Mark boundary faces in a BC's point set (PointRange: 2 endpoints; PointList: npnts explicit
  !> ids), resolving each CGNS NGON_n element id (1-based within the section, offset by
  !> `ngon_estart`) to a Mesh face index via `orig_to_mesh_face`, then to a boundary-local index.
  subroutine mark_patch_faces(patch_id, ptset_type, pnts, ngon_estart, orig_to_mesh_face, n_int_faces, patch_of_bnd)
    integer(ip),       intent(in)    :: patch_id
    integer,           intent(in)    :: ptset_type
    integer(cgsize_t), intent(in)    :: pnts(:)
    integer(cgsize_t), intent(in)    :: ngon_estart
    integer(ip),       intent(in)    :: orig_to_mesh_face(:)
    integer(ip),       intent(in)    :: n_int_faces
    integer(ip),       intent(inout) :: patch_of_bnd(:)

    integer(ip) :: i, gid, local_id, mesh_face, bnd_local
    integer(cgsize_t) :: gid64

    if (ptset_type == CGNS_ENUMV(PointRange)) then
      if (size(pnts) /= 2) error stop 'cgns_read: PointRange BC point set must have exactly 2 entries'
      do gid64 = pnts(1), pnts(2)
        call assign_one(int(gid64, ip))
      end do
    else if (ptset_type == CGNS_ENUMV(PointList)) then
      do i = 1, size(pnts)
        call assign_one(int(pnts(i), ip))
      end do
    else
      error stop 'cgns_read: unsupported BC PointSetType (only PointRange/PointList are supported)'
    end if

  contains
    subroutine assign_one(gid_in)
      integer(ip), intent(in) :: gid_in
      gid = gid_in
      local_id = gid - int(ngon_estart, ip) + 1
      if (local_id < 1 .or. local_id > size(orig_to_mesh_face)) &
        error stop 'cgns_read: BC references a face id outside the NGON_n section range'
      mesh_face = orig_to_mesh_face(local_id)
      bnd_local = mesh_face - n_int_faces
      if (bnd_local < 1) &
        error stop 'cgns_read: BC references an internal (non-boundary) NGON_n face'
      if (patch_of_bnd(bnd_local) /= 0) &
        error stop 'cgns_read: a face is claimed by more than one BC patch'
      patch_of_bnd(bnd_local) = patch_id
    end subroutine assign_one
  end subroutine mark_patch_faces

  !> Regroup the boundary segment [n_int_faces+1, n_faces] of the face CSR arrays so faces of the
  !> same patch are contiguous (stable within each patch), and set patches(:)%start_face/n_faces.
  subroutine regroup_boundary_by_patch(msh, n_int_faces, n_bnd, patch_of_bnd, &
                                        face_offsets, face_vertices, face_owner, face_neigh)
    type(Mesh),  intent(inout) :: msh
    integer(ip), intent(in)    :: n_int_faces, n_bnd
    integer(ip), intent(in)    :: patch_of_bnd(:)
    integer(ip), allocatable, intent(inout) :: face_offsets(:), face_vertices(:)
    integer(ip), allocatable, intent(inout) :: face_owner(:), face_neigh(:)

    integer(ip), allocatable :: counts(:), starts(:), fill(:), new_local_of_old(:)
    integer(ip), allocatable :: new_offsets(:), new_vertices(:), new_owner(:), new_neigh(:)
    integer(ip) :: p, i, np, old_f, new_f, wpos, rpos, nn, n_total

    np = msh%n_patches
    allocate(counts(np), starts(np), fill(np), source=0_ip)
    do i = 1, n_bnd
      counts(patch_of_bnd(i)) = counts(patch_of_bnd(i)) + 1
    end do
    starts(1) = 1
    do p = 2, np
      starts(p) = starts(p-1) + counts(p-1)
    end do
    do p = 1, np
      msh%patches(p)%start_face = starts(p)
      msh%patches(p)%n_faces    = counts(p)
    end do

    allocate(new_local_of_old(n_bnd))
    do i = 1, n_bnd
      p = patch_of_bnd(i)
      new_local_of_old(i) = starts(p) + fill(p)
      fill(p) = fill(p) + 1
    end do

    n_total = size(face_owner)
    allocate(new_offsets(n_total + 1), new_owner(n_total), new_neigh(n_total))
    new_offsets(1) = 1
    do old_f = 1, n_int_faces
      new_offsets(old_f + 1) = new_offsets(old_f) + (face_offsets(old_f+1) - face_offsets(old_f))
      new_owner(old_f) = face_owner(old_f)
      new_neigh(old_f) = face_neigh(old_f)
    end do
    do i = 1, n_bnd
      old_f = n_int_faces + i
      new_f = n_int_faces + new_local_of_old(i)
      new_owner(new_f) = face_owner(old_f)
      new_neigh(new_f) = face_neigh(old_f)
    end do
    ! Offsets for the boundary segment must be computed in final (new_f) order.
    do i = 1, n_bnd
      old_f = n_int_faces + i
      new_f = n_int_faces + new_local_of_old(i)
      new_offsets(new_f + 1) = face_offsets(old_f + 1) - face_offsets(old_f)  ! stash length, fixed up below
    end do
    do new_f = n_int_faces + 1, n_total
      new_offsets(new_f + 1) = new_offsets(new_f) + new_offsets(new_f + 1)
    end do

    allocate(new_vertices(size(face_vertices)))
    do old_f = 1, n_int_faces
      wpos = new_offsets(old_f)
      rpos = face_offsets(old_f)
      nn = face_offsets(old_f+1) - face_offsets(old_f)
      new_vertices(wpos:wpos+nn-1) = face_vertices(rpos:rpos+nn-1)
    end do
    do i = 1, n_bnd
      old_f = n_int_faces + i
      new_f = n_int_faces + new_local_of_old(i)
      wpos = new_offsets(new_f)
      rpos = face_offsets(old_f)
      nn = face_offsets(old_f+1) - face_offsets(old_f)
      new_vertices(wpos:wpos+nn-1) = face_vertices(rpos:rpos+nn-1)
    end do

    call move_alloc(new_offsets, face_offsets)
    call move_alloc(new_vertices, face_vertices)
    call move_alloc(new_owner, face_owner)
    call move_alloc(new_neigh, face_neigh)
  end subroutine regroup_boundary_by_patch

  pure integer(ip) function cgns_bctype_to_id(bocotype) result(bc_id)
    integer, intent(in) :: bocotype
    select case (bocotype)
    case (CGNS_ENUMV(BCWall), CGNS_ENUMV(BCWallViscous), CGNS_ENUMV(BCWallViscousHeatFlux), &
          CGNS_ENUMV(BCWallViscousIsothermal), CGNS_ENUMV(BCWallInviscid))
      bc_id = BC_WALL_MODEL
    case (CGNS_ENUMV(BCInflow), CGNS_ENUMV(BCInflowSubsonic), CGNS_ENUMV(BCInflowSupersonic))
      bc_id = BC_INLET
    case (CGNS_ENUMV(BCOutflow), CGNS_ENUMV(BCOutflowSubsonic), CGNS_ENUMV(BCOutflowSupersonic))
      bc_id = BC_OUTLET
    case (CGNS_ENUMV(BCSymmetryPlane))
      bc_id = BC_SYMMETRY
    case default
      bc_id = BC_ZEROGRADIENT
    end select
  end function cgns_bctype_to_id

  function cgns_bctype_name(bocotype) result(nm)
    integer, intent(in) :: bocotype
    character(len=64) :: nm
    select case (bocotype)
    case (CGNS_ENUMV(BCWall), CGNS_ENUMV(BCWallViscous)); nm = 'wall'
    case (CGNS_ENUMV(BCInflow));  nm = 'inflow'
    case (CGNS_ENUMV(BCOutflow)); nm = 'outflow'
    case (CGNS_ENUMV(BCSymmetryPlane)); nm = 'symmetry'
    case default; nm = 'wall'
    end select
  end function cgns_bctype_name

  subroutine validate_cgns_boundary_connectivity(msh)
    type(Mesh), intent(in) :: msh
    integer(ip) :: p, f_bnd, f_glob, n_assigned

    n_assigned = 0
    do f_bnd = 1, msh%n_bnd_faces
      f_glob = msh%n_int_faces + f_bnd
      if (msh%faces(f_glob)%patch_id < 1 .or. msh%faces(f_glob)%patch_id > msh%n_patches) cycle
      n_assigned = n_assigned + 1
    end do
    if (n_assigned /= msh%n_bnd_faces) &
      error stop 'cgns_read: boundary connectivity invalid (unmapped boundary faces)'

    write(*, '(A,I0,A,I0)') ' CGNS mesh: patches=', msh%n_patches, ', boundary_faces=', msh%n_bnd_faces
    do p = 1, msh%n_patches
      write(*, '(A,I0,A,A,A,I0)') '   patch ', p, ': name=', trim(msh%patches(p)%name), &
          ', bc_id=', msh%patches(p)%bc_id
    end do
  end subroutine validate_cgns_boundary_connectivity

  subroutine cgns_check(ier, context)
    integer,          intent(in) :: ier
    character(len=*), intent(in) :: context
    if (ier /= CG_OK) then
      write(*, '(A,A)') 'CGNS error in ', trim(context)
      error stop 'cgns_read: CGNS library call failed'
    end if
  end subroutine cgns_check

end module mod_cgns_reader
