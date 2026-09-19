! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> OpenFOAM polyMesh reader (constant/polyMesh directory).
module mod_foam_reader
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use iso_fortran_env, only: int8, int32, int64, real32, real64
  implicit none
  private

  public :: foam_read

contains

  subroutine foam_read(polyMesh_dir, msh)
    character(len=*), intent(in)    :: polyMesh_dir
    type(Mesh),       intent(inout) :: msh

    real(rp),    allocatable :: points(:,:)  ! (3, n_points)
    integer(ip), allocatable :: face_pts(:)  ! concatenated node lists
    integer(ip), allocatable :: face_offset(:)  ! start index per face (1-based)
    integer(ip), allocatable :: owner(:)  ! owner cell per face
    integer(ip), allocatable :: neighbour(:)  ! neighbour cell per face (-1 = bnd)
    integer(ip) :: n_points, n_faces, n_int_faces, n_cells
    integer(ip) :: f, c

    call foam_read_points(trim(polyMesh_dir)//'/points', points, n_points)
    call foam_read_faces(trim(polyMesh_dir)//'/faces', face_pts, face_offset, n_faces)
    call foam_read_owner(trim(polyMesh_dir)//'/owner', owner, n_faces, n_cells)
    call foam_read_neighbour(trim(polyMesh_dir)//'/neighbour', neighbour, n_int_faces)

    if (any(owner < 0) .or. any(neighbour < 0)) &
      error stop 'foam_read: negative cell index in owner/neighbour (corrupt polyMesh files)'

    ! Safeguard against underestimating cell count from owner list alone (some faces may only appear in neighbour).
    n_cells = max(maxval(owner), maxval(neighbour)) + 1

    ! OpenFOAM upper-triangular convention: owner(f) < neighbour(f) on internal faces.
    do f = 1, n_int_faces
      if (owner(f) >= neighbour(f)) then
        write(*, '(A,I0,A,I0,A,I0)') 'foam_read: face ', f, ' violates owner<neighbour (', &
                                     owner(f), ' >= ', neighbour(f)
        error stop 'foam_read: OpenFOAM upper-triangular convention violated'
      end if
    end do

    msh%n_faces     = n_faces
    msh%n_int_faces = n_int_faces
    msh%n_bnd_faces = n_faces - n_int_faces
    msh%n_cells     = n_cells
    msh%n_ghost     = 0

    allocate(msh%cells(n_cells))
    allocate(msh%faces(n_faces))
    allocate(msh%global_id(n_cells))

    ! Set face owner/neighbour
    do f = 1, n_faces
      msh%faces(f)%owner = owner(f) + 1  ! convert 0-based to 1-based
      if (f <= n_int_faces) then
        msh%faces(f)%neighbour = neighbour(f) + 1
      else
        msh%faces(f)%neighbour = -1
      end if
    end do

    ! Preserve topology metadata for downstream partitioning/output.
    call build_cell_topology_metadata(msh, points, face_pts, face_offset, owner, neighbour, n_points, n_faces, n_cells)
    allocate(msh%face_offsets(n_faces + 1))
    allocate(msh%face_vertices(size(face_pts)))
    msh%face_offsets = face_offset
    msh%face_vertices = face_pts

    ! Compute face geometry
    do f = 1, n_faces
      call compute_polygon_face_geometry(f, face_pts, face_offset, points, &
                                          msh%faces(f)%centre, &
                                          msh%faces(f)%normal, &
                                          msh%faces(f)%area)
      msh%faces(f)%sf = msh%faces(f)%area * msh%faces(f)%normal
    end do

    ! Compute cell centroids and volumes directly from the polyhedral faces.
    do c = 1, n_cells
      msh%cells(c)%centroid = 0.0_rp
      msh%cells(c)%volume   = 0.0_rp
      msh%global_id(c)      = int(c, lp)
    end do
    call compute_polyhedron_cell_moments(msh, n_cells, n_faces, face_pts, face_offset, points)

    ! Read boundary patches
    call foam_read_boundary(trim(polyMesh_dir)//'/boundary', msh, n_int_faces)
    call validate_foam_boundary_connectivity(msh)

    deallocate(points, face_pts, face_offset, owner, neighbour)
  end subroutine foam_read

  subroutine validate_foam_boundary_connectivity(msh)
    type(Mesh), intent(in) :: msh

    integer(ip) :: p, f_bnd, f_glob, n_assigned, n_periodic

    n_assigned = 0
    do f_bnd = 1, msh%n_bnd_faces
      f_glob = msh%n_int_faces + f_bnd
      if (msh%faces(f_glob)%patch_id < 1 .or. msh%faces(f_glob)%patch_id > msh%n_patches) cycle
      n_assigned = n_assigned + 1
    end do

    if (n_assigned /= msh%n_bnd_faces) then
      error stop 'foam_read: boundary connectivity invalid (some boundary faces are not mapped to a patch)'
    end if

    n_periodic = 0
    do p = 1, msh%n_patches
      if (msh%patches(p)%start_face < 1 .or. msh%patches(p)%n_faces < 0) then
        error stop 'foam_read: boundary connectivity invalid (patch has invalid start_face/n_faces)'
      end if
      if (msh%patches(p)%start_face + msh%patches(p)%n_faces - 1 > msh%n_bnd_faces) then
        error stop 'foam_read: boundary connectivity invalid (patch face range exceeds boundary face count)'
      end if
      if (msh%patches(p)%bc_id == BC_PERIODIC) n_periodic = n_periodic + 1
    end do

    if (mod(n_periodic, 2_ip) /= 0_ip) then
      error stop 'foam_read: periodic boundary connectivity invalid (odd number of periodic patches)'
    end if

    write(*, '(A,I0,A,I0,A,I0)') ' FOAM boundary summary: patches=', msh%n_patches, &
                                  ', boundary_faces=', msh%n_bnd_faces, ', periodic_patches=', n_periodic
    do p = 1, msh%n_patches
      write(*, '(A,I0,A,A,A,A,A,I0,A,I0,A,I0)') '   patch ', p, ': name=', trim(msh%patches(p)%name), &
          ', type=', trim(msh%patches(p)%bc_type), ', bc_id=', msh%patches(p)%bc_id, &
          ', start_face=', msh%patches(p)%start_face, ', n_faces=', msh%patches(p)%n_faces
    end do
  end subroutine validate_foam_boundary_connectivity

  subroutine foam_read_points(filename, points, n_points)
    character(len=*), intent(in)    :: filename
    real(rp), allocatable,intent(out) :: points(:,:)
    integer(ip),          intent(out) :: n_points

    integer(ip) :: iunit, ios, i, j
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line
    logical :: is_binary, is_little_endian
    integer(ip) :: label_bytes, scalar_bytes

    call foam_probe_file(filename, is_binary, is_little_endian, label_bytes, scalar_bytes)
    if (is_binary) then
      call foam_read_points_binary(filename, points, n_points, scalar_bytes, is_little_endian)
      return
    end if

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_read_points: cannot open '//trim(filename)

    call foam_skip_header(iunit)

    ! Read count
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '/' .and. line(2:2) == '/') cycle
      read(line, *) n_points
      exit
    end do

    allocate(points(3, n_points))

    ! Skip opening '('
    do
      read(iunit, '(A)') line
      if (adjustl(trim(line)) == '(') exit
    end do

    do i = 1, n_points
      read(iunit, '(A)') line
      do j = 1, len_trim(line)
        if (line(j:j) == '(' .or. line(j:j) == ')') line(j:j) = ' '
      end do
      read(line, *) points(1,i), points(2,i), points(3,i)
    end do

    close(iunit)
  end subroutine foam_read_points

  subroutine foam_read_faces(filename, face_pts, face_offset, n_faces)
    character(len=*), intent(in)     :: filename
    integer(ip), allocatable, intent(out) :: face_pts(:), face_offset(:)
    integer(ip),              intent(out) :: n_faces

    integer(ip)   :: iunit, ios, i, j, nn, offset
    integer(ip), allocatable :: face_npts(:)
    integer(ip), allocatable :: vals(:)
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line
    logical :: is_binary, is_little_endian
    integer(ip) :: label_bytes, scalar_bytes

    call foam_probe_file(filename, is_binary, is_little_endian, label_bytes, scalar_bytes)
    if (is_binary) then
      call foam_read_faces_binary(filename, face_pts, face_offset, n_faces, label_bytes, is_little_endian)
      return
    end if

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_read_faces: cannot open '//trim(filename)
    call foam_skip_header(iunit)

    ! Read face count
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:2) == '//') cycle
      read(line, *) n_faces
      exit
    end do

    allocate(face_offset(n_faces + 1))
    ! First pass: read number of nodes per face to get total size.
    allocate(face_npts(n_faces))

    do
      read(iunit, '(A)') line
      if (adjustl(trim(line)) == '(') exit
    end do

    do i = 1, n_faces
      read(iunit, '(A)') line
      line = adjustl(line)
      ! Format: n(p0 p1 p2 ...) - strip parens for list-directed read
      do j = 1, len_trim(line)
        if (line(j:j) == '(' .or. line(j:j) == ')') line(j:j) = ' '
      end do
      read(line, *) nn
      face_npts(i) = nn
    end do

    ! Build CSR
    face_offset(1) = 1
    do i = 1, n_faces
      face_offset(i+1) = face_offset(i) + face_npts(i)
    end do
    allocate(face_pts(face_offset(n_faces+1) - 1))

    rewind(iunit)
    call foam_skip_header(iunit)
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:2) == '//') cycle
      read(line, *) j
      exit
    end do
    do
      read(iunit, '(A)') line
      if (adjustl(trim(line)) == '(') exit
    end do

    do i = 1, n_faces
      read(iunit, '(A)') line
      line = adjustl(line)
      do j = 1, len_trim(line)
        if (line(j:j) == '(' .or. line(j:j) == ')') line(j:j) = ' '
      end do
      nn = face_npts(i)
      allocate(vals(nn + 1))
      read(line, *) vals
      if (vals(1) /= nn) then
        error stop 'foam_read_faces: inconsistent node count in faces entry'
      end if
      offset = face_offset(i)
      do j = 1, face_npts(i)
        face_pts(offset + j - 1) = vals(j + 1) + 1  ! 0->1 based
      end do
      deallocate(vals)
    end do

    deallocate(face_npts)
    close(iunit)
  end subroutine foam_read_faces

  subroutine foam_read_owner(filename, owner, n_faces, n_cells)
    character(len=*), intent(in)    :: filename
    integer(ip), allocatable, intent(out) :: owner(:)
    integer(ip), intent(out) :: n_faces, n_cells

    integer(ip) :: iunit, ios, i
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line
    logical :: is_binary, is_little_endian
    integer(ip) :: label_bytes, scalar_bytes

    call foam_probe_file(filename, is_binary, is_little_endian, label_bytes, scalar_bytes)
    if (is_binary) then
      call foam_read_label_list_binary(filename, owner, n_faces, label_bytes, is_little_endian, 'foam_read_owner')
      n_cells = maxval(owner) + 1  ! may be corrected by caller after neighbour is loaded
      return
    end if

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_read_owner: cannot open '//trim(filename)
    call foam_skip_header(iunit)

    ! Parse nCells and nFaces from FoamFile header comment
    n_cells = 0
    n_faces = -1
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line(1:2) == '//') then
        if (index(line, 'nCells') > 0) &
          call foam_parse_int_token(line, 'nCells', n_cells)
        if (index(line, 'nFaces') > 0) &
          call foam_parse_int_token(line, 'nFaces', n_faces)
      end if
      if (len_trim(line) > 0 .and. line(1:1) /= '/' .and. &
          line(1:1) /= ' ') then
        if (verify(trim(line), '0123456789') == 0) then
          read(line, *) n_faces
          exit
        end if
      end if
    end do
    if (n_faces < 0) error stop 'foam_read_owner: face count not found in '//trim(filename)

    allocate(owner(n_faces))
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) error stop 'foam_read_owner: unexpected EOF before "(" in '//trim(filename)
      if (adjustl(trim(line)) == '(') exit
    end do
    do i = 1, n_faces
      read(iunit, *, iostat=ios) owner(i)
      if (ios /= 0) error stop 'foam_read_owner: malformed entry in '//trim(filename)
    end do
    if (n_cells == 0) n_cells = maxval(owner) + 1
    ! n_cells may still underestimate; the caller (foam_read) corrects it after
    ! neighbour is loaded.

    close(iunit)
  end subroutine foam_read_owner

  subroutine foam_read_neighbour(filename, neighbour, n_int_faces)
    character(len=*), intent(in)    :: filename
    integer(ip), allocatable, intent(out) :: neighbour(:)
    integer(ip), intent(out) :: n_int_faces

    integer(ip) :: iunit, ios, i
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line
    logical :: is_binary, is_little_endian
    integer(ip) :: label_bytes, scalar_bytes

    call foam_probe_file(filename, is_binary, is_little_endian, label_bytes, scalar_bytes)
    if (is_binary) then
      call foam_read_label_list_binary(filename, neighbour, n_int_faces, label_bytes, is_little_endian, 'foam_read_neighbour')
      return
    end if

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_read_neighbour: cannot open '//trim(filename)
    call foam_skip_header(iunit)

    n_int_faces = -1
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) > 0 .and. line(1:1) /= '/' .and. &
          verify(trim(line), '0123456789') == 0) then
        read(line, *) n_int_faces
        exit
      end if
    end do
    if (n_int_faces < 0) error stop 'foam_read_neighbour: face count not found in '//trim(filename)

    allocate(neighbour(n_int_faces))
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) error stop 'foam_read_neighbour: unexpected EOF before "(" in '//trim(filename)
      if (adjustl(trim(line)) == '(') exit
    end do
    do i = 1, n_int_faces
      read(iunit, *, iostat=ios) neighbour(i)
      if (ios /= 0) error stop 'foam_read_neighbour: malformed entry in '//trim(filename)
    end do

    close(iunit)
  end subroutine foam_read_neighbour

  subroutine foam_read_boundary(filename, msh, n_int_faces)
    character(len=*), intent(in)    :: filename
    type(Mesh),       intent(inout) :: msh
    integer(ip),      intent(in)    :: n_int_faces

    integer(ip)   :: iunit, ios, np, p, face_start, nfaces
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line, bc_type, patch_name

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_read_boundary: cannot open '//trim(filename)
    call foam_skip_header(iunit)

    ! Read number of patches
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (verify(trim(line),'0123456789') == 0 .and. len_trim(line) > 0) then
        read(line, *) np
        exit
      end if
    end do

    msh%n_patches = np
    allocate(msh%patches(np))

    do
      read(iunit, '(A)') line
      if (adjustl(trim(line)) == '(') exit
    end do

    do p = 1, np
      ! Read patch name
      do
        read(iunit, '(A)') line
        line = adjustl(line)
        if (len_trim(line) > 0 .and. line(1:1) /= '{') exit
      end do
      patch_name = trim(line)

      ! Read opening brace
      do
        read(iunit, '(A)') line
        if (adjustl(trim(line)) == '{') exit
      end do

      bc_type    = 'patch'
      face_start = 0
      nfaces     = 0

      do
        read(iunit, '(A)') line
        line = adjustl(line)
        if (trim(line) == '}') exit
        if (index(line, 'type') > 0) &
          call foam_read_value_str(line, bc_type)
        if (index(line, 'startFace') > 0) &
          call foam_read_value_int(line, face_start)
        if (index(line, 'nFaces') > 0) &
          call foam_read_value_int(line, nfaces)
      end do

      msh%patches(p)%name       = trim(patch_name)
      msh%patches(p)%bc_type    = trim(bc_type)
      msh%patches(p)%start_face = face_start - n_int_faces + 1  ! 1-based local bnd index
      msh%patches(p)%n_faces    = nfaces
      msh%patches(p)%bc_id      = bc_type_to_id(trim(bc_type))

      ! Tag boundary faces with patch ID
      associate(f0 => face_start + 1, fn => face_start + nfaces)
        do ios = f0, fn
          msh%faces(ios)%patch_id = p
        end do
      end associate
    end do
    close(iunit)
  end subroutine foam_read_boundary

  subroutine build_cell_topology_metadata(msh, points, face_pts, face_offset, owner, neighbour, n_points, n_faces, n_cells)
    type(Mesh),       intent(inout) :: msh
    real(rp),         intent(in)    :: points(:,:)
    integer(ip),      intent(in)    :: face_pts(:), face_offset(:)
    integer(ip),      intent(in)    :: owner(:), neighbour(:)
    integer(ip),      intent(in)    :: n_points, n_faces, n_cells

    integer(ip) :: f, c, idx, pos, total_incidence, write_pos, nverts
    integer(ip) :: ntri, nquad, nface
    integer(ip), allocatable :: cell_face_count(:), cell_face_offset(:), cell_face_ids(:)
    integer(ip), allocatable :: fill_count(:), cell_vertex_count(:), vertex_marker(:)

    msh%n_vertices = n_points
    allocate(msh%vertices(3, n_points))
    msh%vertices = points

    allocate(cell_face_count(n_cells), source=0_ip)
    do f = 1, n_faces
      cell_face_count(owner(f) + 1) = cell_face_count(owner(f) + 1) + 1
    end do
    do f = 1, size(neighbour)
      cell_face_count(neighbour(f) + 1) = cell_face_count(neighbour(f) + 1) + 1
    end do

    allocate(cell_face_offset(n_cells + 1))
    cell_face_offset(1) = 1
    do c = 1, n_cells
      cell_face_offset(c + 1) = cell_face_offset(c) + cell_face_count(c)
    end do
    total_incidence = cell_face_offset(n_cells + 1) - 1
    allocate(cell_face_ids(total_incidence), fill_count(n_cells), source=0_ip)

    do f = 1, n_faces
      c = owner(f) + 1
      pos = cell_face_offset(c) + fill_count(c)
      cell_face_ids(pos) = f
      fill_count(c) = fill_count(c) + 1
    end do
    do f = 1, size(neighbour)
      c = neighbour(f) + 1
      pos = cell_face_offset(c) + fill_count(c)
      cell_face_ids(pos) = f
      fill_count(c) = fill_count(c) + 1
    end do

    ! Deallocate working arrays no longer needed for vertex counting
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
      ntri = 0
      nquad = 0
      nface = cell_face_offset(c + 1) - cell_face_offset(c)

      do idx = cell_face_offset(c), cell_face_offset(c + 1) - 1
        f = cell_face_ids(idx)
        select case (face_offset(f + 1) - face_offset(f))
        case (3)
          ntri = ntri + 1
        case (4)
          nquad = nquad + 1
        end select

        do pos = face_offset(f), face_offset(f + 1) - 1
          if (vertex_marker(face_pts(pos)) /= c) then
            vertex_marker(face_pts(pos)) = c
            msh%cell_vertices(write_pos + nverts) = face_pts(pos)
            nverts = nverts + 1
          end if
        end do
      end do

      ! OpenFOAM cells are native polyhedra; keep face-based topology for output.
      msh%cell_types(c) = 42_ip
    end do

    deallocate(cell_face_offset, cell_face_ids, cell_vertex_count, vertex_marker)
  end subroutine build_cell_topology_metadata

  ! Utility: skip FoamFile header block
  subroutine foam_skip_header(iunit)
    integer(ip), intent(in) :: iunit
    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    character(FOAM_LINE_LEN) :: line
    integer(ip)    :: ios, depth

    depth = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) return
      if (index(line, 'FoamFile') > 0) then
        ! Read until matching '}'
        do
          read(iunit, '(A)', iostat=ios) line
          if (ios /= 0) return
          if (index(line, '{') > 0) depth = depth + 1
          if (index(line, '}') > 0) then
            depth = depth - 1
            if (depth == 0) return
          end if
        end do
      end if
    end do
  end subroutine foam_skip_header

  subroutine foam_probe_file(filename, is_binary, is_little_endian, label_bytes, scalar_bytes)
    character(len=*), intent(in)  :: filename
    logical,          intent(out) :: is_binary, is_little_endian
    integer(ip),      intent(out) :: label_bytes, scalar_bytes

    integer(ip), parameter :: FOAM_LINE_LEN = 16384
    integer(ip) :: iunit, ios, depth
    character(FOAM_LINE_LEN) :: line, aline

    is_binary = .false.
    is_little_endian = .true.
    label_bytes = 4
    scalar_bytes = 8

    open(newunit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) error stop 'foam_probe_file: cannot open '//trim(filename)

    depth = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      aline = adjustl(line)
      if (index(aline, 'FoamFile') == 1) then
        do
          read(iunit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          aline = adjustl(line)
          if (index(aline, 'format') == 1) then
            if (index(aline, 'binary') > 0) is_binary = .true.
            if (index(aline, 'ascii') > 0) is_binary = .false.
          end if
          if (index(aline, 'arch') == 1) then
            if (index(aline, 'MSB') > 0) is_little_endian = .false.
            if (index(aline, 'LSB') > 0) is_little_endian = .true.
            if (index(aline, 'label=64') > 0) label_bytes = 8
            if (index(aline, 'label=32') > 0) label_bytes = 4
            if (index(aline, 'scalar=32') > 0) scalar_bytes = 4
            if (index(aline, 'scalar=64') > 0) scalar_bytes = 8
          end if
          if (index(line, '{') > 0) depth = depth + 1
          if (index(line, '}') > 0) then
            depth = depth - 1
            if (depth == 0) exit
          end if
        end do
        exit
      end if
    end do

    close(iunit)
  end subroutine foam_probe_file

  subroutine foam_read_points_binary(filename, points, n_points, scalar_bytes, is_little_endian)
    character(len=*), intent(in) :: filename
    real(rp), allocatable, intent(out) :: points(:,:)
    integer(ip), intent(out) :: n_points
    integer(ip), intent(in) :: scalar_bytes
    logical, intent(in) :: is_little_endian

    integer(int8), allocatable :: bytes(:)
    integer(ip) :: i, pos, list_open_pos, list_close_pos

    call foam_read_raw_bytes(filename, bytes)
    call foam_binary_locate_list(bytes, n_points, list_open_pos, list_close_pos, pos, 'foam_read_points_binary')

    allocate(points(3, n_points))
    do i = 1, n_points
      points(1, i) = foam_read_binary_scalar(bytes, pos, scalar_bytes, is_little_endian)
      points(2, i) = foam_read_binary_scalar(bytes, pos, scalar_bytes, is_little_endian)
      points(3, i) = foam_read_binary_scalar(bytes, pos, scalar_bytes, is_little_endian)
    end do

    deallocate(bytes)
  end subroutine foam_read_points_binary

  subroutine foam_read_faces_binary(filename, face_pts, face_offset, n_faces, label_bytes, is_little_endian)
    character(len=*), intent(in) :: filename
    integer(ip), allocatable, intent(out) :: face_pts(:), face_offset(:)
    integer(ip), intent(out) :: n_faces
    integer(ip), intent(in) :: label_bytes
    logical, intent(in) :: is_little_endian

    integer(int8), allocatable :: bytes(:)
    integer(ip), allocatable :: face_npts(:), compact_off(:)
    integer(ip) :: i, j, nn, pos, pos_data, list_open_pos, list_close_pos
    integer(ip) :: data_start_label, n_compact_pts, n_offsets
    logical :: use_compact
    character(1) :: ch

    call foam_read_raw_bytes(filename, bytes)
    call foam_binary_locate_list(bytes, n_faces, list_open_pos, list_close_pos, pos_data, 'foam_read_faces_binary')

    data_start_label = foam_peek_binary_label(bytes, pos_data, label_bytes, is_little_endian)
    use_compact = (data_start_label == 0)

    if (use_compact) then
      ! OpenFOAM compact face storage: first list is offset table (nFaces+1),
      ! followed by a second list containing concatenated face point labels.
      n_offsets = n_faces
      if (n_offsets < 2) error stop 'foam_read_faces_binary: compact-face offset table too small'

      allocate(compact_off(n_offsets))
      pos = pos_data
      do i = 1, n_offsets
        compact_off(i) = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian)
      end do

      if (compact_off(1) /= 0) then
        error stop 'foam_read_faces_binary: compact-face offset table must start at 0'
      end if
      do i = 2, n_offsets
        if (compact_off(i) < compact_off(i - 1)) then
          error stop 'foam_read_faces_binary: compact-face offsets are not monotonic'
        end if
      end do

      call foam_skip_ws_and_comments_bytes(bytes, pos)
      if (pos > size(bytes)) error stop 'foam_read_faces_binary: unexpected EOF after compact offset table'
      ch = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
      if (ch /= ')') error stop 'foam_read_faces_binary: expected closing parenthesis after compact offset table'
      pos = pos + 1

      call foam_parse_int_bytes(bytes, pos, n_compact_pts, 'foam_read_faces_binary')
      call foam_skip_ws_and_comments_bytes(bytes, pos)
      if (pos > size(bytes)) error stop 'foam_read_faces_binary: unexpected EOF before compact connectivity list'
      ch = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
      if (ch /= '(') error stop 'foam_read_faces_binary: expected opening parenthesis for compact connectivity list'
      pos = pos + 1
      call foam_skip_ws_bytes(bytes, pos)

      n_faces = n_offsets - 1
      allocate(face_offset(n_faces + 1))
      do i = 1, n_faces + 1
        face_offset(i) = compact_off(i) + 1
      end do

      if (n_compact_pts < 0) error stop 'foam_read_faces_binary: negative compact-face payload size'
      if (compact_off(n_faces + 1) /= n_compact_pts) then
        error stop 'foam_read_faces_binary: compact-face offset terminator does not match connectivity count'
      end if

      allocate(face_pts(n_compact_pts))
      do i = 1, n_compact_pts
        face_pts(i) = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian) + 1
      end do

      deallocate(compact_off, bytes)
      return
    end if

    allocate(face_offset(n_faces + 1))
    allocate(face_npts(n_faces))

    pos = pos_data
    do i = 1, n_faces
      nn = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian)
      if (nn < 0) error stop 'foam_read_faces_binary: negative node count in binary faces list'
      face_npts(i) = nn
      do j = 1, nn
        nn = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian)
      end do
    end do

    face_offset(1) = 1
    do i = 1, n_faces
      face_offset(i + 1) = face_offset(i) + face_npts(i)
    end do
    allocate(face_pts(face_offset(n_faces + 1) - 1))

    pos = pos_data
    do i = 1, n_faces
      nn = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian)
      if (nn /= face_npts(i)) then
        error stop 'foam_read_faces_binary: inconsistent node count in binary faces entry'
      end if
      do j = 1, nn
        face_pts(face_offset(i) + j - 1) = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian) + 1
      end do
    end do

    deallocate(face_npts, bytes)
  end subroutine foam_read_faces_binary

  integer(ip) function foam_peek_binary_label(bytes, pos, label_bytes, is_little_endian) result(val)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(in) :: pos
    integer(ip), intent(in) :: label_bytes
    logical, intent(in) :: is_little_endian

    integer(ip) :: tmp_pos

    tmp_pos = pos
    val = foam_read_binary_label(bytes, tmp_pos, label_bytes, is_little_endian)
  end function foam_peek_binary_label

  subroutine foam_read_label_list_binary(filename, vals, n_vals, label_bytes, is_little_endian, context)
    character(len=*), intent(in) :: filename
    integer(ip), allocatable, intent(out) :: vals(:)
    integer(ip), intent(out) :: n_vals
    integer(ip), intent(in) :: label_bytes
    logical, intent(in) :: is_little_endian
    character(len=*), intent(in) :: context

    integer(int8), allocatable :: bytes(:)
    integer(ip) :: i, pos, list_open_pos, list_close_pos

    call foam_read_raw_bytes(filename, bytes)
    call foam_binary_locate_list(bytes, n_vals, list_open_pos, list_close_pos, pos, trim(context))

    allocate(vals(n_vals))
    do i = 1, n_vals
      vals(i) = foam_read_binary_label(bytes, pos, label_bytes, is_little_endian)
    end do

    deallocate(bytes)
  end subroutine foam_read_label_list_binary

  subroutine foam_read_raw_bytes(filename, bytes)
    character(len=*), intent(in) :: filename
    integer(int8), allocatable, intent(out) :: bytes(:)

    integer(ip) :: iunit, ios
    integer(int64) :: nbytes

    open(newunit=iunit, file=trim(filename), status='old', access='stream', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) error stop 'foam_read_raw_bytes: cannot open '//trim(filename)

    inquire(unit=iunit, size=nbytes)
    if (nbytes <= 0_int64) error stop 'foam_read_raw_bytes: empty file '//trim(filename)
    allocate(bytes(int(nbytes, ip)))
    read(iunit, pos=1, iostat=ios) bytes
    if (ios /= 0) error stop 'foam_read_raw_bytes: failed reading '//trim(filename)
    close(iunit)
  end subroutine foam_read_raw_bytes

  subroutine foam_binary_locate_list(bytes, n_entries, list_open_pos, list_close_pos, data_pos, context)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(out) :: n_entries, list_open_pos, list_close_pos, data_pos
    character(len=*), intent(in) :: context

    integer(ip) :: pos, hdr_pos

    hdr_pos = foam_find_token(bytes, 'FoamFile', 1_ip)
    if (hdr_pos <= 0) error stop trim(context)//': could not find FoamFile header'

    pos = hdr_pos
    do
      if (pos > size(bytes)) error stop trim(context)//': malformed FoamFile header'
      if (achar(int(iand(int(bytes(pos), int32), int(z'FF', int32)))) == '}') exit
      pos = pos + 1
    end do
    pos = pos + 1

    call foam_skip_ws_and_comments_bytes(bytes, pos)
    call foam_parse_int_bytes(bytes, pos, n_entries, context)
    call foam_skip_ws_and_comments_bytes(bytes, pos)

    list_open_pos = foam_find_char(bytes, '(', pos)
    if (list_open_pos <= 0) error stop trim(context)//': could not find list opening parenthesis'

    list_close_pos = foam_rfind_char(bytes, ')')
    if (list_close_pos <= list_open_pos) error stop trim(context)//': could not find list closing parenthesis'

    data_pos = list_open_pos + 1
    call foam_skip_ws_bytes(bytes, data_pos)
  end subroutine foam_binary_locate_list

  subroutine foam_skip_ws_bytes(bytes, pos)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    character(1) :: ch

    do while (pos <= size(bytes))
      ch = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
      if (ch /= ' ' .and. ch /= new_line('a') .and. ch /= char(9) .and. ch /= char(13)) exit
      pos = pos + 1
    end do
  end subroutine foam_skip_ws_bytes

  subroutine foam_skip_ws_and_comments_bytes(bytes, pos)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    character(1) :: c1, c2

    do
      call foam_skip_ws_bytes(bytes, pos)
      if (pos + 1 > size(bytes)) return
      c1 = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
      c2 = achar(int(iand(int(bytes(pos + 1), int32), int(z'FF', int32))))
      if (c1 == '/' .and. c2 == '/') then
        do while (pos <= size(bytes))
          c1 = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
          pos = pos + 1
          if (c1 == new_line('a')) exit
        end do
      else
        exit
      end if
    end do
  end subroutine foam_skip_ws_and_comments_bytes

  subroutine foam_parse_int_bytes(bytes, pos, val, context)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    integer(ip), intent(out) :: val
    character(len=*), intent(in) :: context

    integer(int64) :: accum
    integer(ip) :: sign
    character(1) :: ch

    call foam_skip_ws_and_comments_bytes(bytes, pos)
    if (pos > size(bytes)) error stop trim(context)//': expected integer but reached EOF'

    sign = 1
    ch = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
    if (ch == '-') then
      sign = -1
      pos = pos + 1
    end if

    accum = 0_int64
    do while (pos <= size(bytes))
      ch = achar(int(iand(int(bytes(pos), int32), int(z'FF', int32))))
      if (ch < '0' .or. ch > '9') exit
      accum = accum * 10_int64 + int(iachar(ch) - iachar('0'), int64)
      pos = pos + 1
    end do

    if (accum > int(huge(val), int64)) then
      error stop trim(context)//': integer value out of range while parsing binary list header'
    end if
    val = int(sign * accum, ip)
  end subroutine foam_parse_int_bytes

  integer(ip) function foam_find_token(bytes, token, start_pos) result(pos)
    integer(int8), intent(in) :: bytes(:)
    character(len=*), intent(in) :: token
    integer(ip), intent(in) :: start_pos

    integer(ip) :: i, j, n
    character(1) :: bch

    pos = -1
    n = len_trim(token)
    if (n <= 0) return

    do i = max(1_ip, start_pos), size(bytes) - n + 1
      do j = 1, n
        bch = achar(int(iand(int(bytes(i + j - 1), int32), int(z'FF', int32))))
        if (bch /= token(j:j)) exit
      end do
      if (j > n) then
        pos = i
        return
      end if
    end do
  end function foam_find_token

  integer(ip) function foam_find_char(bytes, ch, start_pos) result(pos)
    integer(int8), intent(in) :: bytes(:)
    character(1), intent(in) :: ch
    integer(ip), intent(in) :: start_pos

    integer(ip) :: i
    character(1) :: bch

    pos = -1
    do i = max(1_ip, start_pos), size(bytes)
      bch = achar(int(iand(int(bytes(i), int32), int(z'FF', int32))))
      if (bch == ch) then
        pos = i
        return
      end if
    end do
  end function foam_find_char

  integer(ip) function foam_rfind_char(bytes, ch) result(pos)
    integer(int8), intent(in) :: bytes(:)
    character(1), intent(in) :: ch

    integer(ip) :: i
    character(1) :: bch

    pos = -1
    do i = size(bytes), 1, -1
      bch = achar(int(iand(int(bytes(i), int32), int(z'FF', int32))))
      if (bch == ch) then
        pos = i
        return
      end if
    end do
  end function foam_rfind_char

  integer(ip) function foam_read_binary_label(bytes, pos, label_bytes, is_little_endian) result(val)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    integer(ip), intent(in) :: label_bytes
    logical, intent(in) :: is_little_endian

    integer(int64) :: bits64
    integer(int32) :: bits32

    select case (label_bytes)
    case (4)
      bits32 = foam_read_u32(bytes, pos, is_little_endian)
      val = int(bits32, ip)
    case (8)
      bits64 = foam_read_u64(bytes, pos, is_little_endian)
      if (bits64 < (-int(huge(val), int64) - 1_int64) .or. bits64 > int(huge(val), int64)) then
        error stop 'foam_read_binary_label: 64-bit label out of 32-bit range'
      end if
      val = int(bits64, ip)
    case default
      error stop 'foam_read_binary_label: unsupported label byte width'
    end select
  end function foam_read_binary_label

  real(rp) function foam_read_binary_scalar(bytes, pos, scalar_bytes, is_little_endian) result(val)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    integer(ip), intent(in) :: scalar_bytes
    logical, intent(in) :: is_little_endian

    integer(int64) :: bits64
    integer(int32) :: bits32
    real(real64) :: r64
    real(real32) :: r32

    select case (scalar_bytes)
    case (8)
      bits64 = foam_read_u64(bytes, pos, is_little_endian)
      r64 = transfer(bits64, r64)
      val = real(r64, rp)
    case (4)
      bits32 = foam_read_u32(bytes, pos, is_little_endian)
      r32 = transfer(bits32, r32)
      val = real(r32, rp)
    case default
      error stop 'foam_read_binary_scalar: unsupported scalar byte width'
    end select
  end function foam_read_binary_scalar

  integer(int32) function foam_read_u32(bytes, pos, is_little_endian) result(val)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    logical, intent(in) :: is_little_endian

    integer(int64) :: u
    integer(int64) :: b1, b2, b3, b4

    if (pos + 3 > size(bytes)) error stop 'foam_read_u32: unexpected EOF'

    b1 = iand(int(bytes(pos    ), int64), int(z'FF', int64))
    b2 = iand(int(bytes(pos + 1), int64), int(z'FF', int64))
    b3 = iand(int(bytes(pos + 2), int64), int(z'FF', int64))
    b4 = iand(int(bytes(pos + 3), int64), int(z'FF', int64))

    if (is_little_endian) then
      u = b1 + ishft(b2, 8) + ishft(b3, 16) + ishft(b4, 24)
    else
      u = b4 + ishft(b3, 8) + ishft(b2, 16) + ishft(b1, 24)
    end if

    val = int(u, int32)
    pos = pos + 4
  end function foam_read_u32

  integer(int64) function foam_read_u64(bytes, pos, is_little_endian) result(val)
    integer(int8), intent(in) :: bytes(:)
    integer(ip), intent(inout) :: pos
    logical, intent(in) :: is_little_endian

    integer(int64) :: b(8)

    if (pos + 7 > size(bytes)) error stop 'foam_read_u64: unexpected EOF'

    b(1) = iand(int(bytes(pos    ), int64), int(z'FF', int64))
    b(2) = iand(int(bytes(pos + 1), int64), int(z'FF', int64))
    b(3) = iand(int(bytes(pos + 2), int64), int(z'FF', int64))
    b(4) = iand(int(bytes(pos + 3), int64), int(z'FF', int64))
    b(5) = iand(int(bytes(pos + 4), int64), int(z'FF', int64))
    b(6) = iand(int(bytes(pos + 5), int64), int(z'FF', int64))
    b(7) = iand(int(bytes(pos + 6), int64), int(z'FF', int64))
    b(8) = iand(int(bytes(pos + 7), int64), int(z'FF', int64))

    if (is_little_endian) then
      val = b(1) + ishft(b(2), 8) + ishft(b(3), 16) + ishft(b(4), 24) + &
            ishft(b(5), 32) + ishft(b(6), 40) + ishft(b(7), 48) + ishft(b(8), 56)
    else
      val = b(8) + ishft(b(7), 8) + ishft(b(6), 16) + ishft(b(5), 24) + &
            ishft(b(4), 32) + ishft(b(3), 40) + ishft(b(2), 48) + ishft(b(1), 56)
    end if

    pos = pos + 8
  end function foam_read_u64

  subroutine foam_parse_int_token(line, key, val)
    character(len=*), intent(in)  :: line, key
    integer(ip),      intent(out) :: val
    integer(ip) :: pos
    pos = index(line, trim(key))
    if (pos > 0) read(line(pos+len_trim(key):), *) val
  end subroutine foam_parse_int_token

  subroutine foam_read_value_str(line, val)
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: val
    integer(ip) :: pos
    pos = index(line, ' ')
    if (pos > 0) then
      val = adjustl(line(pos:))
      pos = index(val, ';')
      if (pos > 0) val(pos:) = ' '
      val = adjustl(trim(val))
    end if
  end subroutine foam_read_value_str

  subroutine foam_read_value_int(line, val)
    character(len=*), intent(in)  :: line
    integer(ip),      intent(out) :: val
    integer(ip) :: pos, ios
    character(64) :: tmp
    val = 0
    pos = index(line, ' ')
    if (pos <= 0) return

    tmp = adjustl(line(pos:))
    pos = index(tmp, ';')
    if (pos > 0) tmp(pos:) = ' '
    read(tmp, *, iostat=ios) val
    if (ios /= 0) error stop 'foam_read_value_int: failed to parse integer from line: '//trim(line)
  end subroutine foam_read_value_int

  pure integer(ip) function bc_type_to_id(bctype)
    character(len=*), intent(in) :: bctype
    select case (trim(bctype))
    case ('inlet', 'fixedValue', 'inletOutlet')
      bc_type_to_id = BC_INLET
    case ('outlet', 'inletOutlet_out')
      bc_type_to_id = BC_OUTLET
    case ('zeroGradient')
      bc_type_to_id = BC_ZEROGRADIENT
    case ('wall')
      bc_type_to_id = BC_WALL_MODEL
    case ('fixedTemperature')
      bc_type_to_id = BC_FIXEDTEMP
    case ('movingWall', 'wallVelocity')
      bc_type_to_id = BC_WALL_VELOCITY
    case ('symmetry', 'symmetryPlane')
      bc_type_to_id = BC_SYMMETRY
    case ('empty')
      bc_type_to_id = BC_SLIP
    case ('cyclic', 'cyclicAMI')
      bc_type_to_id = BC_PERIODIC
    case default
      bc_type_to_id = BC_ZEROGRADIENT
    end select
  end function bc_type_to_id

end module mod_foam_reader
