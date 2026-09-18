! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Multi-scalar transport configuration and source helpers.
module mod_scalar
  use mod_kinds
  use mod_parameters, only: SMALL
  implicit none
  private

  integer(ip), parameter, public :: MAX_PASSIVE_SCALARS = 16
  integer(ip), parameter, public :: MAX_SCALAR_SOURCES = 128

  integer(ip), parameter, public :: SCALAR_BC_DIRICHLET = 1_ip
  integer(ip), parameter, public :: SCALAR_BC_NEUMANN = 2_ip
  integer(ip), parameter, public :: SCALAR_BC_ZERO_GRADIENT = 3_ip

  integer(ip), parameter, public :: SCALAR_SRC_POINT = 1_ip
  integer(ip), parameter, public :: SCALAR_SRC_LINE = 2_ip
  integer(ip), parameter, public :: SCALAR_SRC_SURFACE = 3_ip
  integer(ip), parameter, public :: SCALAR_SRC_VOLUME = 4_ip
  integer(ip), parameter, public :: SCALAR_SRC_OBJ = 5_ip

  integer(ip), parameter, public :: SCALAR_SHAPE_SPHERE = 1_ip
  integer(ip), parameter, public :: SCALAR_SHAPE_BOX = 2_ip
  integer(ip), parameter, public :: SCALAR_SHAPE_PLANE = 3_ip

  integer(ip), parameter, public :: SCALAR_TIME_CONSTANT = 1_ip
  integer(ip), parameter, public :: SCALAR_TIME_HARMONIC = 2_ip

  type, public :: ScalarSourceDef
    logical :: enabled = .false.
    integer(ip) :: scalar_id = 1_ip
    integer(ip) :: src_kind = SCALAR_SRC_VOLUME
    integer(ip) :: shape = SCALAR_SHAPE_SPHERE
    integer(ip) :: time_mode = SCALAR_TIME_CONSTANT
    character(len=64) :: name = ''
    character(len=256) :: file = ''
    real(rp) :: magnitude = 0.0_rp
    real(rp) :: omega = 0.0_rp
    real(rp) :: phase = 0.0_rp
    real(rp) :: center(3) = 0.0_rp
    real(rp) :: p0(3) = 0.0_rp
    real(rp) :: p1(3) = 0.0_rp
    real(rp) :: normal(3) = [0.0_rp, 0.0_rp, 1.0_rp]
    real(rp) :: half_thickness = 0.0_rp
    real(rp) :: radius = 0.0_rp
    real(rp) :: box_min(3) = 0.0_rp
    real(rp) :: box_max(3) = 0.0_rp
    real(rp) :: tag_radius = 0.0_rp
    real(rp), allocatable :: tag_points(:,:)  ! (3,n)
  end type ScalarSourceDef

  public :: scalar_bc_mode_from_name
  public :: scalar_source_kind_from_name
  public :: scalar_source_shape_from_name
  public :: scalar_time_mode_from_name
  public :: scalar_source_prepare
  public :: scalar_source_value_at
  public :: source_time_scale

contains

  pure integer(ip) function scalar_bc_mode_from_name(name) result(mode)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(name)))
    case ('dirichlet', 'fixed', 'fixedvalue')
      mode = SCALAR_BC_DIRICHLET
    case ('neumann', 'flux')
      mode = SCALAR_BC_NEUMANN
    case default
      mode = SCALAR_BC_ZERO_GRADIENT
    end select
  end function scalar_bc_mode_from_name

  pure integer(ip) function scalar_source_kind_from_name(name) result(kind)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(name)))
    case ('point')
      kind = SCALAR_SRC_POINT
    case ('line')
      kind = SCALAR_SRC_LINE
    case ('surface', 'plane')
      kind = SCALAR_SRC_SURFACE
    case ('obj', 'stl', 'points')
      kind = SCALAR_SRC_OBJ
    case default
      kind = SCALAR_SRC_VOLUME
    end select
  end function scalar_source_kind_from_name

  pure integer(ip) function scalar_source_shape_from_name(name) result(shape)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(name)))
    case ('box', 'aabb')
      shape = SCALAR_SHAPE_BOX
    case ('plane', 'surface')
      shape = SCALAR_SHAPE_PLANE
    case default
      shape = SCALAR_SHAPE_SPHERE
    end select
  end function scalar_source_shape_from_name

  pure integer(ip) function scalar_time_mode_from_name(name) result(mode)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(name)))
    case ('harmonic', 'sin', 'cos')
      mode = SCALAR_TIME_HARMONIC
    case default
      mode = SCALAR_TIME_CONSTANT
    end select
  end function scalar_time_mode_from_name

  pure real(rp) function source_time_scale(src, t) result(scale)
    type(ScalarSourceDef), intent(in) :: src
    real(rp), intent(in) :: t

    select case (src%time_mode)
    case (SCALAR_TIME_HARMONIC)
      scale = cos(src%omega * t + src%phase)
    case default
      scale = 1.0_rp
    end select
  end function source_time_scale

  subroutine scalar_source_prepare(src)
    type(ScalarSourceDef), intent(inout) :: src

    if (src%src_kind /= SCALAR_SRC_OBJ) return
    if (len_trim(src%file) <= 0) return
    call read_points_from_file(src%file, src%tag_points)
  end subroutine scalar_source_prepare

  pure real(rp) function scalar_source_value_at(src, x, t) result(val)
    type(ScalarSourceDef), intent(in) :: src
    real(rp), intent(in) :: x(3)
    real(rp), intent(in) :: t

    real(rp) :: w, r2, d

    if (.not. src%enabled) then
      val = 0.0_rp
      return
    end if

    w = 0.0_rp
    select case (src%src_kind)
    case (SCALAR_SRC_POINT)
      r2 = sum((x - src%center)**2)
      if (r2 <= max(src%radius, SMALL)**2) w = 1.0_rp
    case (SCALAR_SRC_LINE)
      d = point_segment_distance(x, src%p0, src%p1)
      if (d <= max(src%radius, SMALL)) w = 1.0_rp
    case (SCALAR_SRC_SURFACE)
      d = abs(dot_product(x - src%center, unit3(src%normal)))
      if (d <= max(src%half_thickness, SMALL)) then
        if (src%radius > 0.0_rp) then
          if (sum((x - src%center)**2) <= src%radius**2) w = 1.0_rp
        else
          w = 1.0_rp
        end if
      end if
    case (SCALAR_SRC_OBJ)
      if (allocated(src%tag_points)) then
        if (size(src%tag_points, 2) > 0) then
          d = min_distance_to_points(x, src%tag_points)
          if (d <= max(src%tag_radius, SMALL)) w = 1.0_rp
        end if
      end if
    case default
      select case (src%shape)
      case (SCALAR_SHAPE_BOX)
        if (all(x >= src%box_min) .and. all(x <= src%box_max)) w = 1.0_rp
      case (SCALAR_SHAPE_PLANE)
        d = abs(dot_product(x - src%center, unit3(src%normal)))
        if (d <= max(src%half_thickness, SMALL)) w = 1.0_rp
      case default
        r2 = sum((x - src%center)**2)
        if (r2 <= max(src%radius, SMALL)**2) w = 1.0_rp
      end select
    end select

    val = src%magnitude * source_time_scale(src, t) * w
  end function scalar_source_value_at

  pure real(rp) function point_segment_distance(x, a, b) result(d)
    real(rp), intent(in) :: x(3), a(3), b(3)
    real(rp) :: ab(3), t, proj(3), den

    ab = b - a
    den = dot_product(ab, ab)
    if (den <= SMALL) then
      d = sqrt(sum((x - a)**2))
      return
    end if

    t = dot_product(x - a, ab) / den
    t = max(0.0_rp, min(1.0_rp, t))
    proj = a + t * ab
    d = sqrt(sum((x - proj)**2))
  end function point_segment_distance

  pure real(rp) function min_distance_to_points(x, pts) result(d)
    real(rp), intent(in) :: x(3)
    real(rp), intent(in) :: pts(:,:)
    integer(ip) :: i
    real(rp) :: di

    d = huge(1.0_rp)
    do i = 1, size(pts, 2)
      di = sqrt(sum((x - pts(:, i))**2))
      if (di < d) d = di
    end do
  end function min_distance_to_points

  pure function unit3(v) result(u)
    real(rp), intent(in) :: v(3)
    real(rp) :: u(3)
    real(rp) :: n

    n = sqrt(sum(v**2))
    if (n <= SMALL) then
      u = [0.0_rp, 0.0_rp, 1.0_rp]
    else
      u = v / n
    end if
  end function unit3

  subroutine read_points_from_file(filename, points)
    character(len=*), intent(in) :: filename
    real(rp), allocatable, intent(out) :: points(:,:)

    integer(ip) :: unit_id, ios, n_points, n_alloc
    character(len=1024) :: line
    real(rp), allocatable :: tmp(:,:)
    real(rp) :: x, y, z

    n_points = 0
    n_alloc = 256
    allocate(tmp(3, n_alloc))

    open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      deallocate(tmp)
      allocate(points(3,0))
      return
    end if

    do
      read(unit_id, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call parse_point_line(line, x, y, z, ios)
      if (ios /= 0) cycle
      n_points = n_points + 1
      if (n_points > n_alloc) call grow_points(tmp, n_alloc)
      tmp(:, n_points) = [x, y, z]
    end do
    close(unit_id)

    allocate(points(3, n_points))
    if (n_points > 0) points = tmp(:, 1:n_points)
    deallocate(tmp)
  end subroutine read_points_from_file

  subroutine parse_point_line(line, x, y, z, ios)
    character(len=*), intent(in)  :: line
    real(rp),         intent(out) :: x, y, z
    integer(ip),      intent(out) :: ios

    character(len=1024) :: s

    s = adjustl(trim(line))
    ios = 0

    if (len_trim(s) == 0) then
      ios = 1
      return
    end if
    if (s(1:1) == '#') then
      ios = 1
      return
    end if
    if (index(s, 'vertex ') == 1) then
      read(s(8:), *, iostat=ios) x, y, z
      return
    end if
    if (index(s, 'v ') == 1) then
      read(s(3:), *, iostat=ios) x, y, z
      return
    end if

    read(s, *, iostat=ios) x, y, z
  end subroutine parse_point_line

  subroutine grow_points(buf, n_alloc)
    real(rp), allocatable, intent(inout) :: buf(:,:)
    integer(ip), intent(inout) :: n_alloc

    real(rp), allocatable :: tmp(:,:)
    integer(ip) :: old_n

    old_n = n_alloc
    allocate(tmp(3, old_n))
    tmp = buf(:,1:old_n)
    deallocate(buf)
    n_alloc = 2 * n_alloc
    allocate(buf(3, n_alloc))
    buf(:,1:old_n) = tmp
    deallocate(tmp)
  end subroutine grow_points

end module mod_scalar
