! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Polygon face geometry on concave faces; with a .msh argument, reads it and expects a reader error.
program test_mesh_geometry
  use mod_kinds
  use mod_mesh, only: Mesh, compute_polygon_face_geometry
  use mod_gmsh_reader, only: gmsh_read
  implicit none

  real(rp) :: pts(3,6), centre(3), normal(3), area
  integer(ip) :: fp(6), fo(2), s, k, nfail
  character(len=512) :: arg
  type(Mesh) :: msh

  if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    call gmsh_read(trim(arg), msh)
    print *, 'test_mesh_geometry: reader accepted an open mesh'
    stop 0
  end if

  ! L-shaped planar face: area 3, centroid (5/6, 5/6), whatever vertex the fan starts from.
  pts(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]
  pts(:,2) = [2.0_rp, 0.0_rp, 0.0_rp]
  pts(:,3) = [2.0_rp, 1.0_rp, 0.0_rp]
  pts(:,4) = [1.0_rp, 1.0_rp, 0.0_rp]
  pts(:,5) = [1.0_rp, 2.0_rp, 0.0_rp]
  pts(:,6) = [0.0_rp, 2.0_rp, 0.0_rp]
  fo = [1_ip, 7_ip]
  nfail = 0
  do s = 0, 5
    do k = 1, 6
      fp(k) = mod(k - 1 + s, 6) + 1
    end do
    call compute_polygon_face_geometry(1_ip, fp, fo, pts, centre, normal, area)
    if (abs(centre(1) - 5.0_rp/6.0_rp) > 1.0e-12_rp .or. abs(centre(2) - 5.0_rp/6.0_rp) > 1.0e-12_rp &
        .or. abs(area - 3.0_rp) > 1.0e-12_rp .or. abs(normal(3) - 1.0_rp) > 1.0e-12_rp) then
      nfail = nfail + 1
      print '(A,I0,A,3ES12.4,A,ES12.4)', 'FAIL L-face start ', s, ' centre ', centre, ' area ', area
    end if
  end do

  ! Convex unit square: centre (1/2, 1/2), area 1.
  pts(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]
  pts(:,2) = [1.0_rp, 0.0_rp, 0.0_rp]
  pts(:,3) = [1.0_rp, 1.0_rp, 0.0_rp]
  pts(:,4) = [0.0_rp, 1.0_rp, 0.0_rp]
  fp(1:4) = [1_ip, 2_ip, 3_ip, 4_ip]
  fo = [1_ip, 5_ip]
  call compute_polygon_face_geometry(1_ip, fp, fo, pts, centre, normal, area)
  if (maxval(abs(centre - [0.5_rp, 0.5_rp, 0.0_rp])) > 1.0e-14_rp .or. abs(area - 1.0_rp) > 1.0e-12_rp) then
    nfail = nfail + 1
    print '(A,3ES12.4)', 'FAIL square centre ', centre
  end if

  if (nfail > 0) error stop 'test_mesh_geometry FAILED'
  print *, 'test_mesh_geometry PASSED'
end program test_mesh_geometry
