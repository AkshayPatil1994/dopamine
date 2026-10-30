! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Exact cell projection on a gmsh box (argument 1): adjointness, idempotence, V-orthogonality, divergence removal.
program test_cell_projection
  use mod_kinds
  use mod_mesh, only: Mesh
  use mod_gmsh_reader, only: gmsh_read
  use mod_cell_projection
  use mpi_f08
  implicit none

  type(Mesh) :: msh
  type(cell_proj_t) :: cp
  character(len=512) :: arg
  integer(ip) :: nc, nf, f, c, n_int, n_it, ierr, nfail
  integer(ip), allocatable :: own(:), nbr(:), bcell(:)
  real(rp), allocatable :: sf(:,:), gc(:), bsf(:,:), vol(:), u(:,:), u0(:,:), u1(:,:), v(:,:), chi(:), b(:), bconst(:)
  real(rp) :: lhs, rhs, e0, e1, orth, resid, scale

  call MPI_Init(ierr)
  nfail = 0
  call get_command_argument(1, arg)
  call gmsh_read(trim(arg), msh)
  call msh%compute_geometry()
  nc = msh%n_cells;  n_int = msh%n_int_faces;  nf = n_int
  allocate(own(nf), nbr(nf), sf(3, nf), gc(nf), vol(nc), bsf(3, 0), bconst(nc), bcell(0))
  bconst = 0.0_rp
  do f = 1, nf
    own(f) = msh%faces(f)%owner;  nbr(f) = msh%faces(f)%neighbour
    sf(:, f) = msh%faces(f)%sf;  gc(f) = msh%faces(f)%gc
  end do
  do c = 1, nc
    vol(c) = msh%cells(c)%volume
  end do
  call cp%init(MPI_COMM_WORLD, nc, nc, vol, own, nbr, sf, gc, bcell, bsf)

  allocate(u(3, nc), u0(3, nc), u1(3, nc), v(3, nc), chi(nc), b(nc))
  do c = 1, nc
    u(:, c) = [sin(3.0_rp * c), cos(1.7_rp * c), sin(0.37_rp * c * c)]
  end do

  ! Adjointness: sum_c chi_c (B u)_c = - sum_c V_c u_c . (V^-1 B^T chi)_c with the sign of the outflow convention.
  do c = 1, nc
    chi(c) = cos(0.9_rp * c)
  end do
  call cp%divergence(u, b)
  lhs = dot_product(chi, b)
  call cp%adjoint_gradient(chi, v)
  rhs = 0.0_rp
  do c = 1, nc
    rhs = rhs + vol(c) * dot_product(u(:, c), v(:, c))
  end do
  call check('adjoint identity', abs(lhs - rhs) / max(abs(lhs), 1.0_rp), 1.0e-12_rp)

  u0 = u
  call cp%project(u, bconst, identity_precond, 1.0e-6_rp, 4000_ip, n_it)
  call cp%divergence(u, b)
  scale = maxval(abs(u0)) * maxval(abs(sf))
  call check('divergence removed', maxval(abs(b)) / scale, 1.0e-6_rp)

  ! V-orthogonality: (u0 - Pu) . Pu = 0 and the energy never grows.
  orth = 0.0_rp;  e0 = 0.0_rp;  e1 = 0.0_rp
  do c = 1, nc
    orth = orth + vol(c) * dot_product(u0(:, c) - u(:, c), u(:, c))
    e0 = e0 + vol(c) * dot_product(u0(:, c), u0(:, c))
    e1 = e1 + vol(c) * dot_product(u(:, c), u(:, c))
  end do
  call check('V-orthogonality', abs(orth) / e0, 1.0e-6_rp)
  call check('energy non-increasing', max(e1 / e0 - 1.0_rp, 0.0_rp), 1.0e-12_rp)

  ! Idempotence: a second projection changes nothing.
  u1 = u
  call cp%project(u1, bconst, identity_precond, 1.0e-6_rp, 4000_ip, n_it)
  resid = 0.0_rp
  do c = 1, nc
    resid = max(resid, maxval(abs(u1(:, c) - u(:, c))))
  end do
  call check('idempotence', resid / maxval(abs(u0)), 1.0e-5_rp)

  if (nfail > 0) error stop 'test_cell_projection FAILED'
  print *, 'test_cell_projection PASSED'
  call MPI_Finalize(ierr)

contains

  subroutine identity_precond(r, z)
    real(rp), intent(in)    :: r(:)
    real(rp), intent(inout) :: z(:)
    z = r
  end subroutine identity_precond

  subroutine check(name, got, tol)
    character(len=*), intent(in) :: name
    real(rp), intent(in) :: got, tol
    if (.not. (got <= tol)) then
      nfail = nfail + 1
      print '(A,A,2ES14.6)', 'FAIL ', name, got, tol
    else
      print '(A,A,ES12.4)', 'ok   ', name, got
    end if
  end subroutine check

end program test_cell_projection
