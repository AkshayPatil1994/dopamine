! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Exact V-orthogonal projection of cell velocities onto fields with zero divergence of the interpolated face flux.
module mod_cell_projection
  use mod_kinds
  use mpi_f08
  implicit none
  private

  public :: cell_proj_t
  public :: precond_iface

  abstract interface
    subroutine precond_iface(r, z)
      import :: rp
      real(rp), intent(in)    :: r(:)
      real(rp), intent(inout) :: z(:)
    end subroutine precond_iface
  end interface

  abstract interface
    subroutine halo_s_iface(x)
      import :: rp
      real(rp), intent(inout) :: x(:)
    end subroutine halo_s_iface
    subroutine halo_v_iface(x)
      import :: rp
      real(rp), intent(inout) :: x(:,:)
    end subroutine halo_v_iface
  end interface

  !> Face lists and weights of B = D I (outflow divergence of the gc-interpolated flux) and its V-adjoint.
  type :: cell_proj_t
    integer(ip) :: nc = 0, ntot = 0, nf = 0, nb = 0
    integer(ip), allocatable :: own(:), nbr(:), bcell(:)
    real(rp),    allocatable :: sf(:,:), gc(:), bsf(:,:), vol(:)
    type(MPI_Comm) :: comm
    logical :: parallel = .false.
    procedure(halo_s_iface), pointer, nopass :: halo_s => null()
    procedure(halo_v_iface), pointer, nopass :: halo_v => null()
  contains
    procedure :: init => cp_init
    procedure :: divergence => cp_divergence
    procedure :: adjoint_gradient => cp_adjoint_gradient
    procedure :: project => cp_project
  end type cell_proj_t


contains

  !> own/nbr index cells (nbr may be a ghost), sf is outward from own; bcell/bsf are outlet boundary faces.
  subroutine cp_init(self, comm, nc, ntot, vol, own, nbr, sf, gc, bcell, bsf)
    class(cell_proj_t), intent(inout) :: self
    type(MPI_Comm),     intent(in)    :: comm
    integer(ip),        intent(in)    :: nc, ntot, own(:), nbr(:), bcell(:)
    real(rp),           intent(in)    :: vol(:), sf(:,:), gc(:), bsf(:,:)
    integer(ip) :: nranks, ierr

    self%comm = comm
    call MPI_Comm_size(comm, nranks, ierr)
    self%parallel = (nranks > 1)
    self%nc = nc;  self%ntot = ntot;  self%nf = size(own);  self%nb = size(bcell)
    if (allocated(self%own)) deallocate(self%own, self%nbr, self%bcell, self%sf, self%gc, self%bsf, self%vol)
    allocate(self%own(self%nf), self%nbr(self%nf), self%bcell(self%nb), self%sf(3, self%nf), &
             self%gc(self%nf), self%bsf(3, self%nb), self%vol(nc))
    self%own = own;  self%nbr = nbr;  self%bcell = bcell
    self%sf = sf;  self%gc = gc;  self%bsf = bsf;  self%vol = vol(1:nc)
  end subroutine cp_init

  !> out = outflow sum of (gc u_own + (1-gc) u_nbr).Sf over faces plus u.Sf on outlet faces.
  subroutine cp_divergence(self, u, out)
    class(cell_proj_t), intent(in)  :: self
    real(rp),           intent(in)  :: u(:,:)
    real(rp),           intent(out) :: out(:)
    integer(ip) :: f, o, n
    real(rp)    :: fl

    out(1:self%nc) = 0.0_rp
    do f = 1, self%nf
      o = self%own(f);  n = self%nbr(f)
      fl = sum((self%gc(f) * u(:, o) + (1.0_rp - self%gc(f)) * u(:, n)) * self%sf(:, f))
      out(o) = out(o) + fl
      if (n <= self%nc) out(n) = out(n) - fl
    end do
    do f = 1, self%nb
      o = self%bcell(f)
      out(o) = out(o) + sum(u(:, o) * self%bsf(:, f))
    end do
  end subroutine cp_divergence

  !> v = V^-1 B^T chi, the V-adjoint of the divergence; zero Dirichlet chi on outlet faces.
  subroutine cp_adjoint_gradient(self, chi, v)
    class(cell_proj_t), intent(in)    :: self
    real(rp),           intent(in)    :: chi(:)
    real(rp),           intent(inout) :: v(:,:)
    integer(ip) :: f, o, n, c
    real(rp)    :: d

    v = 0.0_rp
    do f = 1, self%nf
      o = self%own(f);  n = self%nbr(f)
      d = chi(o) - chi(n)
      v(:, o) = v(:, o) + self%gc(f) * d * self%sf(:, f)
      if (n <= self%nc) v(:, n) = v(:, n) + (1.0_rp - self%gc(f)) * d * self%sf(:, f)
    end do
    do f = 1, self%nb
      o = self%bcell(f)
      v(:, o) = v(:, o) + chi(o) * self%bsf(:, f)
    end do
    do c = 1, self%nc
      v(:, c) = v(:, c) / self%vol(c)
    end do
  end subroutine cp_adjoint_gradient

  !> u <- u - V^-1 B^T chi with (B V^-1 B^T) chi = B u by flexible PCG; tol bounds the last relative velocity-energy step.
  subroutine cp_project(self, u, bconst, precond, tol, max_it, n_it)
    class(cell_proj_t), intent(in)    :: self
    real(rp),           intent(inout) :: u(:,:)
    real(rp),           intent(in)    :: bconst(:)  ! fixed boundary outflow per cell (inlets)
    procedure(precond_iface)          :: precond
    real(rp),           intent(in)    :: tol
    integer(ip),        intent(in)    :: max_it
    integer(ip),        intent(out)   :: n_it

    real(rp), allocatable :: b(:), x(:), r(:), z(:), p(:), q(:), rold(:), chi(:), v(:,:)
    real(rp) :: loc(2), glob(2), bn, rn, rz, pq, alpha, beta, en_new, en_tot
    integer(ip) :: it, nc, ntot, c, ierr

    nc = self%nc;  ntot = self%ntot
    allocate(b(ntot), x(ntot), r(ntot), z(ntot), p(ntot), q(ntot), rold(ntot), chi(ntot), v(3, ntot))
    b = 0.0_rp;  x = 0.0_rp;  r = 0.0_rp;  z = 0.0_rp;  p = 0.0_rp;  q = 0.0_rp;  rold = 0.0_rp
    call self%divergence(u, b)
    b(1:nc) = b(1:nc) + bconst(1:nc)
    r = b
    n_it = 0
    loc(1) = maxval(abs(r(1:nc)))
    call MPI_Allreduce(loc(1), bn, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)
    if (bn <= tiny(1.0_rp)) return
    call precond(r(1:nc), z(1:nc))
    p = z
    loc(1) = dot_product(r(1:nc), z(1:nc))
    call MPI_Allreduce(loc(1), rz, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
    en_tot = 0.0_rp
    do it = 1, max_it
      n_it = it
      chi = p
      if (associated(self%halo_s)) call self%halo_s(chi)
      call self%adjoint_gradient(chi, v)
      if (associated(self%halo_v)) call self%halo_v(v)
      call self%divergence(v, q)
      loc(1) = dot_product(p(1:nc), q(1:nc))
      call MPI_Allreduce(loc(1), pq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      if (abs(pq) <= tiny(1.0_rp)) exit
      alpha = rz / pq
      rold(1:nc) = r(1:nc)
      en_new = alpha * alpha * pq
      en_tot = en_tot + en_new
      x(1:nc) = x(1:nc) + alpha * p(1:nc)
      r(1:nc) = r(1:nc) - alpha * q(1:nc)
      loc(1) = maxval(abs(r(1:nc)))
      call MPI_Allreduce(loc(1), rn, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)
      if (rn <= 1.0e-12_rp * bn .or. (it >= 2 .and. en_new <= tol * tol * en_tot)) exit
      z = 0.0_rp
      call precond(r(1:nc), z(1:nc))
      loc(1) = dot_product(r(1:nc) - rold(1:nc), z(1:nc))
      loc(2) = dot_product(r(1:nc), z(1:nc))
      call MPI_Allreduce(loc, glob, 2, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      beta = max(0.0_rp, glob(1) / rz)
      rz = glob(2)
      p(1:nc) = z(1:nc) + beta * p(1:nc)
    end do

    chi = x
    if (associated(self%halo_s)) call self%halo_s(chi)
    call self%adjoint_gradient(chi, v)
    u(:, 1:nc) = u(:, 1:nc) - v(:, 1:nc)
  end subroutine cp_project

end module mod_cell_projection
