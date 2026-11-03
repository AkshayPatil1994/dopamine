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
    logical :: staged = .false.  ! operators are on the device (single rank only)
    procedure(halo_s_iface), pointer, nopass :: halo_s => null()
    procedure(halo_v_iface), pointer, nopass :: halo_v => null()
  contains
    procedure :: init => cp_init
    procedure :: to_device => cp_to_device
    procedure :: leave_device => cp_leave_device
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

  !> Stage the operators on the device; only called from OpenACC builds.
  subroutine cp_to_device(self)
    class(cell_proj_t), intent(inout) :: self
    if (self%staged .or. .not. allocated(self%own)) return
    !$acc enter data copyin(self%own, self%nbr, self%bcell, self%sf, self%gc, self%bsf, self%vol)
    self%staged = .true.
  end subroutine cp_to_device

  subroutine cp_leave_device(self)
    class(cell_proj_t), intent(inout) :: self
    if (.not. self%staged) return
    !$acc exit data delete(self%own, self%nbr, self%bcell, self%sf, self%gc, self%bsf, self%vol)
    self%staged = .false.
  end subroutine cp_leave_device

  !> out = outflow sum of (gc u_own + (1-gc) u_nbr).Sf over faces plus u.Sf on outlet faces.
  subroutine cp_divergence(self, u, out, on_device)
    class(cell_proj_t), intent(in), target :: self
    real(rp),           intent(in)  :: u(:,:)
    real(rp),           intent(inout) :: out(:)
    logical, optional,  intent(in)    :: on_device
    integer(ip) :: f, o, n, c, nc, nf, nb
    real(rp)    :: fl, w
    logical     :: dev
    integer(ip), pointer :: own(:), nbr(:), bcell(:)
    real(rp),    pointer :: sf(:,:), gc(:), bsf(:,:)

    dev = .false.
    if (present(on_device)) dev = on_device
    nc = self%nc;  nf = self%nf;  nb = self%nb
    own => self%own;  nbr => self%nbr;  bcell => self%bcell
    sf => self%sf;  gc => self%gc;  bsf => self%bsf

    !$acc parallel loop if(dev) default(present)
    do c = 1, nc
      out(c) = 0.0_rp
    end do
    !$acc parallel loop if(dev) default(present) private(o, n, fl, w)
    do f = 1, nf
      o = own(f);  n = nbr(f);  w = gc(f)
      fl = (w * u(1, o) + (1.0_rp - w) * u(1, n)) * sf(1, f) &
         + (w * u(2, o) + (1.0_rp - w) * u(2, n)) * sf(2, f) &
         + (w * u(3, o) + (1.0_rp - w) * u(3, n)) * sf(3, f)
      !$acc atomic update
      out(o) = out(o) + fl
      if (n <= nc) then
        !$acc atomic update
        out(n) = out(n) - fl
      end if
    end do
    !$acc parallel loop if(dev) default(present) private(o, fl)
    do f = 1, nb
      o = bcell(f)
      fl = u(1, o) * bsf(1, f) + u(2, o) * bsf(2, f) + u(3, o) * bsf(3, f)
      !$acc atomic update
      out(o) = out(o) + fl
    end do
  end subroutine cp_divergence

  !> v = V^-1 B^T chi, the V-adjoint of the divergence; zero Dirichlet chi on outlet faces.
  subroutine cp_adjoint_gradient(self, chi, v, on_device)
    class(cell_proj_t), intent(in), target :: self
    real(rp),           intent(in)    :: chi(:)
    real(rp),           intent(inout) :: v(:,:)
    logical, optional,  intent(in)    :: on_device
    integer(ip) :: f, o, n, c, nc, nf, nb, ntot
    real(rp)    :: d, w, iv
    logical     :: dev
    integer(ip), pointer :: own(:), nbr(:), bcell(:)
    real(rp),    pointer :: sf(:,:), gc(:), bsf(:,:), vol(:)

    dev = .false.
    if (present(on_device)) dev = on_device
    nc = self%nc;  nf = self%nf;  nb = self%nb;  ntot = self%ntot
    own => self%own;  nbr => self%nbr;  bcell => self%bcell
    sf => self%sf;  gc => self%gc;  bsf => self%bsf;  vol => self%vol

    !$acc parallel loop if(dev) default(present)
    do c = 1, ntot
      v(1, c) = 0.0_rp;  v(2, c) = 0.0_rp;  v(3, c) = 0.0_rp
    end do
    !$acc parallel loop if(dev) default(present) private(o, n, d, w)
    do f = 1, nf
      o = own(f);  n = nbr(f);  w = gc(f)
      d = chi(o) - chi(n)
      !$acc atomic update
      v(1, o) = v(1, o) + w * d * sf(1, f)
      !$acc atomic update
      v(2, o) = v(2, o) + w * d * sf(2, f)
      !$acc atomic update
      v(3, o) = v(3, o) + w * d * sf(3, f)
      if (n <= nc) then
        !$acc atomic update
        v(1, n) = v(1, n) + (1.0_rp - w) * d * sf(1, f)
        !$acc atomic update
        v(2, n) = v(2, n) + (1.0_rp - w) * d * sf(2, f)
        !$acc atomic update
        v(3, n) = v(3, n) + (1.0_rp - w) * d * sf(3, f)
      end if
    end do
    !$acc parallel loop if(dev) default(present) private(o)
    do f = 1, nb
      o = bcell(f)
      !$acc atomic update
      v(1, o) = v(1, o) + chi(o) * bsf(1, f)
      !$acc atomic update
      v(2, o) = v(2, o) + chi(o) * bsf(2, f)
      !$acc atomic update
      v(3, o) = v(3, o) + chi(o) * bsf(3, f)
    end do
    !$acc parallel loop if(dev) default(present) private(iv)
    do c = 1, nc
      iv = 1.0_rp / vol(c)
      v(1, c) = v(1, c) * iv;  v(2, c) = v(2, c) * iv;  v(3, c) = v(3, c) * iv
    end do
  end subroutine cp_adjoint_gradient

  !> u <- u - V^-1 B^T chi with (B V^-1 B^T) chi = B u by flexible PCG; tol bounds the last relative velocity-energy step.
  !> On the device u and bconst are device-resident and precond receives device sections.
  subroutine cp_project(self, u, bconst, precond, tol, max_it, n_it, on_device)
    class(cell_proj_t), intent(in)    :: self
    real(rp),           intent(inout) :: u(:,:)
    real(rp),           intent(in)    :: bconst(:)  ! fixed boundary outflow per cell (inlets)
    procedure(precond_iface)          :: precond
    real(rp),           intent(in)    :: tol
    integer(ip),        intent(in)    :: max_it
    integer(ip),        intent(out)   :: n_it
    logical, optional,  intent(in)    :: on_device  ! u and bconst are device-resident

    real(rp), allocatable :: b(:), x(:), r(:), z(:), p(:), q(:), rold(:), chi(:), v(:,:)
    real(rp) :: loc(2), glob(2), bn, rn, rz, pq, alpha, beta, en_new, en_tot
    real(rp) :: s1, s2, smax
    integer(ip) :: it, nc, ntot, c, ierr
    logical :: dev

    dev = .false.
    if (present(on_device)) dev = on_device .and. self%staged
    nc = self%nc;  ntot = self%ntot
    allocate(b(ntot), x(ntot), r(ntot), z(ntot), p(ntot), q(ntot), rold(ntot), chi(ntot), v(3, ntot))
    !$acc enter data create(b, x, r, z, p, q, rold, chi, v) if(dev)
    !$acc parallel loop if(dev) default(present)
    do c = 1, ntot
      b(c) = 0.0_rp;  x(c) = 0.0_rp;  r(c) = 0.0_rp;  z(c) = 0.0_rp
      p(c) = 0.0_rp;  q(c) = 0.0_rp;  rold(c) = 0.0_rp;  chi(c) = 0.0_rp
    end do
    call self%divergence(u, b, dev)
    smax = 0.0_rp
    !$acc parallel loop if(dev) default(present) reduction(max:smax)
    do c = 1, nc
      b(c) = b(c) + bconst(c)
      r(c) = b(c)
      smax = max(smax, abs(r(c)))
    end do
    n_it = 0
    call MPI_Allreduce(smax, bn, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)
    if (bn > tiny(1.0_rp)) then
      call precond(r(1:nc), z(1:nc))
      s1 = 0.0_rp
      !$acc parallel loop if(dev) default(present) reduction(+:s1)
      do c = 1, nc
        p(c) = z(c)
        s1 = s1 + r(c) * z(c)
      end do
      call MPI_Allreduce(s1, rz, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
      en_tot = 0.0_rp
      do it = 1, max_it
        n_it = it
        !$acc parallel loop if(dev) default(present)
        do c = 1, ntot
          chi(c) = p(c)
        end do
        if (associated(self%halo_s)) call self%halo_s(chi)
        call self%adjoint_gradient(chi, v, dev)
        if (associated(self%halo_v)) call self%halo_v(v)
        call self%divergence(v, q, dev)
        s1 = 0.0_rp
        !$acc parallel loop if(dev) default(present) reduction(+:s1)
        do c = 1, nc
          s1 = s1 + p(c) * q(c)
        end do
        call MPI_Allreduce(s1, pq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
        if (abs(pq) <= tiny(1.0_rp)) exit
        alpha = rz / pq
        en_new = alpha * alpha * pq
        en_tot = en_tot + en_new
        smax = 0.0_rp
        !$acc parallel loop if(dev) default(present) reduction(max:smax)
        do c = 1, nc
          rold(c) = r(c)
          x(c) = x(c) + alpha * p(c)
          r(c) = r(c) - alpha * q(c)
          smax = max(smax, abs(r(c)))
        end do
        call MPI_Allreduce(smax, rn, 1, MPI_DOUBLE_PRECISION, MPI_MAX, self%comm, ierr)
        if (rn <= 1.0e-12_rp * bn .or. (it >= 2 .and. en_new <= tol * tol * en_tot)) exit
        !$acc parallel loop if(dev) default(present)
        do c = 1, ntot
          z(c) = 0.0_rp
        end do
        call precond(r(1:nc), z(1:nc))
        s1 = 0.0_rp;  s2 = 0.0_rp
        !$acc parallel loop if(dev) default(present) reduction(+:s1, s2)
        do c = 1, nc
          s1 = s1 + (r(c) - rold(c)) * z(c)
          s2 = s2 + r(c) * z(c)
        end do
        loc(1) = s1;  loc(2) = s2
        call MPI_Allreduce(loc, glob, 2, MPI_DOUBLE_PRECISION, MPI_SUM, self%comm, ierr)
        beta = max(0.0_rp, glob(1) / rz)
        rz = glob(2)
        !$acc parallel loop if(dev) default(present)
        do c = 1, nc
          p(c) = z(c) + beta * p(c)
        end do
      end do

      !$acc parallel loop if(dev) default(present)
      do c = 1, ntot
        chi(c) = x(c)
      end do
      if (associated(self%halo_s)) call self%halo_s(chi)
      call self%adjoint_gradient(chi, v, dev)
      !$acc parallel loop if(dev) default(present)
      do c = 1, nc
        u(1, c) = u(1, c) - v(1, c)
        u(2, c) = u(2, c) - v(2, c)
        u(3, c) = u(3, c) - v(3, c)
      end do
    end if
    !$acc exit data delete(b, x, r, z, p, q, rold, chi, v) if(dev)
  end subroutine cp_project

end module mod_cell_projection
