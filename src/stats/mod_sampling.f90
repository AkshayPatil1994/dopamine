! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> On-the-fly advanced sampling (planes, lines, curves, surfaces) to NetCDF.
module mod_sampling
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  use mod_fv_solver, only: FVSolver
  use mod_sampling_types
  use mpi_f08
  use netcdf
  implicit none
  private

  type :: SampleObjectRuntime
    logical            :: active = .false.
    character(len=16)  :: kind = ''
    character(len=64)  :: name = ''
    integer(ip)        :: write_frequency = 0
    integer(ip)        :: n_points = 0
    real(rp), allocatable :: points(:,:)  ! (3, n_points)
    integer(ip), allocatable :: owner_cell(:)  ! local owner cell for each point (-1 if non-owner rank)
    integer(ip), allocatable :: owner_rank(:)  ! owning rank for each point

    integer(ip)        :: ncid = -1
    integer(ip)        :: dim_time = -1
    integer(ip)        :: dim_point = -1
    integer(ip)        :: dim_xyz = -1
    integer(ip)        :: dim_vec = -1
    integer(ip)        :: var_time = -1
    integer(ip)        :: var_coords = -1
    integer(ip)        :: var_u = -1
    integer(ip)        :: var_p = -1
    integer(ip)        :: var_t = -1
    integer(ip)        :: var_nu_sgs = -1
    integer(ip)        :: samples_written = 0
  end type SampleObjectRuntime

  type, public :: SamplingCollector
    logical            :: active = .false.
    character(len=256) :: prefix = 'dopamine'
    integer(ip)        :: rank = 0
    type(MPI_Comm)     :: comm
    integer(ip)        :: n_objects = 0
    type(SampleObjectRuntime), allocatable :: objects(:)
  end type SamplingCollector

  public :: sampling_init
  public :: sampling_maybe_sample
  public :: sampling_is_due
  public :: sampling_finalize

contains

  subroutine sampling_init(sampler, prefix, msh, rank, n_planes, planes, n_lines, lines, n_curves, curves, n_surfaces, surfaces, restart_time)
    type(SamplingCollector), intent(inout) :: sampler
    character(len=*),        intent(in)    :: prefix
    type(Mesh),              intent(in)    :: msh
    integer(ip),             intent(in)    :: rank
    integer(ip),             intent(in)    :: n_planes, n_lines, n_curves, n_surfaces
    type(PlaneSampleDef),    intent(in)    :: planes(:)
    type(LineSampleDef),     intent(in)    :: lines(:)
    type(CurveSampleDef),    intent(in)    :: curves(:)
    type(SurfaceSampleDef),  intent(in)    :: surfaces(:)
    real(rp), optional,      intent(in)    :: restart_time  ! continue existing sample files, dropping later records

    integer(ip) :: total, idx, i, j

    sampler%active = .false.
    sampler%prefix = trim(prefix)
    sampler%rank = rank
    sampler%comm = msh%comm

    total = max(0_ip, n_planes) + max(0_ip, n_lines) + max(0_ip, n_curves) + max(0_ip, n_surfaces) + &
            max(0_ip, n_surfaces) * max(0_ip, n_planes)
    if (total <= 0) return

    allocate(sampler%objects(total))
    idx = 0

    do i = 1, n_planes
      if (.not. planes(i)%enabled) cycle
      if (planes(i)%write_frequency <= 0) cycle
      idx = idx + 1
      call init_plane_object(sampler%objects(idx), planes(i))
    end do

    do i = 1, n_lines
      if (.not. lines(i)%enabled) cycle
      if (lines(i)%write_frequency <= 0) cycle
      idx = idx + 1
      call init_line_object(sampler%objects(idx), lines(i))
    end do

    do i = 1, n_curves
      if (.not. curves(i)%enabled) cycle
      if (curves(i)%write_frequency <= 0) cycle
      idx = idx + 1
      call init_file_points_object(sampler%objects(idx), 'curve', curves(i)%name, curves(i)%file, curves(i)%write_frequency)
      if (.not. sampler%objects(idx)%active) idx = idx - 1
    end do

    do i = 1, n_surfaces
      if (.not. surfaces(i)%enabled) cycle
      if (surfaces(i)%write_frequency <= 0) cycle
      idx = idx + 1
      call init_file_points_object(sampler%objects(idx), 'surface', surfaces(i)%name, surfaces(i)%file, surfaces(i)%write_frequency)
      if (.not. sampler%objects(idx)%active) idx = idx - 1
    end do

    do i = 1, n_surfaces
      if (.not. surfaces(i)%enabled) cycle
      if (.not. surfaces(i)%intersect_with_planes) cycle
      do j = 1, n_planes
        if (.not. planes(j)%enabled) cycle
        idx = idx + 1
        call init_surface_plane_intersection_object(sampler%objects(idx), surfaces(i), planes(j))
        if (.not. sampler%objects(idx)%active) idx = idx - 1
      end do
    end do

    sampler%n_objects = idx
    if (sampler%n_objects <= 0) then
      deallocate(sampler%objects)
      return
    end if

    sampler%active = .true.
    do i = 1, sampler%n_objects
      call cache_object_owners(sampler%objects(i), msh, rank, sampler%comm)
    end do

    if (sampler%rank == 0) then
      do i = 1, sampler%n_objects
        call open_object_netcdf(sampler, sampler%objects(i), restart_time)
      end do
    end if
  end subroutine sampling_init

  !> Whether any active sampling object writes at this step.
  logical function sampling_is_due(sampler, step) result(due)
    type(SamplingCollector), intent(in) :: sampler
    integer(ip),             intent(in) :: step
    integer(ip) :: i

    due = .false.
    if (.not. sampler%active) return
    do i = 1, sampler%n_objects
      if (.not. sampler%objects(i)%active) cycle
      if (sampler%objects(i)%write_frequency <= 0) cycle
      if (mod(step, sampler%objects(i)%write_frequency) == 0) then
        due = .true.
        return
      end if
    end do
  end function sampling_is_due

  subroutine sampling_maybe_sample(sampler, solver, time, step)
    type(SamplingCollector), intent(inout) :: sampler
    type(FVSolver),          intent(in)    :: solver
    real(rp),                intent(in)    :: time
    integer(ip),             intent(in)    :: step

    integer(ip) :: i, npts, nc, npts_total, off, ierr
    logical :: has_due
    real(rp), allocatable :: grad_u(:,:,:), grad_p(:,:), grad_t(:,:), grad_nu(:,:), nu_ext(:)
    real(rp), allocatable :: payload_local(:,:), payload_global(:,:)

    has_due = sampling_is_due(sampler, step)
    if (.not. has_due) return

    nc = solver%mesh%n_cells
    allocate(grad_u(3,3,nc), grad_p(3,nc), grad_t(3,nc), grad_nu(3,nc))
    allocate(nu_ext(size(solver%p)))

    nu_ext = 0.0_rp
    nu_ext(1:nc) = solver%nu_sgs
    if (size(nu_ext) > nc) call solver%mesh%exchange_halos_scalar(nu_ext)

    call green_gauss_gradient_vector(solver%U, solver%mesh, grad_u, solver%U_bnd)
    call green_gauss_gradient_scalar(solver%p, solver%mesh, grad_p, solver%p_bnd)
    call green_gauss_gradient_scalar(solver%T, solver%mesh, grad_t, solver%T_bnd)
    call green_gauss_gradient_scalar(nu_ext, solver%mesh, grad_nu)

    npts_total = 0
    do i = 1, sampler%n_objects
      if (.not. sampler%objects(i)%active) cycle
      if (sampler%objects(i)%write_frequency <= 0) cycle
      if (mod(step, sampler%objects(i)%write_frequency) /= 0) cycle
      npts_total = npts_total + max(0_ip, sampler%objects(i)%n_points)
    end do

    if (npts_total > 0) then
      allocate(payload_local(6, npts_total), payload_global(6, npts_total))

      off = 0
      do i = 1, sampler%n_objects
        if (.not. sampler%objects(i)%active) cycle
        if (sampler%objects(i)%write_frequency <= 0) cycle
        if (mod(step, sampler%objects(i)%write_frequency) /= 0) cycle
        npts = sampler%objects(i)%n_points
        if (npts <= 0) cycle
        call sample_points_local(solver, sampler%objects(i), grad_u, grad_p, grad_t, grad_nu, &
                                  payload_local(:, off+1:off+npts))
        off = off + npts
      end do

      call MPI_Allreduce(payload_local, payload_global, 6*npts_total, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          solver%mesh%comm, ierr)

      off = 0
      do i = 1, sampler%n_objects
        if (.not. sampler%objects(i)%active) cycle
        if (sampler%objects(i)%write_frequency <= 0) cycle
        if (mod(step, sampler%objects(i)%write_frequency) /= 0) cycle
        npts = sampler%objects(i)%n_points
        if (npts <= 0) cycle
        if (sampler%rank == 0) then
          call write_object_sample(sampler%objects(i), time, &
                                    payload_global(1:3, off+1:off+npts), payload_global(4, off+1:off+npts), &
                                    payload_global(5, off+1:off+npts), payload_global(6, off+1:off+npts))
        end if
        off = off + npts
      end do

      deallocate(payload_local, payload_global)
    end if

    deallocate(grad_u, grad_p, grad_t, grad_nu, nu_ext)
  end subroutine sampling_maybe_sample

  subroutine sampling_finalize(sampler)
    type(SamplingCollector), intent(inout) :: sampler
    integer(ip) :: i

    if (allocated(sampler%objects)) then
      if (sampler%rank == 0) then
        do i = 1, sampler%n_objects
          if (sampler%objects(i)%ncid >= 0) call check_nc(nf90_close(sampler%objects(i)%ncid), 'nf90_close')
          sampler%objects(i)%ncid = -1
        end do
      end if
      deallocate(sampler%objects)
    end if

    sampler%active = .false.
    sampler%n_objects = 0
  end subroutine sampling_finalize

  subroutine init_plane_object(obj, plane)
    type(SampleObjectRuntime), intent(inout) :: obj
    type(PlaneSampleDef),      intent(in)    :: plane

    integer(ip) :: ix, iy, k, nx, ny
    real(rp) :: ex(3), ey(3), nvec(3), nmag, hmag, tx, ty

    obj%active = .true.
    obj%kind = 'plane'
    obj%name = trim(plane%name)
    obj%write_frequency = plane%write_frequency

    nx = max(2_ip, plane%nx)
    ny = max(2_ip, plane%ny)
    obj%n_points = nx * ny
    allocate(obj%points(3, obj%n_points))

    nvec = plane%normal
    nmag = sqrt(sum(nvec**2))
    if (nmag <= SMALL) nvec = [0.0_rp, 0.0_rp, 1.0_rp]
    nvec = nvec / max(sqrt(sum(nvec**2)), SMALL)

    ex = plane%tangent_hint - dot_product(plane%tangent_hint, nvec) * nvec
    hmag = sqrt(sum(ex**2))
    if (hmag <= SMALL) then
      if (abs(nvec(1)) < 0.9_rp) then
        ex = [1.0_rp, 0.0_rp, 0.0_rp]
      else
        ex = [0.0_rp, 1.0_rp, 0.0_rp]
      end if
      ex = ex - dot_product(ex, nvec) * nvec
      hmag = sqrt(sum(ex**2))
    end if
    ex = ex / max(hmag, SMALL)
    ey = cross3(nvec, ex)

    k = 0
    do iy = 1, ny
      ty = real(iy - 1, rp) / real(max(1_ip, ny - 1), rp)
      ty = (2.0_rp * ty - 1.0_rp) * plane%half_height
      do ix = 1, nx
        tx = real(ix - 1, rp) / real(max(1_ip, nx - 1), rp)
        tx = (2.0_rp * tx - 1.0_rp) * plane%half_width
        k = k + 1
        obj%points(:, k) = plane%origin + tx * ex + ty * ey
      end do
    end do
  end subroutine init_plane_object

  subroutine init_line_object(obj, line)
    type(SampleObjectRuntime), intent(inout) :: obj
    type(LineSampleDef),       intent(in)    :: line

    integer(ip) :: i, npts
    real(rp) :: s

    obj%active = .true.
    obj%kind = 'line'
    obj%name = trim(line%name)
    obj%write_frequency = line%write_frequency

    npts = max(2_ip, line%npts)
    obj%n_points = npts
    allocate(obj%points(3, npts))

    do i = 1, npts
      s = real(i - 1, rp) / real(max(1_ip, npts - 1), rp)
      obj%points(:, i) = (1.0_rp - s) * line%p0 + s * line%p1
    end do
  end subroutine init_line_object

  subroutine init_file_points_object(obj, kind, name, file, write_frequency)
    type(SampleObjectRuntime), intent(inout) :: obj
    character(len=*),          intent(in)    :: kind, name, file
    integer(ip),               intent(in)    :: write_frequency

    obj%active = .true.
    obj%kind = trim(kind)
    obj%name = trim(name)
    obj%write_frequency = write_frequency

    call read_points_from_file(trim(file), obj%points, obj%n_points)
    if (obj%n_points <= 0) obj%active = .false.
  end subroutine init_file_points_object

  subroutine init_surface_plane_intersection_object(obj, surface, plane)
    type(SampleObjectRuntime), intent(inout) :: obj
    type(SurfaceSampleDef),    intent(in)    :: surface
    type(PlaneSampleDef),      intent(in)    :: plane

    real(rp), allocatable :: tris(:,:,:), points(:,:)
    integer(ip) :: ntris, npts
    real(rp) :: decimate_spacing

    obj%active = .false.
    obj%kind = 'curve'
    obj%name = trim(surface%name)//'_x_'//trim(plane%name)
    obj%write_frequency = max(1_ip, surface%write_frequency)

    call read_surface_triangles_from_file(trim(surface%file), tris, ntris)
    if (ntris <= 0) return

    call intersect_surface_with_plane(tris, ntris, plane%origin, plane%normal, points, npts)
    if (allocated(tris)) deallocate(tris)
    if (npts <= 1) then
      if (allocated(points)) deallocate(points)
      return
    end if

    if (surface%sort_intersection_points) call sort_points_greedy(points)

    decimate_spacing = surface%intersection_decimate_spacing
    if (decimate_spacing > 0.0_rp) then
      call decimate_points(points, decimate_spacing)
    end if

    npts = size(points, 2)
    if (npts <= 1) then
      if (allocated(points)) deallocate(points)
      return
    end if

    obj%n_points = npts
    allocate(obj%points(3, npts))
    obj%points = points
    obj%active = .true.
    deallocate(points)
  end subroutine init_surface_plane_intersection_object

  subroutine open_object_netcdf(sampler, obj, restart_time)
    type(SamplingCollector),   intent(in)    :: sampler
    type(SampleObjectRuntime), intent(inout) :: obj
    real(rp), optional,        intent(in)    :: restart_time

    character(len=512) :: fname

    if (.not. obj%active) return
    fname = trim(sampler%prefix)//'_sample_'//trim(obj%kind)//'_'//trim(obj%name)//'.nc'

    if (present(restart_time)) then
      if (restart_time >= 0.0_rp) then  ! negative = fresh run
        if (reopen_object_netcdf(obj, trim(fname), restart_time)) return
      end if
    end if

    call check_nc(nf90_create(trim(fname), nf90_clobber, obj%ncid), 'nf90_create('//trim(fname)//')')
    call check_nc(nf90_def_dim(obj%ncid, 'time', nf90_unlimited, obj%dim_time), 'nf90_def_dim(time)')
    call check_nc(nf90_def_dim(obj%ncid, 'point', obj%n_points, obj%dim_point), 'nf90_def_dim(point)')
    call check_nc(nf90_def_dim(obj%ncid, 'xyz', 3, obj%dim_xyz), 'nf90_def_dim(xyz)')
    call check_nc(nf90_def_dim(obj%ncid, 'vec', 3, obj%dim_vec), 'nf90_def_dim(vec)')

    call check_nc(nf90_def_var(obj%ncid, 'time', nf90_double, (/obj%dim_time/), obj%var_time), 'nf90_def_var(time)')
    call check_nc(nf90_def_var(obj%ncid, 'coords', nf90_double, (/obj%dim_xyz, obj%dim_point/), obj%var_coords), 'nf90_def_var(coords)')
    call check_nc(nf90_def_var(obj%ncid, 'U', nf90_double, (/obj%dim_vec, obj%dim_point, obj%dim_time/), obj%var_u), 'nf90_def_var(U)')
    call check_nc(nf90_def_var(obj%ncid, 'p', nf90_double, (/obj%dim_point, obj%dim_time/), obj%var_p), 'nf90_def_var(p)')
    call check_nc(nf90_def_var(obj%ncid, 'T', nf90_double, (/obj%dim_point, obj%dim_time/), obj%var_t), 'nf90_def_var(T)')
    call check_nc(nf90_def_var(obj%ncid, 'nu_sgs', nf90_double, (/obj%dim_point, obj%dim_time/), obj%var_nu_sgs), 'nf90_def_var(nu_sgs)')

    call check_nc(nf90_put_att(obj%ncid, nf90_global, 'sampler_kind', trim(obj%kind)), 'nf90_put_att(kind)')
    call check_nc(nf90_put_att(obj%ncid, nf90_global, 'sampler_name', trim(obj%name)), 'nf90_put_att(name)')

    call check_nc(nf90_enddef(obj%ncid), 'nf90_enddef')
    call check_nc(nf90_put_var(obj%ncid, obj%var_coords, obj%points), 'nf90_put_var(coords)')
  end subroutine open_object_netcdf

  !> Continue an existing sample file after a restart: keep the records up to restart_time and let the
  !> next samples overwrite the rest. Returns .false. (caller creates a fresh file) when the file is
  !> missing or was written for a different point set.
  logical function reopen_object_netcdf(obj, fname, restart_time) result(ok)
    type(SampleObjectRuntime), intent(inout) :: obj
    character(len=*),          intent(in)    :: fname
    real(rp),                  intent(in)    :: restart_time

    integer(ip) :: ncid, dim_time, dim_point, n_rec, n_pts, keep, stat
    real(rp), allocatable :: times(:)
    logical :: exists

    ok = .false.
    inquire(file=fname, exist=exists)
    if (.not. exists) return
    if (nf90_open(fname, nf90_write, ncid) /= nf90_noerr) return

    stat = nf90_inq_dimid(ncid, 'time', dim_time)
    if (stat == nf90_noerr) stat = nf90_inq_dimid(ncid, 'point', dim_point)
    if (stat == nf90_noerr) stat = nf90_inquire_dimension(ncid, dim_point, len=n_pts)
    if (stat == nf90_noerr) stat = nf90_inquire_dimension(ncid, dim_time, len=n_rec)
    if (stat == nf90_noerr) stat = nf90_inq_varid(ncid, 'time', obj%var_time)
    if (stat == nf90_noerr) stat = nf90_inq_varid(ncid, 'U', obj%var_u)
    if (stat == nf90_noerr) stat = nf90_inq_varid(ncid, 'p', obj%var_p)
    if (stat == nf90_noerr) stat = nf90_inq_varid(ncid, 'T', obj%var_t)
    if (stat == nf90_noerr) stat = nf90_inq_varid(ncid, 'nu_sgs', obj%var_nu_sgs)
    if (stat /= nf90_noerr .or. n_pts /= obj%n_points) then
      write(*, '(A)') ' Sampling: cannot continue '//fname//' (different layout); starting a new file.'
      stat = nf90_close(ncid)
      return
    end if

    keep = n_rec
    if (n_rec > 0) then
      allocate(times(n_rec))
      call check_nc(nf90_get_var(ncid, obj%var_time, times, start=(/1/), count=(/n_rec/)), 'nf90_get_var(time)')
      ! Records are in time order; the tolerance absorbs the sample written exactly at the checkpoint.
      do while (keep > 0)
        if (times(keep) <= restart_time + 1.0e-9_rp * max(abs(restart_time), 1.0_rp)) exit
        keep = keep - 1
      end do
      deallocate(times)
    end if
    if (keep < n_rec) write(*, '(A,I0,A)') ' Sampling: dropping ', n_rec - keep, ' record(s) after the restart in '//fname

    obj%ncid = ncid
    obj%dim_time = dim_time
    obj%dim_point = dim_point
    obj%samples_written = keep
    ok = .true.
  end function reopen_object_netcdf

  subroutine write_object_sample(obj, time, u, p, t, nu_sgs)
    type(SampleObjectRuntime), intent(inout) :: obj
    real(rp),                  intent(in)    :: time
    real(rp),                  intent(in)    :: u(:,:), p(:), t(:), nu_sgs(:)

    integer(ip) :: tidx
    real(rp) :: time_buf(1)

    tidx = obj%samples_written + 1
    time_buf(1) = time

    call check_nc(nf90_put_var(obj%ncid, obj%var_time, time_buf, start=(/tidx/), count=(/1/)), 'nf90_put_var(time)')
    call check_nc(nf90_put_var(obj%ncid, obj%var_u, u, start=(/1, 1, tidx/), count=(/3, obj%n_points, 1/)), 'nf90_put_var(U)')
    call check_nc(nf90_put_var(obj%ncid, obj%var_p, p, start=(/1, tidx/), count=(/obj%n_points, 1/)), 'nf90_put_var(p)')
    call check_nc(nf90_put_var(obj%ncid, obj%var_t, t, start=(/1, tidx/), count=(/obj%n_points, 1/)), 'nf90_put_var(T)')
    call check_nc(nf90_put_var(obj%ncid, obj%var_nu_sgs, nu_sgs, start=(/1, tidx/), count=(/obj%n_points, 1/)), 'nf90_put_var(nu_sgs)')

    obj%samples_written = tidx
  end subroutine write_object_sample

  !> Local (pre-reduction) contribution for one sample object; caller fuses the MPI_Allreduce across objects.
  subroutine sample_points_local(solver, obj, grad_u, grad_p, grad_t, grad_nu, payload_local)
    type(FVSolver),            intent(in)  :: solver
    type(SampleObjectRuntime), intent(in)  :: obj
    real(rp),                  intent(in)  :: grad_u(:,:,:), grad_p(:,:), grad_t(:,:), grad_nu(:,:)
    real(rp),                  intent(out) :: payload_local(:,:)

    integer(ip) :: i, c, npts
    real(rp) :: dx(3)

    npts = obj%n_points
    payload_local(:, 1:npts) = 0.0_rp

    do i = 1, npts
      c = -1
      if (allocated(obj%owner_cell)) c = obj%owner_cell(i)
      if (c >= 1 .and. c <= solver%mesh%n_cells) then
        dx = obj%points(:, i) - solver%mesh%cells(c)%centroid
        payload_local(1, i) = solver%U(1, c) + dot_product(grad_u(:,1,c), dx)
        payload_local(2, i) = solver%U(2, c) + dot_product(grad_u(:,2,c), dx)
        payload_local(3, i) = solver%U(3, c) + dot_product(grad_u(:,3,c), dx)
        payload_local(4, i) = solver%p(c) + dot_product(grad_p(:,c), dx)
        payload_local(5, i) = solver%T(c) + dot_product(grad_t(:,c), dx)
        payload_local(6, i) = solver%nu_sgs(c) + dot_product(grad_nu(:,c), dx)
      end if
    end do
  end subroutine sample_points_local

  subroutine cache_object_owners(obj, msh, rank, comm)
    type(SampleObjectRuntime), intent(inout) :: obj
    type(Mesh),                intent(in)    :: msh
    integer(ip),               intent(in)    :: rank
    type(MPI_Comm),            intent(in)    :: comm

    integer(ip) :: i, c, best_cell, ierr, rank_local, rank_global
    real(rp) :: d2, d2_local, d2_global

    if (.not. obj%active) return
    if (obj%n_points <= 0) return

    if (allocated(obj%owner_cell)) deallocate(obj%owner_cell)
    if (allocated(obj%owner_rank)) deallocate(obj%owner_rank)
    allocate(obj%owner_cell(obj%n_points), obj%owner_rank(obj%n_points))
    obj%owner_cell = -1
    obj%owner_rank = -1

    do i = 1, obj%n_points
      d2_local = huge(1.0_rp)
      best_cell = -1
      do c = 1, msh%n_cells
        d2 = sum((msh%cells(c)%centroid - obj%points(:, i))**2)
        if (d2 < d2_local) then
          d2_local = d2
          best_cell = c
        end if
      end do

      call MPI_Allreduce(d2_local, d2_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)

      rank_local = huge(1_ip)
      if (best_cell >= 1) then
        if (abs(d2_local - d2_global) <= 1.0e-12_rp * (1.0_rp + d2_global)) rank_local = rank
      end if
      call MPI_Allreduce(rank_local, rank_global, 1, MPI_INTEGER, MPI_MIN, comm, ierr)

      obj%owner_rank(i) = rank_global
      if (rank == rank_global) then
        obj%owner_cell(i) = best_cell
      end if
    end do
  end subroutine cache_object_owners

  subroutine read_points_from_file(filename, points, n_points)
    character(len=*),            intent(in)  :: filename
    real(rp), allocatable,       intent(out) :: points(:,:)
    integer(ip),                 intent(out) :: n_points

    integer(ip) :: unit_id, ios, n_alloc
    character(len=1024) :: line
    real(rp), allocatable :: tmp(:,:)
    real(rp) :: x, y, z

    n_points = 0
    n_alloc = 512
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
      if (n_points > n_alloc) call grow_point_buffer(tmp, n_alloc)
      tmp(:, n_points) = [x, y, z]
    end do
    close(unit_id)

    allocate(points(3, n_points))
    if (n_points > 0) points = tmp(:, 1:n_points)
    deallocate(tmp)
  end subroutine read_points_from_file

  subroutine read_surface_triangles_from_file(filename, tris, n_tris)
    character(len=*),            intent(in)  :: filename
    real(rp), allocatable,       intent(out) :: tris(:,:,:)
    integer(ip),                 intent(out) :: n_tris

    character(len=16) :: ext

    ext = lower_ext(filename)
    if (trim(ext) == 'obj') then
      call read_obj_triangles(filename, tris, n_tris)
    else if (trim(ext) == 'stl') then
      call read_stl_triangles(filename, tris, n_tris)
    else
      n_tris = 0
      allocate(tris(3,3,0))
    end if
  end subroutine read_surface_triangles_from_file

  subroutine read_obj_triangles(filename, tris, n_tris)
    character(len=*),            intent(in)  :: filename
    real(rp), allocatable,       intent(out) :: tris(:,:,:)
    integer(ip),                 intent(out) :: n_tris

    integer(ip) :: unit_id, ios, nv_alloc, nt_alloc, nv, nt, nidx, i
    character(len=1024) :: line
    real(rp), allocatable :: verts(:,:)
    integer(ip), allocatable :: tri_conn(:,:)
    integer(ip) :: idx(64)

    nv_alloc = 1024
    nt_alloc = 1024
    nv = 0
    nt = 0
    allocate(verts(3, nv_alloc))
    allocate(tri_conn(3, nt_alloc))

    open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      n_tris = 0
      allocate(tris(3,3,0))
      deallocate(verts, tri_conn)
      return
    end if

    do
      read(unit_id, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(trim(line))
      if (len_trim(line) < 2) cycle
      if (line(1:2) == 'v ') then
        nv = nv + 1
        if (nv > nv_alloc) call grow_points_xyz(verts, nv_alloc)
        read(line(3:), *, iostat=ios) verts(1, nv), verts(2, nv), verts(3, nv)
        if (ios /= 0) nv = nv - 1
      else if (line(1:2) == 'f ') then
        call parse_obj_face_line(line(3:), idx, nidx)
        if (nidx >= 3) then
          do i = 2, nidx - 1
            nt = nt + 1
            if (nt > nt_alloc) call grow_int3(tri_conn, nt_alloc)
            tri_conn(:, nt) = [idx(1), idx(i), idx(i+1)]
          end do
        end if
      end if
    end do
    close(unit_id)

    if (nt <= 0) then
      n_tris = 0
      allocate(tris(3,3,0))
      deallocate(verts, tri_conn)
      return
    end if

    allocate(tris(3,3,nt))
    n_tris = 0
    do i = 1, nt
      if (any(tri_conn(:, i) < 1) .or. any(tri_conn(:, i) > nv)) cycle
      n_tris = n_tris + 1
      tris(:,1,n_tris) = verts(:, tri_conn(1,i))
      tris(:,2,n_tris) = verts(:, tri_conn(2,i))
      tris(:,3,n_tris) = verts(:, tri_conn(3,i))
    end do

    if (n_tris < nt) then
      call shrink_tris(tris, n_tris)
    end if

    deallocate(verts, tri_conn)
  end subroutine read_obj_triangles

  subroutine read_stl_triangles(filename, tris, n_tris)
    character(len=*),            intent(in)  :: filename
    real(rp), allocatable,       intent(out) :: tris(:,:,:)
    integer(ip),                 intent(out) :: n_tris

    integer(ip) :: unit_id, ios, nt_alloc, nt, vcount
    character(len=1024) :: line
    real(rp), allocatable :: tmp(:,:,:)

    nt_alloc = 1024
    nt = 0
    vcount = 0
    allocate(tmp(3,3,nt_alloc))

    open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      n_tris = 0
      allocate(tris(3,3,0))
      deallocate(tmp)
      return
    end if

    do
      read(unit_id, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(trim(line))
      if (index(line, 'vertex ') == 1) then
        if (vcount == 0) then
          nt = nt + 1
          if (nt > nt_alloc) call grow_triangles(tmp, nt_alloc)
        end if
        vcount = vcount + 1
        read(line(8:), *, iostat=ios) tmp(1, vcount, nt), tmp(2, vcount, nt), tmp(3, vcount, nt)
        if (vcount == 3) vcount = 0
      end if
    end do
    close(unit_id)

    n_tris = nt
    allocate(tris(3,3,n_tris))
    if (n_tris > 0) tris = tmp(:,:,1:n_tris)
    deallocate(tmp)
  end subroutine read_stl_triangles

  subroutine intersect_surface_with_plane(tris, n_tris, origin, normal, points, n_points)
    real(rp),             intent(in)  :: tris(:,:,:)
    integer(ip),          intent(in)  :: n_tris
    real(rp),             intent(in)  :: origin(3), normal(3)
    real(rp), allocatable,intent(out) :: points(:,:)
    integer(ip),          intent(out) :: n_points

    integer(ip) :: i, n_alloc, n_used
    real(rp), allocatable :: tmp(:,:)
    logical :: has_seg
    real(rp) :: p1(3), p2(3)

    n_alloc = max(1_ip, 2_ip * n_tris)
    n_used = 0
    allocate(tmp(3, n_alloc))

    do i = 1, n_tris
      call triangle_plane_segment(tris(:,1,i), tris(:,2,i), tris(:,3,i), origin, normal, has_seg, p1, p2)
      if (.not. has_seg) cycle
      n_used = n_used + 1
      if (n_used > n_alloc) call grow_points_xyz(tmp, n_alloc)
      tmp(:, n_used) = p1
      n_used = n_used + 1
      if (n_used > n_alloc) call grow_points_xyz(tmp, n_alloc)
      tmp(:, n_used) = p2
    end do

    if (n_used <= 0) then
      n_points = 0
      allocate(points(3,0))
      deallocate(tmp)
      return
    end if

    call deduplicate_points(tmp(:,1:n_used), 1.0e-9_rp, points, n_points)
    deallocate(tmp)
  end subroutine intersect_surface_with_plane

  subroutine triangle_plane_segment(a, b, c, origin, normal, has_seg, p1, p2)
    real(rp),    intent(in)  :: a(3), b(3), c(3), origin(3), normal(3)
    logical,     intent(out) :: has_seg
    real(rp),    intent(out) :: p1(3), p2(3)

    real(rp), parameter :: eps = 1.0e-12_rp
    real(rp) :: d(3)
    real(rp) :: q(3,6)
    integer(ip) :: n, i

    d(1) = dot_product(a - origin, normal)
    d(2) = dot_product(b - origin, normal)
    d(3) = dot_product(c - origin, normal)

    n = 0
    call edge_intersection(a, b, d(1), d(2), eps, q, n)
    call edge_intersection(b, c, d(2), d(3), eps, q, n)
    call edge_intersection(c, a, d(3), d(1), eps, q, n)

    call unique_points_inplace(q, n, eps)

    if (n < 2) then
      has_seg = .false.
      p1 = 0.0_rp
      p2 = 0.0_rp
      return
    end if

    has_seg = .true.
    p1 = q(:,1)
    p2 = q(:,2)

    do i = 3, n
      if (sum((q(:,i) - p1)**2) > sum((p2 - p1)**2)) p2 = q(:,i)
    end do
  end subroutine triangle_plane_segment

  subroutine edge_intersection(v1, v2, d1, d2, eps, q, n)
    real(rp),    intent(in)    :: v1(3), v2(3), d1, d2, eps
    real(rp),    intent(inout) :: q(:,:)
    integer(ip), intent(inout) :: n
    real(rp) :: t

    if (abs(d1) <= eps .and. abs(d2) <= eps) return

    if (d1 * d2 < -eps*eps) then
      t = d1 / (d1 - d2)
      n = n + 1
      q(:, n) = v1 + t * (v2 - v1)
    else if (abs(d1) <= eps) then
      n = n + 1
      q(:, n) = v1
    else if (abs(d2) <= eps) then
      n = n + 1
      q(:, n) = v2
    end if
  end subroutine edge_intersection

  subroutine unique_points_inplace(points, n, tol)
    real(rp),    intent(inout) :: points(:,:)
    integer(ip), intent(inout) :: n
    real(rp),    intent(in)    :: tol

    integer(ip) :: i, j, k
    logical :: is_dup

    k = 0
    do i = 1, n
      is_dup = .false.
      do j = 1, k
        if (sum((points(:,i) - points(:,j))**2) <= tol*tol) then
          is_dup = .true.
          exit
        end if
      end do
      if (.not. is_dup) then
        k = k + 1
        points(:,k) = points(:,i)
      end if
    end do
    n = k
  end subroutine unique_points_inplace

  subroutine deduplicate_points(points_in, tol, points_out, n_out)
    real(rp),             intent(in)  :: points_in(:,:)
    real(rp),             intent(in)  :: tol
    real(rp), allocatable,intent(out) :: points_out(:,:)
    integer(ip),          intent(out) :: n_out

    integer(ip) :: i, j, k, nin
    logical :: is_dup

    nin = size(points_in, 2)
    if (nin <= 0) then
      n_out = 0
      allocate(points_out(3,0))
      return
    end if

    allocate(points_out(3, nin))
    k = 0
    do i = 1, nin
      is_dup = .false.
      do j = 1, k
        if (sum((points_in(:,i) - points_out(:,j))**2) <= tol*tol) then
          is_dup = .true.
          exit
        end if
      end do
      if (.not. is_dup) then
        k = k + 1
        points_out(:,k) = points_in(:,i)
      end if
    end do

    n_out = k
    if (k < nin) then
      points_out = points_out(:,1:k)
    end if
  end subroutine deduplicate_points

  subroutine sort_points_greedy(points)
    real(rp), intent(inout), allocatable :: points(:,:)

    integer(ip) :: n, i, j, current, next_i, start_i
    logical, allocatable :: used(:)
    real(rp), allocatable :: sorted(:,:)
    real(rp) :: d2, best_d2, key, best_key

    n = size(points, 2)
    if (n <= 2) return

    allocate(used(n), sorted(3,n))
    used = .false.

    start_i = 1
    best_key = huge(1.0_rp)
    do i = 1, n
      key = points(1,i) + points(2,i) + points(3,i)
      if (key < best_key) then
        best_key = key
        start_i = i
      end if
    end do

    current = start_i
    do i = 1, n
      sorted(:,i) = points(:,current)
      used(current) = .true.
      if (i == n) exit
      best_d2 = huge(1.0_rp)
      next_i = -1
      do j = 1, n
        if (used(j)) cycle
        d2 = sum((points(:,j) - points(:,current))**2)
        if (d2 < best_d2) then
          best_d2 = d2
          next_i = j
        end if
      end do
      if (next_i < 1) exit
      current = next_i
    end do

    points = sorted
    deallocate(used, sorted)
  end subroutine sort_points_greedy

  subroutine decimate_points(points, spacing)
    real(rp), intent(inout), allocatable :: points(:,:)
    real(rp), intent(in) :: spacing

    integer(ip) :: n, i, k
    real(rp), allocatable :: out(:,:)

    n = size(points, 2)
    if (n <= 2 .or. spacing <= 0.0_rp) return

    allocate(out(3, n))
    k = 1
    out(:,k) = points(:,1)

    do i = 2, n - 1
      if (sqrt(sum((points(:,i) - out(:,k))**2)) >= spacing) then
        k = k + 1
        out(:,k) = points(:,i)
      end if
    end do

    if (sqrt(sum((points(:,n) - out(:,k))**2)) > 0.0_rp) then
      k = k + 1
      out(:,k) = points(:,n)
    end if

    points = out(:,1:k)
    deallocate(out)
  end subroutine decimate_points

  subroutine parse_obj_face_line(s, idx, nidx)
    character(len=*), intent(in)  :: s
    integer(ip),      intent(out) :: idx(:)
    integer(ip),      intent(out) :: nidx

    integer(ip) :: pos, n, ios
    character(len=256) :: token, t

    nidx = 0
    t = trim(adjustl(s))
    pos = 1

    do
      call next_token(t, pos, token, n)
      if (n <= 0) exit
      nidx = nidx + 1
      if (nidx > size(idx)) then
        nidx = nidx - 1
        exit
      end if
      call parse_obj_index_token(token, idx(nidx), ios)
      if (ios /= 0) nidx = nidx - 1
    end do
  end subroutine parse_obj_face_line

  subroutine next_token(s, pos, token, n)
    character(len=*), intent(in)    :: s
    integer(ip),      intent(inout) :: pos
    character(len=*), intent(out)   :: token
    integer(ip),      intent(out)   :: n

    integer(ip) :: i, j, ls

    token = ''
    n = 0
    ls = len_trim(s)
    if (pos > ls) return

    i = pos
    do while (i <= ls .and. s(i:i) == ' ')
      i = i + 1
    end do
    if (i > ls) then
      pos = ls + 1
      return
    end if

    j = i
    do while (j <= ls .and. s(j:j) /= ' ')
      j = j + 1
    end do

    n = j - i
    if (n > 0) token(1:n) = s(i:j-1)
    pos = j + 1
  end subroutine next_token

  subroutine parse_obj_index_token(token, idx, ios)
    character(len=*), intent(in)  :: token
    integer(ip),      intent(out) :: idx
    integer(ip),      intent(out) :: ios

    integer(ip) :: slash
    character(len=64) :: base

    base = ''
    slash = index(trim(token), '/')
    if (slash > 0) then
      base = token(1:slash-1)
    else
      base = token
    end if

    read(base, *, iostat=ios) idx
  end subroutine parse_obj_index_token

  function lower_ext(filename) result(ext)
    character(len=*), intent(in) :: filename
    character(len=16) :: ext
    character(len=512) :: s
    integer(ip) :: i, dot

    ext = ''
    s = trim(filename)
    dot = 0
    do i = len_trim(s), 1, -1
      if (s(i:i) == '.') then
        dot = i
        exit
      end if
    end do
    if (dot <= 0 .or. dot >= len_trim(s)) return

    ext = s(dot+1:len_trim(s))
    do i = 1, len_trim(ext)
      if (iachar(ext(i:i)) >= iachar('A') .and. iachar(ext(i:i)) <= iachar('Z')) then
        ext(i:i) = achar(iachar(ext(i:i)) + 32)
      end if
    end do
  end function lower_ext

  subroutine grow_points_xyz(buf, n_alloc)
    real(rp), allocatable, intent(inout) :: buf(:,:)
    integer(ip),           intent(inout) :: n_alloc

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
  end subroutine grow_points_xyz

  subroutine grow_triangles(buf, n_alloc)
    real(rp), allocatable, intent(inout) :: buf(:,:,:)
    integer(ip),           intent(inout) :: n_alloc

    real(rp), allocatable :: tmp(:,:,:)
    integer(ip) :: old_n

    old_n = n_alloc
    allocate(tmp(3,3,old_n))
    tmp = buf(:,:,1:old_n)
    deallocate(buf)
    n_alloc = 2 * n_alloc
    allocate(buf(3,3,n_alloc))
    buf(:,:,1:old_n) = tmp
    deallocate(tmp)
  end subroutine grow_triangles

  subroutine grow_int3(buf, n_alloc)
    integer(ip), allocatable, intent(inout) :: buf(:,:)
    integer(ip),              intent(inout) :: n_alloc

    integer(ip), allocatable :: tmp(:,:)
    integer(ip) :: old_n

    old_n = n_alloc
    allocate(tmp(3, old_n))
    tmp = buf(:,1:old_n)
    deallocate(buf)
    n_alloc = 2 * n_alloc
    allocate(buf(3, n_alloc))
    buf(:,1:old_n) = tmp
    deallocate(tmp)
  end subroutine grow_int3

  subroutine shrink_tris(tris, n_keep)
    real(rp), allocatable, intent(inout) :: tris(:,:,:)
    integer(ip), intent(in) :: n_keep

    real(rp), allocatable :: tmp(:,:,:)

    allocate(tmp(3,3,n_keep))
    if (n_keep > 0) tmp = tris(:,:,1:n_keep)
    deallocate(tris)
    allocate(tris(3,3,n_keep))
    if (n_keep > 0) tris = tmp
    deallocate(tmp)
  end subroutine shrink_tris

  subroutine grow_point_buffer(buf, n_alloc)
    real(rp), allocatable, intent(inout) :: buf(:,:)
    integer(ip),           intent(inout) :: n_alloc

    call grow_points_xyz(buf, n_alloc)
  end subroutine grow_point_buffer

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

  pure function cross3(a, b) result(c)
    real(rp), intent(in) :: a(3), b(3)
    real(rp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  subroutine check_nc(status, where)
    integer(ip),       intent(in) :: status
    character(len=*),  intent(in) :: where

    if (status /= nf90_noerr) then
      write(*, '(A,1X,A,1X,A)') 'ERROR(NetCDF):', trim(where), trim(nf90_strerror(status))
      error stop 'NetCDF failure in sampling module'
    end if
  end subroutine check_nc

end module mod_sampling
