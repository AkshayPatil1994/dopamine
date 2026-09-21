! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Weather data I/O and Davies nudging for large-scale atmospheric coupling
module mod_abl_weather
  use mod_kinds
  use mod_abl_types
  use netcdf
  implicit none
  private

  real(rp), parameter :: SMALL = 1.0e-16_rp

  public :: nudging_load, nudging_apply_rhs

contains

  !> Load weather profile data from file (netCDF or CSV)
  !> CSV format expected: one row per (time, level) pair
  subroutine nudging_load(data, filename)
    type(NudgingData), intent(inout) :: data
    character(len=*),  intent(in)    :: filename

    data%loaded = .false.
    if (len_trim(filename) == 0) then
      write(*, '(A)') 'ERROR: nudging_file is empty'
      return
    end if

    ! Detect file type by extension; each loader sets data%loaded only on success.
    if (index(filename, '.nc') > 0) then
      call nudging_load_netcdf(data, filename)
    else
      call nudging_load_csv(data, filename)
    end if

  end subroutine nudging_load

  !> Load from CSV file (simple implementation)
  subroutine nudging_load_csv(data, filename)
    type(NudgingData), intent(inout) :: data
    character(len=*),  intent(in)    :: filename

    integer(ip) :: nt, nz, it, iz, n_rows, k
    real(rp), allocatable :: row(:,:)  ! (6, n_rows): t, z, u, v, w, T

    call csv_read_rows(filename, row, n_rows)
    if (n_rows <= 0) return

    ! Levels per time block = leading rows sharing the first time value.
    nz = 0
    do k = 1, n_rows
      if (abs(row(1, k) - row(1, 1)) > 1.0e-9_rp * max(abs(row(1, 1)), 1.0_rp)) exit
      nz = nz + 1
    end do
    if (nz < 1) nz = 1
    nt = n_rows / nz

    if (nt < 1 .or. nt * nz /= n_rows) then
      write(*, '(A,I0,A,I0,A,I0)') 'ERROR: weather file is not a rectangular grid: ', &
        n_rows, ' rows, ', nz, ' levels, remainder ', n_rows - nt * nz
      return
    end if

    data%n_times  = nt
    data%n_levels = nz
    allocate(data%times(nt), source=0.0_rp)
    allocate(data%z_levels(nz), source=0.0_rp)
    allocate(data%U_profile(3, nz, nt), source=0.0_rp)
    allocate(data%T_profile(nz, nt), source=0.0_rp)

    do k = 1, n_rows
      it = (k - 1) / nz + 1
      iz = mod(k - 1, nz) + 1
      data%times(it)          = row(1, k)
      data%z_levels(iz)       = row(2, k)
      data%U_profile(:, iz, it) = row(3:5, k)
      data%T_profile(iz, it)  = row(6, k)
    end do

    deallocate(row)
    data%loaded = .true.
  end subroutine nudging_load_csv

  !> Read all non-comment rows of a weather CSV into row(6, n_rows); n_rows = 0 on failure.
  subroutine csv_read_rows(filename, row, n_rows)
    character(len=*),      intent(in)  :: filename
    real(rp), allocatable, intent(out) :: row(:,:)
    integer(ip),           intent(out) :: n_rows

    integer(ip) :: iunit, ios, k
    character(len=512) :: line
    real(rp) :: vals(6)

    n_rows = 0
    open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') 'ERROR: cannot open weather file: ' // trim(filename)
      return
    end if

    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) vals
      if (ios == 0) n_rows = n_rows + 1
    end do

    if (n_rows == 0) then
      write(*, '(A)') 'ERROR: no valid data in weather file'
      close(iunit)
      return
    end if

    allocate(row(6, n_rows), source=0.0_rp)
    rewind(iunit)
    k = 0
    do
      read(iunit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) vals
      if (ios /= 0) cycle
      k = k + 1
      if (k > n_rows) exit
      row(:, k) = vals
    end do
    close(iunit)
  end subroutine csv_read_rows

  !> Load from netCDF file
  !> Expects file to contain: time(:), z(:), U(time,z), V(time,z), W(time,z), T(time,z)
  subroutine nudging_load_netcdf(data, filename)
    type(NudgingData), intent(inout) :: data
    character(len=*),  intent(in)    :: filename

    integer(ip) :: ncid, status, dimid_time, dimid_z
    integer(ip) :: nt, nz, stat_close
    real(rp), allocatable :: times_tmp(:), z_tmp(:), U_tmp(:,:,:), T_tmp(:,:)

    ! Attempt to open netCDF file
    status = nf90_open(trim(filename), nf90_nowrite, ncid)
    if (status /= nf90_noerr) then
      write(*, '(A)') 'WARNING: Cannot open as netCDF: ' // trim(filename)
      write(*, '(A)') '         Attempting to read as CSV format instead...'
      call nudging_load_csv(data, filename)
      if (.not. data%loaded) then
        write(*, '(A)') 'ERROR: Unable to load weather file in either format'
      end if
      return
    end if

    ! Inquire dimension IDs
    status = nf90_inq_dimid(ncid, 'time', dimid_time)
    if (status /= nf90_noerr) then
      write(*, '(A)') 'ERROR: netCDF file missing "time" dimension'
      stat_close = nf90_close(ncid)
      return
    end if

    status = nf90_inq_dimid(ncid, 'z', dimid_z)
    if (status /= nf90_noerr) then
      write(*, '(A)') 'ERROR: netCDF file missing "z" dimension'
      stat_close = nf90_close(ncid)
      return
    end if

    ! Get dimension sizes
    status = nf90_inquire_dimension(ncid, dimid_time, len=nt)
    if (status /= nf90_noerr) nt = 0

    status = nf90_inquire_dimension(ncid, dimid_z, len=nz)
    if (status /= nf90_noerr) nz = 0

    if (nt < 1 .or. nz < 1) then
      write(*, '(A)') 'ERROR: Invalid netCDF dimensions'
      stat_close = nf90_close(ncid)
      return
    end if

    ! Allocate temporary arrays
    allocate(times_tmp(nt), source=0.0_rp)
    allocate(z_tmp(nz), source=0.0_rp)
    allocate(U_tmp(3, nz, nt), source=0.0_rp)
    allocate(T_tmp(nz, nt), source=0.0_rp)

    ! Every variable is required: a silent zero would nudge toward U=0 / T=0 K.
    call nc_read_required(ncid, 'time', filename, vec1=times_tmp)
    call nc_read_required(ncid, 'z',    filename, vec1=z_tmp)
    call nc_read_required(ncid, 'U',    filename, mat=U_tmp(1, :, :))
    call nc_read_required(ncid, 'V',    filename, mat=U_tmp(2, :, :))
    call nc_read_required(ncid, 'W',    filename, mat=U_tmp(3, :, :))
    call nc_read_required(ncid, 'T',    filename, mat=T_tmp)

    ! Close file
    stat_close = nf90_close(ncid)

    ! Store in data structure
    data%n_times = nt
    data%n_levels = nz
    allocate(data%times(nt), source=times_tmp)
    allocate(data%z_levels(nz), source=z_tmp)
    allocate(data%U_profile(3, nz, nt), source=U_tmp)
    allocate(data%T_profile(nz, nt), source=T_tmp)

    deallocate(times_tmp, z_tmp, U_tmp, T_tmp)

    write(*, '(A,I0,A,I0,A)') 'Loaded netCDF weather data: ', nt, ' times, ', nz, ' levels'

    data%loaded = .true.
  end subroutine nudging_load_netcdf

  subroutine nc_read_required(ncid, name, filename, vec1, mat)
    integer(ip),      intent(in)  :: ncid
    character(len=*), intent(in)  :: name, filename
    real(rp), intent(out), optional :: vec1(:)
    real(rp), intent(out), optional :: mat(:,:)
    integer(ip) :: st, varid

    st = nf90_inq_varid(ncid, name, varid)
    if (st == nf90_noerr) then
      if (present(vec1)) st = nf90_get_var(ncid, varid, vec1)
      if (present(mat))  st = nf90_get_var(ncid, varid, mat)
    end if
    if (st /= nf90_noerr) then
      write(*, '(A)') 'ERROR: weather file '//trim(filename)//' variable "'//name//'": '//trim(nf90_strerror(st))
      error stop 'nudging_load_netcdf: required variable missing or unreadable'
    end if
  end subroutine nc_read_required

  !> Apply Davies nudging: add relaxation term to the explicit RHS
  !> For each cell with z > z_start, add source:
  subroutine nudging_apply_rhs(data, U_cells, T_cells, z_cells, nc, time, tau, z_start, &
                                rhs_U, rhs_T, z_full, horizontal_only, dt)
    type(NudgingData), intent(in)    :: data
    real(rp),          intent(in)    :: U_cells(:,:)  ! (3, nc+ng)
    real(rp),          intent(in)    :: T_cells(:)  ! (nc+ng)
    real(rp),          intent(in)    :: z_cells(:)  ! (nc) - cell z-coordinates
    integer(ip),       intent(in)    :: nc
    real(rp),          intent(in)    :: time, tau, z_start
    real(rp),          intent(inout) :: rhs_U(:,:)  ! (3, nc)
    real(rp),          intent(inout) :: rhs_T(:)  ! (nc)
    real(rp),          intent(in), optional :: z_full  ! height of full strength; <= z_start = no ramp
    logical,           intent(in), optional :: horizontal_only  ! leave w unforced
    real(rp),          intent(in), optional :: dt  ! used to bound the explicit relaxation

    integer(ip) :: c
    real(rp)    :: U_ref(3), T_ref, nudge_coeff, ramp, z_top, coeff_max
    logical     :: horiz

    if (.not. data%loaded .or. data%n_times < 1) return

    z_top = z_start
    if (present(z_full)) z_top = z_full
    horiz = .true.
    if (present(horizontal_only)) horiz = horizontal_only

    ! Explicit relaxation is stable only for dt/tau below order one; clamp when dt is known.
    nudge_coeff = 1.0_rp / max(tau, SMALL)
    if (present(dt)) then
      if (dt > SMALL) then
        coeff_max = 0.5_rp / dt
        nudge_coeff = min(nudge_coeff, coeff_max)
      end if
    end if

    do c = 1, nc
      if (z_cells(c) < z_start) cycle

      ! Ramp from zero at z_start to full strength at z_full to avoid a step in the forcing.
      if (z_top > z_start) then
        ramp = min(1.0_rp, (z_cells(c) - z_start) / (z_top - z_start))
      else
        ramp = 1.0_rp
      end if

      call nudging_interpolate(data, time, z_cells(c), U_ref, T_ref)

      rhs_U(1, c) = rhs_U(1, c) - ramp * nudge_coeff * (U_cells(1, c) - U_ref(1))
      rhs_U(2, c) = rhs_U(2, c) - ramp * nudge_coeff * (U_cells(2, c) - U_ref(2))
      ! Nudging w would damp resolved vertical motion, so it is off by default.
      if (.not. horiz) &
        rhs_U(3, c) = rhs_U(3, c) - ramp * nudge_coeff * (U_cells(3, c) - U_ref(3))
      rhs_T(c)    = rhs_T(c)    - ramp * nudge_coeff * (T_cells(c) - T_ref)
    end do

  end subroutine nudging_apply_rhs

  !> Bilinear interpolation in (time, height) space
  !> Returns reference U(3) and T at current time and height
  subroutine nudging_interpolate(data, time, z, U_ref, T_ref)
    type(NudgingData), intent(in)  :: data
    real(rp),          intent(in)  :: time, z
    real(rp),          intent(out) :: U_ref(3), T_ref

    integer(ip) :: it1, it2, iz1, iz2
    real(rp)    :: wt, wz, U_11(3), U_12(3), U_21(3), U_22(3), T_11, T_12, T_21, T_22

    U_ref = 0.0_rp
    T_ref = 0.0_rp

    if (.not. data%loaded) return
    if (.not. allocated(data%times) .or. .not. allocated(data%z_levels)) return

    ! Find surrounding time indices
    call interp_find_indices(data%times, time, it1, it2, wt)

    ! Find surrounding height indices
    call interp_find_indices(data%z_levels, z, iz1, iz2, wz)

    ! Bilinear interpolation
    U_11 = data%U_profile(:, iz1, it1)
    U_12 = data%U_profile(:, iz1, it2)
    U_21 = data%U_profile(:, iz2, it1)
    U_22 = data%U_profile(:, iz2, it2)

    U_ref = (1.0_rp - wt) * (1.0_rp - wz) * U_11 &
          + wt * (1.0_rp - wz) * U_12 &
          + (1.0_rp - wt) * wz * U_21 &
          + wt * wz * U_22

    T_11 = data%T_profile(iz1, it1)
    T_12 = data%T_profile(iz1, it2)
    T_21 = data%T_profile(iz2, it1)
    T_22 = data%T_profile(iz2, it2)

    T_ref = (1.0_rp - wt) * (1.0_rp - wz) * T_11 &
          + wt * (1.0_rp - wz) * T_12 &
          + (1.0_rp - wt) * wz * T_21 &
          + wt * wz * T_22

  end subroutine nudging_interpolate

  !> Find surrounding indices and weight for 1D linear interpolation
  subroutine interp_find_indices(arr, val, i1, i2, weight)
    real(rp), intent(in)  :: arr(:), val
    integer(ip), intent(out) :: i1, i2
    real(rp), intent(out) :: weight

    integer(ip) :: n, i
    real(rp)    :: denom

    n = size(arr)
    i1 = 1
    i2 = n

    ! Find bracketing indices
    do i = 1, n - 1
      if (val >= arr(i) .and. val <= arr(i+1)) then
        i1 = i
        i2 = i + 1
        exit
      end if
    end do

    ! Compute weight in [0, 1]
    denom = arr(i2) - arr(i1)
    if (abs(denom) < SMALL) then
      weight = 0.0_rp
    else
      weight = (val - arr(i1)) / denom
      weight = max(0.0_rp, min(1.0_rp, weight))
    end if

  end subroutine interp_find_indices

end module mod_abl_weather
