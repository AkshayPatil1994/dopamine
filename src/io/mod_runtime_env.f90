! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Runtime environment configuration helpers.
module mod_runtime_env
  use iso_c_binding, only: c_char, c_int, c_null_char
  use mod_kinds, only: ip
  implicit none
  private

  public :: configure_hdf5_file_locking

  interface
    integer(c_int) function c_setenv(name, value, overwrite) bind(C, name='setenv')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: name(*)
      character(kind=c_char), intent(in) :: value(*)
      integer(c_int), value :: overwrite
    end function c_setenv
  end interface

contains

  !> Set HDF5 locking behavior for this process before NetCDF/HDF5 files are opened.
  subroutine configure_hdf5_file_locking(use_locking, rank)
    logical,     intent(in) :: use_locking
    integer(ip), intent(in) :: rank

    character(kind=c_char, len=*), parameter :: env_name = 'HDF5_USE_FILE_LOCKING'//c_null_char
    character(kind=c_char, len=*), parameter :: env_true = 'TRUE'//c_null_char
    character(kind=c_char, len=*), parameter :: env_false = 'FALSE'//c_null_char
    integer(c_int) :: ierr

    if (use_locking) then
      ierr = c_setenv(env_name, env_true, 1_c_int)
    else
      ierr = c_setenv(env_name, env_false, 1_c_int)
    end if

    if (rank == 0) then
      if (ierr /= 0_c_int) then
        write(*, '(A)') 'WARNING: failed to set HDF5_USE_FILE_LOCKING environment variable.'
      else if (use_locking) then
        write(*, '(A)') ' HDF5_USE_FILE_LOCKING=TRUE (strict locking enabled)'
      else
        write(*, '(A)') ' HDF5_USE_FILE_LOCKING=FALSE (concurrent read-friendly mode)'
      end if
    end if
  end subroutine configure_hdf5_file_locking

end module mod_runtime_env
