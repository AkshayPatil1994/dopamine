! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Precision kinds for dopamine
module mod_kinds
  use iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private

  integer, parameter, public :: rp = real64  ! working real precision
  integer, parameter, public :: sp = real32  ! single precision
  integer, parameter, public :: ip = int32  ! default integer
  integer, parameter, public :: lp = int64  ! long integer (global cell IDs)
  integer, parameter, public :: mp = kind(1)  ! MPI argument kind: counts, ranks, tags, ierror

end module mod_kinds
