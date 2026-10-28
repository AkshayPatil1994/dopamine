! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Named profiler ranges (NVTX); no-ops unless DOPAMINE_NVTX is defined.
module mod_profiling
  use iso_c_binding, only: c_char, c_null_char
  implicit none
  private

  public :: prof_push, prof_pop

#ifdef DOPAMINE_NVTX
  interface
    subroutine dopamine_nvtx_push(name) bind(C)
      import :: c_char
      character(kind=c_char), intent(in) :: name(*)
    end subroutine
    subroutine dopamine_nvtx_pop() bind(C)
    end subroutine
  end interface
#endif

contains

  subroutine prof_push(name)
    character(*), intent(in) :: name
#ifdef DOPAMINE_NVTX
    call dopamine_nvtx_push(name//c_null_char)
#else
    if (.false.) print *, name
#endif
  end subroutine prof_push

  subroutine prof_pop()
#ifdef DOPAMINE_NVTX
    call dopamine_nvtx_pop()
#endif
  end subroutine prof_pop

end module mod_profiling
