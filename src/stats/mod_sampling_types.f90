! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Shared data types for on-the-fly sampling definitions.
module mod_sampling_types
  use mod_kinds
  implicit none
  private

  integer(ip), parameter, public :: MAX_SAMPLE_OBJECTS = 16_ip

  type, public :: PlaneSampleDef
    logical              :: enabled = .false.
    character(len=64)    :: name = ''
    real(rp)             :: origin(3) = [0.0_rp, 0.0_rp, 0.0_rp]
    real(rp)             :: normal(3) = [0.0_rp, 0.0_rp, 1.0_rp]
    real(rp)             :: tangent_hint(3) = [1.0_rp, 0.0_rp, 0.0_rp]
    real(rp)             :: half_width = 0.5_rp
    real(rp)             :: half_height = 0.5_rp
    integer(ip)          :: nx = 64_ip
    integer(ip)          :: ny = 64_ip
    integer(ip)          :: write_frequency = 0_ip
  end type PlaneSampleDef

  type, public :: LineSampleDef
    logical              :: enabled = .false.
    character(len=64)    :: name = ''
    real(rp)             :: p0(3) = [0.0_rp, 0.0_rp, 0.0_rp]
    real(rp)             :: p1(3) = [1.0_rp, 0.0_rp, 0.0_rp]
    integer(ip)          :: npts = 128_ip
    integer(ip)          :: write_frequency = 0_ip
  end type LineSampleDef

  type, public :: CurveSampleDef
    logical              :: enabled = .false.
    character(len=64)    :: name = ''
    character(len=256)   :: file = ''
    integer(ip)          :: write_frequency = 0_ip
  end type CurveSampleDef

  type, public :: SurfaceSampleDef
    logical              :: enabled = .false.
    character(len=64)    :: name = ''
    character(len=256)   :: file = ''
    integer(ip)          :: write_frequency = 0_ip
    logical              :: intersect_with_planes = .true.
    real(rp)             :: intersection_decimate_spacing = 0.0_rp
    logical              :: sort_intersection_points = .true.
  end type SurfaceSampleDef

end module mod_sampling_types
