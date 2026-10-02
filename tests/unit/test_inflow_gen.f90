! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> The offline generator: statistics against the targets, independence of segments, temporal
!> resolution (frame-to-frame continuity), turnover-based length, grid-file round trip.
program test_inflow_gen
  use mod_kinds
  use mod_abl_types
  use mod_inflow_db
  use mod_inflow_gen
  implicit none

  integer :: n_fail = 0
  type(ABLConfig) :: cfg
  type(InflowGenParams) :: gp
  type(InflowDB) :: db
  real(rp), allocatable :: p1(:,:,:), p2(:,:,:), z(:)
  real(rp) :: a, b, c12, c11, c22, lag1, l11, l22, W, y0
  integer(ip) :: it, k, ny
  logical :: ok

  cfg%u_ref = 1.0_rp; cfg%z_ref = 1.0_rp; cfg%z0_m = 0.001_rp; cfg%wind_dir_deg = 270.0_rp
  cfg%sem_sigma = 0.15_rp; cfg%sem_L_exponent = 0.3_rp
  cfg%sem_mean_law = 'power'; cfg%sem_alpha = 0.15_rp; cfg%sem_Iu_ref = 0.1_rp

  write(*,'(A)') 'inflow generator'
  gp%file = 'unit_gen.nc'; gp%W = 1.0_rp; gp%y0 = 0.03125_rp; gp%ny = 16; gp%nz = 16
  gp%z_min = 0.03125_rp; gp%z_max = 0.96875_rp; gp%n_segments = 3; gp%n_turnovers = 2.0_rp
  gp%verify = .true.
  call inflow_gen_generate(cfg, gp, ok)
  call check(ok, 'plane statistics within 15% of the targets (rms, u''w'')')

  call inflow_db_open(db, trim(gp%file))
  call check(db%h%nseg == 3 .and. db%h%ny == 16 .and. db%h%nz == 16, 'grid and segment count')
  call check(abs(db%h%T_turnover - 0.15_rp / 0.1_rp) < 1.0e-9_rp, 'turnover time = Lx / sigma_u at z_ref = 1.5 s')
  call check(inflow_db_segment_time(db%h) >= 2.0_rp * db%h%T_turnover - db%h%dt_p, 'segment length is two turnovers')
  call check(abs(db%h%n_in(1) - 1.0_rp) < 1.0e-12_rp .and. abs(db%h%t1(2) - 1.0_rp) < 1.0e-12_rp, &
             'plane frame: flow along +x, span along +y')

  ! Segments are independent realisations.
  allocate(p1(16, 16, 3), p2(16, 16, 3))
  c12 = 0.0_rp; c11 = 0.0_rp; c22 = 0.0_rp
  lag1 = 0.0_rp; l11 = 0.0_rp; l22 = 0.0_rp
  do it = 1, db%h%nt - 1
    call inflow_db_read_frame(db, 1_ip, it, p1)
    call inflow_db_read_frame(db, 2_ip, it, p2)
    c12 = c12 + sum(p1(:, :, 1) * p2(:, :, 1)); c11 = c11 + sum(p1(:, :, 1)**2); c22 = c22 + sum(p2(:, :, 1)**2)
    ! frame-to-frame continuity within segment 1
    call inflow_db_read_frame(db, 1_ip, it + 1_ip, p2)
    lag1 = lag1 + sum(p1(:, :, 1) * p2(:, :, 1)); l11 = l11 + sum(p1(:, :, 1)**2); l22 = l22 + sum(p2(:, :, 1)**2)
  end do
  write(*,'(A,F7.3,A,F7.4)') '    corr(seg1,seg2) = ', c12 / sqrt(c11 * c22), &
    '   lag-1 frame correlation = ', lag1 / sqrt(l11 * l22)
  call check(abs(c12 / sqrt(c11 * c22)) < 0.05_rp, 'segments uncorrelated (|rho| < 0.05)')
  call check(lag1 / sqrt(l11 * l22) > 0.95_rp, 'frame spacing resolves the signal (lag-1 correlation > 0.95)')
  call inflow_db_close(db)

  ! Grid file round trip, and generation on a graded z grid taken from it.
  allocate(z(12))
  do k = 1, 12
    z(k) = 0.02_rp * 1.35_rp**(k - 1)
  end do
  call inflow_write_grid_file('unit_grid.txt', 1.0_rp, 0.0625_rp, 8_ip, [1.0_rp, 0.0_rp, 0.0_rp], z)
  block
    real(rp), allocatable :: zr(:)
    real(rp) :: nn(3)
    call inflow_read_grid_file('unit_grid.txt', W, y0, ny, nn, zr)
    call check(abs(W - 1.0_rp) + abs(y0 - 0.0625_rp) < 1.0e-14_rp .and. ny == 8 .and. size(zr) == 12 .and. &
               maxval(abs(zr - z)) < 1.0e-14_rp .and. abs(nn(1) - 1.0_rp) < 1.0e-14_rp, 'grid file round trip')
  end block
  gp%grid_file = 'unit_grid.txt'; gp%file = 'unit_gen2.nc'; gp%n_segments = 2; gp%n_turnovers = 1.0_rp; gp%verify = .false.
  cfg%n_eddies = 0
  call inflow_gen_generate(cfg, gp, ok)
  call inflow_db_open(db, 'unit_gen2.nc')
  call check(db%h%ny == 8 .and. db%h%nz == 12 .and. maxval(abs(db%h%z - z)) < 1.0e-14_rp .and. &
             abs(inflow_db_node_y(db%h, 1_ip) - 0.0625_rp) < 1.0e-14_rp, 'database generated on the exported (graded) grid')
  call inflow_db_close(db)
  call execute_command_line('rm -f unit_gen.nc unit_gen2.nc unit_grid.txt')

  if (n_fail > 0) then
    write(*,'(A,I0,A)') 'test_inflow_gen: ', n_fail, ' FAILED'
    error stop 1
  end if
  write(*,'(A)') 'test_inflow_gen: all passed'

contains
  subroutine check(cond, msg)
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (cond) then
      write(*,'(A,A)') '  ok    ', msg
    else
      write(*,'(A,A)') '  FAIL  ', msg
      n_fail = n_fail + 1
    end if
  end subroutine check
end program test_inflow_gen
