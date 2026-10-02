! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Round trip of the inflow-plane database file: header, coordinates and frames.
program test_inflow_db
  use mod_kinds
  use mod_inflow_db
  implicit none

  integer :: n_fail = 0
  character(len=*), parameter :: fn = 'inflow_db_test.nc'
  type(InflowDB) :: db
  type(InflowDBHeader) :: h
  real(rp), allocatable :: plane(:,:,:), back(:,:,:)
  integer(ip) :: is, it, j, k, c
  real(rp) :: err

  write(*,'(A)') 'inflow database round trip'
  h%ny = 8; h%nz = 5; h%nt = 6; h%nseg = 2
  h%W = 4.0_rp; h%y0 = 0.25_rp; h%dt_p = 0.05_rp
  h%frame = reshape([1.0_rp,0.0_rp,0.0_rp, 0.0_rp,1.0_rp,0.0_rp, 0.0_rp,0.0_rp,1.0_rp], [3,3])
  h%n_in = [1.0_rp, 0.0_rp, 0.0_rp]; h%t1 = [0.0_rp, 1.0_rp, 0.0_rp]
  h%wind_dir_deg = 270.0_rp; h%u_ref = 10.0_rp; h%z_ref = 5.0_rp; h%T_turnover = 3.5_rp
  h%seed = 42; h%generator = 'unit-test'
  allocate(h%z(5), h%U_mean(5), h%R_target(6, 5))
  h%z = [0.1_rp, 0.3_rp, 0.8_rp, 1.9_rp, 4.0_rp]    ! graded
  h%U_mean = 2.0_rp * h%z
  h%R_target = 0.0_rp; h%R_target(1, :) = 1.0_rp + h%z

  allocate(plane(h%ny, h%nz, 3), back(h%ny, h%nz, 3))
  call inflow_db_create(db, fn, h)
  do is = 1, 2
    do it = 1, 6
      do c = 1, 3
        do k = 1, 5
          do j = 1, 8
            plane(j, k, c) = sin(0.7_rp * j + 1.3_rp * k) * c + 0.1_rp * it + 10.0_rp * is
          end do
        end do
      end do
      call inflow_db_write_frame(db, is, it, plane)
    end do
  end do
  call inflow_db_close(db)

  call inflow_db_open(db, fn)
  call check(db%h%ny == 8 .and. db%h%nz == 5 .and. db%h%nt == 6 .and. db%h%nseg == 2, 'dimensions')
  call check(abs(db%h%W - 4.0_rp) + abs(db%h%y0 - 0.25_rp) + abs(db%h%dt_p - 0.05_rp) < 1.0e-14_rp, 'W, y0, dt_p')
  call check(maxval(abs(db%h%frame - h%frame)) == 0.0_rp .and. maxval(abs(db%h%n_in - h%n_in)) == 0.0_rp, 'frame, n_in')
  call check(abs(db%h%wind_dir_deg - 270.0_rp) + abs(db%h%T_turnover - 3.5_rp) < 1.0e-14_rp .and. &
             db%h%seed == 42 .and. trim(db%h%generator) == 'unit-test', 'scalars and generator')
  call check(maxval(abs(db%h%z - h%z)) == 0.0_rp .and. maxval(abs(db%h%R_target - h%R_target)) == 0.0_rp, 'z levels, targets')
  call check(abs(inflow_db_node_y(db%h, 3_ip) - 1.25_rp) < 1.0e-14_rp, 'node coordinates y0 + j W/ny')
  call check(abs(inflow_db_segment_time(db%h) - 0.25_rp) < 1.0e-14_rp, 'segment time = (nt-1) dt_p')

  err = 0.0_rp
  do is = 1, 2
    do it = 1, 6
      call inflow_db_read_frame(db, is, it, back)
      do c = 1, 3
        do k = 1, 5
          do j = 1, 8
            err = max(err, abs(back(j, k, c) - (sin(0.7_rp * j + 1.3_rp * k) * c + 0.1_rp * it + 10.0_rp * is)))
          end do
        end do
      end do
    end do
  end do
  call check(err < 2.0e-6_rp * 40.0_rp, 'frames read back to float32 precision')
  call inflow_db_close(db)
  call execute_command_line('rm -f '//fn)

  if (n_fail > 0) then
    write(*,'(A,I0,A)') 'test_inflow_db: ', n_fail, ' FAILED'
    error stop 1
  end if
  write(*,'(A)') 'test_inflow_db: all passed'

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
end program test_inflow_db
