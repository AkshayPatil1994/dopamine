! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Offline inflow-plane generator: dopamine-inflow-gen <input>.
program dopamine_inflow_gen
  use mod_kinds
  use mod_input_parser
  use mod_inflow_gen
  implicit none

  character(256)    :: input_file
  type(InputConfig) :: config
  logical           :: ok

  if (command_argument_count() < 1) then
    write(*,'(A)') 'Usage: dopamine-inflow-gen <input_file>'
    stop
  end if
  call get_command_argument(1, input_file)
  call parse_input_file(trim(input_file), config)

  write(*,'(A)') '================================================='
  write(*,'(A)') ' DOPAMINE inflow-plane generator'
  write(*,'(A)') '================================================='
  call inflow_gen_generate(config%abl_config, config%inflow_gen, ok)
  if (.not. ok) error stop 1
end program dopamine_inflow_gen
