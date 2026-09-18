! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Offline mesh decomposition: dopamine-decompose <input> --nranks N [--chunk-size M] [--output-dir DIR].

program dopamine_decompose
  use mod_kinds
  use mod_input_parser
  use mod_streaming_decompose, only: decompose_streaming
  use mpi_f08
  implicit none

  ! -- MPI
  integer(ip) :: rank_mpi, nranks_mpi
  integer(mp) :: ierr

  ! -- Timing
  integer(lp) :: t_prog_start, t_prog_end, clock_rate_main

  ! -- Input/config
  character(256) :: input_file, output_dir, chunk_str, mem_str, thread_str
  integer(ip) :: nranks_target, chunk_size, n_threads
  real(rp) :: max_memory_gb
  logical :: weight_by_volume
  type(InputConfig) :: config
  integer(ip) :: i, ios

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_mpi, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks_mpi, ierr)
  call system_clock(t_prog_start, clock_rate_main)

  ! Parse command-line arguments
  input_file = ''
  nranks_target = -1
  output_dir = '.'
  chunk_size = 2_ip * 1000_ip * 1000_ip  ! default 2M cells per chunk
  max_memory_gb = 0.0_rp
  weight_by_volume = .false.
  n_threads = 0

  i = 1
  do while (i <= command_argument_count())
    block
      character(len=256) :: arg
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--nranks')
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *, iostat=ios) nranks_target
        if (ios /= 0) then
          if (rank_mpi == 0) write(*, '(A)') 'ERROR: --nranks requires integer argument'
          call MPI_Finalize(ierr)
          stop 1
        end if
      case ('--output-dir')
        i = i + 1
        call get_command_argument(i, output_dir)
      case ('--chunk-size')
        i = i + 1
        call get_command_argument(i, chunk_str)
        read(chunk_str, *, iostat=ios) chunk_size
        if (ios /= 0) chunk_size = 2_ip * 1000_ip * 1000_ip
      case ('--weight-by-volume')
        weight_by_volume = .true.
      case ('--threads')
        i = i + 1
        call get_command_argument(i, thread_str)
        read(thread_str, *, iostat=ios) n_threads
        if (ios /= 0) n_threads = 0
      case ('--max-memory-gb')
        i = i + 1
        call get_command_argument(i, mem_str)
        read(mem_str, *, iostat=ios) max_memory_gb
        if (ios /= 0) max_memory_gb = 0.0_rp
      case default
        if (input_file == '') then
          input_file = trim(arg)
        else
          if (rank_mpi == 0) write(*, '(A)') 'WARNING: unexpected argument: '//trim(arg)
        end if
      end select
    end block
    i = i + 1
  end do

  ! -- Validate arguments
  if (input_file == '' .or. nranks_target < 1) then
    if (rank_mpi == 0) then
      write(*, '(A)') 'Usage: dopamine-decompose <input_file> --nranks N [--chunk-size M] [--output-dir DIR] [--max-memory-gb G]'
      write(*, '(A)') ''
      write(*, '(A)') 'Example:'
      write(*, '(A)') '  dopamine-decompose mesh.input --nranks 128 --chunk-size 2000000 --output-dir ./parts --max-memory-gb 48'
      write(*, '(A)') ''
      write(*, '(A)') 'Memory-efficient streaming decomposition:'
      write(*, '(A)') '  Reads mesh in chunks of --chunk-size cells (default: 2M)'
      write(*, '(A)') '  Assigns to ranks using a space-filling curve, equal cells per rank'
      write(*, '(A)') '  --threads N          writer threads (default: all cores)'
      write(*, '(A)') '  --weight-by-volume   balance total volume instead of cell count'
      write(*, '(A)') '  Writes local mesh files: processorN/mesh.bin'
      write(*, '(A)') ''
      write(*, '(A)') 'Scales to 200M+ cells on workstations with 64GB RAM'
    end if
    call MPI_Finalize(ierr)
    stop 1
  end if

  if (rank_mpi == 0) then
    write(*, '(A)') '================================================='
    write(*, '(A)') ' DOPAMINE -- Offline Mesh Decomposer'
    write(*, '(A)') '================================================='
    write(*, '(A,A,A,I0)') ' Input : ', trim(input_file), '  |  Target ranks: ', nranks_target
    write(*, '(A,A)') ' Output: ', trim(output_dir)
    write(*, '(A)') '-------------------------------------------------'
  end if

  ! Parse input file
  call parse_input_file(trim(input_file), config)

  ! Perform streaming decomposition (rank 0 only)
  if (rank_mpi == 0) then
    call decompose_streaming(trim(input_file), config, nranks_target, &
                            trim(output_dir), chunk_size, max_memory_gb, &
                            weight_by_volume=weight_by_volume, n_threads=n_threads)

    call system_clock(t_prog_end)
    write(*, '(A)') ' '
    write(*, '(A)') ' [6/6]  Done.'
    write(*, '(A,F8.2,A)') ' Total time: ', &
        real(t_prog_end - t_prog_start, rp)/real(clock_rate_main, rp), ' s'
    write(*, '(A)') '-------------------------------------------------'
    write(*, '(A)') ' Decomposition complete.'
  end if

  call MPI_Finalize(ierr)

end program dopamine_decompose
