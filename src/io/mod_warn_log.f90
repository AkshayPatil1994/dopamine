! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Aggregated solver warning log: rank-local counters reduced on an interval, appended to a lazily-opened text file.
module mod_warn_log
  use mod_kinds
  use mpi_f08
  implicit none
  private

  public :: warn_log_init, warn_log_bump, warn_log_flush, warn_log_message, warn_log_finalize
  public :: warn_log_total, warn_log_filename, warn_log_set_step_start
  public :: WARN_WALL_UTAU_CAP, WARN_TVD_FALLBACK, WARN_CFL_CLAMP, WARN_PSOLVER_RESIDUAL

  integer(ip), parameter, public :: NWARN = 8

  integer(ip), parameter :: WARN_WALL_UTAU_CAP    = 1_ip
  integer(ip), parameter :: WARN_TVD_FALLBACK     = 2_ip
  integer(ip), parameter :: WARN_CFL_CLAMP        = 3_ip
  integer(ip), parameter :: WARN_PSOLVER_RESIDUAL = 4_ip

  character(*), parameter :: TAG(NWARN) = [ character(12) :: &
    '[wall-model]', '[advection] ', '[cfl]       ', '[psolver]   ', &
    '[spare5]    ', '[spare6]    ', '[spare7]    ', '[spare8]    ' ]

  integer(ip), save :: cnt(NWARN) = 0_ip           ! rank-local events since last flush
  real(rp),    save :: aux(NWARN) = 0.0_rp         ! rank-local worst-case magnitude since last flush
  integer(ip), save :: n_total = 0_ip              ! warning blocks written (rank 0)
  integer(ip), save :: step_lo = 1_ip              ! first step of the current accumulation window

  type(MPI_Comm), save :: log_comm
  integer(ip),  save :: log_rank = 0_ip
  integer(ip),  save :: log_nranks = 1_ip
  integer(ip),  save :: flush_every = 100_ip
  integer(ip),  save :: unit_no = -1_ip
  logical,      save :: is_open = .false.
  logical,      save :: is_init = .false.
  character(256), save :: log_file = 'solver_warnings.txt'
  character(256), save :: session_tag = ''

contains

  !> Store MPI handles and settings; does not touch the filesystem.
  subroutine warn_log_init(comm, rank, nranks, filename, interval, session, step_start)
    type(MPI_Comm), intent(in) :: comm
    integer(ip),    intent(in) :: rank, nranks, interval
    character(*),   intent(in) :: filename, session
    integer(ip),    intent(in), optional :: step_start
    log_comm    = comm
    log_rank    = rank
    log_nranks  = nranks
    flush_every = max(1_ip, interval)
    log_file    = trim(filename)
    session_tag = trim(session)
    cnt = 0_ip
    aux = 0.0_rp
    n_total = 0_ip
    step_lo = 1_ip
    if (present(step_start)) step_lo = step_start
    is_init = .true.
  end subroutine warn_log_init

  !> Accumulate n events on slot id, tracking the worst magnitude. Hot path: no I/O, no MPI.
  subroutine warn_log_bump(id, n, aux_val)
    integer(ip), intent(in) :: id
    integer(ip), intent(in), optional :: n
    real(rp),    intent(in), optional :: aux_val
    if (id < 1_ip .or. id > NWARN) return
    if (present(n)) then
      cnt(id) = cnt(id) + n
    else
      cnt(id) = cnt(id) + 1_ip
    end if
    if (present(aux_val)) aux(id) = max(aux(id), aux_val)
  end subroutine warn_log_bump

  !> Reduce counters across ranks and append a block if anything fired; called on the flush interval.
  subroutine warn_log_flush(step, time, force)
    integer(ip), intent(in) :: step
    real(rp),    intent(in) :: time
    logical,     intent(in), optional :: force
    integer(ip) :: cnt_g(NWARN), i, ierr
    real(rp)    :: aux_g(NWARN)
    logical     :: do_now

    if (.not. is_init) return
    do_now = (mod(step, flush_every) == 0_ip)
    if (present(force)) do_now = do_now .or. force
    if (.not. do_now) return

    call MPI_Reduce(cnt, cnt_g, NWARN, MPI_INTEGER, MPI_SUM, 0, log_comm, ierr)
    call MPI_Reduce(aux, aux_g, NWARN, MPI_DOUBLE_PRECISION, MPI_MAX, 0, log_comm, ierr)

    if (log_rank == 0_ip .and. any(cnt_g > 0_ip)) then
      call ensure_open()
      if (is_open) then
        do i = 1, NWARN
          if (cnt_g(i) > 0_ip) call write_line(step, time, i, cnt_g(i), aux_g(i))
        end do
        n_total = n_total + 1_ip
        flush(unit_no)
      end if
    end if

    cnt = 0_ip
    aux = 0.0_rp
    step_lo = step + 1_ip
  end subroutine warn_log_flush

  !> Append a one-shot warning from rank 0 (setup checks, solver failures).
  subroutine warn_log_message(tag_in, text, step, time)
    character(*), intent(in) :: tag_in, text
    integer(ip),  intent(in), optional :: step
    real(rp),     intent(in), optional :: time
    if (.not. is_init .or. log_rank /= 0_ip) return
    call ensure_open()
    if (.not. is_open) return
    if (present(step) .and. present(time)) then
      write(unit_no, '(A,I8,A,ES11.4,2A,1X,A)') 'step ', step, '  t=', time, '  ', tag_in, trim(text)
    else
      write(unit_no, '(2A,1X,A)') 'setup                     ', tag_in, trim(text)
    end if
    n_total = n_total + 1_ip
    flush(unit_no)
  end subroutine warn_log_message

  !> Final reduction and close.
  subroutine warn_log_finalize(step, time)
    integer(ip), intent(in) :: step
    real(rp),    intent(in) :: time
    if (.not. is_init) return
    call warn_log_flush(step, time, force=.true.)
    if (is_open) then
      close(unit_no)
      is_open = .false.
    end if
  end subroutine warn_log_finalize

  !> Set the first step of the accumulation window, once the restart step is known.
  subroutine warn_log_set_step_start(step)
    integer(ip), intent(in) :: step
    step_lo = step
  end subroutine warn_log_set_step_start

  integer(ip) function warn_log_total() result(n)
    n = n_total
  end function warn_log_total

  function warn_log_filename() result(f)
    character(:), allocatable :: f
    f = trim(log_file)
  end function warn_log_filename

  !> Open the log on first use and stamp a session header. Never called unless a warning fired.
  subroutine ensure_open()
    integer(ip) :: ios
    character(8)  :: d
    character(10) :: t
    if (is_open) return
    open(newunit=unit_no, file=trim(log_file), status='unknown', position='append', &
         action='write', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A)') 'WARNING: cannot open warning log ', trim(log_file)
      return
    end if
    is_open = .true.
    call date_and_time(date=d, time=t)
    write(unit_no, '(A)') '# ---------------------------------------------------------------'
    write(unit_no, '(9A)') '# dopamine solver warnings -- session ', &
      d(1:4), '-', d(5:6), '-', d(7:8), ' ', t(1:2)//':'//t(3:4)//':'//t(5:6), ''
    write(unit_no, '(A,I0,2A)') '# ranks=', log_nranks, '  ', trim(session_tag)
  end subroutine ensure_open

  !> Emit one aggregated warning line for slot id over the current step window.
  subroutine write_line(step, time, id, n, worst)
    integer(ip), intent(in) :: step, id, n
    real(rp),    intent(in) :: time, worst
    character(256) :: msg

    select case (id)
    case (WARN_WALL_UTAU_CAP)
      write(msg, '(A,I0,A,I0,A,I0,A)') 'matching point at y+ < 2: ', n, ' face-events (all ranks, steps ', &
        step_lo, '-', step, '); mesh too fine for a wall model'
    case (WARN_TVD_FALLBACK)
      write(msg, '(A,I0,A,I0,A,I0,A)') 'TVD far-upwind stencil missing: ', n, &
        ' face-flux events (all ranks, steps ', step_lo, '-', step, ') used 1st-order upwind'
    case (WARN_CFL_CLAMP)
      write(msg, '(A,I0,A,I0,A,I0,A,F0.3)') 'dt clamped to CFL-stable value on ', n, &
        ' steps (', step_lo, '-', step, '); worst requested/stable = ', worst
    case (WARN_PSOLVER_RESIDUAL)
      write(msg, '(A,I0,A,I0,A,I0,A,ES11.4)') 'pressure residual above atol floor on ', n, &
        ' solves (steps ', step_lo, '-', step, '); worst ||r||_2 = ', worst
    case default
      write(msg, '(A,I0,A)') 'unnamed warning slot: ', n, ' events'
    end select

    write(unit_no, '(A,I8,A,ES11.4,2A,1X,A)') 'step ', step, '  t=', time, '  ', TAG(id), trim(msg)
  end subroutine write_line

end module mod_warn_log
