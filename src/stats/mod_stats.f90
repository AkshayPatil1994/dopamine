! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> On-the-fly run statistics: CFL, KE, dissipation, divU; screen + CSV output.
module mod_stats
  use mod_kinds
  use mod_parameters
  use mod_mesh
  use mod_gradient
  use mod_fv_solver, only: FVSolver
  use mpi_f08
  implicit none
  private

  public :: StatsCollector
  public :: stats_init
  public :: stats_sample
  public :: stats_finalize

  type, public :: StatsCollector
    integer(ip)   :: freq     = 10  ! sample every this many steps
    character(256):: prefix   = 'dopamine'  ! output file prefix
    logical       :: header_written = .false.
    logical       :: file_opened = .false.
    integer(ip)   :: file_unit = -1
    real(rp)      :: E0_ref = -1.0_rp
  end type StatsCollector

contains

  !> Initialise stats collector: open CSV on rank 0; append to existing file on restart, else write header.
  subroutine stats_init(stats, prefix, freq, rank, restart, restart_step)
    type(StatsCollector), intent(inout) :: stats
    character(len=*),     intent(in)    :: prefix
    integer(ip),          intent(in)    :: freq
    integer(ip),          intent(in)    :: rank
    logical, optional,    intent(in)    :: restart
    integer(ip), optional, intent(in)   :: restart_step  ! rows after this step are dropped from an existing CSV
    logical    :: do_restart
    integer(ip) :: ios

    stats%freq   = max(1_ip, freq)
    stats%prefix = trim(prefix)

    do_restart = .false.
    if (present(restart)) do_restart = restart

    if (rank == 0) then
      if (do_restart) then
        ! Drop the rows of the abandoned continuation and recover E0 before appending.
        if (present(restart_step)) call stats_prepare_restart(stats, restart_step)
        ! Try to append to an existing file; fall back to creating a new one.
        open(newunit=stats%file_unit, &
             file=trim(stats%prefix)//'_stats.csv', &
             status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) then
          ! File absent - create it and write the header.
          open(newunit=stats%file_unit, &
               file=trim(stats%prefix)//'_stats.csv', &
               status='replace', action='write', iostat=ios)
          if (ios /= 0) then
            write(*, '(A,I0)') 'WARNING: could not open stats CSV, iostat=', ios
            stats%file_unit = -1
            stats%file_opened = .false.
            return
          end if
          write(stats%file_unit, '(A)') &
            'step,time,E,eps,E_over_E0,eps_ref,Umax,max_abs_divU,volavg_abs_divU,' // &
            'max_abs_divU_cell,volavg_abs_divU_cell,' // &
            'Ub,tauw_mean,Cf_bulk,nu_sgs_mean,nu_sgs_over_nu,dt'
          flush(stats%file_unit)
        end if
      else
        open(newunit=stats%file_unit, &
             file=trim(stats%prefix)//'_stats.csv', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
          write(*, '(A,I0)') 'WARNING: could not open stats CSV, iostat=', ios
          stats%file_unit = -1
          stats%file_opened = .false.
          return
        end if
        ! CSV header
        write(stats%file_unit, '(A)') &
          'step,time,E,eps,E_over_E0,eps_ref,Umax,max_abs_divU,volavg_abs_divU,' // &
          'max_abs_divU_cell,volavg_abs_divU_cell,' // &
          'Ub,tauw_mean,Cf_bulk,nu_sgs_mean,nu_sgs_over_nu,dt'
        flush(stats%file_unit)
      end if
      stats%file_opened = .true.
    end if
  end subroutine stats_init

  !> Rewrite an existing stats CSV without the rows past restart_step and recover the E0 that
  !> E_over_E0 was normalised with (E / E_over_E0 of the first usable row).
  subroutine stats_prepare_restart(stats, restart_step)
    type(StatsCollector), intent(inout) :: stats
    integer(ip),          intent(in)    :: restart_step

    character(len=1024), allocatable :: lines(:), keep(:)
    character(len=1024) :: line
    character(len=512)  :: fname
    integer(ip) :: iu, ios, n, nkeep, i, step_r
    real(rp)    :: t_r, e_r, eps_r, eoe0_r
    logical     :: exists

    fname = trim(stats%prefix)//'_stats.csv'
    inquire(file=trim(fname), exist=exists)
    if (.not. exists) return
    open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    allocate(lines(1024))
    n = 0
    do
      read(iu, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (n == size(lines)) then
        allocate(keep(2 * n))
        keep(1:n) = lines
        call move_alloc(keep, lines)
      end if
      n = n + 1
      lines(n) = line
    end do
    close(iu)
    if (n < 1) return

    allocate(keep(n))
    nkeep = 1
    keep(1) = lines(1)  ! header
    do i = 2, n
      read(lines(i), *, iostat=ios) step_r, t_r, e_r, eps_r, eoe0_r
      if (ios /= 0) cycle
      if (step_r > restart_step) cycle
      nkeep = nkeep + 1
      keep(nkeep) = lines(i)
      if (stats%E0_ref <= SMALL .and. eoe0_r > SMALL .and. e_r > SMALL) stats%E0_ref = e_r / eoe0_r
    end do

    if (nkeep < n) then
      write(*, '(A,I0,A,I0,A)') ' Stats: dropping ', n - nkeep, ' row(s) after restart step ', restart_step, '.'
      open(newunit=iu, file=trim(fname), status='replace', action='write', iostat=ios)
      if (ios == 0) then
        do i = 1, nkeep
          write(iu, '(A)') trim(keep(i))
        end do
        close(iu)
      end if
    end if
  end subroutine stats_prepare_restart

  !> Sample statistics and write to screen + CSV.
  subroutine stats_sample(stats, solver, time, step, rank, E0, dt, output_written)
    type(StatsCollector), intent(inout) :: stats
    type(FVSolver),       intent(inout) :: solver
    real(rp),             intent(in)    :: time
    integer(ip),          intent(in)    :: step
    integer(ip),          intent(in)    :: rank
    real(rp), optional,   intent(in)    :: E0
    real(rp), optional,   intent(in)    :: dt
    logical,  optional,   intent(in)    :: output_written

    real(rp) :: E, eps, umax, max_div, abs_div, cfl
    real(rp) :: max_div_cell, abs_div_cell
    real(rp) :: ubulk, tauw_mean, cf_bulk, nu_sgs_mean, sgs_to_mol
    real(rp) :: E_norm, eps_ref, E0_val, dt_val
    character(len=16) :: output_tag

    ! One cell/boundary pass builds every partial sum, settled by a single SUM- and MAX-Allreduce (was 6 separate collectives).
    call compute_stats_sums(solver, E, eps, umax, ubulk, tauw_mean, nu_sgs_mean)
    call solver%divergence_metrics(max_div, abs_div, max_div_cell, abs_div_cell)

    dt_val = solver%params%dt
    if (present(dt)) dt_val = dt
    cfl = solver%compute_cfl_rate() * dt_val

    cf_bulk = tauw_mean / max(0.5_rp * ubulk**2, SMALL)
    sgs_to_mol = nu_sgs_mean / max(solver%params%nu, SMALL)

    ! Derived normalised quantities
    E0_val  = stats%E0_ref
    if (present(E0)) E0_val = E0
    if (E0_val <= SMALL) then
      E0_val = E
      stats%E0_ref = E0_val
    end if
    if (E0_val > SMALL) then
      E_norm  = E / E0_val
    else
      E_norm  = 1.0_rp
    end if
    ! Reference dissipation from initial energy and a characteristic time scale
    ! eps_ref = -dE/dt analytic (only meaningful for TGV); set to eps itself here
    eps_ref = eps
    output_tag = ' '
    if (present(output_written)) then
      if (output_written) output_tag = right_header_col('3D output ok', 16)
    end if

    if (rank /= 0) return

    ! ---- screen output -------------------------------------------------------
    if (.not. stats%header_written) then
      write(*, '(A)') ''
      write(*, '(A)') repeat('-', 104)
      write(*, '(A)') right_header_col('step', 8) // right_header_col('dt', 10) // &
                      right_header_col('time', 14) // right_header_col('CFL', 14) // &
                      right_header_col('max|U|', 14) // right_header_col('max|divU|', 14) // &
                      right_header_col('<|divU|>V', 14) // right_header_col('output', 16)
      write(*, '(A)') repeat('-', 104)
      stats%header_written = .true.
    end if

    if (present(dt)) then
      write(*, '(I8,ES10.3,5ES14.6,A)') step, dt, time, cfl, umax, max_div, abs_div, output_tag
    else
      write(*, '(I8,A10,5ES14.6,A)') step, '-', time, cfl, umax, max_div, abs_div, output_tag
    end if

    ! ---- CSV output ----------------------------------------------------------
    if (stats%file_opened .and. stats%file_unit /= -1) then
      if (present(dt)) then
        write(stats%file_unit, '(I0,16(",",ES18.9))') &
          step, time, E, eps, E_norm, eps_ref, umax, max_div, abs_div, &
          max_div_cell, abs_div_cell, &
          ubulk, tauw_mean, cf_bulk, nu_sgs_mean, sgs_to_mol, dt
      else
        write(stats%file_unit, '(I0,16(",",ES18.9))') &
          step, time, E, eps, E_norm, eps_ref, umax, max_div, abs_div, &
          max_div_cell, abs_div_cell, &
          ubulk, tauw_mean, cf_bulk, nu_sgs_mean, sgs_to_mol, 0.0_rp
      end if
      flush(stats%file_unit)
    end if

  end subroutine stats_sample

  !> Finalize stats output resources.
  subroutine stats_finalize(stats, rank)
    type(StatsCollector), intent(inout) :: stats
    integer(ip),          intent(in)    :: rank

    if (rank == 0 .and. stats%file_opened .and. stats%file_unit /= -1) then
      close(stats%file_unit)
      stats%file_unit = -1
      stats%file_opened = .false.
    end if
  end subroutine stats_finalize

  ! Worker routine (MPI-collective so must be called by ALL ranks)

  !> Compute KE, dissipation, max|U|, bulk |U|, mean wall shear, and mean nu_sgs with one SUM- and one MAX-Allreduce.
  subroutine compute_stats_sums(solver, E, eps, umax, ubulk, tauw_mean, nu_sgs_mean)
    type(FVSolver), intent(inout) :: solver
    real(rp), intent(out) :: E, eps, umax, ubulk, tauw_mean, nu_sgs_mean

    integer(ip) :: c, p, f, f_bnd, ierr
    real(rp)    :: vol_c, umag, nu_eff
    real(rp)    :: S11, S22, S33, S12, S13, S23, SijSij
    real(rp)    :: sums_local(9), sums_global(9)  ! 1:E_num 2:V_tot 3:eps_num 4-6:momentum 7:sgs_num 8:tau_num 9:tau_area
    real(rp)    :: umax_local

    associate(msh => solver%mesh, U => solver%U, grad_U => solver%grad_U, nu_sgs => solver%nu_sgs)

      sums_local = 0.0_rp
      umax_local = 0.0_rp
      do c = 1, msh%n_cells
        vol_c = msh%cells(c)%volume

        sums_local(1) = sums_local(1) + 0.5_rp * (U(1,c)**2 + U(2,c)**2 + U(3,c)**2) * vol_c
        sums_local(2) = sums_local(2) + vol_c

        ! Diagonal Sij (= symmetric part of velocity gradient) and off-diagonal Sij = 1/2(dUi/dxj + dUj/dxi)
        S11 = grad_U(1, 1, c); S22 = grad_U(2, 2, c); S33 = grad_U(3, 3, c)
        S12 = 0.5_rp * (grad_U(2, 1, c) + grad_U(1, 2, c))
        S13 = 0.5_rp * (grad_U(3, 1, c) + grad_U(1, 3, c))
        S23 = 0.5_rp * (grad_U(3, 2, c) + grad_U(2, 3, c))
        SijSij = S11**2 + S22**2 + S33**2 + 2.0_rp * (S12**2 + S13**2 + S23**2)
        nu_eff = solver%params%nu + nu_sgs(c)
        sums_local(3) = sums_local(3) + nu_eff * SijSij * vol_c

        sums_local(4:6) = sums_local(4:6) + U(:, c) * vol_c
        sums_local(7)   = sums_local(7)   + nu_sgs(c) * vol_c

        umag = sqrt(U(1,c)**2 + U(2,c)**2 + U(3,c)**2)
        if (umag > umax_local) umax_local = umag
      end do

      do p = 1, msh%n_patches
        if (msh%patches(p)%bc_id /= BC_WALL_MODEL) cycle
        do f = 1, msh%patches(p)%n_faces
          f_bnd = msh%patches(p)%start_face + f - 1
          sums_local(8) = sums_local(8) + sqrt(sum(solver%tau_w_bnd(:, f_bnd)**2)) * &
                          msh%faces(msh%n_int_faces + f_bnd)%area
          sums_local(9) = sums_local(9) + msh%faces(msh%n_int_faces + f_bnd)%area
        end do
      end do

      call MPI_Allreduce(sums_local, sums_global, 9, MPI_DOUBLE_PRECISION, MPI_SUM, msh%comm, ierr)
      call MPI_Allreduce(umax_local, umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, msh%comm, ierr)

    end associate

    if (sums_global(2) > SMALL) then
      E           = sums_global(1) / sums_global(2)
      eps         = 2.0_rp * sums_global(3) / sums_global(2)
      ubulk       = sqrt(sum((sums_global(4:6) / sums_global(2))**2))
      nu_sgs_mean = sums_global(7) / sums_global(2)
    else
      E = 0.0_rp; eps = 0.0_rp; ubulk = 0.0_rp; nu_sgs_mean = 0.0_rp
    end if

    if (sums_global(9) > SMALL) then
      tauw_mean = sums_global(8) / sums_global(9)
    else
      tauw_mean = 0.0_rp
    end if
  end subroutine compute_stats_sums

  pure function right_header_col(label, width) result(col)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: width
    character(len=width) :: col
    integer :: n

    col = ' '
    n = min(len_trim(label), width)
    if (n > 0) col(width-n+1:width) = trim(label(1:n))
  end function right_header_col

end module mod_stats
