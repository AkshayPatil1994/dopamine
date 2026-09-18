! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Per-rank binary checkpoint serialization/deserialization of FVSolver state.
module mod_checkpoint_io
  use mod_kinds
  implicit none
  private
  public :: write_checkpoint, read_checkpoint

contains

  !> Write checkpoint file for a single rank
  !> File contains: header, time, global_id, fields (U,p,T,C)
  subroutine write_checkpoint(solver, msh, checkpoint_dir, step, rank)
    ! Import types locally to avoid module-level dependency issues
    use mod_mesh,    only: Mesh
    use mod_fv_solver, only: FVSolver
    type(FVSolver),       intent(in) :: solver
    type(Mesh),           intent(in) :: msh
    character(len=*),     intent(in) :: checkpoint_dir
    integer(ip),          intent(in) :: step
    integer(ip),          intent(in) :: rank

    integer(ip) :: iunit, ios, nc, ns, i
    character(len=512) :: filename, dir_cmd
    character(len=16)  :: rank_str, step_str
    integer(ip) :: header(8)
    integer(ip) :: has_C
    logical :: dir_exists

    nc = msh%n_cells
    ns = solver%params%n_scalars

    ! Determine which optional fields to write
    has_C = 0;  if (ns > 0) has_C = 1

    ! Pack header
    header(1) = nc
    header(2) = ns
    header(3) = solver%step
    header(4) = solver%step
    header(5) = has_C
    header(6) = 0  ! reserved
    header(7) = 0  ! reserved
    header(8) = 0  ! reserved

    ! Create processorR directory if needed (each rank creates its own to avoid race condition)
    write(rank_str, '(I0)') rank
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)
    inquire(file=trim(filename), exist=dir_exists)
    if (.not. dir_exists) then
      write(dir_cmd, '(A)') 'mkdir -p '//trim(filename)
      call execute_command_line(trim(dir_cmd), wait=.true., exitstat=ios)
      if (ios /= 0) error stop 'write_checkpoint: mkdir failed'
    end if

    ! Construct checkpoint filename: checkpoint_XXXXXXXX.bin (8-digit zero-padded step)
    write(step_str, '(I8.8)') step
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)//'/checkpoint_'//trim(step_str)//'.bin'

    ! Open file for writing
    open(newunit=iunit, file=trim(filename), form='unformatted', access='sequential', &
         action='write', status='replace', iostat=ios)
    if (ios /= 0) error stop 'write_checkpoint: cannot open '//trim(filename)

    ! Record 1: header
    write(iunit) header

    ! Record 2: time
    write(iunit) solver%time

    ! Record 3: global_id (owned cells only, 1..nc)
    if (nc > 0) then
      if (allocated(msh%global_id)) then
        write(iunit) msh%global_id(1:nc)
      else
        ! Fallback: write 1-based indices if global_id not allocated
        write(iunit) [(int(i, lp), i=1,nc)]
      end if
    end if

    ! Record 4-7: Primary fields (owned cells only)
    if (nc > 0) then
      write(iunit) solver%U(:, 1:nc)
      write(iunit) solver%p(1:nc)
      write(iunit) solver%T(1:nc)
      write(iunit) solver%nu_sgs(1:nc)
    end if

    ! Record 8: Passive scalars (if present)
    if (has_C == 1 .and. nc > 0) then
      write(iunit) solver%C(1:ns, 1:nc)
    end if

    close(iunit)

    ! Rank 0 writes checkpoint_latest manifest
    if (rank == 0) then
      filename = trim(checkpoint_dir)//'/checkpoint_latest'
      open(newunit=iunit, file=trim(filename), form='formatted', &
           action='write', status='replace', iostat=ios)
      if (ios /= 0) error stop 'write_checkpoint: cannot open checkpoint_latest'
      write(iunit, '(I0)') step
      close(iunit)
    end if

  end subroutine write_checkpoint


  !> Read checkpoint file for a single rank
  !> Populates FVSolver fields from checkpoint_XXXXXXXX.bin
  subroutine read_checkpoint(solver, msh, checkpoint_dir, step, rank)
    ! Import types locally to avoid module-level dependency issues
    use mod_mesh,    only: Mesh
    use mod_fv_solver, only: FVSolver
    type(FVSolver),       intent(inout) :: solver
    type(Mesh),           intent(in) :: msh
    character(len=*),     intent(in) :: checkpoint_dir
    integer(ip),          intent(in) :: step
    integer(ip),          intent(in) :: rank

    integer(ip) :: iunit, ios, nc, ns, nc_ckpt, ns_ckpt, step_ckpt
    character(len=512) :: filename
    character(len=16)  :: rank_str, step_str
    integer(ip) :: header(8)
    integer(lp), allocatable :: global_id_ckpt(:)

    nc = msh%n_cells
    ns = solver%params%n_scalars

    ! Construct checkpoint filename
    write(rank_str, '(I0)') rank
    write(step_str, '(I8.8)') step
    filename = trim(checkpoint_dir)//'/processor'//trim(rank_str)//'/checkpoint_'//trim(step_str)//'.bin'

    ! Open file for reading
    open(newunit=iunit, file=trim(filename), form='unformatted', access='sequential', &
         action='read', status='old', iostat=ios)
    if (ios /= 0) error stop 'read_checkpoint: cannot open '//trim(filename)

    ! Record 1: header
    read(iunit) header
    nc_ckpt = header(1)
    ns_ckpt = header(2)
    step_ckpt = header(4)

    ! Validate dimensions
    if (nc_ckpt /= nc) then
      write(*, '(A,I0,A,I0)') 'ERROR: checkpoint n_cells=', nc_ckpt, &
                               ' does not match mesh n_cells=', nc
      error stop 'read_checkpoint: dimension mismatch (n_cells)'
    end if
    if (ns_ckpt /= ns) then
      write(*, '(A,I0,A,I0)') 'ERROR: checkpoint n_scalars=', ns_ckpt, &
                               ' does not match solver n_scalars=', ns
      error stop 'read_checkpoint: dimension mismatch (n_scalars)'
    end if

    ! Record 2: time
    read(iunit) solver%time

    ! Record 3: global_id (just for validation/tracking; not used in read)
    if (nc > 0) then
      allocate(global_id_ckpt(nc))
      read(iunit) global_id_ckpt
      deallocate(global_id_ckpt)
    end if

    ! Record 4-7: Primary fields
    if (nc > 0) then
      read(iunit) solver%U(1:3, 1:nc)
      read(iunit) solver%p(1:nc)
      read(iunit) solver%T(1:nc)
      read(iunit) solver%nu_sgs(1:nc)
    end if

    ! Record 8: Passive scalars
    if (ns_ckpt > 0 .and. nc > 0) then
      read(iunit) solver%C(1:ns, 1:nc)
    end if

    ! Restore step counter
    solver%step = step_ckpt

    close(iunit)

  end subroutine read_checkpoint

end module mod_checkpoint_io
