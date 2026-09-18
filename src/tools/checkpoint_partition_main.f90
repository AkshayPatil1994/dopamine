! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> dopamine-checkpoint-partition: Repartition checkpoint files
program dopamine_checkpoint_partition
  use mod_kinds
  use mod_mesh
  use mod_mesh_io
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  ! -- MPI
  integer(ip)    :: rank_mpi, nranks_mpi, ierr

  ! -- Input/config
  character(256) :: old_dir, new_dir
  integer(ip) :: old_nranks, new_nranks, step_target
  integer(ip) :: i, ios

  ! -- Mesh
  type(Mesh) :: tmp_mesh

  ! -- Global checkpoint data (gathered on rank 0)
  integer(lp) :: nc_global
  real(rp), allocatable :: U_global(:,:), p_global(:), T_global(:), nu_sgs_global(:)
  real(rp), allocatable :: C_global(:,:)
  real(rp), allocatable :: H_U_n1_global(:,:), H_U_n2_global(:,:)
  real(rp), allocatable :: H_T_n1_global(:), H_T_n2_global(:)
  real(rp), allocatable :: H_C_n1_global(:,:), H_C_n2_global(:,:)
  real(rp) :: time_ckpt
  integer(ip) :: step_ckpt, history_ckpt, n_scalars_ckpt

  ! -- Book-keeping
  integer(ip) :: r, iunit, nc_r, ns_r, hist_r, nc_new
  character(len=16) :: rank_str, step_str
  character(len=512) :: old_ckpt_file, new_ckpt_file, old_dir_r, new_dir_r
  character(len=512) :: checkpoint_latest_file
  integer(ip) :: header(8)
  integer(lp), allocatable :: global_id_r(:)
  integer(ip), allocatable :: gid_filled(:)
  real(rp), allocatable :: U_r(:,:), p_r(:), T_r(:), nu_sgs_r(:)
  real(rp), allocatable :: C_r(:,:)
  real(rp), allocatable :: H_U_n1_r(:,:), H_U_n2_r(:,:)
  real(rp), allocatable :: H_T_n1_r(:), H_T_n2_r(:)
  real(rp), allocatable :: H_C_n1_r(:,:), H_C_n2_r(:,:)
  logical :: file_exists
  integer(ip) :: c

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_mpi, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks_mpi, ierr)

  ! Parse command-line arguments
  old_dir = ''
  new_dir = ''
  old_nranks = -1
  new_nranks = -1
  step_target = -1

  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, rank_str)
    select case (trim(rank_str))
    case ('--old-dir')
      i = i + 1
      call get_command_argument(i, old_dir)
    case ('--new-dir')
      i = i + 1
      call get_command_argument(i, new_dir)
    case ('--old-nranks')
      i = i + 1
      call get_command_argument(i, rank_str)
      read(rank_str, *) old_nranks
    case ('--new-nranks')
      i = i + 1
      call get_command_argument(i, rank_str)
      read(rank_str, *) new_nranks
    case ('--step')
      i = i + 1
      call get_command_argument(i, rank_str)
      read(rank_str, *) step_target
    case default
      if (rank_mpi == 0) write(*, '(A)') 'Unknown argument: '//trim(rank_str)
    end select
    i = i + 1
  end do

  ! Validate arguments
  if (rank_mpi == 0) then
    if (len_trim(old_dir) == 0) then
      write(*, '(A)') 'ERROR: --old-dir not specified'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (len_trim(new_dir) == 0) then
      write(*, '(A)') 'ERROR: --new-dir not specified'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (old_nranks < 1) then
      write(*, '(A)') 'ERROR: --old-nranks must be >= 1'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (new_nranks < 1) then
      write(*, '(A)') 'ERROR: --new-nranks must be >= 1'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
  end if

  ! If step not specified, read checkpoint_latest
  if (rank_mpi == 0 .and. step_target < 0) then
    checkpoint_latest_file = trim(old_dir)//'/checkpoint_latest'
    inquire(file=trim(checkpoint_latest_file), exist=file_exists)
    if (file_exists) then
      open(newunit=iunit, file=trim(checkpoint_latest_file), form='formatted', &
           action='read', iostat=ios)
      if (ios == 0) then
        read(iunit, *, iostat=ios) step_target
        close(iunit)
      end if
    end if
    if (step_target < 0) then
      write(*, '(A)') 'ERROR: cannot determine checkpoint step (no --step and no checkpoint_latest)'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    write(*, '(A,I0)') ' Auto-detected step from checkpoint_latest: ', step_target
  end if

  ! Only rank 0 does the work
  if (rank_mpi /= 0) then
    call MPI_Finalize(ierr)
    stop
  end if

  write(*, '(A)') '================================================='
  write(*, '(A)') ' dopamine-checkpoint-partition'
  write(*, '(A)') '================================================='
  write(*, '(A,A)') ' Old dir     : ', trim(old_dir)
  write(*, '(A,I0)') ' Old nranks  : ', old_nranks
  write(*, '(A,A)') ' New dir     : ', trim(new_dir)
  write(*, '(A,I0)') ' New nranks  : ', new_nranks
  write(*, '(A,I0)') ' Step        : ', step_target
  write(*, '(A)') '-------------------------------------------------'

  ! Read all old-rank checkpoints into global arrays
  if (rank_mpi == 0) write(*, '(A)') ' Step 1: Reading old checkpoints...'

  nc_global = 0
  n_scalars_ckpt = 0
  step_ckpt = step_target
  history_ckpt = 0
  time_ckpt = 0.0_rp

  do r = 0, old_nranks - 1
    write(rank_str, '(I0)') r
    write(step_str, '(I8.8)') step_target
    old_dir_r = trim(old_dir)//'/processor'//trim(rank_str)
    old_ckpt_file = trim(old_dir_r)//'/checkpoint_'//trim(step_str)//'.bin'

    open(newunit=iunit, file=trim(old_ckpt_file), form='unformatted', access='sequential', &
         action='read', status='old', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') 'ERROR: cannot open '//trim(old_ckpt_file)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Read header
    read(iunit) header
    nc_r = header(1)
    ns_r = header(2)
    hist_r = header(3)

    ! Validate step and history consistency across ranks
    if (r == 0) then
      step_ckpt = header(4)
    else
      if (header(4) /= step_ckpt) then
        write(*, '(A,I0,A,I0,A,I0)') 'ERROR: rank ', r, ' has step=', header(4), &
                                      ' but rank 0 had step=', step_ckpt
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
      if (hist_r /= history_ckpt) then
        write(*, '(A,I0,A,I0,A,I0)') 'ERROR: rank ', r, ' has history=', hist_r, &
                                      ' but rank 0 had history=', history_ckpt
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
    end if

    ! On first rank, allocate global arrays and read mesh for nc_global
    if (r == 0) then
      ! Read mesh.bin to get nc_global
      old_dir_r = trim(old_dir)//'/processor0'
      call read_local_mesh_file(tmp_mesh, trim(old_dir_r)//'/mesh.bin')
      nc_global = tmp_mesh%n_cells_global
      n_scalars_ckpt = ns_r
      history_ckpt = hist_r

      if (nc_global <= 0) then
        write(*, '(A)') 'ERROR: invalid nc_global from mesh.bin'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

      allocate(U_global(3, nc_global), source=0.0_rp)
      allocate(p_global(nc_global), source=0.0_rp)
      allocate(T_global(nc_global), source=0.0_rp)
      allocate(nu_sgs_global(nc_global), source=0.0_rp)
      allocate(gid_filled(nc_global), source=0_ip)
      if (n_scalars_ckpt > 0) allocate(C_global(n_scalars_ckpt, nc_global), source=0.0_rp)
      if (history_ckpt > 0) then
        allocate(H_U_n1_global(3, nc_global), source=0.0_rp)
        allocate(H_U_n2_global(3, nc_global), source=0.0_rp)
        allocate(H_T_n1_global(nc_global), source=0.0_rp)
        allocate(H_T_n2_global(nc_global), source=0.0_rp)
        if (n_scalars_ckpt > 0) then
          allocate(H_C_n1_global(n_scalars_ckpt, nc_global), source=0.0_rp)
          allocate(H_C_n2_global(n_scalars_ckpt, nc_global), source=0.0_rp)
        end if
      end if
      call tmp_mesh%destroy()
      write(*, '(A,I0,A)') '   Allocating global arrays for ', int(nc_global), ' cells'
    end if

    ! Validate
    if (ns_r /= n_scalars_ckpt) then
      write(*, '(A,I0,A,I0,A,I0)') 'ERROR: rank ', r, ' has n_scalars=', ns_r, &
                                   ' but expected ', n_scalars_ckpt
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Read time (same for all ranks, so just overwrite)
    read(iunit) time_ckpt
    if (.not. ieee_is_finite(time_ckpt)) then
      write(*, '(A,I0,A,E14.6)') 'ERROR: rank ', r, ' has non-finite time_ckpt=', time_ckpt
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Read global_id and fields
    allocate(global_id_r(nc_r))
    read(iunit) global_id_r

    allocate(U_r(3, nc_r))
    read(iunit) U_r
    allocate(p_r(nc_r))
    read(iunit) p_r
    allocate(T_r(nc_r))
    read(iunit) T_r
    allocate(nu_sgs_r(nc_r))
    read(iunit) nu_sgs_r

    if (n_scalars_ckpt > 0) then
      allocate(C_r(n_scalars_ckpt, nc_r))
      read(iunit) C_r
    end if

    if (history_ckpt > 0) then
      allocate(H_U_n1_r(3, nc_r))
      read(iunit) H_U_n1_r
      allocate(H_U_n2_r(3, nc_r))
      read(iunit) H_U_n2_r
      allocate(H_T_n1_r(nc_r))
      read(iunit) H_T_n1_r
      allocate(H_T_n2_r(nc_r))
      read(iunit) H_T_n2_r
      if (n_scalars_ckpt > 0) then
        allocate(H_C_n1_r(n_scalars_ckpt, nc_r))
        read(iunit) H_C_n1_r
        allocate(H_C_n2_r(n_scalars_ckpt, nc_r))
        read(iunit) H_C_n2_r
      end if
    end if

    close(iunit)

    ! Scatter into global arrays using global_id
    do c = 1, nc_r
      if (global_id_r(c) < 1 .or. global_id_r(c) > nc_global) then
        write(*, '(A,I0,A,I0,A,I0)') 'ERROR: rank ', r, ' cell ', c, &
                                      ' has invalid global_id=', int(global_id_r(c))
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
      if (.not. all(ieee_is_finite(U_r(:, c)))) then
        write(*, '(A,I0,A,I0)') 'ERROR: NaN/Inf in U at rank ', r, ' cell ', c
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
      U_global(:, global_id_r(c)) = U_r(:, c)
      p_global(global_id_r(c)) = p_r(c)
      T_global(global_id_r(c)) = T_r(c)
      nu_sgs_global(global_id_r(c)) = nu_sgs_r(c)
      gid_filled(global_id_r(c)) = gid_filled(global_id_r(c)) + 1_ip
      if (n_scalars_ckpt > 0) C_global(:, global_id_r(c)) = C_r(:, c)
      if (history_ckpt > 0) then
        H_U_n1_global(:, global_id_r(c)) = H_U_n1_r(:, c)
        H_U_n2_global(:, global_id_r(c)) = H_U_n2_r(:, c)
        H_T_n1_global(global_id_r(c)) = H_T_n1_r(c)
        H_T_n2_global(global_id_r(c)) = H_T_n2_r(c)
        if (n_scalars_ckpt > 0) then
          H_C_n1_global(:, global_id_r(c)) = H_C_n1_r(:, c)
          H_C_n2_global(:, global_id_r(c)) = H_C_n2_r(:, c)
        end if
      end if
    end do

    ! Clean up rank-local arrays
    deallocate(global_id_r, U_r, p_r, T_r, nu_sgs_r)
    if (n_scalars_ckpt > 0) deallocate(C_r)
    if (history_ckpt > 0) then
      deallocate(H_U_n1_r, H_U_n2_r, H_T_n1_r, H_T_n2_r)
      if (n_scalars_ckpt > 0) deallocate(H_C_n1_r, H_C_n2_r)
    end if

    if (mod(r, max(1, old_nranks/10)) == 0 .or. r == old_nranks - 1) &
      write(*, '(A,I0,A,I0)') '   Rank ', r, ': scattered ', nc_r, ' cells'
  end do

  write(*, '(A)') ' Step 1 complete'

  ! Verify bijection: every global cell filled exactly once.
  if (any(gid_filled /= 1_ip)) then
    write(*, '(A,I0,A,I0,A)') 'ERROR: global_id coverage broken (missing=', &
         count(gid_filled == 0_ip), ', duplicates=', count(gid_filled > 1_ip), ')'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  deallocate(gid_filled)

  ! Write new-rank checkpoints
  if (rank_mpi == 0) write(*, '(A)') ' Step 2: Writing new checkpoints...'

  do r = 0, new_nranks - 1
    write(rank_str, '(I0)') r
    write(step_str, '(I8.8)') step_target
    new_dir_r = trim(new_dir)//'/processor'//trim(rank_str)

    ! Read mesh.bin to get new rank's global_id array
    call read_local_mesh_file(tmp_mesh, trim(new_dir_r)//'/mesh.bin')
    nc_new = tmp_mesh%n_cells
    if (.not. allocated(tmp_mesh%global_id) .and. tmp_mesh%n_cells > 0) then
      write(*, '(A,I0)') 'ERROR: global_id not allocated for new rank ', r
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Create checkpoint file
    new_ckpt_file = trim(new_dir_r)//'/checkpoint_'//trim(step_str)//'.bin'
    open(newunit=iunit, file=trim(new_ckpt_file), form='unformatted', access='sequential', &
         action='write', status='replace', iostat=ios)
    if (ios /= 0) then
      write(*, '(A)') 'ERROR: cannot open '//trim(new_ckpt_file)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! Write header
    header(1) = nc_new
    header(2) = n_scalars_ckpt
    header(3) = history_ckpt
    header(4) = step_ckpt
    header(5) = 0; if (history_ckpt > 0) header(5) = 1
    header(6) = 0; if (n_scalars_ckpt > 0) header(6) = 1
    header(7) = 0
    header(8) = 0
    write(iunit) header

    ! Write time
    write(iunit) time_ckpt

    ! Write global_id
    if (nc_new > 0) then
      write(iunit) tmp_mesh%global_id(1:nc_new)
    end if

    ! Write primary fields (gather from global arrays)
    allocate(U_r(3, nc_new), p_r(nc_new), T_r(nc_new), nu_sgs_r(nc_new))
    if (nc_new > 0) then
      do c = 1, nc_new
        if (tmp_mesh%global_id(c) < 1 .or. tmp_mesh%global_id(c) > nc_global) then
          write(*, '(A,I0,A,I0,A,I0)') 'ERROR: new rank ', r, ' cell ', c, &
                                        ' has out-of-range global_id=', int(tmp_mesh%global_id(c))
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        U_r(:, c) = U_global(:, tmp_mesh%global_id(c))
        p_r(c) = p_global(tmp_mesh%global_id(c))
        T_r(c) = T_global(tmp_mesh%global_id(c))
        nu_sgs_r(c) = nu_sgs_global(tmp_mesh%global_id(c))
      end do
      write(iunit) U_r
      write(iunit) p_r
      write(iunit) T_r
      write(iunit) nu_sgs_r
    end if

    ! Write scalars
    if (n_scalars_ckpt > 0 .and. nc_new > 0) then
      allocate(C_r(n_scalars_ckpt, nc_new))
      do c = 1, nc_new
        C_r(:, c) = C_global(:, tmp_mesh%global_id(c))
      end do
      write(iunit) C_r
      deallocate(C_r)
    end if

    ! Write history
    if (history_ckpt > 0 .and. nc_new > 0) then
      allocate(H_U_n1_r(3, nc_new), H_U_n2_r(3, nc_new))
      allocate(H_T_n1_r(nc_new), H_T_n2_r(nc_new))
      do c = 1, nc_new
        H_U_n1_r(:, c) = H_U_n1_global(:, tmp_mesh%global_id(c))
        H_U_n2_r(:, c) = H_U_n2_global(:, tmp_mesh%global_id(c))
        H_T_n1_r(c) = H_T_n1_global(tmp_mesh%global_id(c))
        H_T_n2_r(c) = H_T_n2_global(tmp_mesh%global_id(c))
      end do
      write(iunit) H_U_n1_r
      write(iunit) H_U_n2_r
      write(iunit) H_T_n1_r
      write(iunit) H_T_n2_r
      deallocate(H_U_n1_r, H_U_n2_r, H_T_n1_r, H_T_n2_r)

      if (n_scalars_ckpt > 0) then
        allocate(H_C_n1_r(n_scalars_ckpt, nc_new))
        allocate(H_C_n2_r(n_scalars_ckpt, nc_new))
        do c = 1, nc_new
          H_C_n1_r(:, c) = H_C_n1_global(:, tmp_mesh%global_id(c))
          H_C_n2_r(:, c) = H_C_n2_global(:, tmp_mesh%global_id(c))
        end do
        write(iunit) H_C_n1_r
        write(iunit) H_C_n2_r
        deallocate(H_C_n1_r, H_C_n2_r)
      end if
    end if

    deallocate(U_r, p_r, T_r, nu_sgs_r)
    close(iunit)
    call tmp_mesh%destroy()

    if (mod(r, max(1, new_nranks/10)) == 0 .or. r == new_nranks - 1) &
      write(*, '(A,I0,A,I0)') '   Rank ', r, ': wrote ', nc_new, ' cells'
  end do

  ! Write checkpoint_latest
  checkpoint_latest_file = trim(new_dir)//'/checkpoint_latest'
  open(newunit=iunit, file=trim(checkpoint_latest_file), form='formatted', &
       action='write', status='replace', iostat=ios)
  if (ios /= 0) then
    write(*, '(A)') 'WARNING: cannot write checkpoint_latest'
  else
    write(iunit, '(I0)') step_ckpt
    close(iunit)
  end if

  write(*, '(A)') ' Step 2 complete'

  ! Clean up global arrays
  deallocate(U_global, p_global, T_global, nu_sgs_global)
  if (n_scalars_ckpt > 0) deallocate(C_global)
  if (history_ckpt > 0) then
    deallocate(H_U_n1_global, H_U_n2_global, H_T_n1_global, H_T_n2_global)
    if (n_scalars_ckpt > 0) deallocate(H_C_n1_global, H_C_n2_global)
  end if

  write(*, '(A)') '================================================='
  write(*, '(A)') ' Repartitioning complete'
  write(*, '(A)') '================================================='

  call MPI_Finalize(ierr)

end program dopamine_checkpoint_partition
