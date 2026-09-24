! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Rank -> GPU binding. One MPI rank per device: node-local rank % ndev, shared by the OpenACC
!> runtime and PETSc's CUDA backend so a rank's kernels and its pressure solve use the same GPU.
module mod_device
  use iso_c_binding, only: c_int
  use mod_kinds, only: ip
  use mpi_f08
#ifdef _OPENACC
  use openacc
#endif
  implicit none
  private

  public :: device_bind

  !> Device this rank is bound to (-1 = host-only build or no device found).
  integer(ip), public, protected :: device_id = -1
  !> Rank and size within the node (MPI_COMM_TYPE_SHARED split of the world communicator).
  integer(ip), public, protected :: device_local_rank = 0
  integer(ip), public, protected :: device_local_size = 1

  interface
    subroutine fvm_petsc_set_device_local_rank(local_rank) bind(C)
      import :: c_int
      integer(c_int), value :: local_rank
    end subroutine
  end interface

contains

  !> Bind this rank to a GPU. Call right after MPI_Init, before any !$acc region and before
  !> PETSc is initialised (PETSc reads its device choice once, at PetscInitialize).
  subroutine device_bind(comm)
    type(MPI_Comm), intent(in) :: comm

    type(MPI_Comm) :: node_comm
    integer(ip) :: rank, nranks, ndev, ierr, max_local_size
    integer(ip), allocatable :: all_dev(:), all_local(:)

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, node_comm, ierr)
    call MPI_Comm_rank(node_comm, device_local_rank, ierr)
    call MPI_Comm_size(node_comm, device_local_size, ierr)
    call MPI_Comm_free(node_comm, ierr)

    call fvm_petsc_set_device_local_rank(int(device_local_rank, c_int))

    ndev = 0
#ifdef _OPENACC
    ndev = acc_get_num_devices(acc_device_nvidia)
    if (ndev > 0) then
      device_id = mod(device_local_rank, ndev)
      call acc_set_device_num(device_id, acc_device_nvidia)
      call acc_init(acc_device_nvidia)
    end if
#endif

    allocate(all_dev(nranks), all_local(nranks))
    call MPI_Gather(device_id, 1, MPI_INTEGER, all_dev, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Gather(device_local_rank, 1, MPI_INTEGER, all_local, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Allreduce(device_local_size, max_local_size, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

    if (rank == 0) then
#ifdef _OPENACC
      write(*, '(A,I0,A,I0,A)') ' GPU binding  : ', ndev, ' device(s) per node, ', &
                                max_local_size, ' rank(s) per node'
      if (nranks <= 16) then
        block
          integer(ip) :: r
          do r = 1, nranks
            write(*, '(A,I0,A,I0,A,I0)') '   rank ', r-1, ' (local ', all_local(r), ') -> device ', all_dev(r)
          end do
        end block
      end if
      if (ndev > 0 .and. max_local_size > ndev) &
        write(*, '(A,I0,A,I0,A)') ' WARNING: ', max_local_size, ' ranks per node share ', ndev, &
          ' device(s); GPUs are oversubscribed (expected one rank per GPU)'
      if (ndev == 0) write(*, '(A)') ' WARNING: OpenACC build but no NVIDIA device found; kernels run on host'
#endif
    end if
    deallocate(all_dev, all_local)
  end subroutine device_bind

end module mod_device
