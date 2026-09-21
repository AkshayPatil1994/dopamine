! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Boundary condition configuration and patch setup
module mod_bc_setup
  use mod_kinds
  use mod_parameters
  use mod_scalar
  use mod_mesh
  use mod_bc, only: BCConfig
  use mod_input_parser, only: InputConfig
  implicit none
  private

  public :: build_bc_configs

contains

  subroutine build_bc_configs(config, msh, bc_configs_out)
    type(InputConfig), intent(in)  :: config
    type(Mesh),        intent(inout)  :: msh
    type(BCConfig), allocatable, intent(out) :: bc_configs_out(:)

    integer(ip) :: p, q, n_missing
    logical :: matched
    logical, allocatable :: bc_used(:)

    allocate(bc_configs_out(msh%n_patches))
    allocate(bc_used(config%n_bcs), source=.false.)

    bc_configs_out%bc_id = BC_ZEROGRADIENT
    bc_configs_out%T_val = config%params%T_ref
    if (config%n_scalars > 0) then
      do p = 1, msh%n_patches
        allocate(bc_configs_out(p)%scalar_bc_mode(config%n_scalars), source=SCALAR_BC_ZERO_GRADIENT)
        allocate(bc_configs_out(p)%scalar_value(config%n_scalars), source=0.0_rp)
        allocate(bc_configs_out(p)%scalar_flux(config%n_scalars), source=0.0_rp)
        allocate(bc_configs_out(p)%scalar_flux_omega(config%n_scalars), source=0.0_rp)
        allocate(bc_configs_out(p)%scalar_flux_phase(config%n_scalars), source=0.0_rp)
      end do
    end if

    do p = 1, msh%n_patches
      matched = .false.
      do q = 1, config%n_bcs
        if (trim(msh%patches(p)%name) == trim(config%bc_patch_names(q))) then
          bc_configs_out(p) = config%bc_configs(q)
          msh%patches(p)%bc_id = bc_configs_out(p)%bc_id
          matched = .true.
          bc_used(q) = .true.
          exit
        end if
      end do

      if (.not. matched .and. msh%patches(p)%bc_id /= 0) then
        bc_configs_out(p)%bc_id = msh%patches(p)%bc_id
      end if
      if (.not. matched .and. msh%rank == 0) &
        write(*, '(A,A,A,I0,A)') ' WARNING: mesh patch "', trim(msh%patches(p)%name), &
          '" is not listed in &dopamine_bc; using default bc_id=', &
          merge(msh%patches(p)%bc_id, BC_ZEROGRADIENT, msh%patches(p)%bc_id /= 0), ' (open boundary if 9)'

      msh%patches(p)%bc_id = bc_configs_out(p)%bc_id
    end do

    block
      use mpi_f08
      integer(ip), allocatable :: nf_local(:), nf_global(:)
      integer :: ierr
      allocate(nf_local(msh%n_patches), nf_global(msh%n_patches))
      do p = 1, msh%n_patches
        nf_local(p) = msh%patches(p)%n_faces
      end do
      if (msh%nranks > 1) then
        call MPI_Allreduce(nf_local, nf_global, int(msh%n_patches), MPI_INTEGER, MPI_SUM, msh%comm, ierr)
      else
        nf_global = nf_local
      end if
      if (msh%rank == 0) then
        do p = 1, msh%n_patches
          if (nf_global(p) == 0) then
            write(*, '(A,A,A)') '   WARNING: patch "', trim(msh%patches(p)%name), '" has zero faces globally'
          end if
        end do
      end if
      deallocate(nf_local, nf_global)
    end block

    if (msh%rank == 0) then
      n_missing = 0
      do q = 1, config%n_bcs
        if (.not. bc_used(q)) then
          n_missing = n_missing + 1
          write(*, '(A,A)') 'WARNING: input BC patch not found in mesh: ', trim(config%bc_patch_names(q))
        end if
      end do
      if (n_missing > 0) then
        error stop 'build_bc_configs: one or more input BC patch names do not exist in the mesh'
      end if
    end if

    call build_periodic_pairs(msh, bc_configs_out)
    block
      use mpi_f08
      integer(ip) :: n_local_pair, n_global_pair, n_local_nbr, n_global_nbr
      integer :: ierr
      n_local_pair = 0
      n_local_nbr  = 0
      if (allocated(msh%periodic_pair))    n_local_pair = count(msh%periodic_pair > 0_ip)
      if (allocated(msh%periodic_neighbor)) n_local_nbr = count(msh%periodic_neighbor > 0_ip)
      if (msh%nranks > 1) then
        call MPI_Reduce(n_local_pair, n_global_pair, 1, MPI_INTEGER, MPI_SUM, 0, msh%comm, ierr)
        call MPI_Reduce(n_local_nbr,  n_global_nbr,  1, MPI_INTEGER, MPI_SUM, 0, msh%comm, ierr)
      else
        n_global_pair = n_local_pair
        n_global_nbr  = n_local_nbr
      end if
      if (msh%rank == 0 .and. n_global_pair > 0) then
        write(*, '(2X,I0,A,I0,A)') msh%n_patches, ' patches  |  ', n_global_pair, ' periodic face pairs'
      end if
    end block
    deallocate(bc_used)
  end subroutine build_bc_configs

  subroutine build_periodic_pairs(msh, bc_configs)
    type(Mesh),     intent(inout) :: msh
    type(BCConfig), intent(in)    :: bc_configs(:)

    integer(ip), allocatable :: periodic_patches(:), map(:), candidate_map(:)
    logical, allocatable :: patch_used(:)
    integer(ip) :: p, q, i, np, best_q
    real(rp) :: best_err, err
    logical :: ok

    ! In MPI runs the partitioner builds periodic_pair/periodic_neighbor.
    ! Skip geometric matching on already-partitioned meshes (n_cells_global > 0)
    if (msh%nranks > 1 .and. msh%n_cells_global > 0_lp) then
      if (.not. allocated(msh%periodic_pair) .and. msh%n_bnd_faces > 0) then
        allocate(msh%periodic_pair(msh%n_bnd_faces), source=0_ip)
      end if
      if (.not. allocated(msh%periodic_pair_dcf) .and. msh%n_bnd_faces > 0) then
        allocate(msh%periodic_pair_dcf(msh%n_bnd_faces), source=0.0_rp)
      end if
      call finalize_periodic_partner_metrics(msh)
      call validate_periodic_mapping_mpi(msh, bc_configs)
      return
    end if

    if (allocated(msh%periodic_pair)) deallocate(msh%periodic_pair)
    if (allocated(msh%periodic_pair_dcf)) deallocate(msh%periodic_pair_dcf)
    if (msh%n_bnd_faces <= 0) return
    allocate(msh%periodic_pair(msh%n_bnd_faces), source=0_ip)
    allocate(msh%periodic_pair_dcf(msh%n_bnd_faces), source=0.0_rp)

    np = count([(bc_configs(i)%bc_id == BC_PERIODIC, i=1,msh%n_patches)])
    if (np == 0) return

    allocate(periodic_patches(np))
    allocate(patch_used(msh%n_patches), source=.false.)
    q = 0
    do p = 1, msh%n_patches
      if (bc_configs(p)%bc_id == BC_PERIODIC) then
        q = q + 1
        periodic_patches(q) = p
      end if
    end do

    do i = 1, np
      p = periodic_patches(i)
      if (patch_used(p)) cycle

      best_q = 0
      best_err = huge(1.0_rp)
      do q = i + 1, np
        if (patch_used(periodic_patches(q))) cycle
        call try_pair_patches(msh, p, periodic_patches(q), ok, err, candidate_map)
        if (ok .and. err < best_err) then
          ! Found a better match; save it
          if (allocated(map)) deallocate(map)
          allocate(map(size(candidate_map)))
          map = candidate_map
          best_err = err
          best_q = periodic_patches(q)
        end if
        if (allocated(candidate_map)) deallocate(candidate_map)
      end do

      if (best_q == 0) then
        if (allocated(map)) deallocate(map)
        error stop 'build_periodic_pairs: unable to find periodic counterpart for patch '//trim(msh%patches(p)%name)
      end if

      ! Use map from best match found in search loop
      call store_patch_pairing(msh, p, best_q, map)
      patch_used(p) = .true.
      patch_used(best_q) = .true.
      deallocate(map)
    end do

    call finalize_periodic_partner_metrics(msh)

    deallocate(periodic_patches, patch_used)
  end subroutine build_periodic_pairs

  subroutine try_pair_patches(msh, p, q, ok, max_err, map)
    type(Mesh), intent(in) :: msh
    integer(ip), intent(in) :: p, q
    logical, intent(out) :: ok
    real(rp), intent(out) :: max_err
    integer(ip), allocatable, intent(out) :: map(:)

    integer(ip) :: n, i, fbp, fbq, fgp, fgq
    real(rp), allocatable :: cp(:,:), cq(:,:), np(:,:), nq(:,:)
    real(rp) :: shift(3), avg_dot, scale, tol, dist
    real(rp), allocatable :: key_p(:), key_q(:)
    integer(ip), allocatable :: ord_p(:), ord_q(:)
    ! Irrational multipliers: minimise probability of key collisions for any mesh geometry.
    real(rp), parameter :: KC1 = 3.14159265358979_rp  ! pi
    real(rp), parameter :: KC2 = 2.71828182845904_rp  ! e
    real(rp), parameter :: KC3 = 1.41421356237310_rp  ! sqrt(2)

    ok = .false.
    max_err = huge(1.0_rp)

    n = msh%patches(p)%n_faces
    if (n == 0) return
    if (n /= msh%patches(q)%n_faces) return

    allocate(cp(3,n), cq(3,n), np(3,n), nq(3,n))
    cp = 0.0_rp; cq = 0.0_rp; np = 0.0_rp; nq = 0.0_rp

    do i = 1, n
      fbp = msh%patches(p)%start_face + i - 1
      fgp = msh%n_int_faces + fbp
      cp(:,i) = msh%faces(fgp)%centre
      np(:,i) = msh%faces(fgp)%normal

      fbq = msh%patches(q)%start_face + i - 1
      fgq = msh%n_int_faces + fbq
      cq(:,i) = msh%faces(fgq)%centre
      nq(:,i) = msh%faces(fgq)%normal
    end do

    shift = sum(cq, dim=2) / real(n, rp) - sum(cp, dim=2) / real(n, rp)
    avg_dot = dot_product(sum(np, dim=2), sum(nq, dim=2))
    scale = max(1.0_rp, maxval(abs(cp)), maxval(abs(cq)))
    tol = 1.0e-5_rp * scale

    ! Sort both face sets by the same scalar key; paired by rank -> O(n log n) matching.
    allocate(key_p(n), key_q(n), ord_p(n), ord_q(n))
    do i = 1, n
      key_p(i) = KC1*(cp(1,i)+shift(1)) + KC2*(cp(2,i)+shift(2)) + KC3*(cp(3,i)+shift(3))
      key_q(i) = KC1*cq(1,i)            + KC2*cq(2,i)            + KC3*cq(3,i)
      ord_p(i) = i
      ord_q(i) = i
    end do
    call argsort_rp(key_p, ord_p)
    call argsort_rp(key_q, ord_q)

    allocate(map(n))
    max_err = 0.0_rp
    do i = 1, n
      map(ord_p(i)) = ord_q(i)
      dist = sqrt(sum((cq(:, ord_q(i)) - (cp(:, ord_p(i)) + shift))**2))
      max_err = max(max_err, dist)
    end do

    deallocate(cp, cq, np, nq, key_p, key_q, ord_p, ord_q)

    if (max_err <= tol .and. avg_dot < 0.0_rp) ok = .true.
  end subroutine try_pair_patches

  subroutine store_patch_pairing(msh, p, q, map)
    type(Mesh), intent(inout) :: msh
    integer(ip), intent(in) :: p, q
    integer(ip), intent(in) :: map(:)

    integer(ip) :: i, fbp, fbq

    do i = 1, size(map)
      fbp = msh%patches(p)%start_face + i - 1
      fbq = msh%patches(q)%start_face + map(i) - 1
      msh%periodic_pair(fbp) = fbq
      msh%periodic_pair(fbq) = fbp
      if (allocated(msh%periodic_pair_dcf)) then
        msh%periodic_pair_dcf(fbp) = msh%faces(msh%n_int_faces + fbq)%d_cf
        msh%periodic_pair_dcf(fbq) = msh%faces(msh%n_int_faces + fbp)%d_cf
      end if
    end do
  end subroutine store_patch_pairing

  subroutine finalize_periodic_partner_metrics(msh)
    type(Mesh), intent(inout) :: msh

    integer(ip) :: f_bnd, pair_bnd, pair_face, i

    if (msh%n_bnd_faces <= 0) return
    if (.not. allocated(msh%periodic_pair_dcf)) allocate(msh%periodic_pair_dcf(msh%n_bnd_faces), source=0.0_rp)

    do f_bnd = 1, msh%n_bnd_faces
      if (msh%periodic_pair_dcf(f_bnd) > 0.0_rp) cycle
      if (allocated(msh%periodic_pair)) then
        if (f_bnd <= size(msh%periodic_pair)) then
          pair_bnd = msh%periodic_pair(f_bnd)
          if (pair_bnd >= 1 .and. pair_bnd <= msh%n_bnd_faces) then
            pair_face = msh%n_int_faces + pair_bnd
            if (pair_face >= 1 .and. pair_face <= size(msh%faces)) then
              msh%periodic_pair_dcf(f_bnd) = msh%faces(pair_face)%d_cf
            end if
          end if
        end if
      end if
      if (msh%periodic_pair_dcf(f_bnd) <= 0.0_rp) then
        i = msh%n_int_faces + f_bnd
        if (i >= 1 .and. i <= size(msh%faces)) msh%periodic_pair_dcf(f_bnd) = msh%faces(i)%d_cf
      end if
    end do
  end subroutine finalize_periodic_partner_metrics

  subroutine validate_periodic_mapping_mpi(msh, bc_configs)    
    use mpi_f08

    type(Mesh),     intent(in) :: msh
    type(BCConfig), intent(in) :: bc_configs(:)

    integer(ip) :: q, f, f_bnd, f_glob, o, nb, pair_bnd, pair_face, ierr
    integer(ip) :: local_faces, local_resolved, local_unresolved, global_faces, global_resolved, global_unresolved
    integer(ip) :: local_missing_metric, global_missing_metric

    local_faces = 0
    local_resolved = 0
    local_unresolved = 0
    local_missing_metric = 0

    do q = 1, msh%n_patches
      if (bc_configs(q)%bc_id /= BC_PERIODIC) cycle
      do f = 1, msh%patches(q)%n_faces
        f_bnd = msh%patches(q)%start_face + f - 1
        f_glob = msh%n_int_faces + f_bnd
        if (f_glob < 1 .or. f_glob > size(msh%faces)) cycle
        o = msh%faces(f_glob)%owner
        if (o < 1) cycle
        local_faces = local_faces + 1

        nb = 0
        if (allocated(msh%periodic_neighbor)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_neighbor)) nb = msh%periodic_neighbor(f_bnd)
        end if
        if (nb <= 0 .and. allocated(msh%periodic_pair)) then
          if (f_bnd >= 1 .and. f_bnd <= size(msh%periodic_pair)) then
            pair_bnd = msh%periodic_pair(f_bnd)
            if (pair_bnd >= 1 .and. pair_bnd <= msh%n_bnd_faces) then
              pair_face = msh%n_int_faces + pair_bnd
              if (pair_face >= 1 .and. pair_face <= size(msh%faces)) nb = msh%faces(pair_face)%owner
            end if
          end if
        end if
        ! For preprocessed distributed meshes the decomposer embeds the periodic partner
        ! as a ghost cell in faces(f_glob)%neighbour; use it as the final fallback.
        if (nb <= 0) nb = msh%faces(f_glob)%neighbour

        if (nb >= 1 .and. nb <= msh%n_cells + msh%n_ghost .and. nb /= o) then
          local_resolved = local_resolved + 1
          if (.not. allocated(msh%periodic_pair_dcf)) then
            local_missing_metric = local_missing_metric + 1
          else if (f_bnd < 1 .or. f_bnd > size(msh%periodic_pair_dcf)) then
            local_missing_metric = local_missing_metric + 1
          else if (msh%periodic_pair_dcf(f_bnd) <= SMALL) then
            local_missing_metric = local_missing_metric + 1
          end if
        else
          local_unresolved = local_unresolved + 1
        end if
      end do
    end do

    call MPI_Allreduce(local_faces, global_faces, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_resolved, global_resolved, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_unresolved, global_unresolved, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)
    call MPI_Allreduce(local_missing_metric, global_missing_metric, 1, MPI_INTEGER, MPI_SUM, msh%comm, ierr)

    if (global_faces > 0 .and. global_unresolved > 0) then
      error stop 'validate_periodic_mapping_mpi: unresolved periodic neighbor mapping on distributed mesh'
    end if
    if (global_faces > 0 .and. global_missing_metric > 0) then
      error stop 'validate_periodic_mapping_mpi: periodic partner metric (periodic_pair_dcf) is missing'
    end if
  end subroutine validate_periodic_mapping_mpi

  !> Argsort real(rp) array in ascending order; on entry idx(i)=i, on exit idx is sorted permutation.
  subroutine argsort_rp(arr, idx)
    real(rp),    intent(in)    :: arr(:)
    integer(ip), intent(inout) :: idx(:)
    call qs_rp(arr, idx, 1_ip, int(size(arr), ip))
  end subroutine argsort_rp

  recursive subroutine qs_rp(arr, idx, lo, hi)
    real(rp),    intent(in)    :: arr(:)
    integer(ip), intent(inout) :: idx(:)
    integer(ip), intent(in)    :: lo, hi
    integer(ip) :: i, j, tmp
    real(rp)    :: pivot
    if (lo >= hi) return
    pivot = arr(idx((lo + hi) / 2))
    i = lo
    j = hi
    do while (i <= j)
      do while (arr(idx(i)) < pivot)
        i = i + 1
      end do
      do while (arr(idx(j)) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
        i = i + 1
        j = j - 1
      end if
    end do
    call qs_rp(arr, idx, lo, j)
    call qs_rp(arr, idx, i, hi)
  end subroutine qs_rp

end module mod_bc_setup
