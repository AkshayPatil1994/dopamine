! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Binary mesh I/O module
module mod_mesh_io
  use mod_kinds,  only: rp, ip, lp
  use mod_mesh,   only: Mesh
  implicit none
  private
  public :: write_local_mesh_file, read_local_mesh_file

contains

  !> Write a local mesh to a Fortran unformatted binary file.
  !> Record layout mirrors send_local_mesh protocol from mod_partitioner.f90.
  subroutine write_local_mesh_file(msh, filename)
    type(Mesh),       intent(in) :: msh
    character(len=*), intent(in) :: filename

    integer(ip) :: iunit, ios, nc, ng, nf, nif, nbf, np, nv, nh, i, format_version
    integer(ip) :: header(13), total_cell_verts, total_face_verts
    real(rp),    allocatable :: centroid(:,:), volume(:), area(:), &
                                centre(:,:), normal_a(:,:), sf(:,:), dcf(:), gc(:)
    integer(ip), allocatable :: nfaces_a(:), offsets(:), &
                                owner(:), nei(:), patch_id(:), &
                                bcid(:), start(:), pnf(:)
    character(len=128), allocatable :: names(:)
    character(len=64),  allocatable :: types(:)

    nc  = msh%n_cells;  ng  = msh%n_ghost
    nf  = msh%n_faces;  nif = msh%n_int_faces;  nbf = msh%n_bnd_faces
    np  = msh%n_patches; nv  = msh%n_vertices
    nh  = 0;  if (allocated(msh%halos)) nh = size(msh%halos)
    ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()
    ! must never be evaluated in the same expression as the allocated() check.
    total_cell_verts = 0
    if (allocated(msh%cell_offsets)) then
      if (size(msh%cell_offsets) > nc) total_cell_verts = msh%cell_offsets(nc + 1) - 1
    end if
    total_face_verts = 0
    if (allocated(msh%face_offsets)) then
      if (size(msh%face_offsets) > nf) total_face_verts = msh%face_offsets(nf + 1) - 1
    end if

    header(1)  = nc;  header(2)  = ng;   header(3)  = nf
    header(4)  = nif; header(5)  = nbf;  header(6)  = np
    header(7)  = nv;  header(8)  = int(msh%n_cells_global, ip)
    header(9)  = 0;   if (allocated(msh%vertices)) header(9) = 1
    header(10) = 0;   if (allocated(msh%face_offsets)) header(10) = 1
    header(11) = nh;  header(12) = total_cell_verts;  header(13) = total_face_verts

    ! Format marker in header(13): v1 uses non-negative total_face_verts.
    ! v2 stores header(13)=-(total_face_verts+1) and includes periodic_pair_dcf.
    format_version = 2_ip
    if (format_version >= 2_ip) header(13) = - (total_face_verts + 1_ip)

    open(newunit=iunit, file=trim(filename), form='unformatted', access='stream', &
         action='write', status='replace', iostat=ios)
    if (ios /= 0) error stop 'write_local_mesh_file: cannot open '//trim(filename)

    ! Record 1: header
    write(iunit) header

    ! Records 2-5: cell data (if nc+ng > 0)
    if (nc + ng > 0) then
      allocate(centroid(3,nc+ng), volume(nc+ng), nfaces_a(nc+ng), offsets(nc+ng))
      do i = 1, nc+ng
        centroid(:,i) = msh%cells(i)%centroid
        volume(i) = msh%cells(i)%volume
        nfaces_a(i) = msh%cells(i)%n_faces
        offsets(i) = msh%cells(i)%face_offset
      end do
      write(iunit) centroid;  write(iunit) volume
      write(iunit) nfaces_a;  write(iunit) offsets
      deallocate(centroid, volume, nfaces_a, offsets)
    end if

    ! Records 6-14: face data (if nf > 0)
    if (nf > 0) then
      allocate(owner(nf), nei(nf), patch_id(nf), area(nf), &
               centre(3,nf), normal_a(3,nf), sf(3,nf), dcf(nf), gc(nf))
      do i = 1, nf
        owner(i) = msh%faces(i)%owner; nei(i) = msh%faces(i)%neighbour
        patch_id(i) = msh%faces(i)%patch_id; area(i) = msh%faces(i)%area
        centre(:,i) = msh%faces(i)%centre; normal_a(:,i) = msh%faces(i)%normal
        sf(:,i) = msh%faces(i)%sf; dcf(i) = msh%faces(i)%d_cf; gc(i) = msh%faces(i)%gc
      end do
      write(iunit) owner;  write(iunit) nei;  write(iunit) patch_id
      write(iunit) area;   write(iunit) centre; write(iunit) normal_a
      write(iunit) sf;     write(iunit) dcf;    write(iunit) gc
      deallocate(owner, nei, patch_id, area, centre, normal_a, sf, dcf, gc)
    end if

    ! Records 15-19: patch data (if np > 0)
    if (np > 0) then
      allocate(bcid(np), start(np), pnf(np), names(np), types(np))
      do i = 1, np
        bcid(i) = msh%patches(i)%bc_id; start(i) = msh%patches(i)%start_face
        pnf(i) = msh%patches(i)%n_faces
        names(i) = msh%patches(i)%name; types(i) = msh%patches(i)%bc_type
      end do
      write(iunit) bcid;   write(iunit) start;  write(iunit) pnf
      write(iunit) names;  write(iunit) types
      deallocate(bcid, start, pnf, names, types)
    end if

    ! Records 20-21: periodic pairing (if nbf > 0)
    if (nbf > 0) then
      if (allocated(msh%periodic_pair)) then
        write(iunit) msh%periodic_pair(1:nbf)
      else
        write(iunit) [(0_ip, i=1,nbf)]
      end if
      if (allocated(msh%periodic_neighbor)) then
        write(iunit) msh%periodic_neighbor(1:nbf)
      else
        write(iunit) [(0_ip, i=1,nbf)]
      end if
      if (allocated(msh%periodic_pair_dcf)) then
        write(iunit) msh%periodic_pair_dcf(1:nbf)
      else if (allocated(msh%faces)) then
        write(iunit) [(msh%faces(msh%n_int_faces + i)%d_cf, i=1,nbf)]
      else
        write(iunit) [(0.0_rp, i=1,nbf)]
      end if
    end if

    ! Record 23: global_id (if nc+ng > 0)
    if (nc + ng > 0) then
      if (allocated(msh%global_id)) then
        write(iunit) msh%global_id(1:nc+ng)
      else
        write(iunit) [(int(i, lp), i=1,nc+ng)]
      end if
    end if

    ! Records 24-27: vertex topology (if has_vertices==1)
    if (header(9) == 1 .and. nv > 0) then
      write(iunit) msh%vertices(1:3, 1:nv)
      write(iunit) msh%cell_offsets(1:nc+1)
      write(iunit) msh%cell_vertices(1:total_cell_verts)
      write(iunit) msh%cell_types(1:nc)
    end if

    ! Records 28-29: face topology (if has_face_offsets==1)
    if (header(10) == 1 .and. nf > 0) then
      write(iunit) msh%face_offsets(1:nf+1)
      write(iunit) msh%face_vertices(1:total_face_verts)
    end if

    ! Record 30: cell_faces size (always written)
    if (allocated(msh%cell_faces)) then
      write(iunit) int(size(msh%cell_faces), ip)
      write(iunit) msh%cell_faces
    else
      write(iunit) 0_ip
    end if

    ! Records 32+: halos
    if (nh > 0) then
      do i = 1, nh
        write(iunit) [msh%halos(i)%rank, msh%halos(i)%n_send, msh%halos(i)%n_recv]
        if (msh%halos(i)%n_send > 0) write(iunit) msh%halos(i)%send_cells
        if (msh%halos(i)%n_recv > 0) write(iunit) msh%halos(i)%recv_cells
      end do
    end if

    ! Trailing record (absent in older files): vertex-adjacent owned/ghost pairs for the LS stencil.
    write(iunit) msh%n_vx
    if (msh%n_vx > 0) then
      write(iunit) msh%vx_owned(1:msh%n_vx)
      write(iunit) msh%vx_ghost(1:msh%n_vx)
    end if

    close(iunit)
  end subroutine write_local_mesh_file


  !> Read a local mesh from a Fortran unformatted binary file.
  !> Mirrors recv_local_mesh protocol, replacing MPI_Recv with file reads.
  subroutine read_local_mesh_file(msh, filename)
    type(Mesh),       intent(inout) :: msh
    character(len=*), intent(in)    :: filename

    integer(ip) :: iunit, ios, nc, ng, nf, nif, nbf, np, nv, nh
    integer(ip) :: header(13), total_cell_verts, total_face_verts, ncf, format_version
    integer(ip) :: i, rank_halo, n_send, n_recv, n_vx_read
    real(rp),    allocatable :: centroid(:,:), volume(:), area(:), &
                                centre(:,:), normal_a(:,:), sf(:,:), dcf(:), gc(:)
    integer(ip), allocatable :: nfaces_a(:), offsets(:), &
                                owner(:), nei(:), patch_id(:), &
                                bcid(:), start(:), pnf(:)
    character(len=128), allocatable :: names(:)
    character(len=64),  allocatable :: types(:)
    integer(lp), allocatable :: global_id(:)

    open(newunit=iunit, file=trim(filename), form='unformatted', access='stream', &
         action='read', status='old', iostat=ios)
    if (ios /= 0) error stop 'read_local_mesh_file: cannot open '//trim(filename)

    ! Record 1: header
    read(iunit) header

    nc  = header(1);  ng  = header(2);   nf  = header(3)
    nif = header(4);  nbf = header(5);   np  = header(6)
    nv  = header(7)
    msh%n_cells = nc;  msh%n_ghost = ng
    msh%n_faces = nf;  msh%n_int_faces = nif;  msh%n_bnd_faces = nbf
    msh%n_patches = np; msh%n_vertices = nv
    msh%n_cells_global = int(header(8), lp)
    nh = header(11)
    total_cell_verts = header(12)
    if (header(13) < 0_ip) then
      format_version = 2_ip
      total_face_verts = -header(13) - 1_ip
    else
      format_version = 1_ip
      total_face_verts = header(13)
    end if

    ! Allocate cell and face arrays
    if (.not. allocated(msh%cells)) allocate(msh%cells(nc+ng))
    if (.not. allocated(msh%faces)) allocate(msh%faces(nf))
    if (.not. allocated(msh%patches)) allocate(msh%patches(np))

    ! Records 2-5: cell data (if nc+ng > 0)
    if (nc + ng > 0) then
      allocate(centroid(3,nc+ng), volume(nc+ng), nfaces_a(nc+ng), offsets(nc+ng))
      read(iunit) centroid;  read(iunit) volume
      read(iunit) nfaces_a;  read(iunit) offsets
      do i = 1, nc+ng
        msh%cells(i)%centroid = centroid(:,i)
        msh%cells(i)%volume = volume(i)
        msh%cells(i)%n_faces = nfaces_a(i)
        msh%cells(i)%face_offset = offsets(i)
      end do
      deallocate(centroid, volume, nfaces_a, offsets)
    end if

    ! Records 6-14: face data (if nf > 0)
    if (nf > 0) then
      allocate(owner(nf), nei(nf), patch_id(nf), area(nf), &
               centre(3,nf), normal_a(3,nf), sf(3,nf), dcf(nf), gc(nf))
      read(iunit) owner;  read(iunit) nei;  read(iunit) patch_id
      read(iunit) area;   read(iunit) centre; read(iunit) normal_a
      read(iunit) sf;     read(iunit) dcf;    read(iunit) gc
      do i = 1, nf
        msh%faces(i)%owner = owner(i); msh%faces(i)%neighbour = nei(i)
        msh%faces(i)%patch_id = patch_id(i); msh%faces(i)%area = area(i)
        msh%faces(i)%centre = centre(:,i); msh%faces(i)%normal = normal_a(:,i)
        msh%faces(i)%sf = sf(:,i); msh%faces(i)%d_cf = dcf(i); msh%faces(i)%gc = gc(i)
      end do
      deallocate(owner, nei, patch_id, area, centre, normal_a, sf, dcf, gc)
    end if

    ! Records 15-19: patch data (if np > 0)
    if (np > 0) then
      allocate(bcid(np), start(np), pnf(np), names(np), types(np))
      read(iunit) bcid;   read(iunit) start;  read(iunit) pnf
      read(iunit) names;  read(iunit) types
      do i = 1, np
        msh%patches(i)%bc_id = bcid(i); msh%patches(i)%start_face = start(i)
        msh%patches(i)%n_faces = pnf(i)
        msh%patches(i)%name = names(i); msh%patches(i)%bc_type = types(i)
      end do
      deallocate(bcid, start, pnf, names, types)
    end if

    ! Records 20-22: periodic pairing (if nbf > 0)
    if (nbf > 0) then
      allocate(msh%periodic_pair(nbf), msh%periodic_neighbor(nbf), msh%periodic_pair_dcf(nbf))
      read(iunit) msh%periodic_pair
      read(iunit) msh%periodic_neighbor
      if (format_version >= 2_ip) then
        read(iunit) msh%periodic_pair_dcf
      else
        do i = 1, nbf
          if (msh%periodic_pair(i) >= 1 .and. msh%periodic_pair(i) <= nbf) then
            msh%periodic_pair_dcf(i) = msh%faces(msh%n_int_faces + msh%periodic_pair(i))%d_cf
          else
            msh%periodic_pair_dcf(i) = msh%faces(msh%n_int_faces + i)%d_cf
          end if
        end do
      end if
    end if

    ! Record 23: global_id (if nc+ng > 0)
    if (nc + ng > 0) then
      allocate(global_id(nc+ng))
      read(iunit) global_id
      if (.not. allocated(msh%global_id)) allocate(msh%global_id(nc+ng))
      msh%global_id = global_id
      deallocate(global_id)
    end if

    ! Records 24-27: vertex topology (if header(9)==1)
    if (header(9) == 1 .and. nv > 0) then
      allocate(msh%vertices(3,nv))
      allocate(msh%cell_offsets(nc+1))
      allocate(msh%cell_vertices(total_cell_verts))
      allocate(msh%cell_types(nc))
      read(iunit) msh%vertices
      read(iunit) msh%cell_offsets
      read(iunit) msh%cell_vertices
      read(iunit) msh%cell_types
    end if

    ! Records 28-29: face topology (if header(10)==1)
    if (header(10) == 1 .and. nf > 0) then
      allocate(msh%face_offsets(nf+1))
      allocate(msh%face_vertices(total_face_verts))
      read(iunit) msh%face_offsets
      read(iunit) msh%face_vertices
    end if

    ! Record 30: cell_faces
    read(iunit) ncf
    if (ncf > 0) then
      allocate(msh%cell_faces(ncf))
      read(iunit) msh%cell_faces
    end if

    ! Records 32+: halos
    if (nh > 0) then
      allocate(msh%halos(nh))
      do i = 1, nh
        read(iunit) rank_halo, n_send, n_recv
        msh%halos(i)%rank = rank_halo
        msh%halos(i)%n_send = n_send
        msh%halos(i)%n_recv = n_recv
        if (n_send > 0) then
          allocate(msh%halos(i)%send_cells(n_send))
          read(iunit) msh%halos(i)%send_cells
        end if
        if (n_recv > 0) then
          allocate(msh%halos(i)%recv_cells(n_recv))
          read(iunit) msh%halos(i)%recv_cells
        end if
      end do
    end if

    ! Trailing vertex-pair record; files written before it existed end here.
    msh%n_vx = 0
    read(iunit, iostat=ios) n_vx_read
    if (ios == 0 .and. n_vx_read > 0) then
      allocate(msh%vx_owned(n_vx_read), msh%vx_ghost(n_vx_read))
      read(iunit) msh%vx_owned
      read(iunit) msh%vx_ghost
      msh%n_vx = n_vx_read
    end if

    close(iunit)
  end subroutine read_local_mesh_file

end module mod_mesh_io
