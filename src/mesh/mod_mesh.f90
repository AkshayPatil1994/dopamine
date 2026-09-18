! SPDX-FileCopyrightText: 2026 DOPAMINE contributors
! SPDX-License-Identifier: AGPL-3.0-only
!> Core mesh data structures: cell-centred unstructured FV layout.
module mod_mesh
  use mod_kinds
  use mod_parameters
  use mpi_f08
  implicit none
  private
  public :: mesh_assign_compact_ids
  public :: mesh_ls_degenerate_count
  public :: refine_blend_weight

  ! Relative determinant below which an LS stencil is treated as rank-deficient.
  real(rp), parameter :: LS_DET_TOL  = 1.0e-12_rp
  real(rp), parameter :: LS_REG_FRAC = 1.0e-6_rp

  !> Cells whose LS gradient stencil needed regularisation; non-zero means sliver or coplanar cells.
  integer(ip), save :: ls_degenerate_count = 0_ip

  !> Face skewness above which a face is treated as a refinement interface; 2:1 gives 0.26, 1.5-graded gives 0.10.
  real(rp), parameter, public :: REFINE_SKEW_TOL = 0.15_rp

  ! Skewness-correction blend range: full correction below SKEW_BLEND_LO, none above SKEW_BLEND_HI.
  ! Straddles REFINE_SKEW_TOL so ordinary (non-refine) mesh skew keeps the full O(h^2) correction.
  real(rp), parameter, public :: SKEW_BLEND_LO = 0.15_rp
  real(rp), parameter, public :: SKEW_BLEND_HI = 0.25_rp

  integer(ip), parameter :: HALO_TAG_SCALAR = 1101_ip
  integer(ip), parameter :: HALO_TAG_VECTOR = 1102_ip
  integer(ip), parameter :: HALO_TAG_MULTI  = 1103_ip
  integer(ip), parameter :: HALO_TAG_STATE  = 1104_ip

  ! Face data type
  type, public :: FaceData
    integer(ip) :: owner       = -1  ! owning cell index (local)
    integer(ip) :: neighbour   = -1  ! neighbouring cell index (-1 = boundary)
    integer(ip) :: patch_id    = -1  ! boundary patch index (-1 = internal)
    real(rp)    :: area        = 0.0_rp
    real(rp)    :: centre(3)   = 0.0_rp
    real(rp)    :: normal(3)   = 0.0_rp  ! unit outward normal (owner -> neighbour)
    real(rp)    :: sf(3)       = 0.0_rp  ! area-weighted normal = area * normal
    real(rp)    :: d_cf        = 0.0_rp  ! |owner_centroid -> neighbour_centroid|
    real(rp)    :: gc          = 0.0_rp  ! geometric interpolation weight (owner side)
    real(rp)    :: skew        = 0.0_rp  ! |face centre - centroid midpoint| / |d_on|
    logical     :: is_refine   = .false.  ! refinement-interface face (skew above threshold)
    real(rp)    :: eps(3)      = 0.0_rp  ! skewness offset: centre - (owner + (1-gc)*d_on); interior faces only
  end type FaceData

  ! Cell data type
  type, public :: CellData
    real(rp)    :: centroid(3) = 0.0_rp
    real(rp)    :: volume      = 0.0_rp
    integer(ip) :: n_faces     = 0
    ! Faces are stored in a compressed adjacency list (see Mesh%cell_faces)
    integer(ip) :: face_offset = 0  ! index into cell_faces array
  end type CellData

  ! Boundary patch
  type, public :: PatchData
    character(len=128) :: name      = ''
    character(len=64)  :: bc_type   = ''
    integer(ip)        :: bc_id     = 0  ! BC_* constant
    integer(ip)        :: start_face = 0  ! first face index in patch
    integer(ip)        :: n_faces    = 0
  end type PatchData

  ! MPI halo descriptor
  type, public :: HaloData
    integer(ip)              :: rank  ! remote MPI rank
    integer(ip), allocatable :: send_cells(:)  ! local cells to send
    integer(ip), allocatable :: recv_cells(:)  ! local ghost cell indices
    integer(ip)              :: n_send = 0
    integer(ip)              :: n_recv = 0
  end type HaloData

  ! Main mesh type
  type, public :: Mesh
    ! -- dimensions
    integer(ip) :: n_cells       = 0  ! local interior cells
    integer(ip) :: n_ghost       = 0  ! ghost (halo) cells
    integer(ip) :: n_faces       = 0  ! all faces (internal + boundary)
    integer(ip) :: n_int_faces   = 0  ! internal faces only
    integer(ip) :: n_bnd_faces   = 0  ! boundary faces
    integer(ip) :: n_patches     = 0
    integer(lp) :: n_cells_global = 0  ! global cell count
    integer(ip) :: n_vertices    = 0  ! optional vertex count (if available)

    ! -- geometry arrays
    type(CellData), allocatable :: cells(:)  ! size n_cells + n_ghost
    type(FaceData), allocatable :: faces(:)  ! size n_faces
    type(PatchData),allocatable :: patches(:)  ! size n_patches
    integer(ip), allocatable     :: periodic_pair(:)  ! boundary-face pairing map (size n_bnd_faces)
    integer(ip), allocatable     :: periodic_neighbor(:)  ! paired owner cell index (local/ghost) per boundary face
    real(rp),    allocatable     :: periodic_pair_dcf(:)  ! paired-face owner->face distance per boundary face

    ! Optional topology arrays (used by geometry-aware output).
    ! cell_types store VTK codes (e.g. 10=tet,12=hex,13=wedge,14=pyramid,42=polyhedron).
    real(rp),    allocatable :: vertices(:,:)  ! (3, n_vertices)
    integer(ip), allocatable :: cell_offsets(:)  ! size n_cells+1, 1-based CSR offsets
    integer(ip), allocatable :: cell_vertices(:)  ! concatenated vertex indices (1-based)
    integer(ip), allocatable :: cell_types(:)  ! size n_cells
    integer(ip), allocatable :: face_offsets(:)  ! size n_faces+1, 1-based CSR offsets
    integer(ip), allocatable :: face_vertices(:)  ! concatenated face vertex indices (1-based)

    ! -- compressed cell->face adjacency
    integer(ip), allocatable :: cell_faces(:)  ! concatenated face lists

    ! -- least-squares gradient stencil weights (dx,dy,dz per neighbour)
    ! stored as: ls_weights(3, i) for each entry in ls_stencil
    integer(ip), allocatable :: ls_stencil_offset(:)  ! size n_cells+1
    integer(ip), allocatable :: ls_stencil(:)  ! neighbour cell indices
    real(rp),    allocatable :: ls_weights(:,:)  ! (3, size(ls_stencil))
    integer(ip), allocatable :: ls_stencil_face(:)  ! face index per stencil entry; -1 for non-boundary

    ! -- global cell IDs (for output / partitioning)
    integer(lp), allocatable :: global_id(:)  ! size n_cells + n_ghost
    integer(lp), allocatable :: compact_id(:)  ! compact global IDs (size n_cells + n_ghost)

    ! -- MPI
    type(MPI_Comm)             :: comm
    integer(ip)                :: rank = 0
    integer(ip)                :: nranks = 1
    type(HaloData), allocatable :: halos(:)  ! one per neighbour rank

    ! -- cached halo communication layout and reusable pack buffers
    logical                    :: halo_layout_ready = .false.
    logical                    :: halo_exchange_active = .false.
    integer(ip)                :: halo_send_cells_total = 0
    integer(ip)                :: halo_recv_cells_total = 0
    integer(ip), allocatable   :: halo_send_off_cells(:)
    integer(ip), allocatable   :: halo_recv_off_cells(:)
    type(MPI_Request), allocatable :: halo_reqs(:)
    real(rp), allocatable      :: halo_sbuf_scalar(:), halo_rbuf_scalar(:)
    real(rp), allocatable      :: halo_sbuf_vec(:), halo_rbuf_vec(:)
    real(rp), allocatable      :: halo_sbuf_scalars(:), halo_rbuf_scalars(:)

    ! -- F9: face coloring for thread-safe parallel face accumulation.
    ! face_color_perm(:) is a permutation of 1..n_faces grouped by color so
    integer(ip)              :: n_colors = 0
    integer(ip), allocatable :: color_offsets(:)  ! size n_colors+1
    integer(ip), allocatable :: face_color_perm(:)  ! size n_faces

    ! Flat SoA mirror of face/cell geometry for !$acc kernels that can't map AoS faces(:)/cells(:); built once, static after setup.
    logical                  :: soa_ready = .false.
    integer(ip), allocatable :: face_owner_soa(:)      ! size n_faces
    integer(ip), allocatable :: face_neighbour_soa(:)  ! size n_faces
    real(rp),    allocatable :: face_sf_soa(:,:)       ! (3, n_faces)
    real(rp),    allocatable :: face_gc_soa(:)         ! size n_faces
    real(rp),    allocatable :: face_eps_soa(:,:)      ! (3, n_faces)
    real(rp),    allocatable :: face_skew_soa(:)        ! size n_faces
    real(rp),    allocatable :: cell_volume_soa(:)     ! size n_cells+n_ghost
    ! Additional geometry needed by viscous/advection flux kernels 
    real(rp),    allocatable :: cell_centroid_soa(:,:)     ! (3, n_cells+n_ghost)
    real(rp),    allocatable :: face_normal_soa(:,:)       ! (3, n_faces) unit outward normal
    real(rp),    allocatable :: face_area_soa(:)           ! size n_faces
    real(rp),    allocatable :: face_centre_soa(:,:)       ! (3, n_faces)
    real(rp),    allocatable :: face_dcf_soa(:)            ! size n_faces
    integer(ip), allocatable :: face_patch_id_soa(:)       ! size n_faces
    integer(ip), allocatable :: patch_bc_id_soa(:)         ! size max(n_patches,1)
    integer(ip), allocatable :: periodic_neighbor_soa(:)   ! size max(n_bnd_faces,1)
    integer(ip), allocatable :: periodic_pair_soa(:)       ! size max(n_bnd_faces,1)
    real(rp),    allocatable :: periodic_pair_dcf_soa(:)   ! size max(n_bnd_faces,1)
    real(rp),    allocatable :: periodic_gc_soa(:)         ! owner weight for periodic-pair interpolation, size max(n_bnd_faces,1)
    ! Flat CSR mirror of cells(:)%face_offset for kernels that walk a cell's face list via cell_faces(:).
    integer(ip), allocatable :: cell_face_offset_soa(:)    ! size n_cells+1

  contains
    procedure :: build_connectivity        => mesh_build_connectivity
    procedure :: compute_geometry          => mesh_compute_geometry
    procedure :: tag_refinement_faces      => mesh_tag_refinement_faces
    procedure :: build_ls_stencil          => mesh_build_ls_stencil
    procedure :: build_face_coloring       => mesh_build_face_coloring
    procedure :: build_soa_mirror          => mesh_build_soa_mirror
    procedure :: exchange_halos_scalar     => mesh_exchange_halos_scalar
    procedure :: begin_halos_scalar        => mesh_begin_halos_scalar
    procedure :: end_halos_scalar          => mesh_end_halos_scalar
    procedure :: exchange_halos_vec        => mesh_exchange_halos_vec
    procedure :: begin_halos_vec           => mesh_begin_halos_vec
    procedure :: end_halos_vec             => mesh_end_halos_vec
    procedure :: exchange_halos_scalars    => mesh_exchange_halos_scalars
    procedure :: begin_halos_scalars       => mesh_begin_halos_scalars
    procedure :: end_halos_scalars         => mesh_end_halos_scalars
    procedure :: exchange_halos_state      => mesh_exchange_halos_state
    procedure :: begin_halos_state         => mesh_begin_halos_state
    procedure :: end_halos_state           => mesh_end_halos_state
    procedure :: destroy                   => mesh_destroy
  end type Mesh

contains

  !> Build cell-to-face and face-to-cell connectivity from raw face owner/
  !> neighbour arrays (already populated by mesh readers).
  subroutine mesh_build_connectivity(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: i, c, f, cnt
    integer(ip), allocatable :: degree(:)

    if (allocated(self%cell_faces)) deallocate(self%cell_faces)  ! safe to call twice (e.g. after a preprocessed-mesh read that already carried it)

    ! Count faces per cell
    allocate(degree(self%n_cells + self%n_ghost), source=0)
    do f = 1, self%n_faces
      c = self%faces(f)%owner
      if (c >= 1) degree(c) = degree(c) + 1
      c = self%faces(f)%neighbour
      if (c >= 1) degree(c) = degree(c) + 1
    end do

    ! Build offsets (CSR style)
    allocate(self%cell_faces(sum(degree)))
    do i = 1, self%n_cells + self%n_ghost
      self%cells(i)%n_faces     = degree(i)
      if (i == 1) then
        self%cells(i)%face_offset = 1
      else
        self%cells(i)%face_offset = self%cells(i-1)%face_offset + degree(i-1)
      end if
    end do

    ! Fill cell_faces
    degree = 0
    do f = 1, self%n_faces
      c = self%faces(f)%owner
      if (c >= 1) then
        cnt = self%cells(c)%face_offset + degree(c)
        self%cell_faces(cnt) = f
        degree(c) = degree(c) + 1
      end if
      c = self%faces(f)%neighbour
      if (c >= 1) then
        cnt = self%cells(c)%face_offset + degree(c)
        self%cell_faces(cnt) = f
        degree(c) = degree(c) + 1
      end if
    end do
    deallocate(degree)
  end subroutine mesh_build_connectivity

  !> F9: build a greedy face coloring for thread-safe accumulation.
  !> Two faces may share a color iff they have no common cell (owner/nb).
  subroutine mesh_build_face_coloring(self)
    class(Mesh), intent(inout) :: self
    integer(ip), allocatable :: face_color(:), counts(:)
    logical, allocatable :: cell_color_used(:,:)
    integer(ip) :: f, o, nb, col, ntot, max_deg, max_colors, c
    logical :: ok

    if (self%n_faces <= 0) return
    if (allocated(self%face_color_perm)) return  ! already built

    ntot = self%n_cells + self%n_ghost
    if (ntot <= 0) return

    ! Upper bound on chromatic number ~= max cell degree + small slack.
    max_deg = 1
    if (self%n_cells > 0 .and. allocated(self%cells)) then
      max_deg = max(1_ip, maxval(self%cells(1:self%n_cells)%n_faces))
    end if
    max_colors = max_deg + 8_ip  ! safety slack for irregular boundary stencils

    allocate(face_color(self%n_faces), source=0_ip)
    allocate(cell_color_used(max_colors, ntot), source=.false.)

    do f = 1, self%n_faces
      o  = self%faces(f)%owner
      nb = self%faces(f)%neighbour
      do col = 1, max_colors
        ok = .true.
        if (o >= 1 .and. o <= ntot) then
          if (cell_color_used(col, o)) ok = .false.
        end if
        if (ok .and. nb >= 1 .and. nb <= ntot) then
          if (cell_color_used(col, nb)) ok = .false.
        end if
        if (ok) exit
      end do
      if (col > max_colors) then
        ! Should not happen with the chosen slack; fall back to a fresh color
        ! by growing the bitmap.
        block
          logical, allocatable :: tmp(:,:)
          integer(ip) :: new_size
          new_size = max_colors * 2_ip
          allocate(tmp(new_size, ntot), source=.false.)
          tmp(1:max_colors, :) = cell_color_used
          call move_alloc(tmp, cell_color_used)
          col = max_colors + 1_ip
          max_colors = new_size
        end block
      end if
      face_color(f) = col
      if (o  >= 1 .and. o  <= ntot) cell_color_used(col, o)  = .true.
      if (nb >= 1 .and. nb <= ntot) cell_color_used(col, nb) = .true.
    end do

    self%n_colors = maxval(face_color)
    if (self%n_colors <= 0) then
      deallocate(face_color, cell_color_used)
      return
    end if

    allocate(counts(self%n_colors), source=0_ip)
    do f = 1, self%n_faces
      counts(face_color(f)) = counts(face_color(f)) + 1
    end do

    allocate(self%color_offsets(self%n_colors + 1))
    self%color_offsets(1) = 1
    do c = 1, self%n_colors
      self%color_offsets(c+1) = self%color_offsets(c) + counts(c)
    end do

    allocate(self%face_color_perm(self%n_faces))
    counts = 0_ip
    do f = 1, self%n_faces
      col = face_color(f)
      self%face_color_perm(self%color_offsets(col) + counts(col)) = f
      counts(col) = counts(col) + 1_ip
    end do

    deallocate(face_color, cell_color_used, counts)
  end subroutine mesh_build_face_coloring

  !> Flatten faces(:)/cells(:) into plain arrays so !$acc kernels can map them without deep-copying the AoS derived types.
  subroutine mesh_build_soa_mirror(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: f, c, ntot

    if (self%soa_ready) return
    if (self%n_faces <= 0) return

    allocate(self%face_owner_soa(self%n_faces))
    allocate(self%face_neighbour_soa(self%n_faces))
    allocate(self%face_sf_soa(3, self%n_faces))
    allocate(self%face_gc_soa(self%n_faces))
    allocate(self%face_eps_soa(3, self%n_faces))
    allocate(self%face_skew_soa(self%n_faces))

    do f = 1, self%n_faces
      self%face_owner_soa(f)     = self%faces(f)%owner
      self%face_neighbour_soa(f) = self%faces(f)%neighbour
      self%face_sf_soa(:, f)     = self%faces(f)%sf
      self%face_gc_soa(f)        = self%faces(f)%gc
      self%face_eps_soa(:, f)    = self%faces(f)%eps
      self%face_skew_soa(f)      = self%faces(f)%skew
    end do

    ntot = self%n_cells + self%n_ghost
    allocate(self%cell_volume_soa(ntot))
    allocate(self%cell_centroid_soa(3, ntot))
    do c = 1, ntot
      self%cell_volume_soa(c) = self%cells(c)%volume
      self%cell_centroid_soa(:, c) = self%cells(c)%centroid
    end do

    allocate(self%face_normal_soa(3, self%n_faces))
    allocate(self%face_area_soa(self%n_faces))
    allocate(self%face_centre_soa(3, self%n_faces))
    allocate(self%face_dcf_soa(self%n_faces))
    allocate(self%face_patch_id_soa(self%n_faces))
    do f = 1, self%n_faces
      self%face_normal_soa(:, f) = self%faces(f)%normal
      self%face_area_soa(f)      = self%faces(f)%area
      self%face_centre_soa(:, f) = self%faces(f)%centre
      self%face_dcf_soa(f)       = self%faces(f)%d_cf
      self%face_patch_id_soa(f)  = self%faces(f)%patch_id
    end do

    allocate(self%patch_bc_id_soa(max(self%n_patches, 1_ip)))
    self%patch_bc_id_soa = 0_ip
    if (self%n_patches >= 1 .and. allocated(self%patches)) then
      do c = 1, self%n_patches
        self%patch_bc_id_soa(c) = self%patches(c)%bc_id
      end do
    end if

    allocate(self%periodic_neighbor_soa(max(self%n_bnd_faces, 1_ip)))
    allocate(self%periodic_pair_soa(max(self%n_bnd_faces, 1_ip)))
    allocate(self%periodic_pair_dcf_soa(max(self%n_bnd_faces, 1_ip)))
    self%periodic_neighbor_soa = -1_ip
    self%periodic_pair_soa     = -1_ip
    self%periodic_pair_dcf_soa = 0.0_rp
    if (allocated(self%periodic_neighbor)) &
      self%periodic_neighbor_soa(1:size(self%periodic_neighbor)) = self%periodic_neighbor
    if (allocated(self%periodic_pair)) &
      self%periodic_pair_soa(1:size(self%periodic_pair)) = self%periodic_pair
    if (allocated(self%periodic_pair_dcf)) &
      self%periodic_pair_dcf_soa(1:size(self%periodic_pair_dcf)) = self%periodic_pair_dcf

    ! Owner weight for periodic-pair interpolation: gc_per = d_cf_nb/(d_cf_o+d_cf_nb), matching the
    ! interior-face gc convention. faces(:)%gc is hard-coded to 1 on boundary/periodic faces (see
    ! mesh_compute_geometry) so consumers must read this array instead for BC_PERIODIC faces.
    allocate(self%periodic_gc_soa(max(self%n_bnd_faces, 1_ip)))
    self%periodic_gc_soa = 0.5_rp
    do f = 1, self%n_bnd_faces
      if (self%periodic_pair_dcf_soa(f) > 0.0_rp) then
        self%periodic_gc_soa(f) = self%periodic_pair_dcf_soa(f) / &
          (self%face_dcf_soa(self%n_int_faces + f) + self%periodic_pair_dcf_soa(f) + SMALL)
      end if
    end do

    allocate(self%cell_face_offset_soa(self%n_cells + 1))
    if (self%n_cells >= 1) then
      do c = 1, self%n_cells
        self%cell_face_offset_soa(c) = self%cells(c)%face_offset
      end do
      self%cell_face_offset_soa(self%n_cells + 1) = self%cells(self%n_cells)%face_offset + self%cells(self%n_cells)%n_faces
    else
      self%cell_face_offset_soa(1) = 1
    end if

    self%soa_ready = .true.
  end subroutine mesh_build_soa_mirror

  !> Compute face area vectors and geometric interpolation weights.
  !> Cell centroids and volumes must have been set by the reader.
  subroutine mesh_compute_geometry(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: f
    real(rp)    :: d_oc(3), d_cn(3), d_on(3), d_on2, d_len, w

    do f = 1, self%n_faces
      associate(face => self%faces(f))
        ! Area-weighted face normal
        face%sf = face%area * face%normal

        ! Distance owner -> neighbour (or owner -> face centre for boundary)
        if (face%neighbour >= 1) then
          d_on   = self%cells(face%neighbour)%centroid - &
                   self%cells(face%owner)%centroid
          d_on2  = sum(d_on**2)
          d_len  = sqrt(d_on2)
          face%d_cf = max(d_len, 1.0e-30_rp)  ! true geometric distance, zero-protected

          ! gc = 1 - (d_oc.d_on)/|d_on|^2; SMALL in denominator only as zero guard.
          d_oc = face%centre - self%cells(face%owner)%centroid
          w    = dot_product(d_oc, d_on) / (d_on2 + SMALL)
          face%gc = max(0.0_rp, min(1.0_rp, 1.0_rp - w))  ! owner weight

          ! Skewness offset for the O(h^2) gc-interpolation correction (Kim & Choi Eq. 26);
          ! gc is the owner weight, so the gc-interpolation point is x_o + (1-gc)*d_on.
          face%eps = face%centre - (self%cells(face%owner)%centroid + (1.0_rp - face%gc) * d_on)
        else
          ! Boundary: gc = 1 (use owner value directly)
          d_cn      = face%centre - self%cells(face%owner)%centroid
          face%d_cf = max(sqrt(sum(d_cn**2)), 1.0e-30_rp)
          face%gc   = 1.0_rp
          face%eps  = 0.0_rp
        end if
      end associate
    end do

    call mesh_tag_refinement_faces(self)
  end subroutine mesh_compute_geometry

  !> Tag refinement-interface faces from serialized geometry; the preprocessed path skips compute_geometry.
  subroutine mesh_tag_refinement_faces(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: f, o, nb
    real(rp)    :: eps_mid(3)

    do f = 1, self%n_faces
      o  = self%faces(f)%owner
      nb = self%faces(f)%neighbour
      if (nb < 1 .or. o < 1 .or. o > size(self%cells) .or. nb > size(self%cells)) then
        self%faces(f)%skew = 0.0_rp
        self%faces(f)%is_refine = .false.
        cycle
      end if
      eps_mid = self%faces(f)%centre - 0.5_rp * (self%cells(o)%centroid + self%cells(nb)%centroid)
      self%faces(f)%skew = sqrt(sum(eps_mid**2)) / max(self%faces(f)%d_cf, 1.0e-30_rp)
      self%faces(f)%is_refine = (self%faces(f)%skew > REFINE_SKEW_TOL)
    end do
  end subroutine mesh_tag_refinement_faces

  !> Smoothstep blend weight for the eps.grad skewness correction: 1 below SKEW_BLEND_LO, 0 above SKEW_BLEND_HI.
  pure function refine_blend_weight(skew) result(w)
    !$acc routine seq
    real(rp), intent(in) :: skew
    real(rp) :: t, w

    t = (skew - SKEW_BLEND_LO) / (SKEW_BLEND_HI - SKEW_BLEND_LO)
    t = max(0.0_rp, min(1.0_rp, t))
    w = 1.0_rp - t * t * (3.0_rp - 2.0_rp * t)
  end function refine_blend_weight

  !> Build least-squares gradient stencil weights for each interior cell.
  !> Uses Green-Gauss neighbour stencil with inverse-distance weighting.
  subroutine mesh_build_ls_stencil(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: i, f, nb, idx, offset, nstenc, f_bnd, pair_bnd, pair_face, patch_id
    integer(ip) :: max_stencil, jj, c, k, kk, v, cnb
    integer(ip), allocatable :: nbcount(:), tmp_nb(:,:)
    logical,     allocatable :: tmp_is_periodic(:,:)
    integer(ip), allocatable :: tmp_face(:,:)
    logical                  :: is_per, found, use_vertex_ext
    real(rp)    :: dx, dy, dz, w, ww, d_per_o, d_per_nb
    real(rp)    :: AtA(3,3), AtAinv(3,3), rhs(3)
    real(rp)    :: det, ata_scale
    ! Vertex-to-cell CSR (built if cell_offsets/cell_vertices are available)
    integer(ip), allocatable :: voff(:), vcell(:), vptr(:)

    ! Allocate tmp stencil arrays: 4*max_faces covers face + edge neighbours (hex: 6+12=18).
    allocate(nbcount(self%n_cells), source=0)
    max_stencil = max(32_ip, 4_ip * max(1_ip, maxval(self%cells(1:self%n_cells)%n_faces)))
    allocate(tmp_nb(max_stencil, self%n_cells))
    allocate(tmp_is_periodic(max_stencil, self%n_cells), source=.false.)
    allocate(tmp_face(max_stencil, self%n_cells), source=-1_ip)

    ! First pass: collect face-sharing neighbours (as before)
    do i = 1, self%n_cells
      offset = self%cells(i)%face_offset
      nstenc = self%cells(i)%n_faces
      do idx = offset, offset + nstenc - 1
        f  = self%cell_faces(idx)
        nb = self%faces(f)%owner
        if (nb == i) nb = self%faces(f)%neighbour
        is_per = .false.
        if (nb < 1) then
          ! Periodic boundary: use paired owner cell as stencil neighbour.
          f_bnd = f - self%n_int_faces
          pair_bnd = 0
          patch_id = self%faces(f)%patch_id
          if (patch_id >= 1 .and. patch_id <= self%n_patches) then
            if (self%patches(patch_id)%bc_id == BC_PERIODIC) then
              if (allocated(self%periodic_neighbor)) then
                if (f_bnd >= 1 .and. f_bnd <= size(self%periodic_neighbor)) then
                  nb = self%periodic_neighbor(f_bnd)
                end if
              end if
              if (nb < 1 .and. allocated(self%periodic_pair)) then
                if (f_bnd >= 1 .and. f_bnd <= size(self%periodic_pair)) then
                  pair_bnd = self%periodic_pair(f_bnd)
                  if (pair_bnd > 0) then
                    pair_face = self%n_int_faces + pair_bnd
                    if (pair_face >= 1 .and. pair_face <= size(self%faces)) nb = self%faces(pair_face)%owner
                  end if
                end if
              end if
              if (nb >= 1) is_per = .true.
            end if
          end if
          if (nb < 1) nb = -1  ! physical boundary: face-centre pseudo-neighbour
        end if
        nbcount(i) = nbcount(i) + 1
        tmp_nb(nbcount(i), i) = nb
        tmp_is_periodic(nbcount(i), i) = is_per
        if (nb == -1) tmp_face(nbcount(i), i) = f
      end do
    end do

    ! Optional second pass: add vertex-adjacent cells to improve LS conditioning near walls.
    ! Split into nested ifs: .and. does not short-circuit in Fortran, so size(self%cell_offsets)
    ! must never be evaluated in the same expression as the allocated() check.
    use_vertex_ext = .false.
    if (allocated(self%cell_offsets) .and. allocated(self%cell_vertices) &
        .and. self%n_vertices > 0) then
      if (size(self%cell_offsets) >= self%n_cells + 1) use_vertex_ext = .true.
    end if
    if (use_vertex_ext) then
      ! Build vertex-to-cell CSR
      allocate(voff(self%n_vertices + 1), source=0)
      do c = 1, self%n_cells
        do k = self%cell_offsets(c), self%cell_offsets(c+1) - 1
          v = self%cell_vertices(k)
          if (v >= 1 .and. v <= self%n_vertices) voff(v+1) = voff(v+1) + 1
        end do
      end do
      voff(1) = 1
      do v = 1, self%n_vertices
        voff(v+1) = voff(v+1) + voff(v)
      end do
      allocate(vcell(voff(self%n_vertices+1) - 1), source=0_ip)
      allocate(vptr(self%n_vertices+1))
      vptr = voff
      do c = 1, self%n_cells
        do k = self%cell_offsets(c), self%cell_offsets(c+1) - 1
          v = self%cell_vertices(k)
          if (v >= 1 .and. v <= self%n_vertices) then
            vcell(vptr(v)) = c
            vptr(v) = vptr(v) + 1
          end if
        end do
      end do
      deallocate(vptr)
      ! For each cell, add vertex-sharing cells not already in the face stencil
      do i = 1, self%n_cells
        if (i+1 > size(self%cell_offsets)) cycle
        do k = self%cell_offsets(i), self%cell_offsets(i+1) - 1
          v = self%cell_vertices(k)
          if (v < 1 .or. v > self%n_vertices) cycle
          do kk = voff(v), voff(v+1) - 1
            cnb = vcell(kk)
            if (cnb < 1 .or. cnb == i .or. cnb > self%n_cells) cycle
            ! Check not already in current stencil
            found = .false.
            do jj = 1, nbcount(i)
              if (tmp_nb(jj, i) == cnb) then; found = .true.; exit; end if
            end do
            if (found) cycle
            if (nbcount(i) < max_stencil) then
              nbcount(i) = nbcount(i) + 1
              tmp_nb(nbcount(i), i)         = cnb
              tmp_is_periodic(nbcount(i), i) = .false.
              tmp_face(nbcount(i), i)        = -1_ip
            end if
          end do
        end do
      end do
      deallocate(voff, vcell)
    end if

    ! Build CSR structure
    allocate(self%ls_stencil_offset(self%n_cells + 1))
    self%ls_stencil_offset(1) = 1
    do i = 1, self%n_cells
      self%ls_stencil_offset(i+1) = self%ls_stencil_offset(i) + nbcount(i)
    end do
    nstenc = self%ls_stencil_offset(self%n_cells + 1) - 1
    allocate(self%ls_stencil(nstenc))
    allocate(self%ls_weights(3, nstenc))
    allocate(self%ls_stencil_face(nstenc), source=-1_ip)

    ! Second pass: compute inverse-distance LS weights via 3x3 normal equations
    do i = 1, self%n_cells
      offset = self%ls_stencil_offset(i)
      AtA    = 0.0_rp

      do idx = 1, nbcount(i)
        nb = tmp_nb(idx, i)
        self%ls_stencil(offset + idx - 1) = nb
        self%ls_stencil_face(offset + idx - 1) = tmp_face(idx, i)

        if (nb >= 1 .and. tmp_is_periodic(idx, i)) then
          ! Periodic neighbour: use (d_cf_o + d_cf_nb)*n_hat (exact for any spacing; 2*d_cf approximation is uniform-grid only)
          f       = self%cell_faces(self%cells(i)%face_offset + idx - 1)
          d_per_o = self%faces(f)%d_cf
          d_per_nb = d_per_o  ! default: assume symmetric
          f_bnd = f - self%n_int_faces
          ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()/
          ! the array element must never be evaluated in the same expression as allocated().
          if (allocated(self%periodic_pair_dcf)) then
            if (f_bnd >= 1 .and. f_bnd <= size(self%periodic_pair_dcf)) then
              if (self%periodic_pair_dcf(f_bnd) > 0.0_rp) d_per_nb = self%periodic_pair_dcf(f_bnd)
            end if
          end if
          dx = (d_per_o + d_per_nb) * self%faces(f)%normal(1)
          dy = (d_per_o + d_per_nb) * self%faces(f)%normal(2)
          dz = (d_per_o + d_per_nb) * self%faces(f)%normal(3)
        else if (nb >= 1) then
          dx = self%cells(nb)%centroid(1) - self%cells(i)%centroid(1)
          dy = self%cells(nb)%centroid(2) - self%cells(i)%centroid(2)
          dz = self%cells(nb)%centroid(3) - self%cells(i)%centroid(3)
        else
          ! Boundary: use face centre
          f  = self%cell_faces(self%cells(i)%face_offset + idx - 1)
          dx = self%faces(f)%centre(1) - self%cells(i)%centroid(1)
          dy = self%faces(f)%centre(2) - self%cells(i)%centroid(2)
          dz = self%faces(f)%centre(3) - self%cells(i)%centroid(3)
        end if

        w  = 1.0_rp / (sqrt(dx**2 + dy**2 + dz**2) + SMALL)
        ww = w * w
        ! Accumulate AtA = sum_j w^2 * [dx,dy,dz]^T [dx,dy,dz]
        AtA(1,1) = AtA(1,1) + ww*dx*dx; AtA(1,2) = AtA(1,2) + ww*dx*dy
        AtA(1,3) = AtA(1,3) + ww*dx*dz; AtA(2,2) = AtA(2,2) + ww*dy*dy
        AtA(2,3) = AtA(2,3) + ww*dy*dz; AtA(3,3) = AtA(3,3) + ww*dz*dz
      end do
      AtA(2,1) = AtA(1,2); AtA(3,1) = AtA(1,3); AtA(3,2) = AtA(2,3)

      ! Tikhonov-regularise ill-conditioned stencils (sliver or near-coplanar cells) before inverting.
      ata_scale = (AtA(1,1) + AtA(2,2) + AtA(3,3)) / 3.0_rp
      call invert3x3(AtA, AtAinv, det)
      if (ata_scale <= 0.0_rp .or. abs(det) < LS_DET_TOL * ata_scale**3) then
        do k = 1, 3
          AtA(k,k) = AtA(k,k) + LS_REG_FRAC * max(ata_scale, SMALL)
        end do
        call invert3x3(AtA, AtAinv, det)
        ls_degenerate_count = ls_degenerate_count + 1_ip
      end if

      ! Compute weight row for each neighbour: w_j = AtAinv * w^2 * [dx,dy,dz]
      do idx = 1, nbcount(i)
        nb = tmp_nb(idx, i)
        if (nb >= 1 .and. tmp_is_periodic(idx, i)) then
          f       = self%cell_faces(self%cells(i)%face_offset + idx - 1)
          d_per_o = self%faces(f)%d_cf
          d_per_nb = d_per_o
          f_bnd = f - self%n_int_faces
          ! Split into nested ifs: .and. does not short-circuit in Fortran, so size()/
          ! the array element must never be evaluated in the same expression as allocated().
          if (allocated(self%periodic_pair_dcf)) then
            if (f_bnd >= 1 .and. f_bnd <= size(self%periodic_pair_dcf)) then
              if (self%periodic_pair_dcf(f_bnd) > 0.0_rp) d_per_nb = self%periodic_pair_dcf(f_bnd)
            end if
          end if
          dx = (d_per_o + d_per_nb) * self%faces(f)%normal(1)
          dy = (d_per_o + d_per_nb) * self%faces(f)%normal(2)
          dz = (d_per_o + d_per_nb) * self%faces(f)%normal(3)
        else if (nb >= 1) then
          dx = self%cells(nb)%centroid(1) - self%cells(i)%centroid(1)
          dy = self%cells(nb)%centroid(2) - self%cells(i)%centroid(2)
          dz = self%cells(nb)%centroid(3) - self%cells(i)%centroid(3)
        else
          f  = self%cell_faces(self%cells(i)%face_offset + idx - 1)
          dx = self%faces(f)%centre(1) - self%cells(i)%centroid(1)
          dy = self%faces(f)%centre(2) - self%cells(i)%centroid(2)
          dz = self%faces(f)%centre(3) - self%cells(i)%centroid(3)
        end if
        w = 1.0_rp / (sqrt(dx**2 + dy**2 + dz**2) + SMALL)
        rhs = w*w * [dx, dy, dz]
        self%ls_weights(:, offset + idx - 1) = matmul(AtAinv, rhs)
      end do
    end do

    deallocate(nbcount, tmp_nb, tmp_is_periodic, tmp_face)
  end subroutine mesh_build_ls_stencil

  integer(ip) function mesh_ls_degenerate_count() result(n)
    n = ls_degenerate_count
  end function mesh_ls_degenerate_count

  !> MPI halo exchange for a scalar field (n_cells + n_ghost).
  subroutine mesh_exchange_halos_scalar(self, phi)
    class(Mesh),  intent(inout) :: self
    real(rp),     intent(inout) :: phi(:)
    call mesh_begin_halos_scalar(self, phi)
    call mesh_end_halos_scalar(self, phi)
  end subroutine mesh_exchange_halos_scalar

  subroutine mesh_begin_halos_scalar(self, phi)
    class(Mesh),  intent(inout) :: self
    real(rp),     intent(inout) :: phi(:)

    integer(ip) :: h, i, nh
    integer(ip) :: nreq, tag, send_total, recv_total

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (self%halo_exchange_active) error stop 'mesh_begin_halos_scalar: halo exchange already active'

    call ensure_halo_layout(self)

    nh = size(self%halos)
    if (nh <= 0) return

    send_total = self%halo_send_cells_total
    recv_total = self%halo_recv_cells_total
    call ensure_real_buffer_size(self%halo_sbuf_scalar, send_total)
    call ensure_real_buffer_size(self%halo_rbuf_scalar, recv_total)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_send
          self%halo_sbuf_scalar(self%halo_send_off_cells(h) + i) = phi(halo%send_cells(i))
        end do
      end associate
    end do

    if (2*nh > size(self%halo_reqs)) error stop 'halo_reqs array too small - internal error'
    nreq = 0
    tag = HALO_TAG_SCALAR
    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Irecv(self%halo_rbuf_scalar(self%halo_recv_off_cells(h) + 1), halo%n_recv, MPI_DOUBLE_PRECISION, &
                       halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Isend(self%halo_sbuf_scalar(self%halo_send_off_cells(h) + 1), halo%n_send, MPI_DOUBLE_PRECISION, &
                       halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    self%halo_exchange_active = .true.
  end subroutine mesh_begin_halos_scalar

  subroutine mesh_end_halos_scalar(self, phi)
    class(Mesh),  intent(inout) :: self
    real(rp),     intent(inout) :: phi(:)

    integer(ip) :: h, i, nh, nreq

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (.not. self%halo_exchange_active) return

    nh = size(self%halos)
    if (nh <= 0) then
      self%halo_exchange_active = .false.
      return
    end if

    nreq = 2 * nh
    if (nreq > 0) call MPI_Waitall(nreq, self%halo_reqs(1:nreq), MPI_STATUSES_IGNORE)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_recv
          phi(halo%recv_cells(i)) = self%halo_rbuf_scalar(self%halo_recv_off_cells(h) + i)
        end do
      end associate
    end do

    ! See mesh_end_halos_vec: keeps a device-resident field's ghost tail in sync with the
    ! host scatter above; if_present makes this a no-op when phi isn't on device (required,
    ! not automatic -- update on absent data errors without it).
    !$acc update device(phi(self%n_cells+1:self%n_cells+self%n_ghost)) if(self%n_ghost > 0) if_present

    self%halo_exchange_active = .false.
  end subroutine mesh_end_halos_scalar

  !> MPI halo exchange for a 3-vector field.
  subroutine mesh_exchange_halos_vec(self, U)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! U(3, n_cells+n_ghost)
    call mesh_begin_halos_vec(self, U)
    call mesh_end_halos_vec(self, U)
  end subroutine mesh_exchange_halos_vec

  subroutine mesh_begin_halos_vec(self, U)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! U(3, n_cells+n_ghost)

    integer(ip) :: h, i, nh
    integer(ip) :: nreq, tag, send_total, recv_total

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (self%halo_exchange_active) error stop 'mesh_begin_halos_vec: halo exchange already active'

    call ensure_halo_layout(self)

    nh = size(self%halos)
    if (nh <= 0) return

    send_total = 3 * self%halo_send_cells_total
    recv_total = 3 * self%halo_recv_cells_total
    call ensure_real_buffer_size(self%halo_sbuf_vec, send_total)
    call ensure_real_buffer_size(self%halo_rbuf_vec, recv_total)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_send
          self%halo_sbuf_vec(3*self%halo_send_off_cells(h) + 3*(i-1) + 1:3*self%halo_send_off_cells(h) + 3*i) = &
            U(:, halo%send_cells(i))
        end do
      end associate
    end do

    if (2*nh > size(self%halo_reqs)) error stop 'halo_reqs array too small - internal error'
    nreq = 0
    tag = HALO_TAG_VECTOR
    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Irecv(self%halo_rbuf_vec(3*self%halo_recv_off_cells(h) + 1), 3*halo%n_recv, MPI_DOUBLE_PRECISION, &
                       halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Isend(self%halo_sbuf_vec(3*self%halo_send_off_cells(h) + 1), 3*halo%n_send, MPI_DOUBLE_PRECISION, &
                       halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    self%halo_exchange_active = .true.
  end subroutine mesh_begin_halos_vec

  subroutine mesh_end_halos_vec(self, U)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! U(3, n_cells+n_ghost)

    integer(ip) :: h, i, nh, nreq

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (.not. self%halo_exchange_active) return

    nh = size(self%halos)
    if (nh <= 0) then
      self%halo_exchange_active = .false.
      return
    end if

    nreq = 2 * nh
    if (nreq > 0) call MPI_Waitall(nreq, self%halo_reqs(1:nreq), MPI_STATUSES_IGNORE)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_recv
          U(:, halo%recv_cells(i)) = self%halo_rbuf_vec(3*self%halo_recv_off_cells(h) + 3*(i-1) + 1: &
                                 3*self%halo_recv_off_cells(h) + 3*i)
        end do
      end associate
    end do

    ! Ghost cells occupy the contiguous tail [n_cells+1, n_cells+n_ghost] of every field shaped
    ! like U, so one range update covers whatever this call's halos touched. Host U above is
    ! always authoritative (BC application and the predictor correction are host-only), so this
    ! is purely pushing that forward to keep a device-resident U in sync for the numerics kernels
    ! that read ghost cells next (gradient/advection/viscous). if_present is required for safety
    ! when U isn't on device yet (e.g. before advance_rk3w's data region opens) -- without it,
    ! update device on absent data is a hard runtime error, not a no-op.
    !$acc update device(U(:, self%n_cells+1:self%n_cells+self%n_ghost)) if(self%n_ghost > 0) if_present

    self%halo_exchange_active = .false.
  end subroutine mesh_end_halos_vec

  !> MPI halo exchange for multiple scalar fields packed in one communication.
  subroutine mesh_exchange_halos_scalars(self, phi)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: phi(:,:)  ! (n_fields, n_cells+n_ghost)

    integer(ip) :: h, i, n_fields, nh
    integer(ip) :: nreq, tag, send_total, recv_total

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (self%halo_exchange_active) error stop 'mesh_exchange_halos_scalars: halo exchange already active'
    call ensure_halo_layout(self)

    n_fields = size(phi, 1)
    if (n_fields <= 0) return

    nh = size(self%halos)
    if (nh <= 0) return

    send_total = n_fields * self%halo_send_cells_total
    recv_total = n_fields * self%halo_recv_cells_total
    call ensure_real_buffer_size(self%halo_sbuf_scalars, send_total)
    call ensure_real_buffer_size(self%halo_rbuf_scalars, recv_total)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_send
          self%halo_sbuf_scalars(n_fields*self%halo_send_off_cells(h) + n_fields*(i-1) + 1: &
                                 n_fields*self%halo_send_off_cells(h) + n_fields*i) = phi(:, halo%send_cells(i))
        end do
      end associate
    end do

    if (2*nh > size(self%halo_reqs)) error stop 'halo_reqs array too small - internal error'
    nreq = 0
    tag = HALO_TAG_MULTI

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Irecv(self%halo_rbuf_scalars(n_fields*self%halo_recv_off_cells(h) + 1), &
                       n_fields*halo%n_recv, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Isend(self%halo_sbuf_scalars(n_fields*self%halo_send_off_cells(h) + 1), &
                       n_fields*halo%n_send, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    if (nreq > 0) call MPI_Waitall(nreq, self%halo_reqs(1:nreq), MPI_STATUSES_IGNORE)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_recv
          phi(:, halo%recv_cells(i)) = self%halo_rbuf_scalars(n_fields*self%halo_recv_off_cells(h) + n_fields*(i-1) + 1: &
                                  n_fields*self%halo_recv_off_cells(h) + n_fields*i)
        end do
      end associate
    end do
  end subroutine mesh_exchange_halos_scalars

  subroutine mesh_begin_halos_scalars(self, phi)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: phi(:,:)  ! (n_fields, n_cells+n_ghost)

    integer(ip) :: h, i, n_fields, nh
    integer(ip) :: nreq, tag, send_total, recv_total

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (self%halo_exchange_active) error stop 'mesh_begin_halos_scalars: halo exchange already active'

    call ensure_halo_layout(self)

    n_fields = size(phi, 1)
    if (n_fields <= 0) return

    nh = size(self%halos)
    if (nh <= 0) return

    send_total = n_fields * self%halo_send_cells_total
    recv_total = n_fields * self%halo_recv_cells_total
    call ensure_real_buffer_size(self%halo_sbuf_scalars, send_total)
    call ensure_real_buffer_size(self%halo_rbuf_scalars, recv_total)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_send
          self%halo_sbuf_scalars(n_fields*self%halo_send_off_cells(h) + n_fields*(i-1) + 1: &
                                 n_fields*self%halo_send_off_cells(h) + n_fields*i) = phi(:, halo%send_cells(i))
        end do
      end associate
    end do

    if (2*nh > size(self%halo_reqs)) error stop 'halo_reqs array too small - internal error'
    nreq = 0
    tag = HALO_TAG_MULTI

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Irecv(self%halo_rbuf_scalars(n_fields*self%halo_recv_off_cells(h) + 1), &
                       n_fields*halo%n_recv, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Isend(self%halo_sbuf_scalars(n_fields*self%halo_send_off_cells(h) + 1), &
                       n_fields*halo%n_send, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    self%halo_exchange_active = .true.
  end subroutine mesh_begin_halos_scalars

  subroutine mesh_end_halos_scalars(self, phi)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: phi(:,:)  ! (n_fields, n_cells+n_ghost)

    integer(ip) :: h, i, n_fields, nh, nreq

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (.not. self%halo_exchange_active) return

    n_fields = size(phi, 1)
    nh = size(self%halos)
    if (n_fields <= 0 .or. nh <= 0) then
      self%halo_exchange_active = .false.
      return
    end if

    nreq = 2 * nh
    if (nreq > 0) call MPI_Waitall(nreq, self%halo_reqs(1:nreq), MPI_STATUSES_IGNORE)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_recv
          phi(:, halo%recv_cells(i)) = self%halo_rbuf_scalars(n_fields*self%halo_recv_off_cells(h) + n_fields*(i-1) + 1: &
                                  n_fields*self%halo_recv_off_cells(h) + n_fields*i)
        end do
      end associate
    end do

    ! See mesh_end_halos_vec: keeps a device-resident field's ghost tail in sync with the
    ! host scatter above; if_present is required so this is a no-op when phi isn't on device.
    !$acc update device(phi(:, self%n_cells+1:self%n_cells+self%n_ghost)) if(self%n_ghost > 0) if_present

    self%halo_exchange_active = .false.
  end subroutine mesh_end_halos_scalars

  !> Batched halo exchange for solver state: velocity, pressure, temperature,
  !> and optional passive scalars.
  subroutine mesh_exchange_halos_state(self, U, p, T, C)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! (3, n_cells+n_ghost)
    real(rp),    intent(inout) :: p(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout) :: T(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout), optional :: C(:,:)  ! (n_scalars, n_cells+n_ghost)
    call mesh_begin_halos_state(self, U, p, T, C)
    call mesh_end_halos_state(self, U, p, T, C)
  end subroutine mesh_exchange_halos_state

  subroutine mesh_begin_halos_state(self, U, p, T, C)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! (3, n_cells+n_ghost)
    real(rp),    intent(inout) :: p(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout) :: T(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout), optional :: C(:,:)  ! (n_scalars, n_cells+n_ghost)

    integer(ip) :: h, i, nh, n_scalars, n_fields, cidx, base
    integer(ip) :: nreq, tag, send_total, recv_total

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (self%halo_exchange_active) error stop 'mesh_begin_halos_state: halo exchange already active'
    call ensure_halo_layout(self)

    nh = size(self%halos)
    if (nh <= 0) return

    n_scalars = 0
    if (present(C)) n_scalars = size(C, 1)
    n_fields = 5 + n_scalars

    send_total = n_fields * self%halo_send_cells_total
    recv_total = n_fields * self%halo_recv_cells_total
    call ensure_real_buffer_size(self%halo_sbuf_scalars, send_total)
    call ensure_real_buffer_size(self%halo_rbuf_scalars, recv_total)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_send
          cidx = halo%send_cells(i)
          base = n_fields*self%halo_send_off_cells(h) + n_fields*(i-1)
          self%halo_sbuf_scalars(base + 1:base + 3) = U(:, cidx)
          self%halo_sbuf_scalars(base + 4) = p(cidx)
          self%halo_sbuf_scalars(base + 5) = T(cidx)
          if (n_scalars > 0) self%halo_sbuf_scalars(base + 6:base + 5 + n_scalars) = C(:, cidx)
        end do
      end associate
    end do

    if (2*nh > size(self%halo_reqs)) error stop 'halo_reqs array too small - internal error'
    nreq = 0
    tag = HALO_TAG_STATE
    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Irecv(self%halo_rbuf_scalars(n_fields*self%halo_recv_off_cells(h) + 1), &
                       n_fields*halo%n_recv, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    do h = 1, nh
      associate(halo => self%halos(h))
        nreq = nreq + 1
        call MPI_Isend(self%halo_sbuf_scalars(n_fields*self%halo_send_off_cells(h) + 1), &
                       n_fields*halo%n_send, MPI_DOUBLE_PRECISION, halo%rank, tag, self%comm, self%halo_reqs(nreq))
      end associate
    end do

    self%halo_exchange_active = .true.
  end subroutine mesh_begin_halos_state

  subroutine mesh_end_halos_state(self, U, p, T, C)
    class(Mesh), intent(inout) :: self
    real(rp),    intent(inout) :: U(:,:)  ! (3, n_cells+n_ghost)
    real(rp),    intent(inout) :: p(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout) :: T(:)  ! (n_cells+n_ghost)
    real(rp),    intent(inout), optional :: C(:,:)  ! (n_scalars, n_cells+n_ghost)

    integer(ip) :: h, i, nh, n_scalars, n_fields, cidx, base, nreq

    if (self%nranks == 1 .or. .not. allocated(self%halos)) return
    if (.not. self%halo_exchange_active) return

    nh = size(self%halos)
    if (nh <= 0) then
      self%halo_exchange_active = .false.
      return
    end if

    n_scalars = 0
    if (present(C)) n_scalars = size(C, 1)
    n_fields = 5 + n_scalars
    nreq = 2 * nh
    if (nreq > 0) call MPI_Waitall(nreq, self%halo_reqs(1:nreq), MPI_STATUSES_IGNORE)

    do h = 1, nh
      associate(halo => self%halos(h))
        do i = 1, halo%n_recv
          cidx = halo%recv_cells(i)
          base = n_fields*self%halo_recv_off_cells(h) + n_fields*(i-1)
          U(:, cidx) = self%halo_rbuf_scalars(base + 1:base + 3)
          p(cidx) = self%halo_rbuf_scalars(base + 4)
          T(cidx) = self%halo_rbuf_scalars(base + 5)
          if (n_scalars > 0) C(:, cidx) = self%halo_rbuf_scalars(base + 6:base + 5 + n_scalars)
        end do
      end associate
    end do

    ! See mesh_end_halos_vec: keeps device-resident fields' ghost tails in sync with the host
    ! scatter above; if_present is required so this is a no-op for arrays not on device.
    !$acc update device(U(:, self%n_cells+1:self%n_cells+self%n_ghost), &
    !$acc      p(self%n_cells+1:self%n_cells+self%n_ghost), &
    !$acc      T(self%n_cells+1:self%n_cells+self%n_ghost)) if(self%n_ghost > 0) if_present
    ! C is optional: referencing it in an !$acc clause when absent is invalid regardless of the
    ! clause's own runtime if(), since data clauses are resolved statically -- guard with a plain
    ! Fortran present() check instead (see mod_gradient.f90 for the same constraint on U_bnd).
    if (present(C)) then
      !$acc update device(C(:, self%n_cells+1:self%n_cells+self%n_ghost)) if(self%n_ghost > 0 .and. n_scalars > 0) if_present
    end if

    self%halo_exchange_active = .false.
  end subroutine mesh_end_halos_state

  subroutine mesh_destroy(self)
    class(Mesh), intent(inout) :: self
    if (allocated(self%cells))              deallocate(self%cells)
    if (allocated(self%faces))              deallocate(self%faces)
    if (allocated(self%patches))            deallocate(self%patches)
    if (allocated(self%periodic_pair))      deallocate(self%periodic_pair)
    if (allocated(self%periodic_neighbor))  deallocate(self%periodic_neighbor)
    if (allocated(self%periodic_pair_dcf))  deallocate(self%periodic_pair_dcf)
    if (allocated(self%cell_faces))         deallocate(self%cell_faces)
    if (allocated(self%vertices))           deallocate(self%vertices)
    if (allocated(self%cell_offsets))       deallocate(self%cell_offsets)
    if (allocated(self%cell_vertices))      deallocate(self%cell_vertices)
    if (allocated(self%cell_types))         deallocate(self%cell_types)
    if (allocated(self%face_offsets))       deallocate(self%face_offsets)
    if (allocated(self%face_vertices))      deallocate(self%face_vertices)
    if (allocated(self%ls_stencil_offset))  deallocate(self%ls_stencil_offset)
    if (allocated(self%ls_stencil))         deallocate(self%ls_stencil)
    if (allocated(self%ls_weights))         deallocate(self%ls_weights)
    if (allocated(self%ls_stencil_face))    deallocate(self%ls_stencil_face)
    if (allocated(self%global_id))          deallocate(self%global_id)
    if (allocated(self%compact_id))         deallocate(self%compact_id)
    if (allocated(self%halos))              deallocate(self%halos)
    if (allocated(self%halo_send_off_cells)) deallocate(self%halo_send_off_cells)
    if (allocated(self%halo_recv_off_cells)) deallocate(self%halo_recv_off_cells)
    if (allocated(self%halo_reqs))           deallocate(self%halo_reqs)
    if (allocated(self%halo_sbuf_scalar))    deallocate(self%halo_sbuf_scalar)
    if (allocated(self%halo_rbuf_scalar))    deallocate(self%halo_rbuf_scalar)
    if (allocated(self%halo_sbuf_vec))       deallocate(self%halo_sbuf_vec)
    if (allocated(self%halo_rbuf_vec))       deallocate(self%halo_rbuf_vec)
    if (allocated(self%halo_sbuf_scalars))   deallocate(self%halo_sbuf_scalars)
    if (allocated(self%halo_rbuf_scalars))   deallocate(self%halo_rbuf_scalars)
    self%halo_layout_ready = .false.
    self%halo_exchange_active = .false.
    self%halo_send_cells_total = 0
    self%halo_recv_cells_total = 0
  end subroutine mesh_destroy

  subroutine ensure_halo_layout(self)
    class(Mesh), intent(inout) :: self
    integer(ip) :: h, nh

    if (self%halo_layout_ready) return
    if (.not. allocated(self%halos)) return

    nh = size(self%halos)
    if (nh <= 0) then
      self%halo_layout_ready = .true.
      return
    end if

    allocate(self%halo_send_off_cells(nh), self%halo_recv_off_cells(nh))
    self%halo_send_cells_total = 0
    self%halo_recv_cells_total = 0
    do h = 1, nh
      self%halo_send_off_cells(h) = self%halo_send_cells_total
      self%halo_recv_off_cells(h) = self%halo_recv_cells_total
      self%halo_send_cells_total = self%halo_send_cells_total + self%halos(h)%n_send
      self%halo_recv_cells_total = self%halo_recv_cells_total + self%halos(h)%n_recv
    end do

    allocate(self%halo_reqs(max(1_ip, 2*nh)))
    self%halo_layout_ready = .true.
  end subroutine ensure_halo_layout

  subroutine ensure_real_buffer_size(buf, needed)
    real(rp), allocatable, intent(inout) :: buf(:)
    integer(ip),           intent(in)    :: needed
    integer(ip) :: nalloc

    nalloc = max(1_ip, needed)
    if (.not. allocated(buf)) then
      allocate(buf(nalloc))
      return
    end if
    if (size(buf) < nalloc) then
      deallocate(buf)
      allocate(buf(nalloc))
    end if
  end subroutine ensure_real_buffer_size

  ! Helper: invert a symmetric 3x3 matrix via cofactors
  pure subroutine invert3x3(A, Ainv, det)
    real(rp), intent(in)  :: A(3,3)
    real(rp), intent(out) :: Ainv(3,3), det
    real(rp) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, idet

    c11 =  A(2,2)*A(3,3) - A(2,3)*A(3,2)
    c12 = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    c13 =  A(2,1)*A(3,2) - A(2,2)*A(3,1)
    c21 = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    c22 =  A(1,1)*A(3,3) - A(1,3)*A(3,1)
    c23 = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    c31 =  A(1,2)*A(2,3) - A(1,3)*A(2,2)
    c32 = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    c33 =  A(1,1)*A(2,2) - A(1,2)*A(2,1)

    det = A(1,1)*c11 + A(1,2)*c12 + A(1,3)*c13
    ! Sign-preserving guard: det + SMALL can still be near zero for negative det.
    idet = 1.0_rp / sign(max(abs(det), SMALL), det)

    Ainv(1,1) = c11*idet; Ainv(1,2) = c21*idet; Ainv(1,3) = c31*idet
    Ainv(2,1) = c12*idet; Ainv(2,2) = c22*idet; Ainv(2,3) = c32*idet
    Ainv(3,1) = c13*idet; Ainv(3,2) = c23*idet; Ainv(3,3) = c33*idet
  end subroutine invert3x3

  !> Assign compact global IDs and exchange for ghost cells
  subroutine mesh_assign_compact_ids(self)
    use mpi_f08
    class(Mesh), intent(inout) :: self
    integer(ip) :: i, h, ierr, nh
    integer(lp) :: ownership_start, n_cells_lp
    integer(lp), allocatable :: sendbuf(:), recvbuf(:)
    integer(ip) :: tag
    if (allocated(self%compact_id)) deallocate(self%compact_id)
    allocate(self%compact_id(self%n_cells + self%n_ghost))

    ! Scan an int64 copy: MPI_INTEGER8 on a 4-byte n_cells would read past it into n_ghost
    n_cells_lp = int(self%n_cells, lp)
    call MPI_Scan(n_cells_lp, ownership_start, 1, MPI_INTEGER8, MPI_SUM, self%comm, ierr)
    ownership_start = ownership_start - n_cells_lp + 1
    do i = 1, self%n_cells
      self%compact_id(i) = ownership_start + i - 1
    end do

    ! Exchange compact IDs for ghost cells
    nh = 0
    if (allocated(self%halos)) nh = size(self%halos)
    if (nh > 0) then
      allocate(sendbuf(maxval(self%halos(:)%n_send)))
      allocate(recvbuf(maxval(self%halos(:)%n_recv)))
      tag = 9001
      do h = 1, nh
        block
          type(MPI_Request) :: req_recv, req_send
          type(MPI_Request) :: req_pair(2)
          call MPI_Irecv(recvbuf, self%halos(h)%n_recv, MPI_INTEGER8, &
                         self%halos(h)%rank, tag, self%comm, req_recv, ierr)
          do i = 1, self%halos(h)%n_send
            sendbuf(i) = self%compact_id(self%halos(h)%send_cells(i))
          end do
          call MPI_Isend(sendbuf, self%halos(h)%n_send, MPI_INTEGER8, &
                         self%halos(h)%rank, tag, self%comm, req_send, ierr)
          req_pair(1) = req_recv
          req_pair(2) = req_send
          call MPI_Waitall(2, req_pair, MPI_STATUSES_IGNORE, ierr)
        end block
        do i = 1, self%halos(h)%n_recv
          self%compact_id(self%halos(h)%recv_cells(i)) = recvbuf(i)
        end do
      end do
      deallocate(sendbuf)
      deallocate(recvbuf)
    end if
  end subroutine mesh_assign_compact_ids

end module mod_mesh
