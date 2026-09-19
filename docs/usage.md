---
title: Usage
nav_order: 3
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# DOPAMINE Usage Guide

## Build

### Standard build

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DDOPAMINE_ENABLE_PETSC=ON \
  -DDOPAMINE_PETSC_INSTALL_DIR=/path/to/petsc \
  -DDOPAMINE_PETSC_ARCH_DIR=/path/to/petsc-arch
make -j
```

**Cleaning:** `make clean` only removes DOPAMINE's own build artifacts (object
files, `libdopamine_lib.a`, the executables) — it leaves the fetched/built
PETSc and hypre install alone, since that takes minutes to rebuild. Use
`make cleanall` to also wipe the external libraries (`external/` by default,
or `DOPAMINE_EXTERNAL_DIR` if set) when you actually want a from-scratch
external rebuild — e.g. after changing `DOPAMINE_PETSC_WITH_CUDA` or another
option that changes how PETSc itself is configured.

### Build options

`cmake -D<option>=<value>`:

| Option | Default | Effect |
|---|---|---|
| `DOPAMINE_ENABLE_OPENMP` | `ON` | OpenMP threading in the Fortran loops |
| `DOPAMINE_ENABLE_PETSC` | `ON` | PETSc pressure backend when available |
| `DOPAMINE_FETCH_PETSC` | `ON` | Download and build PETSc automatically |
| `DOPAMINE_PETSC_WITH_HYPRE` | `ON` | Build PETSc with hypre (BoomerAMG) |
| `DOPAMINE_PETSC_WITH_CUDA` | `OFF` | Build PETSc with CUDA support |
| `DOPAMINE_PETSC_USE_GPU_AWARE_MPI` | `ON` | Request GPU-aware MPI in CUDA builds |
| `DOPAMINE_ENABLE_OPENACC` | `OFF` | OpenACC GPU offload; requires NVHPC/`nvfortran` |
| `DOPAMINE_OPENACC_GPU_ARCH` | *(empty)* | Optional `-gpu=` target, e.g. `cc86` |
| `DOPAMINE_FETCH_NETCDF` | `ON` | Auto-build NetCDF C/Fortran when missing |
| `DOPAMINE_BUILD_TESTS` | `ON` | Build the CTest unit tests |
| `DOPAMINE_ENABLE_LTO` | `OFF` | Interprocedural/link-time optimization |
| `DOPAMINE_EXTERNAL_DIR` | `<source>/external` | Where fetched dependencies are built/installed |

### Offline / HPC (no internet) build

If your machine does not have internet access, download and build dependencies
on another machine, then transfer them.

| Library         | Official Repository / Download Link                                                                 |
|-----------------|---------------------------------------------------------------------------------------------------|
| CMake           | https://github.com/Kitware/CMake/releases                                                         |
| MPI             | [OpenMPI](https://github.com/open-mpi/ompi) or [MPICH](https://github.com/pmodels/mpich)          |
| PETSc           | https://gitlab.com/petsc/petsc/-/releases                                                         |
| NetCDF-C        | https://github.com/Unidata/netcdf-c/releases                                                      |
| NetCDF-Fortran  | https://github.com/Unidata/netcdf-fortran/releases                                                |
| GCC (Fortran/C) | https://gcc.gnu.org/releases.html                                                                 |

Set environment variables so the build system can find your local installations:

```sh
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
export NETCDF_DIR=/path/to/netcdf
export NETCDFF_DIR=/path/to/netcdf-fortran
export PATH=/path/to/cmake/bin:$PATH
export PATH=/path/to/mpi/bin:$PATH
export LD_LIBRARY_PATH=/path/to/your/libs:$LD_LIBRARY_PATH
```

Configure and build as above, using `-DCMAKE_PREFIX_PATH` if needed.

#### PETSc standalone build

Standalone PETSc build instructions to link (important: `--download-hypre=1`
and `--download-fblaslapack=1` are needed for DOPAMINE):

```sh
git clone -b release https://gitlab.com/petsc/petsc.git && cd petsc

mkdir petsc-install

./configure \
  --prefix=$HOME/petsc/petsc-install \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --download-hypre=1 \
  --download-fblaslapack=1 \
  --with-debugging=0 \
  COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3"

make all
make install
```

### GPU (CUDA) build

Single-GPU support for the pressure solve (the dominant cost — roughly 83% of
step time on the CPU path) works as follows: PETSc's pressure matrix/vectors
run as `MATAIJCUSPARSE`/`VECCUDA`, and BoomerAMG (hypre) runs its setup and
solve on-device. In this build nothing outside the pressure solve runs on
GPU; the explicit kernels are offloaded separately by the OpenACC build
described below. Validated on a single-GPU workstation (RTX 3060, CUDA 12.0);
multi-GPU/multi-node has not been exercised. See
`docs/external-lib-guidance/PETSc_GPU_porting.md` (local, gitignored) for
implementation notes.

#### Requirements

- An NVIDIA GPU and `nvcc` on `PATH`.
- **A CUDA-aware MPI**, if you want `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=1`
  (device-to-device transfers without host staging). Most distro-packaged MPI
  (Ubuntu/Debian `libopenmpi-dev`, etc.) is **not** built with CUDA support —
  check with `ompi_info --parsable | grep cuda`; if `mca:btl:smcuda` isn't
  listed, build your own:

  ```bash
  curl -fsSLo openmpi.tar.gz https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz
  tar xzf openmpi.tar.gz && cd openmpi-4.1.6
  ./configure --prefix="$HOME/opt/openmpi-cuda-4.1.6" --with-cuda=/usr CC=gcc CXX=g++ FC=gfortran
  make -j && make install
  ```

  (`--with-cuda=/usr` works when CUDA is the distro `nvidia-cuda-toolkit`
  package, where headers/libs sit under the normal `/usr/include` and
  `/usr/lib/<arch>` search paths rather than a unified `/usr/local/cuda`
  tree — check `find /usr -iname cuda_runtime.h` if `configure` can't find
  it.)

  Without a CUDA-aware MPI, run with `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0`
  (PETSc's CUDA Mat/Vec objects still work; halo/scatter traffic is staged
  through the host instead of going device-to-device).

#### Configure and build

Use a separate build directory and external-libs prefix so this doesn't
disturb an existing CPU build sharing the same source tree
(`DOPAMINE_EXTERNAL_DIR` defaults to `<source>/external` regardless of build
directory):

```bash
cmake -S . -B build-gpu \
  -DCMAKE_PREFIX_PATH="$HOME/opt/openmpi-cuda-4.1.6" \
  -DDOPAMINE_EXTERNAL_DIR="$(pwd)/external-gpu" \
  -DDOPAMINE_ENABLE_PETSC=ON \
  -DDOPAMINE_FETCH_PETSC=ON \
  -DDOPAMINE_PETSC_WITH_HYPRE=ON \
  -DDOPAMINE_PETSC_WITH_CUDA=ON \
  -DDOPAMINE_PETSC_USE_GPU_AWARE_MPI=ON
cmake --build build-gpu --target dopamine -j
```

Omit `-DCMAKE_PREFIX_PATH` (or point it at the system MPI) if you're running
with `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0`.

`DOPAMINE_PETSC_WITH_CUDA=OFF` (the default) forces `--with-cuda=0` on
PETSc's own configure — without that, PETSc auto-probes for `nvcc` and
silently builds CUDA support anyway on any machine that has the toolkit
installed, which then fails to link (no `-lcudart`/`-lstdc++` added to the
executables). If you're switching an existing tree between CUDA and non-CUDA
(or changed any other PETSc-configure-affecting option), run `make cleanall`
first — a plain `make clean` deliberately leaves the external PETSc/hypre
install alone (see "Cleaning" above), so switching `DOPAMINE_PETSC_WITH_CUDA`
without a full external rebuild would otherwise link against the previous
build's PETSc.

### OpenACC (nvfortran) build

Offloading the explicit kernels — advection, viscous, SGS, gradients — with
OpenACC needs `nvfortran`, since gfortran doesn't support OpenACC. `!$acc`
directives cover the SGS kernels (Vreman, Smagorinsky, dynamic Smagorinsky),
the gradient kernels (Green-Gauss and least-squares, scalar and vector), the
viscous fluxes (momentum and scalar diffusion), and the advection kernels
(CDS, QUICK/TVD reconstruction, and KEP), all driven through the mesh's flat
SoA mirror. See the "OpenACC (nvfortran) Run" section below for the kernel
list, and `docs/external-lib-guidance/PETSc_GPU_porting.md` (local,
gitignored) for implementation notes.

> **Known quirk**: a bare `cmake -S . -B build-acc26` (or `cmake .` from
> inside the build directory) on an *already-configured* OpenACC build tree
> has been observed to silently reset the whole toolchain — compiler,
> `DOPAMINE_ENABLE_OPENACC`, `DOPAMINE_PETSC_WITH_CUDA` — back to plain
> gfortran/CPU defaults, with no error (root cause not fully understood; see
> the porting-plan doc). The build still "succeeds" afterwards, just silently
> without GPU support, so check `grep CMAKE_Fortran_COMPILER: build-acc26/CMakeCache.txt`
> resolves to `nvfortran`, not `/usr/bin/f95`, if you're ever unsure. Always
> reconfigure with the **full** `cmake -S . -B build-acc26 ...` flag list
> below (never a bare `cmake .`); if the cache still looks wrong afterwards,
> delete the build directory (`rm -rf build-acc26`, **not** `external-acc26`
> — that's the expensive PETSc/hypre/CUDA build, kept separate on purpose)
> and reconfigure from scratch.

#### Requirements

- NVHPC (NVIDIA HPC SDK) providing `nvfortran`. A **user-local, no-root
  install** works fine and sidesteps apt entirely:

  ```bash
  wget https://developer.download.nvidia.com/hpc-sdk/26.5/nvhpc_2026_265_Linux_x86_64_cuda_13.2.tar.gz
  tar xpzf nvhpc_2026_265_Linux_x86_64_cuda_13.2.tar.gz
  cd nvhpc_2026_265_Linux_x86_64_cuda_13.2
  NVHPC_SILENT=true NVHPC_INSTALL_TYPE=single NVHPC_INSTALL_DIR="$HOME/nvhpc_26" ./install
  ```

  (~6.85 GB download; check the [NVIDIA HPC SDK downloads page](https://developer.nvidia.com/hpc-sdk-downloads)
  for the current version/URL, since these change.) Use a recent version —
  NVHPC 23.3 hit two internal compiler errors (ICEs) that a newer release
  (26.5) doesn't reproduce; if you're stuck on an old NVHPC install, expect
  to need the per-file `-O1 -Mnovect` / `-Mnoopenmp` workarounds already in
  `CMakeLists.txt` to still apply, and possibly more.

- **A Fortran-compiler-matched MPI.** Fortran `.mod` files are
  compiler-specific (not just ABI — the format itself differs), so
  gfortran-built MPI (the system one, or the custom one built for the CUDA
  GPU path above) cannot be used with `nvfortran`'s `mpi_f08`. NVHPC bundles
  its own CUDA-aware MPI (HPC-X) specifically for this — use it, don't try to
  reuse an existing MPI build:

  ```bash
  export PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/hpcx/bin:$HOME/nvhpc_26/Linux_x86_64/26.5/compilers/bin:$PATH"
  ```

#### Configure and build

Separate build directory and external-libs prefix again, same reasoning as
the CUDA build above:

```bash
NVBIN="$HOME/nvhpc_26/Linux_x86_64/26.5/compilers/bin"
MPIBIN="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/hpcx/bin"
export PATH="$MPIBIN:$NVBIN:$PATH"

cmake -S . -B build-acc26 \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_C_COMPILER=nvc \
  -DCMAKE_CXX_COMPILER=nvc++ \
  -DDOPAMINE_PETSC_CUDAC="$NVBIN/nvcc" \
  -DDOPAMINE_EXTERNAL_DIR="$(pwd)/external-acc26" \
  -DDOPAMINE_ENABLE_PETSC=ON \
  -DDOPAMINE_FETCH_PETSC=ON \
  -DDOPAMINE_PETSC_WITH_HYPRE=ON \
  -DDOPAMINE_PETSC_WITH_CUDA=ON \
  -DDOPAMINE_PETSC_USE_GPU_AWARE_MPI=ON \
  -DDOPAMINE_ENABLE_OPENACC=ON \
  -DDOPAMINE_OPENACC_GPU_ARCH=cc86
cmake --build build-acc26 --target dopamine -j
```

(`cc86` targets an Ampere GPU, e.g. RTX 3060 — swap for your card's compute
capability, or drop `-DDOPAMINE_OPENACC_GPU_ARCH` entirely to let `nvfortran`
pick a default.)

`-DDOPAMINE_PETSC_CUDAC="$NVBIN/nvcc"` matters: point PETSc at NVHPC's own
bundled `nvcc`, not the system one — mixing CUDA toolkit versions between
PETSc's CUDA compiler and the project's C++ compiler produces ambiguous
operator-overload errors deep in PETSc's CUDA complex-number code (two
different Thrust header versions visible at once).

`-DDOPAMINE_PETSC_WITH_CUDA=ON` + `-DDOPAMINE_ENABLE_OPENACC=ON` together
means the final `dopamine` link pulls in both PETSc's CUDA/cupm objects (for
the GPU pressure solve) and the OpenACC runtime (for the offloaded kernels
above) — two independent GPU paths sharing one binary. If you only want one
of them, drop the other's flag (`DOPAMINE_PETSC_WITH_CUDA=OFF` for an
OpenACC-only build with a CPU pressure solve, or `DOPAMINE_ENABLE_OPENACC=OFF`
for the reverse); either combination builds and links cleanly.

CMake automatically adds NVHPC's `math_libs/lib64` (derived from the Fortran
compiler's own path) to the link search when `DOPAMINE_PETSC_WITH_CUDA=ON`
under NVHPC — this is what makes `cublas`/`cusparse`/`cusolver` symbols
resolve at final link. `petsc.pc`'s own `Libs.private` lists those as bare
`-lcublas`-style entries with no matching `-L` for where NVHPC actually keeps
them (it ships them alongside `nvfortran`, not next to `nvcc` in a
conventional CUDA-toolkit layout), so without this the link fails with
`undefined reference to cublasDtrsv_v2` (or `cublasDgemv_v2`,
`cublasGetPointerMode_v2`, etc.) — a chain of dozens of near-identical
messages from `libpetsc.a`'s `cd_cupm.o`/`blas_cyclic_cupm.o`. If you still
hit that error, your NVHPC install layout differs from the expected
`<root>/compilers/bin/nvfortran` + `<root>/math_libs/lib64`; CMake prints a
`WARNING` naming the directory it looked for — point
`DOPAMINE_PETSC_CUDAC`/`CMAKE_Fortran_COMPILER` at the right install, or add
the correct `math_libs/lib64` path yourself via
`-DCMAKE_EXE_LINKER_FLAGS=-L/path/to/math_libs/lib64`.

## Run

Single rank:

```bash
./build/dopamine examples/taylor_green/taylor_green.input
```

MPI run:

```bash
mpirun -np 8 ./build/dopamine examples/channel_flow_retau395/channel_retau395.input
```

Hybrid MPI + OpenMP run:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./build/dopamine examples/channel_flow_retau395/channel_retau395.input
```

## Checkpoint and Restart

**Enable checkpoint writing:**

```ini
&dopamine_checkpoint
  checkpoint_interval = 5000    ! Write every 5000 steps
/
```

Generates `output_prefix/processorR/checkpoint_XXXXXXXX.bin` for each rank
plus a manifest file `checkpoint_latest`.

**Restart from checkpoint (same rank count):**

```ini
&dopamine_checkpoint
  restart_enabled = .true.
  restart_dir     = 'output1'
  restart_step    = -1          ! Auto-read checkpoint_latest
/
```

**Repartition checkpoints to a different rank count** (e.g., 4 → 6 ranks):

```bash
# Step 1: Pre-decompose mesh for 6 ranks
mpirun -np 1 ./dopamine-decompose my_input.input --nranks 6 --output-dir mesh_6ranks

# Step 2: Repartition checkpoint from old decomposition to new
dopamine-checkpoint-partition --old-dir output_4ranks --old-nranks 4 \
                              --new-dir mesh_6ranks --new-nranks 6 \
                              --step 10000

# Step 3: Create new input file pointing to 6-rank mesh and checkpoint
# (Update mesh_file, mesh_type='preprocessed', restart_dir, etc.)

# Step 4: Run 6-rank simulation from repartitioned checkpoint
mpirun -np 6 ./dopamine input_6ranks.input
```

## GPU (CUDA) Run

In this build only the pressure solve (PETSc + hypre BoomerAMG) runs on GPU;
the explicit kernels are offloaded by the separate OpenACC build below. See
the "GPU (CUDA) build" section above for build instructions and MPI
requirements. Once built:

```bash
export LD_LIBRARY_PATH="$HOME/opt/openmpi-cuda-4.1.6/lib:$LD_LIBRARY_PATH"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=1 \
  mpirun -np 1 ./build-gpu/dopamine examples/channel_flow_retau395/channel_retau395.input
```

Set `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` if your MPI isn't CUDA-aware (most
distro packages aren't — see "GPU (CUDA) build" above); PETSc aborts at
startup rather than silently falling back if you leave it on without one.
Compare the `pressure`
line in the performance snapshot against a `DOPAMINE_PETSC_USE_GPU=0` run on
the same binary to see the actual speedup, since that's the step this currently
accelerates.

## OpenACC (nvfortran) Run

`!$acc` directives cover:

- `vreman_sgs_viscosity`, `smagorinsky_sgs` (`mod_vreman.f90`) and
  `dynamic_smagorinsky_sgs` (`mod_dynamic_smagorinsky.f90`).
- `green_gauss_gradient_scalar`/`vector` and `ls_gradient_scalar`/`vector`
  (`mod_gradient.f90`).
- `viscous_rhs_vector`/`viscous_rhs_scalar` (`mod_viscous.f90`, momentum and
  scalar diffusion).
- The advection kernels (`mod_advection.f90`): the CDS, reconstruction
  (QUICK/van Leer/Koren TVD) and KEP scalar/vector SoA paths.

All of these read geometry through the mesh's flat SoA mirror (cell centroid
and volume, face owner/neighbour, normal/area/centre/`d_cf`, interpolation and
skewness factors, least-squares stencils, patch BC ids and periodic pairing
maps), which `main.f90` stages on the device once with a single
`!$acc enter data` so each kernel's own `copyin` finds it already resident. The
mesh is static, so this never needs refreshing. Boundary-condition application
and halo exchange remain host-side. On a single rank the solver keeps `U`,
`grad_U`, `nu_sgs` and the scalars device-resident across the SGS/gradient/
viscous block of each RK3 stage; under MPI (`nranks > 1`) that residency is
disabled and each kernel transfers its own data, since the halo exchange
between stages runs on the host. See the "OpenACC (nvfortran) build" section
above for the build recipe (`build-acc26/`).

```bash
export PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/hpcx/bin:$HOME/nvhpc_26/Linux_x86_64/26.5/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/cuda/13.2/targets/x86_64-linux/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/math_libs/lib64:$LD_LIBRARY_PATH"
export OPAL_PREFIX="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0 mpirun -np 1 ./build-acc26/dopamine examples/channel_flow_retau395/channel_retau395.input
```

The OpenACC path is validated on a single-GPU workstation against a
`channel_flow_retau395`-style case (73728 cells, `-np 1`): `max|divU|` holds at
machine precision (~1e-11) every logged step with `p_resid_max` ~1.4e-9, the
same correctness signature as the CUDA-only run in the "GPU (CUDA) build"
section above. See `OPAL_PREFIX` and
`DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` above — both needed on this HPC-X install
(see the paragraphs above for why). If it fails to even start with a
driver-version error, see the "OpenACC (nvfortran) build" section above for
the fallback (pair `nvfortran` 26.5 with an older `nvcc` via
`DOPAMINE_PETSC_CUDAC`).

## Choosing a Mesh Path

`mesh_type` selects between two quite different startup paths. The choice is
about where the global mesh has to fit, not about accuracy — both produce the
same partitions.

| | Direct read (`gmsh`, `openfoam`, `icem`, `cgns`) | Preprocessed |
|---|---|---|
| Preparation | none | one offline `dopamine-decompose` run |
| Global mesh held by | rank 0, entirely | nobody at solver time |
| Startup | re-partitioned every run | each rank reads its own file |
| Rank count | change freely | fixed at decomposition time |
| Practical range | up to ~10M cells | above that |

**Direct read** is the default and the right choice while the mesh still fits
on one rank. Rank 0 reads the whole mesh, partitions it and scatters; nothing to
prepare and you can vary `-np` between runs.

```fortran
&dopamine_domain
  mesh_file = 'channel.msh'
  mesh_type = 'gmsh'
/
```

`mesh_type = 'icem'` (alias `'fluent'`) reads a Fluent-format `.msh` file, as
exported by Ansys ICEM CFD — including polyhedral meshes. Both ASCII and
binary variants are supported. Boundary patch names/types come from the
mesh's own named zones (or a numeric Fluent BC-type fallback if unnamed) and
are matched against `&dopamine_bc` patch names exactly as with `gmsh`/
`openfoam`.

```fortran
&dopamine_domain
  mesh_file = 'channel.msh'
  mesh_type = 'icem'
/
```

`mesh_type = 'cgns'` reads a CGNS (CFD General Notation System) file via the
CGNS reference library (HDF5 and legacy ADF files). Two layouts are supported,
both single-zone: polyhedral `NGON_n`/`NFACE_n` sections with named
`ZoneBC_t`/`BC_t` patches at `FaceCenter` grid location, and standard
`TETRA_4`/`HEXA_8`/`PENTA_6`/`PYRA_5` volume sections with `TRI_3`/`QUAD_4`
boundary sections (patches take the boundary section names, BC types come from
the matching `BC_t`). `Family_t`-based boundary naming is not yet supported.
Requires the CGNS library, which is fetched and built automatically
(`DOPAMINE_FETCH_CGNS=ON`, the default) alongside the HDF5 dependency the
output writer already requires.

```fortran
&dopamine_domain
  mesh_file = 'channel.cgns'
  mesh_type = 'cgns'
/
```

**Preprocessed** removes the rank-0 ceiling. Decompose once, then every rank
reads only its own `processorN/mesh.bin`:

```bash
mpirun -np 1 ./build/dopamine-decompose channel.input --nranks 128 \
  --output-dir ./channel_128 --threads 16
```

```fortran
&dopamine_domain
  mesh_file = './channel_128'
  mesh_type = 'preprocessed'
/
```

Run it with exactly the rank count it was decomposed for
(`mpirun -np 128`); changing `-np` means decomposing again.

Switch to preprocessed when rank 0 runs out of memory during startup, when the
mesh is beyond roughly 10M cells, or when you will launch the same mesh many
times at a fixed rank count and want the startup cost paid once.

### Decomposer options

| Flag | Meaning |
|---|---|
| `--nranks N` | partitions to write (required) |
| `--output-dir DIR` | destination for `processorN/mesh.bin` |
| `--threads N` | writer threads; defaults to all cores |
| `--max-memory-gb G` | cap peak memory by staging face geometry to disk |
| `--chunk-size M` | face scan chunk, default 2M |
| `--weight-by-volume` | balance total volume per rank instead of cell count |

Partitions carry an equal number of cells by default, because solver cost scales
with cell count. `--weight-by-volume` restores the older behaviour of equalising
total volume, which on a wall-graded channel or ABL mesh leaves the near-wall
ranks holding several times more cells than average — roughly 3.5x at 1024
ranks on a tanh-graded mesh. Leave it off unless you specifically want it.

The decomposer itself runs on a single core (`-np 1`); `--threads` parallelises
the write phase only. It prints the resulting balance, which is the number to
check:

```
Cells/rank  min=130918  max=131204  imbalance=1.001
Load balance: equal cell count
```

Note that the direct-read path has always balanced by cell count, so this
setting only affects meshes that go through the decomposer.

## Supported Mesh Formats

- **Gmsh v4 (.msh):** tet4, hex8, prism6, pyr5, tet10, hex20. Use
  `$PhysicalNames` for boundaries. Generate with
  `gmsh -3 -format msh4 geometry.geo`.
- **OpenFOAM polyMesh:** Reads `constant/polyMesh/{points,faces,owner,
  neighbour,boundary}`. All OpenFOAM-supported cell types.

## Input Structure

The solver reads these namelists:

- `&dopamine_domain`
- `&dopamine_numerics`
- `&dopamine_physics`
- `&dopamine_abl` (optional, for atmospheric boundary layer)
- `&dopamine_bc`
- `&dopamine_initial_condition`
- `&dopamine_output`
- `&dopamine_checkpoint` (optional)
- `&dopamine_stats` (optional)
- `&dopamine_PetscParams` (optional)
- `&dopamine_scalars` (optional)
- `&dopamine_scalar_sources` (optional)
- `&dopamine_sampling_planes` (optional)
- `&dopamine_sampling_lines` (optional)
- `&dopamine_sampling_curves` (optional)
- `&dopamine_sampling_surfaces` (optional)

See `reference.input` for the complete option list.

## Boundary Condition Types

`patch_type(p)` in `&dopamine_bc` selects the treatment for patch `p`. The
accepted names (aliases separated by `/`) are:

| `patch_type` | Treatment |
|---|---|
| `velocity` / `inlet` / `fixedValue` | Dirichlet velocity inlet |
| `outlet` / `fixedPressure` / `outflow` | Outflow with pressure anchoring |
| `zeroGradient` | Homogeneous Neumann on all fields |
| `noSlip` | Stationary no-slip wall |
| `movingWall` / `wallVelocity` | No-slip wall with prescribed wall velocity |
| `wallModel` / `wall` | Wall-modelled wall (log law or MOST, per `wall_model`) |
| `slip` | Free-slip wall |
| `symmetry` | Symmetry plane |
| `periodic` / `cyclic` | Periodic pairing with the matching patch |
| `fixedTemperature` | Dirichlet temperature only; this type sets no velocity BC |
| `ablInlet` | ABL inflow profile, optionally with SEM fluctuations |

An unrecognised name falls back to `zeroGradient` rather than aborting, so
check the startup log if a patch does not behave as expected.

## Output and Visualization

For `prefix = 'case'`, the main outputs are:

- `case.nc`: merged unstructured mesh + time series cell fields.
- `case.xmf`: temporal XDMF sidecar for ParaView.
- `case_stats.csv`: runtime scalar statistics (when enabled).

If advanced sampling is enabled, additional files are emitted:

- `case_sample_plane_*.nc`
- `case_sample_line_*.nc`
- `case_sample_curve_*.nc`
- `case_sample_surface_*.nc`

Load `case.xmf` in ParaView for mesh and cell fields. For decomposed MPI
output, apply these filters in order:

1. `Merge Blocks` (work on merged grid)
2. `Clean to Grid`
3. `Cell Data to Point Data`
4. Color by velocity magnitude (`|U|`)

## Advection Scheme Names

Use the following names in `advection_scheme` (momentum) and `temperature_scheme` (temperature):

- `CENTRAL` or `CDS`: linear central interpolation (second order on uniform orthogonal regions).
- `KEP` or `SKEW`: midpoint momentum convective state for lower dissipation and better resolved kinetic-energy behavior in LES-oriented runs.
- `QUICK`: geometry-aware quadratic upwind reconstruction evaluated in projected streamline coordinates. On aligned uniform meshes it reduces to classical QUICK [@leonard1979].
- `VANLEER_TVD`: bounded MUSCL reconstruction with van Leer limiter using geometry-aware projected slopes [@vanleer1974; @sweby1984].
- `KOREN_TVD`: bounded MUSCL reconstruction with Koren limiter using geometry-aware projected slopes [@koren1993; @sweby1984].

`advection_scheme` controls momentum; `temperature_scheme` controls the temperature (T) field independently. Both accept the same set of scheme names. Passive scalar advection is controlled separately by `scalar_tvd_scheme` in `&dopamine_scalars`.

## Recommended Numerics by Problem Type

- Transitional/turbulent wall-bounded flows:
  - `advection_scheme = 'KEP'`
  - `temperature_scheme = 'VANLEER_TVD'`
  - `time_integrator = 'RK3W'` [@williamson1980]
  - `adapt_dt = .true.` with conservative `cfl_max`
- Smooth periodic benchmark (Taylor-Green):
  - `KEP + RK3W` for robust long-time kinetic-energy behavior
- Thermally coupled buoyant flow:
  - `advection_scheme = 'KEP'`, `temperature_scheme = 'VANLEER_TVD'`, `time_integrator = 'RK3W'`
  - Monitor both momentum and scalar residual metrics

## Atmospheric Boundary Layer (ABL) Module

DOPAMINE includes an optional ABL module for realistic wind-driven and
thermally-stratified urban/environmental flows. Three mechanisms work together:

### Enabling ABL Features

The module has three independent pieces, and you can enable any combination:

| Feature | Flag | What it needs |
|---|---|---|
| Wall treatment | `wall_model = 'loglaw' \| 'most' \| 'none'` | a `wallModel` patch |
| Synthetic-eddy inflow | `use_sem = .true.` | an `ablInlet` patch |
| Davies nudging | `use_nudging = .true.` | a profile file |

None of them requires the others. SEM inflow on a smooth-wall channel is
`use_sem = .true.` with `wall_model = 'loglaw'`; a stratified urban case with no
synthetic inflow is `wall_model = 'most'` with `use_sem = .false.`.

`use_abl = .true.` is kept as a legacy alias meaning `wall_model = 'most'` when
`wall_model` is not given. It no longer gates SEM or nudging — previously it did,
so `use_sem = .true.` alone was silently inert.

In the input file, set the following in `&dopamine_abl`:

```fortran
&dopamine_abl
  wall_model = 'most'                 ! 'loglaw' (default) | 'most' | 'none'
  u_ref = 8.0                         ! Reference wind speed [m/s] at z_ref
  z_ref = 100.0                       ! Reference height [m]
  z0_m = 0.1                          ! Default momentum roughness [m]; per-patch override available
  z0_h = 0.01                         ! Default thermal roughness [m]
  T_surface = 300.0                   ! Surface temperature [K]
  T_ref = 300.0                       ! Reference/background temperature [K]
  L_obukhov = 0.0                     ! Initial Obukhov length [m] (ignored when compute_L=.true.)
  compute_L = .true.                  ! Compute Obukhov length dynamically each step
  wind_dir_deg = 270.0                ! Wind direction (meteorological convention)
  use_sem = .true.                    ! Enable Synthetic Eddy inflow
  n_eddies = 500                      ! Number of synthetic eddies
  sem_sigma = 0.0                     ! Eddy length scale [m]; 0 = auto (0.1*z_ref)
  use_nudging = .false.               ! Enable large-scale relaxation
  nudging_file = ''                   ! Path to weather profile (CSV/netCDF)
  nudging_tau = 1800.0                ! Relaxation timescale [s]
  nudging_z_start = 0.0               ! Height above which nudging applies [m]
/
```

### MOST Wall Model

With `wall_model = 'most'`, `wallModel` patches use the Monin-Obukhov
Similarity Theory (MOST) wall model instead of the neutral log-law model.
The model:

- Computes friction velocity (`u_*`) and temperature scale (`T_*`) dynamically
- Applies Businger-Dyer stability functions for stable, neutral, and unstable regimes
- Clips Obukhov length to $$[-5000, 5000]$$ m for numerical stability
- Couples heat flux to temperature field via implicit source

The Obukhov length $$L = u_*^2 T_{ref} / (\kappa g T_*)$$ is tracked **per boundary
face**, so each face responds to its own local stratification: $$L < 0$$ unstable
(warm surface), $$L > 0$$ stable. Each face reuses its converged value as the next
step's initial guess. The domain-mean `L_obukhov` reported by the solver is a
diagnostic; it seeds the very first call and is otherwise not used in the solve.

The MOST model is fully compatible with existing wall-model BC patches.
No mesh modifications required.

Temperature is transported only when `boussinesq = .true.`. With `boussinesq = .false.`
the run is neutral and isothermal: `T` stays at `T_ref`, so MOST reduces to the log law.

### Per-Patch Surface Roughness

`z0_m` and `z0_h` in `&dopamine_abl` set the domain-wide default. To vary roughness
across surface types — water, vegetation, terrain — set `patch_z0_m` and
`patch_z0_h` per patch in `&dopamine_bc`. These apply to `wallModel` patches under
MOST and to the ABL inlet profile:

```fortran
&dopamine_bc
  n_patches = 4
  patch_name(3) = 'lake'
  patch_type(3) = 'wallModel'
  patch_z0_m(3) = 2.0e-4               ! Open water
  patch_z0_h(3) = 2.0e-5

  patch_name(4) = 'forest'
  patch_type(4) = 'wallModel'
  patch_z0_m(4) = 1.0                  ! Forest canopy
  patch_z0_h(4) = -1.0                 ! Derive from z0_m via kB^-1
/
```

Conventions:

- `patch_z0_m = 0` (the default) inherits `z0_m` from `&dopamine_abl`; existing
  input files are unaffected.
- `patch_z0_h = 0` inherits `z0_h`; a **negative** value derives the thermal
  roughness from that patch's `z0_m` via $$\mathrm{k}B^{-1} = \ln(z_{0m}/z_{0h}) = 2$$,
  which is the usual approximation over vegetated land.

Representative `z0_m` values [m]:

| Surface            | `z0_m`      |
|--------------------|-------------|
| Open water         | 1e-4 – 1e-3 |
| Short grass        | 0.01 – 0.03 |
| Crops              | 0.05 – 0.25 |
| Forest             | 0.5 – 1.5   |
| Urban              | 0.5 – 2.0   |

Roughness spans four orders of magnitude across these classes, so it materially
changes surface drag — a single domain value over mixed terrain is usually the
largest error in an ABL setup. Note that open water is not strictly a constant
`z0_m` (the Charnock relation makes it depend on $$u_*$$); a fixed small value is a
reasonable approximation for moderate winds.

Roughness varying *within* a single patch (a land-cover raster) is not supported;
resolve distinct surfaces as separate mesh patches.

### Synthetic Eddy Method (SEM) for Turbulent Inflow

When `use_sem = .true.`:

- SEM generates $$N_e$$ eddies (default 500) in a recycle box upstream of inlet
- Eddies are convected downstream and create anisotropic turbulent fluctuations
- Target Reynolds stress anisotropy matches ABL log-layer profiles
- Eddy length scale $$\sigma$$ defaults to $$0.1 \times$$ `z_ref` (override with `sem_sigma`)

**Usage notes:**
- Apply SEM only to inlet patches; other BCs are unaffected
- SEM injects velocity fluctuations only; it does not require `boussinesq = .true.`
- Recycle box depth is estimated from mesh bounds automatically
- Eddy recycling is not enabled in the current SEM path

#### Prescribed Inflow Profiles

By default SEM builds its own mean profile (neutral log law from `u_ref`,
`z_ref`, `z0_m`) and its own stresses (Panofsky–Dutton ratios
$$\sigma_u/u_* = 2.4$$, $$\sigma_v/u_* = 1.9$$, $$\sigma_w/u_* = 1.25$$, constant with
height, no shear stress). To drive the inflow from measured, wind-tunnel or
WRF-derived data instead, point `sem_profile_file` at a CSV:

```fortran
&dopamine_abl
  use_sem = .true.
  sem_profile_file = 'inflow_profiles.csv'
/
```

Two layouts are accepted, chosen automatically from the column count.
**Six columns**, wind-aligned — the common case when you have a target profile
and turbulence intensities:

```
# z      U     uu     vv     ww     uw
10.0   5.20  1.020  0.560  0.290  -0.185
20.0   6.05  0.910  0.500  0.262  -0.170
...
```

`U` is a speed along `wind_dir_deg`, and the stresses are streamwise, lateral
and vertical; both are rotated into the wind direction together.

**Ten columns**, global axes — for data already in domain coordinates:

```
# z     U     V     W     uu     vv     ww     uv     uw     vw
10.0  5.20  0.31  0.00  1.020  0.560  0.290  0.041  -0.185  0.012
```

Here `wind_dir_deg` is *not* applied; the file defines the direction.

Conventions:

- Stresses are kinematic, m²/s² (the solver carries no density).
- Rows need not be sorted; they are ordered by height on read.
- Values are interpolated linearly in `z` and **clamped** outside the table, so
  extend the table to at least the domain top or the upper region inherits the
  last row.
- `uw` is the momentum flux and is the component that matters most. Supplying it
  is the main reason to use this path: the analytic default sets $$\overline{u'w'}=0$$,
  so the shear stress has to regenerate over fetch.
- A non-positive-definite tensor is clipped rather than aborting, so a noisy
  measured table will still run — check the profile if the realised stresses
  look low.

The prescribed mean **replaces** the analytic log law at `ablInlet` faces, so
`u_ref`, `z_ref` and `z0_m` no longer affect the inlet once a profile is loaded.
They still set the MOST wall model, which is a separate consumer.

The mean profile is imposed exactly. The Reynolds stresses are matched in the
statistical sense: fluctuations are built by applying the Cholesky factor of the
target tensor to a unit-variance eddy field, so the realised stresses approach
the target as the eddy count rises and the sample lengthens. Expect a few
percent scatter with `n_eddies = 1000`, and note the pressure projection adjusts
the field immediately downstream of the inlet plane.

### Davies Nudging (Large-Scale Coupling)

When `use_nudging = .true.`:

1. Load a weather profile file: CSV or netCDF format
   - **CSV:** Columns are `[time(s), z(m), u(m/s), v(m/s), w(m/s), T(K)]`, one row per point
   - **netCDF:** Variables `times(:)`, `z_levels(:)`, `U_profile(3,z,t)`, `T_profile(z,t)`

2. At each RHS computation, the solver relaxes toward these profiles:
   - **Below** `nudging_z_start`: free LES evolution
   - **Above** `nudging_z_start`: Newtonian relaxation with timescale `nudging_tau`

3. Typical use case:
   - Download WRF or ERA5 profile for your domain
   - Set `nudging_tau = 1800` s (weather model timescale)
   - Apply nudging above ABL height to couple LES to synoptic flow

### Example ABL Setup for Urban Canyon LES

```fortran
&dopamine_abl
  wall_model = 'most'
  u_ref = 10.0                        ! Urban reference wind
  z_ref = 250.0                       ! Above building height
  z0_m = 1.0                          ! Rough urban surface (default; override per patch)
  z0_h = 0.1                          ! Rough heat transfer
  T_surface = 310.0                   ! Warm surface (summer day)
  T_ref = 295.0                       ! Background temperature
  compute_L = .true.                  ! Track stability dynamically
  use_sem = .true.
  n_eddies = 1000                     ! More eddies for finer mesh
  use_nudging = .true.
  nudging_file = 'wrf_profile.csv'
  nudging_tau = 3600.0                ! 1 hour relaxation
  nudging_z_start = 100.0             ! Above buildings
/
```

## Pressure Projection and Linear Solve

DOPAMINE uses a projection method [@chorin1968; @kim1985] with a PETSc-backed pressure linear system path [@balay2021petsc].

Practical settings:

- Keep `p_tol` strict enough for long-time runs (`1e-8` is a robust default in many LES cases).
- Increase `p_max_iter` if pressure residuals stagnate before tolerance.
- In parallel runs, monitor both pressure residual trend and `max|divU|`.
- With a CUDA-enabled PETSc this linear solve can run on GPU — see "GPU (CUDA)
  Run" above. Custom `-pc_hypre_boomeramg_interp_type`/`coarsen_type` overrides
  via `&dopamine_PetscParams` should stick to hypre's GPU-supported values
  (`ext+i`/`PMIS`) when GPU mode is active; host-only choices segfault against a
  device-resident matrix.

## Long-Time Stability Checklist

1. Use the current executable (`build/dopamine`).
2. Check `n_ranks` in log header to confirm intended MPI launch.
3. For wall-model cases, monitor:
   - `max|U|`
   - `E`
   - `max|divU|`
   If only `max|U|` grows while global metrics remain bounded, investigate localized outliers.
4. Keep `p_tol` and `p_max_iter` strict enough for projection quality.
5. Use `solver_stats_frequency` consistently across cases.

## Minimal Numerics Example

```fortran
&dopamine_numerics
  advection_scheme = 'KEP'
  time_integrator = 'RK3W'
  dt = 5.0e-4
  cfl_max = 0.5
  adapt_dt = .true.
  p_tol = 1.0e-8
  p_max_iter = 1000
/
```

## Troubleshooting & Tips

- If CMake cannot find a library, check `*_DIR` variables and `CMAKE_PREFIX_PATH`.
- For runtime missing symbols, ensure `LD_LIBRARY_PATH` includes all required directories.
- Use strict `p_tol` (e.g., `1e-8`) for projection quality in LES.
- Monitor `max|divU|` and pressure residuals for stability.
- Use `solver_stats_frequency` to control statistics print interval.
