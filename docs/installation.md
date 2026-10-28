---
title: Installation
nav_order: 2
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Installation

How to build DOPAMINE on a workstation, an HPC cluster without internet access, and on GPUs.
After building, continue with [Usage](usage.md) to run a case, or the [Examples](examples.md).

## Prerequisites

| Requirement | Needed for | Notes |
|---|---|---|
| Fortran 2008 compiler | everything | `gfortran` by default; NVHPC `nvfortran` for OpenACC builds |
| C and C++ compilers | tools and shims | C11 and C++17 |
| MPI (Fortran, C, C++) | everything | OpenMPI or MPICH |
| CMake >= 3.22 | build | |
| PETSc (+ hypre) | pressure solve | Downloaded and built automatically when `DOPAMINE_FETCH_PETSC=ON` (needs internet) |
| NetCDF C/Fortran (+ HDF5) | output | Built automatically when missing (`DOPAMINE_FETCH_NETCDF=ON`) |
| OpenMP | optional | Threading in Fortran loops (`DOPAMINE_ENABLE_OPENMP`) |
| `gmsh` | meshes, regression tests | Used to generate and convert meshes |
| ParaView | optional | Visualisation of `.xmf` output |

## Verify the build

Unit tests are built on demand; regression cases need `gmsh` on `PATH`.

```sh
cmake --build build -j
ctest --test-dir build -L regression --output-on-failure
cmake --build build --target test_amd test_inflow_db test_inflow_gen test_inflow_planes test_sem_dfsem
ctest --test-dir build -L unit --output-on-failure
```

Test layout and how to add tests: [Development](development.md#testing).

## Troubleshooting

- CMake cannot find a library: check the matching `*_DIR` variable and `CMAKE_PREFIX_PATH`.
- Missing symbols at run time: make sure `LD_LIBRARY_PATH` contains every required library directory.
- A fetched dependency is stale after changing a PETSc option: run `make cleanall` and rebuild.
- No internet on the build machine: follow the offline build below.
- Illegal-instruction crash on compute nodes: the build tunes for the build machine; configure with `-DDOPAMINE_NATIVE_ARCH=OFF` or build on the target architecture.


## Standard build

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

## Build options

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
| `DOPAMINE_NATIVE_ARCH` | `ON` | Tune for the build machine (`-march=native`); `OFF` for portable binaries |
| `DOPAMINE_ENABLE_LTO` | `OFF` | Interprocedural/link-time optimization |
| `DOPAMINE_EXTERNAL_DIR` | `<source>/external` | Where fetched dependencies are built/installed |

## Offline / HPC (no internet) build

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

### PETSc standalone build

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

## GPU (CUDA) build

Single-GPU support for the pressure solve (the dominant cost — roughly 83% of
step time on the CPU path) works as follows: PETSc's pressure matrix/vectors
run as `MATAIJCUSPARSE`/`VECCUDA`, and BoomerAMG (hypre) runs its setup and
solve on-device. In this build nothing outside the pressure solve runs on
GPU; the explicit kernels are offloaded separately by the OpenACC build
described below. Validated on a single-GPU workstation (RTX 3060, CUDA 12.0)
and, in the combined OpenACC + CUDA build, on two RTX A6000 (one rank per GPU,
see "Multi-GPU" below); multi-node has not been exercised.

### Requirements

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

### Configure and build

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

## OpenACC (nvfortran) build

Offloading the explicit kernels — advection, viscous, SGS, gradients — with
OpenACC needs `nvfortran`, since gfortran doesn't support OpenACC. `!$acc`
directives cover the SGS kernels (Vreman, Smagorinsky, dynamic Smagorinsky),
the gradient kernels (Green-Gauss and least-squares, scalar and vector), the
viscous fluxes (momentum and scalar diffusion), and the advection kernels
(CDS, QUICK/TVD reconstruction, and KEP), all driven through the mesh's flat
SoA mirror. See the "OpenACC (nvfortran) Run" section below for the kernel
list.

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

### Requirements

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

### Configure and build

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

