<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

## DOPAMINE - Discrete Operator-based Parallel Algorithm for Modeling Incompressible Navier–Stokes Equations

**DOPAMINE** is an unstructured, cell-centered Large Eddy Simulation (LES) code for incompressible flows, written in modern Fortran with a PETSc-backed pressure solver and NetCDF/XMF output.

---

## Physics


Solves the incompressible **Boussinesq Navier-Stokes equations** in addition to multiple scalar quantities given by

$$
\partial_t U_i + \partial_j (U_j U_i) = -\frac{1}{\rho}\partial_i p + \partial_j \left[ (\nu + \nu_t) (\partial_j U_i + \partial_i U_j) \right] - \beta (T - T_r) g_i + F_i^{ext} ... ... ... ... (1)
$$

$$
\partial_j U_j = 0 ... ... ... ... (2) 
$$

$$
\partial_t T + \partial_j (U_j T) = \partial_j \left[ \left(\alpha + \frac{\nu_t}{Pr_t}\right) \partial_j T \right] ... ... ... ... (3) 
$$

$$
\partial_t C + \partial_j (U_j C) = \partial_j \left[ D_C \partial_j C \right] ... ... ... ... (4)
$$


where $U_i$ is the velocity vector, $t$ is time, $x_i$ is the space coordinate, $\rho$ is the density of the fluid, $p$ is the pressure, $\nu$ is the kinematic viscosity, $\nu_t$ is the turbulent/eddy viscosity, $\beta$ is the thermal expansion coefficient, $T$ is the temperature field, $T_r$ is the reference temperature field, $g_i$ is the gravitational acceleration vector, $F_i^{ext}$ can represent any additional external forcing, $\alpha$ is the thermal diffusivity constant, $Pr_t$ is the turbulent Prandtl number, $C$ is a generic passive scalar, and $D_C$ its (possibly spatially varying) diffusivity. In


Additional passive scalars can be enabled through `&dopamine_scalars`.
Each scalar uses the same finite-volume transport path with configurable
advection (`VANLEER_TVD`, `KOREN_TVD`, `QUICK`, `CDS`/`CENTRAL`) and
scalar-specific diffusion/closure coefficients.

$F_i^{ext}$ can include constant forcing, rotating-frame Coriolis forcing (`-2 Ω × U`), and harmonic forcing (`A cos(ωt+φ) d`).


### Subgrid-scale models
Selected via `sgs_model` in `&dopamine_numerics`.

**Vreman (2004)** (default) — zero near walls without a damping function:

$$
\nu_t = C_v \cdot \Delta^2 \cdot \sqrt{\frac{B_\beta}{\alpha_{ij}\alpha_{ij}}}
$$

where $\alpha_{ij} = \partial U_j / \partial x_i$, $B_\beta = \beta_{11}\beta_{22} - \beta_{12}^2 + \beta_{11}\beta_{33} - \beta_{13}^2 + \beta_{22}\beta_{33} - \beta_{23}^2$, $\beta_{ij} = \alpha_{ki}\alpha_{kj}$.

**Constant-coefficient Smagorinsky** (`SMAGORINSKY`) — $\nu_t = (C_s \Delta)^2 |S|$ with $C_s$ fixed at 0.1.

**Dynamic Smagorinsky** (`DYNAMIC_SMAG`) — Germano-Lilly dynamic procedure: $C_s^2$ is computed locally every step from a compact test-filtered Germano identity and clipped to $\geq 0$. The mesh is unstructured, so there is no homogeneous direction to plane-average over; a one-ring volume-weighted test filter and one-ring averaging of the Lilly numerator/denominator are used instead, which is noisier than the classic channel-flow formulation.

### Wall model options

**Neutral log-law (default):** Equilibrium log-law wall model (high-Re):

$$
u^+ = \begin{cases}
  (1/\kappa) \ln(y^+) + B & (y^+ > 11.5) \\
  y^+ & (y^+ \leq 11.5)
\end{cases}
$$

Friction velocity $u_\tau$ is found via Newton iteration; wall shear stress $\tau_w = u_\tau^2$ is applied as a Neumann BC on the momentum equations.

#### Atmospheric Modelling (Mainly for Urban Micro-climate simulations)

**Monin-Obukhov Similarity Theory (MOST):** Advanced wall model for thermal coupling and atmospheric flows:
- Coupled momentum and heat flux via Newton iteration
- Dynamic Obukhov length (stability-aware corrections)
- Businger-Dyer functions for stable/neutral/unstable regimes
- Realistic ABL surface layer parameterization
- Enabled via `&dopamine_abl` namelist

### Advection schemes
- **CENTRAL / CDS**: linear central face interpolation; lowest numerical dissipation of the currently implemented options and the preferred choice for WMLES.
- **KEP / SKEW**: momentum uses a midpoint face state with projected face mass flux (Ham-Iaccarino style collocated FV convective form) to reduce numerical dissipation and better preserve resolved kinetic energy on interior faces.
- **QUICK**: geometry-aware quadratic upwind reconstruction in projected streamline coordinates; it reduces to classical QUICK on aligned uniform meshes.
- **KOREN_TVD**: bounded Koren-TVD MUSCL reconstruction using geometry-aware projected slopes, retained as a robust nonlinear option for under-resolved transients.
- **VANLEER_TVD**: bounded van-Leer-TVD MUSCL reconstruction using geometry-aware projected slopes.

Momentum advection is selected via `advection_scheme`, temperature via `temperature_scheme`, and passive scalars via `scalar_tvd_scheme` in `&dopamine_scalars`.

### Viscous scheme
Central differencing (2nd order) with over-relaxed non-orthogonality correction (Jasak 1996).

---

## Quick Start

### 1. Build (Standard)

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DDOPAMINE_ENABLE_PETSC=ON \
  -DDOPAMINE_PETSC_INSTALL_DIR=/path/to/petsc \
  -DDOPAMINE_PETSC_ARCH_DIR=/path/to/petsc-arch
make -j
```

**Cleaning:** `make clean` only removes dopamine's own build artifacts (object
files, `libdopamine_lib.a`, the executables) — it leaves the fetched/built
PETSc and hypre install alone, since that takes minutes to rebuild. Use
`make cleanall` to also wipe the external libraries (`external/` by default,
or `DOPAMINE_EXTERNAL_DIR` if set) when you actually want a from-scratch
external rebuild — e.g. after changing `DOPAMINE_PETSC_WITH_CUDA` or another
option that changes how PETSc itself is configured.

**Build options** (`cmake -D<option>=<value>`):

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

### 2. Offline Mesh Decomposition (for large meshes)

For meshes >10M cells, pre-decompose using the `dopamine-decompose` tool:

```sh
# For large meshes, set a memory budget to enable staged face geometry I/O
mpirun -np 1 ./dopamine-decompose input.input --nranks 128 \
  --output-dir ./mesh_decomposed --max-memory-gb 48

# After decomposition, run solver with preprocessed mesh:
# In input.input:
#   mesh_file = './mesh_decomposed'
#   mesh_type = 'preprocessed'
```

Options:

| Flag | Meaning |
|---|---|
| `--nranks N` | number of partitions to write (required) |
| `--output-dir DIR` | destination for `processorN/mesh.bin` |
| `--threads N` | writer threads; defaults to all cores |
| `--max-memory-gb G` | cap peak memory by staging face geometry to disk |
| `--chunk-size M` | face scan chunk, default 2M |
| `--weight-by-volume` | balance total volume per rank instead of cell count |

**Note:** Decomposition memory scales with mesh size and topology. Use `--max-memory-gb` to constrain peak memory; the decomposer stages face geometry to disk when needed.

**Load balancing:** partitions carry an equal number of cells by default, since solver cost scales with cell count. `--weight-by-volume` restores the old behaviour of equalising total volume per rank — on a wall-graded mesh that imbalances cell counts by roughly 3.5x at 1024 ranks, so prefer the default unless you have a specific reason.

The decomposer runs on a single core (`-np 1`); `--threads` parallelises the write phase only.

### 3. Run Solver

```sh
mpirun -np 16 ./dopamine examples/channel_flow_retau395/channel_retau395.input
```

### 4. Visualize

Open `{prefix}.xmf` in ParaView for mesh and field visualization. See below for output details.

### 5. Checkpoint and Restart

**Enable checkpoint writing:**
```ini
&dopamine_checkpoint
  checkpoint_interval = 5000    ! Write every 5000 steps
/
```
Generates: `output_prefix/processorR/checkpoint_XXXXXXXX.bin` for each rank + manifest file `checkpoint_latest`

**Restart from checkpoint (same rank count):**
```ini
&dopamine_checkpoint
  restart_enabled = .true.
  restart_dir     = 'output1'
  restart_step    = -1          ! Auto-read checkpoint_latest
/
```

**Repartition checkpoints to different rank count** (e.g., 4 → 6 ranks):
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

---

## How to Compile (Offline/HPC/No Internet)

If your machine does not have internet access, download and build dependencies on another machine, then transfer them.

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
Standalone PETSc build instructions to link (Important to note that `--download-hypre=1` & `--download-fblaslapack=1` are needed for dopamine)

```
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
---

## GPU (CUDA) Build

Single-GPU support for the pressure solve (the dominant cost — roughly 83% of step
time on the CPU path) works as follows: PETSc's pressure
matrix/vectors run as `MATAIJCUSPARSE`/`VECCUDA`, and BoomerAMG (hypre) runs its
setup and solve on-device. In this build nothing outside the pressure solve runs
on GPU; the explicit kernels are offloaded separately by the OpenACC build
described below. Validated on a single-GPU workstation (RTX 3060, CUDA 12.0); multi-GPU/
multi-node has not been exercised. See
`docs/external-lib-guidance/PETSc_GPU_porting.md` (local, gitignored) for
implementation notes.

### Requirements

- An NVIDIA GPU and `nvcc` on `PATH`.
- **A CUDA-aware MPI**, if you want `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=1` (device-
  to-device transfers without host staging). Most distro-packaged MPI (Ubuntu/
  Debian `libopenmpi-dev`, etc.) is **not** built with CUDA support — check with
  `ompi_info --parsable | grep cuda`; if `mca:btl:smcuda` isn't listed, build
  your own:

  ```bash
  curl -fsSLo openmpi.tar.gz https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz
  tar xzf openmpi.tar.gz && cd openmpi-4.1.6
  ./configure --prefix="$HOME/opt/openmpi-cuda-4.1.6" --with-cuda=/usr CC=gcc CXX=g++ FC=gfortran
  make -j && make install
  ```

  (`--with-cuda=/usr` works when CUDA is the distro `nvidia-cuda-toolkit` package,
  where headers/libs sit under the normal `/usr/include` and `/usr/lib/<arch>`
  search paths rather than a unified `/usr/local/cuda` tree — check
  `find /usr -iname cuda_runtime.h` if `configure` can't find it.)

  Without a CUDA-aware MPI, run with `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0`
  (PETSc's CUDA Mat/Vec objects still work; halo/scatter traffic is staged
  through the host instead of going device-to-device).

### Configure and build

Use a separate build directory and external-libs prefix so this doesn't disturb
an existing CPU build sharing the same source tree (`DOPAMINE_EXTERNAL_DIR`
defaults to `<source>/external` regardless of build directory):

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

`DOPAMINE_PETSC_WITH_CUDA=OFF` (the default) forces `--with-cuda=0` on PETSc's
own configure — without that, PETSc auto-probes for `nvcc` and silently builds
CUDA support anyway on any machine that has the toolkit installed, which then
fails to link (no `-lcudart`/`-lstdc++` added to the executables). If you're
switching an existing tree between CUDA and non-CUDA (or changed any other
PETSc-configure-affecting option), run `make cleanall` first — a plain
`make clean` deliberately leaves the external PETSc/hypre install alone (see
"Cleaning" above), so switching `DOPAMINE_PETSC_WITH_CUDA` without a full
external rebuild would otherwise link against the previous build's PETSc.

### Run

```bash
export LD_LIBRARY_PATH="$HOME/opt/openmpi-cuda-4.1.6/lib:$LD_LIBRARY_PATH"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=1 \
  mpirun -np 1 ./build-gpu/dopamine your_case.input
```

| Env var | Default | Effect |
|---|---|---|
| `DOPAMINE_PETSC_USE_GPU` | on (if PETSc has CUDA) | `0` forces the CPU PETSc objects even in a CUDA-enabled build |
| `DOPAMINE_PETSC_USE_GPU_AWARE_MPI` | on (if PETSc has CUDA) | `0` disables `-use_gpu_aware_mpi`; needed unless your MPI is CUDA-aware (see above) |

Both toggles are independent of the CMake flags of the same name: the CMake
flags control what PETSc/hypre are *built* with; the env vars control what's
actually *used* at runtime, and default to "on" whenever the build supports it
— so a CUDA-enabled build with a non-CUDA-aware MPI needs
`DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` explicitly, or PETSc aborts at startup
with "your MPI is not GPU-aware".

If you tune BoomerAMG yourself via `&dopamine_PetscParams` (`petsc_options_tmp`
in `reference.input`), avoid overriding `-pc_hypre_boomeramg_interp_type` to a
host-only value (e.g. `FF1`, `FF`, `classical`) when running on GPU — hypre's
setup segfaults if a host-only interpolation routine is applied to a
device-resident matrix. Stick to `ext+i` (the default the GPU path already
sets) or another value from hypre's GPU-supported set.

---

## OpenACC (nvfortran) Build

Offloading the explicit kernels — advection, viscous, SGS, gradients — with
OpenACC needs `nvfortran`, since gfortran doesn't support OpenACC. `!$acc`
directives cover the SGS kernels (Vreman, Smagorinsky, dynamic Smagorinsky),
the gradient kernels (Green-Gauss and least-squares, scalar and vector), the
viscous fluxes (momentum and scalar diffusion), and the advection kernels (CDS,
QUICK/TVD reconstruction, and KEP), all driven through the mesh's flat SoA
mirror. See `docs/usage.md`'s "OpenACC (nvfortran) Run" section for the kernel
list, and `docs/external-lib-guidance/PETSc_GPU_porting.md` (local, gitignored)
for implementation notes.

> **Known quirk**: a bare `cmake -S . -B build-acc26` (or `cmake .` from
> inside the build directory) on an *already-configured* OpenACC build tree
> has been observed to silently reset the whole toolchain — compiler,
> `DOPAMINE_ENABLE_OPENACC`, `DOPAMINE_PETSC_WITH_CUDA` — back to plain
> gfortran/CPU defaults, with no error (root cause not fully understood; see
> the porting-plan doc). The build still "succeeds" afterwards, just silently
> without GPU support, so check `grep CMAKE_Fortran_COMPILER: build-acc26/CMakeCache.txt`
> resolves to `nvfortran`, not `/usr/bin/f95`, if you're ever unsure. Always
> reconfigure with the **full** `cmake -S . -B build-acc26 ...` flag list below
> (never a bare `cmake .`); if the cache still looks wrong afterwards, delete
> the build directory (`rm -rf build-acc26`, **not** `external-acc26` — that's
> the expensive PETSc/hypre/CUDA build, kept separate on purpose) and
> reconfigure from scratch.

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
  (26.5) doesn't reproduce; if you're stuck on an old NVHPC install, expect to
  need the per-file `-O1 -Mnovect` / `-Mnoopenmp` workarounds already in
  `CMakeLists.txt` to still apply, and possibly more.

- **A Fortran-compiler-matched MPI.** Fortran `.mod` files are compiler-specific
  (not just ABI — the format itself differs), so gfortran-built MPI (the
  system one, or the custom one built for the CUDA GPU path above) cannot be
  used with `nvfortran`'s `mpi_f08`. NVHPC bundles its own CUDA-aware MPI
  (HPC-X) specifically for this — use it, don't try to reuse an existing MPI
  build:

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

### Run

```bash
export LD_LIBRARY_PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/cuda/13.2/targets/x86_64-linux/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/math_libs/lib64:$LD_LIBRARY_PATH"
export OPAL_PREFIX="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0 mpirun -np 1 ./build-acc26/dopamine your_case.input
```

`LD_LIBRARY_PATH` here is for the *runtime* loader (`ld.so`) resolving the
same CUDA/MPI shared libraries at process startup — a separate concern from
the link-time `-L` search path CMake now handles automatically above. Both
are needed; neither substitutes for the other.

`OPAL_PREFIX` avoids a cosmetic-but-alarming-looking startup message:
NVIDIA's prebuilt HPC-X tarball has its own build-time install prefix
(something like `/proj/libraries/nv/26.5/.../hpcx-2.50/<build-id>/ompi5`)
baked into the binaries for locating PMIx/OMPI's help-text files, which
doesn't match wherever you actually extracted it — so at every `mpirun`
startup it prints `Sorry! ... couldn't open the help file:
.../help-pcompress.txt: No such file or directory` (harmless; PMIx just
couldn't find its own error-explanation text, which has nothing to do with
DOPAMINE). Pointing `OPAL_PREFIX` at the real HPC-X install root — `comm_libs/<cuda-ver>/hpcx/hpcx-<ver>/ompi5` — fixes the path resolution and
silences it.

`DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` matters even though the build was
configured with `-DDOPAMINE_PETSC_USE_GPU_AWARE_MPI=ON`: that CMake flag only
controls what PETSc/hypre are *built* with (`--with-cuda-aware-mpi=1`); at
*runtime* PETSc still probes whether the MPI actually in use is GPU-aware,
and aborts ("your MPI is not GPU-aware") if HPC-X's UCX/PML stack doesn't
report as such — which it may not, depending on how HPC-X negotiates
transports on your hardware. Harmless to disable for single-GPU/single-rank
runs (no inter-rank traffic to accelerate); for genuine multi-rank multi-GPU
runs, only drop this once you've confirmed
`ompi_info --parsable | grep -i cuda` and your UCX transport are actually
wired up, and expect to pass `mpirun ... --mca pml ucx` alongside it.

**Runtime-validated (2026-08-18)** on a single-GPU workstation: built clean, ran a
`channel_flow_retau395`-style case end to end with `max|divU|` at machine precision every logged
step — see `docs/usage.md`'s "OpenACC (nvfortran) Run" section for the numbers. The CUDA
13-runtime-vs-13.0-driver mismatch some machines may hit didn't materialize here, but if you do
hit a driver-insufficient-version error at startup, try building against an older CUDA toolkit's
`nvcc` for `DOPAMINE_PETSC_CUDAC` while keeping `nvfortran` 26.5 as the Fortran/OpenACC compiler.

---

## Input File Structure

DOPAMINE uses Fortran namelists. See `reference.input` for a full template. Key blocks:

- `&dopamine_domain` (mesh file/type)
- `&dopamine_numerics` (schemes, time integration, SGS model, tolerances)
- `&dopamine_physics` (viscosity, buoyancy, forcing)
- `&dopamine_bc` (per-patch boundary conditions)
- `&dopamine_initial_condition` (initial field setup)
- `&dopamine_abl` (optional atmospheric boundary layer)
- `&dopamine_output` (output prefix, intervals)
- `&dopamine_checkpoint` (optional checkpoint/restart)
- `&dopamine_stats` (optional runtime statistics cadence)
- `&dopamine_PetscParams` (optional PETSc option overrides)
- `&dopamine_scalars` (optional passive scalars)
- `&dopamine_scalar_sources` (optional sources)
- `&dopamine_sampling_planes`/`lines`/`curves`/`surfaces` (optional sampling)

Example:

```fortran
&dopamine_domain
    mesh_file = 'channel.msh'
    mesh_type = 'gmsh'
/
&dopamine_numerics
    advection_scheme    = 'KEP'
    temperature_scheme  = 'VANLEER_TVD'
    time_integrator     = 'RK3W'
    dt                  = 5.0e-4
    t_end               = 200.0
    cfl_max             = 0.5
    adapt_dt            = .true.
    p_tol               = 1.0e-8
    p_max_iter          = 500
/
```

---

## Output Files

- `{prefix}.nc` — merged unstructured mesh + time-series cell fields (`p`, `U`, `T`, scalars)
- `{prefix}.xmf` — temporal XDMF sidecar for ParaView
- `{prefix}_stats.csv` — runtime statistics (if enabled)
- `{prefix}_sample_{kind}_{name}.nc` — optional sampling outputs

Open `{prefix}.xmf` in ParaView. For decomposed MPI output, use:
1. `Merge Blocks` (work on merged grid)
2. `Clean to Grid`
3. `Cell Data to Point Data`
4. Color by velocity magnitude (`|U|`)

---

## Supported Mesh Formats

- **Gmsh v4 (.msh):** tet4, hex8, prism6, pyr5, tet10, hex20. Use `$PhysicalNames` for boundaries. Generate with `gmsh -3 -format msh4 geometry.geo`.
- **OpenFOAM polyMesh:** Reads `constant/polyMesh/{points,faces,owner,neighbour,boundary}`. All OpenFOAM-supported cell types.

---

## Example Cases

See `examples/` for ready-to-run input files:

- Channel flow (`channel_flow_retau395/`)
- Lid-driven cavity (`lid_driven_cavity/`)
- Taylor-Green vortex (`taylor_green/`)
- Thermal convection (`thermal_convection/`)

---

## Troubleshooting & Tips

- If CMake cannot find a library, check `*_DIR` variables and `CMAKE_PREFIX_PATH`.
- For runtime missing symbols, ensure `LD_LIBRARY_PATH` includes all required directories.
- Use strict `p_tol` (e.g., `1e-8`) for projection quality in LES.
- Monitor `max|divU|` and pressure residuals for stability.
- Use `solver_stats_frequency` to control statistics print interval.

---

## Documentation

The `docs/` folder is published as a [GitHub Pages site](https://akshaypatil1994.github.io/dopamine/)
(Settings → Pages → source `/docs`) and also builds standalone as a single
PDF/HTML manual via `bash docs/build.sh` (see [docs/README.md](docs/README.md)).

- [docs/theory.md](docs/theory.md): Governing equations, discretization, and numerical implementation
- [docs/usage.md](docs/usage.md): Usage guide (build, run, input, output, examples)
- [docs/development.md](docs/development.md): Code layout, build system, testing, style, contributing
- [reference.input](reference.input): Complete input template

---

## Acknowledgments

I started writing this code as a way to learn a bit more about the finite-volume discretisation from scratch. Initially, it began as an exercise in python which quickly turned to Fortran-2008 standard. A lot of other exisiting codes have inspired including but not limited to **CaNS**, **uDALES**, and **OpenFOAM** that shaped the structure and design of **Dopamine**. A large part of the code took inspiration and foundational knowledge provided from the Coastal Ocean Modelling course (Prof. Oliver Fringer [https://web.stanford.edu/~fringer/]) at Stanford during my PhD, and to the classic textbook by Ferziger, Peric, and Street, (4th Edition) which shaped much of the solver's philosophy and implementation. Also I would like to thank my post-doc advisor Clara García-Sánchez for inspiring many of the ideas implemented here and the wonderful support of the 3DGeoinformation research group.

## Contributing

If you find a bug, see a feature that is interesting, please open a pull request or start a discussion. 

## License

AGPL-3.0-only (see `LICENSE`).

