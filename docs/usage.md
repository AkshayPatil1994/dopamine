---
title: Usage
nav_order: 3
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Usage

## Build and mesh preparation

- Building the solver, dependencies, GPU builds and the offline build: [Installation](installation.md).
- Mesh formats, patch names, direct read vs. pre-decomposed meshes, and `dopamine-decompose`: [Pre-processing](preprocessing.md).

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

Generates `<prefix>/processorR/checkpoint_XXXXXXXX.bin` for each rank
plus a manifest file `checkpoint_latest`, where `<prefix>` is `prefix` from `&dopamine_output`.

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
mpirun -np 1 ./build/dopamine-decompose my_input.input --nranks 6 --output-dir mesh_6ranks

# Step 2: Repartition checkpoint from old decomposition to new
./build/dopamine-checkpoint-partition --old-dir output_4ranks --old-nranks 4 \
                              --new-dir mesh_6ranks --new-nranks 6 \
                              --step 10000

# Step 3: Create new input file pointing to 6-rank mesh and checkpoint
# (Update mesh_file, mesh_type='preprocessed', restart_dir, etc.)

# Step 4: Run 6-rank simulation from repartitioned checkpoint
mpirun -np 6 ./build/dopamine input_6ranks.input
```

## GPU (CUDA) Run

In this build only the pressure solve (PETSc + hypre BoomerAMG) runs on GPU;
the explicit kernels are offloaded by the separate OpenACC build below. See
[Installation: GPU (CUDA) build](installation.md#gpu-cuda-build) for build instructions and MPI
requirements. Once built:

```bash
export LD_LIBRARY_PATH="$HOME/opt/openmpi-cuda-4.1.6/lib:$LD_LIBRARY_PATH"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=1 \
  mpirun -np 1 ./build-gpu/dopamine examples/channel_flow_retau395/channel_retau395.input
```

Set `DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` if your MPI isn't CUDA-aware (most
distro packages aren't — see [Installation](installation.md#gpu-cuda-build)); PETSc aborts at
startup rather than silently falling back if you leave it on without one.
Compare the `pressure`
line in the performance snapshot against a `DOPAMINE_PETSC_USE_GPU=0` run on
the same binary to see the actual speedup, since that's the step this currently
accelerates.

## OpenACC (nvfortran) Run

`!$acc` directives cover the whole time step:

- the SGS models (`vreman_sgs_viscosity`, `smagorinsky_sgs`,
  `dynamic_smagorinsky_sgs` incl. its test filter) and the log-law wall model;
- the gradients (Green-Gauss and least-squares, scalar and vector), the
  viscous fluxes, and the advection kernels (CDS, QUICK/van Leer/Koren TVD, KEP,
  scalar and vector);
- the boundary-condition application (velocity, pressure, temperature,
  passive scalars), the face-flux assembly, the pressure-RHS and the
  velocity/flux projection, the RK3 arithmetic, the time-step and bulk-velocity
  reductions, and the halo pack/unpack.

All of these read geometry through the mesh's flat SoA mirror (cell centroid
and volume, face owner/neighbour, normal/area/centre/`d_cf`, interpolation and
skewness factors, least-squares stencils, patch BC ids and periodic pairing
maps), which `main.f90` stages on the device once with a single
`!$acc enter data`; the static projection operator (face indices, coefficients,
non-orthogonal caches) is staged by `solver_init`. `solver%to_device` then moves
every field and work array to the device once, before the time loop, and the
device copy stays authoritative: `solver%sync_host` refreshes the host copy
only on steps that write output, statistics, samples or checkpoints. Still on
the host, with narrow updates around them: scalar sources, nudging, SEM
inflow, and the PETSc interface when the backend cannot take device pointers.
The MOST wall model and `cell_projection = 'EXACT'` run on the device; `'EXACT'`
is single-rank on the GPU and falls back to `'LS'` with a warning on more ranks.
Its PCG stops at a loose tolerance, so GPU results differ from the CPU at about
1e-5 (use `REGRESSION_RTOL=1e-4` for the `*_exact` regression cases). A
single-GPU timing run is `tests/bench/run_bench.sh <N> <workdir>`. See [Installation: OpenACC build](installation.md#openacc-nvfortran-build)
for the build recipe (`build-acc26/`).

```bash
export PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/hpcx/bin:$HOME/nvhpc_26/Linux_x86_64/26.5/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/cuda/13.2/targets/x86_64-linux/lib:$HOME/nvhpc_26/Linux_x86_64/26.5/math_libs/lib64:$LD_LIBRARY_PATH"
export OPAL_PREFIX="$HOME/nvhpc_26/Linux_x86_64/26.5/comm_libs/13.2/hpcx/hpcx-2.50/ompi5"
DOPAMINE_PETSC_USE_GPU=1 DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0 mpirun -np 1 ./build-acc26/dopamine examples/channel_flow_retau395/channel_retau395.input
```

The OpenACC path is validated on a single-GPU workstation against a
`channel_flow_retau395`-style case (73728 cells, `-np 1`): `max|divU|` holds at
machine precision (~1e-11) every logged step with `p_resid_max` ~1.4e-9, the
same correctness signature as the CUDA-only run in the
[GPU (CUDA) build](installation.md#gpu-cuda-build). `OPAL_PREFIX` and
`DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0` are both needed on this HPC-X install
(see the paragraphs above for why). If it fails to even start with a
driver-version error, see [Installation: OpenACC build](installation.md#openacc-nvfortran-build) for
the fallback (pair `nvfortran` 26.5 with an older `nvcc` via
`DOPAMINE_PETSC_CUDAC`).


### Multi-GPU (one MPI rank per GPU)

Run one rank per GPU: `mpirun -np 2 ./build-acc/dopamine case.input` on a two-GPU
node. At startup the code binds each rank to `local_rank % ndevices` (for both
OpenACC and PETSc) and prints the mapping and whether the MPI is CUDA-aware. All
solver fields stay on the device for the whole run — the host copy is refreshed
only on steps that write output, statistics, samples or checkpoints — and halos
are packed on the GPU and sent as device buffers when the MPI is CUDA-aware
(HPC-X from NVHPC is; `DOPAMINE_GPU_AWARE_MPI=0` forces host staging of the
packed buffers).

- Two GPUs only pay off above a few million cells per GPU: the explicit part
  scales well, but at ~1M cells per GPU every kernel and AMG level is
  launch/latency bound.
- **The pressure solve does not scale with hypre's parallel BoomerAMG on this
  class of node.** On two PCIe-connected RTX A6000 (no NVLink) the parallel
  V-cycle splits every SpMV into a diagonal and an off-diagonal block and each
  launch costs as much as the whole-domain one on a single GPU, so pressure gets
  *slower* than on one GPU (TGV 192^3: 1.46 s/step on one GPU, 2.01 s on two).
  Block-Jacobi with a local AMG per GPU (additive Schwarz with overlap) avoids the
  inter-GPU traffic inside the cycle; add to the input file:

  ```fortran
  &dopamine_PetscParams
    n_petsc_options = 12
    petsc_options_tmp(1)  = '-pc_type asm'
    petsc_options_tmp(2)  = '-pc_asm_overlap 16'
    petsc_options_tmp(3)  = '-pc_asm_type basic'
    petsc_options_tmp(4)  = '-sub_ksp_type preonly'
    petsc_options_tmp(5)  = '-sub_pc_type hypre'
    petsc_options_tmp(6)  = '-sub_pc_hypre_type boomeramg'
    petsc_options_tmp(7)  = '-sub_pc_hypre_boomeramg_strong_threshold 0.25'
    petsc_options_tmp(8)  = '-sub_pc_hypre_boomeramg_coarsen_type PMIS'
    petsc_options_tmp(9)  = '-sub_pc_hypre_boomeramg_interp_type ext+i'
    petsc_options_tmp(10) = '-sub_pc_hypre_boomeramg_relax_type_all l1scaled-Jacobi'
    petsc_options_tmp(11) = '-sub_pc_hypre_boomeramg_agg_nl 0'
    petsc_options_tmp(12) = '-sub_pc_hypre_boomeramg_P_max 4'
  /
  ```

  The `-sub_pc_hypre_boomeramg_*` lines are required: they are the GPU-safe
  settings the backend applies to the unprefixed `-pc_hypre_boomeramg_*`
  options, and hypre's host-only defaults crash on a device matrix. Overlap trades
  iterations for setup cost (TGV 128^3, 2 GPUs: overlap 2/8/32 → 46/34/21 CG
  iterations per solve; step 0.73/0.67/0.63 s versus 0.82 s with the default).
  This is specific to few, well-connected GPUs; with many ranks stay with the
  default parallel BoomerAMG.

## Mesh Path and Formats

See [Pre-processing](preprocessing.md#choosing-a-mesh-path) for choosing between direct read and a pre-decomposed mesh, and for the supported mesh formats.

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
- `&dopamine_inflow`, `&dopamine_inflow_gen`, `&dopamine_inflow_record` (optional, see the inflow-database sections)

Every parameter, with its default, units and allowed values, is listed in the
[Input Reference](input-reference.md). `reference.input` is an annotated template.

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

Names are case-insensitive. An unrecognised name stops the run. A mesh patch that
is not listed in `&dopamine_bc` takes the default type with a warning, so check the startup log.

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

- `CENTRAL` or `CDS`: linear central interpolation (second order on uniform orthogonal regions). It is the default but does not conserve kinetic energy.
- `KEP` or `SKEW`: midpoint momentum convective state for lower dissipation and better resolved kinetic-energy behavior in LES-oriented runs.
- `QUICK`: geometry-aware quadratic upwind reconstruction evaluated in projected streamline coordinates. On aligned uniform meshes it reduces to classical QUICK [@leonard1979].
- `VANLEER_TVD`: bounded MUSCL reconstruction with van Leer limiter using geometry-aware projected slopes [@vanleer1974; @sweby1984].
- `KOREN_TVD`: bounded MUSCL reconstruction with Koren limiter using geometry-aware projected slopes [@koren1993; @sweby1984].

`advection_scheme` controls momentum; `temperature_scheme` controls the temperature (T) field independently. `KEP`/`SKEW` is a momentum scheme only, so `temperature_scheme` accepts `CENTRAL`/`CDS`, `QUICK`, `VANLEER_TVD` and `KOREN_TVD`. Passive scalar advection is controlled separately by `scalar_tvd_scheme` in `&dopamine_scalars` (same four names as `temperature_scheme`). An unknown name stops the run and lists the allowed ones. The default `advection_scheme` is `CENTRAL`; set `KEP` for LES.

## Recommended Numerics by Problem Type

The pressure-velocity coupling defaults (`rhie_chow_dual = .true.`, `cell_projection = 'EXACT'`) need no change for
any of the cases below. On a uniform orthogonal mesh (the Taylor-Green, cavity and channel examples)
the `'LS'` projection is already energy-stable ([Theory](theory.md)), so
`cell_projection = 'LS'` is a cheaper choice there; see
[Pressure Projection and Linear Solve](#pressure-projection-and-linear-solve).

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
  sem_sigma = 0.0                     ! Eddy length scale [m]; 0 = auto (0.1*min inlet span)
  sem_type = 'dfsem'                  ! 'dfsem' (divergence-free eddies, default) or 'jarrin' (original tent kernel)
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
- `sem_type = 'dfsem'` (default) uses the divergence-free eddies of Poletto et al. (2013): ellipsoidal eddies in the principal frame of the Reynolds tensor, so the full tensor including shear is reproduced; `sem_type = 'jarrin'` keeps the original tent-kernel SEM
- Eddies are convected downstream and create anisotropic turbulent fluctuations
- Target Reynolds stress anisotropy matches ABL log-layer profiles
- `sem_sigma` is the longest eddy half-extent (default $$0.1 \times$$ the smaller in-plane extent of the inlet); with `dfsem` the other axes are smaller ($$\sigma\sqrt{\lambda_i/\lambda_{ref}}$$, floored at $$0.2\sigma$$), so resolve the *smallest* half-extent with at least 8 cells

**Usage notes:**
- Apply SEM only to inlet patches; other BCs are unaffected
- SEM injects velocity fluctuations only; it does not require `boussinesq = .true.`
- Recycle box depth is estimated from mesh bounds automatically
- Eddies leaving the box are re-seeded upstream; across a periodic in-plane direction they wrap, so the fluctuation is periodic there (`sem_sigma` must then be below half the periodic width; the run stops otherwise)
- `jarrin` is not divergence-free and gives a large inlet pressure spike; `dfsem` does not, and the inlet is evaluated at every Runge-Kutta stage time so the first cell sees no stale normal flux. The remaining first-cell pressure is the field's own nonlinear pressure plus a discretisation term that falls like $$(\Delta x/\sigma)^2$$: resolve the smallest eddy half-extent with at least 8 cells (a startup message reports it and warns below 5). See *Measured behaviour* in the theory notes and `tests/sem_pressure/`
- SEM is evaluated on the host each step; the face fluctuation is then uploaded, so GPU runs and any rank count use the same replicated eddy state (checkpoints store eddy positions and signs; frames and amplitudes are re-derived)

#### Inputs: mean profile and turbulence intensities

Three ways to specify the inflow, in order of precedence.

**1. A table** (`sem_profile_file`, `sem_profile_format`). Rows need not be sorted; values are
interpolated linearly in `z` and clamped outside the table, so extend it past the domain top.

| `sem_profile_format` | columns | frame |
|---|---|---|
| `'ti'` | `z U Iu Iv Iw [rho_uw]` | wind-aligned |
| `'stress'` (default), 6 col | `z U uu vv ww uw` | wind-aligned |
| `'stress'`, 10 col | `z U V W uu vv ww uv uw vw` | global axes |

The `ti` layout is the wind-tunnel one:
$$I_u=\sigma_u/U$$ etc. are turbulence intensities relative to the local mean speed, and
`rho_uw` is the correlation coefficient $$\overline{u'w'}/(\sigma_u\sigma_w)$$, typically
$$-0.1$$ to $$-0.35$$ in the surface layer.

```
# z     U      Iu    Iv     Iw    rho_uw
 0.5   4.0   0.30  0.22   0.15   -0.30
10.0   7.1   0.19  0.14   0.095  -0.30
60.0  12.0   0.10  0.075  0.05   -0.30
```

A missing `rho_uw` uses `sem_rho_uw`. With a table the wind direction comes from `wind_dir_deg`
(`ti` and 6-column `stress`) or from the file itself (10-column `stress`, where the axes are global).

**2. Namelist laws** (no file; set `sem_Iu_ref > 0`):

```fortran
&dopamine_abl
  use_sem = .true.
  u_ref = 10.0    z_ref = 100.0    z0_m = 0.03
  sem_mean_law = 'power'       ! or 'log' (default)
  sem_alpha    = 0.2           ! U = u_ref (z/z_ref)^alpha
  sem_Iu_ref   = 0.15          ! Iu at z_ref; decays like ln(z_ref/z0)/ln(z/z0) below it
  sem_Iv_ratio = 0.75          ! Iv/Iu
  sem_Iw_ratio = 0.5           ! Iw/Iu
  sem_rho_uw   = -0.3          ! u'w' correlation coefficient
/
```

Intensity is evaluated no lower than $$10 z_0$$, where the log law stops being meaningful.

**3. Panofsky–Dutton defaults** (neither of the above): neutral log law and
$$\sigma_{u,v,w}/u_* = 2.4,\,1.9,\,1.25$$ with $$\overline{u'w'}=-u_*^2$$. The momentum flux was zero in
earlier versions, so the shear stress had to regenerate over fetch; it is now present at the inlet.

**Eddy size.** With `dfsem` the size comes from one parameter, `sem_sigma` (the longest eddy half-extent;
0 = 0.1 * smaller inlet span); the other axes follow the principal stresses of the target tensor at each eddy's
height, so there is no separate length-scale input and no height-dependent scale.

`sem_seed` changes the random sequence (0 = built-in seed) for ensemble runs; every rank
uses the same sequence, so the field is identical for any rank count.

**Choosing `n_eddies`.** Fluctuation statistics are unbiased for any count, but the field
becomes intermittent when too few eddies overlap. The eddies are sorted into spatial bins (width at least the eddy support, rebuilt whenever they move), so a face or plane
point only visits the eddies of its own and the neighbouring bins: the cost is $$\mathcal{O}(N_e + N_{faces}\,n_{local})$$ per
evaluation, with $$n_{local}$$ the number of eddies overlapping a point, independent of `n_eddies`. This holds for the online SEM and `dopamine-inflow-gen`.

**Mean profile.** A table or `sem_mean_law = 'power'` **replaces** the analytic log law at
`ablInlet` faces, so `u_ref`, `z_ref` and `z0_m` then no longer set the inlet mean. They
still set the MOST wall model, which is a separate consumer.

**What is and is not exact.** The mean profile is imposed exactly. The fluctuation is
divergence-free as an analytic field in 3D; the stresses are matched in the statistical
sense and approach the target as sampling lengthens.
The per-step uniform flux correction at the inlet remains, and the pressure projection still
adjusts the field in the first cells downstream of the inlet plane. The heights are taken as
`z` of the face centre, so terrain-following meshes need ground-relative heights in the table.

The binned evaluation is checked against the all-eddy sum for both kernels (periodic, open and two-bin inlets, look-ahead). The kernel is checked by `tests/unit/test_sem_dfsem.f90`
(`cmake --build build --target test_sem_dfsem && ctest --test-dir build -L unit`): eigen-solver, single-eddy
divergence, quadrature and Monte Carlo second moments.

### Inflow-plane databases (precursor and successor)

Instead of running the eddies online, the inflow can be **pre-computed** into a file of
spanwise-periodic fluctuation planes and replayed by the solver. The same file format is written
by two producers and read by one consumer:

```
 dopamine-inflow-gen (SEM) ──┐
                             ├─►  planes.nc  ──►  successor run: inflow_source = 'planes'
 precursor run (recorder) ───┘
```

The database holds the **fluctuation only** (mean removed). The successor takes its mean from the
usual `&dopamine_abl` mean-profile options (`sem_profile_file`, `sem_mean_law`).
The inlet must be periodic in the spanwise direction; the run stops otherwise.

File conventions (NetCDF-4, several independent segments):

- The plane is spanwise periodic and uniform: $$y_j = y_0 + jW/n_y$$, $$j = 0 \ldots n_y-1$$.
- Heights $$z$$ ascend from the ground and may be graded.
- Values are the fluctuation in the profile frame; component 1 is streamwise.
- Frame $$i$$ of a segment is at $$t = (i-1)\,\Delta t_p$$, so a segment spans $$(n_t-1)\Delta t_p$$.
- The plane position along the normal is implicit: the database is the inlet.

#### 1. Generate planes from the SEM

`dopamine-inflow-gen input` reads the SEM inflow description from `&dopamine_abl` (exactly as for
the online SEM) plus `&dopamine_inflow_gen`:

```fortran
&dopamine_inflow_gen
  gen_file        = 'planes.nc'
  gen_grid_file   = 'inlet_grid.txt'  ! plane grid exported by the solver (see below); or give
  ! gen_W = 1.0  gen_y0 = 0.0  gen_ny = 64  gen_nz = 48  gen_z_min = 0.0  gen_z_max = 200.0
  gen_n_turnovers = 2.0               ! length of each segment in eddy turnover times
  gen_n_segments  = 4                 ! independent realisations, replayed in random order
  gen_n_sigma     = 3.0               ! eddies per sigma^3 of box volume (sets the eddy count)
  gen_dt          = 0.0               ! frame spacing; 0 = 6 frames per shortest eddy passage
  gen_seed        = 1
/
```

- The eddy turnover time is $$T_e = \sigma/\sigma_u(z_{ref})$$ (`sem_sigma` over the rms streamwise velocity), i.e. it follows from `u_ref`
  and the turbulence intensity. `gen_T_turnover` overrides it.
- The SEM runs in a box exactly one period wide with minimum-image distances, so the planes are
  periodic across the span. The period must exceed twice the largest eddy half-extent (the run
  stops with a message otherwise).
- The generator prints the number of nodes that span the smallest eddy and verifies the written
  planes against the targets (rms of each component and $$\overline{u'w'}$$ within 15%).

**Generating on the successor's own lattice.** Run the solver once with
`&dopamine_inflow inflow_export_grid = 'inlet_grid.txt'`. It writes the inlet's spanwise nodes and
(graded) heights and stops; pass the file as `gen_grid_file`. The inlet must be uniformly spaced
across the period. With this, the stored nodes *are* the inlet faces, so the solver copies values
rather than interpolating them: no spatial filtering at all.

#### 2. Record planes from a precursor

```fortran
&dopamine_inflow_record
  record_enabled = .true.
  record_file    = 'planes.nc'
  record_point   = 100.0 0.0 0.0     ! a point on the plane
  record_normal  = 1.0 0.0 0.0       ! flow direction through the plane (horizontal)
  record_dt      = 0.02              ! frame spacing
  record_t_start = 600.0             ! discard the spin-up before this time
  record_segment_time = 60.0         ! length of one segment
  record_n_segments   = 4
/
```

The recorder samples the **layer of cell centres** on that plane (the solver's own cell values,
nothing interpolated). The layer must form a regular lattice, uniform across a spanwise-periodic
domain, otherwise the run stops. Frames are placed at `record_t_start + n*record_dt` by
interpolating linearly between solver steps (adaptive time stepping is fine; keep the solver step
below `record_dt`). Once the segments are complete the mean is removed, the measured stresses are
stored, and the replay line is removed (below). The precursor must be statistically stationary
over the recording: use `record_t_start` to skip the spin-up. Its mesh should be at least as fine
as the successor's inlet; a finer successor receives smoothed data and re-grows the small scales
over a development fetch.

The database stores the flow direction; set `wind_dir_deg` in the successor to the value the
recorder prints (270 for flow along +x).

#### 3. Replay in the successor

```fortran
&dopamine_inflow
  inflow_source  = 'planes'
  inflow_file    = 'planes.nc'
  inflow_seed    = 1               ! seeds the segment order and spanwise shifts
  inflow_t_offset = 0.0            ! start offset into the replay sequence [s]
  inflow_crossfade_time = 0.0      ! 0 = 10% of a segment
/
```

Replay is a pure function of time and seed, so a restart needs no extra state and gives the same
inflow. Period $$P = T_{seg} - T_x$$ uses a randomly chosen segment (never the same twice in a
row) shifted spanwise by a whole number of nodes (golden-ratio sequence); consecutive periods
cross-fade over $$T_x$$ with weights $$\sin\theta,\cos\theta$$, which keeps the variance for
decorrelated segments and does not disturb the in-plane divergence (the weights are uniform in space).

**Why independent segments are needed.** A spanwise shift moves every structure sideways but cannot
change the spanwise-mean ($$k_y=0$$) part of the plane. Replaying $$K$$ segments in random order
still leaves their *average* waveform repeating every period, a spectral line at $$1/P$$ holding
about $$1/K$$ of the spanwise-mean variance. Both producers therefore subtract, at every frame
index, the mean over segments of the spanwise-mean component and rescale it by
$$\sqrt{K/(K-1)}$$ (`gen_decorrelate = .true.` by default). The unit test measures the line
before (about 160 times the neighbouring spectral bins) and after (none). Use at least 4 segments.

**Interpolation to the faces:**

| direction | method | notes |
|---|---|---|
| span | exact Fourier (trigonometric) | spans are periodic and uniform; a shift is whole nodes, so it costs nothing |
| height | 6-point Lagrange | works on graded levels; a matching level is copied |
| time | Catmull-Rom cubic over four frames | needs frames well resolved: see the printed frame-to-frame correlation |

At start-up rank 0 reports how many faces coincide with database nodes (copied, not filtered),
the frame-to-frame correlation, and, when faces are interpolated, a vertical resolution check
that warns if the stored heights under-resolve the data. A database whose heights do not cover
the inlet, whose period differs from the domain width, whose flow direction differs from the inlet
normal, or whose wind direction differs from `wind_dir_deg`, stops the run.

The chain (export grid, generate, replay on 1 and 2 ranks, stop checks, restart, recorder against
the solver's own output) is exercised by `ctest -L inflow` (`tests/inflow/chain.sh`).

### Davies Nudging (Large-Scale Coupling)

When `use_nudging = .true.`:

1. Load a weather profile file: CSV or netCDF format
   - **CSV:** Columns are `[time(s), z(m), u(m/s), v(m/s), w(m/s), T(K)]`, one row per point
   - **netCDF:** Variables `times(:)`, `z_levels(:)`, `U_profile(3,z,t)`, `T_profile(z,t)`

2. At each RHS computation, the solver relaxes toward these profiles:
   - **Below** `nudging_z_start`: free LES evolution
   - **Above** `nudging_z_start`: Newtonian relaxation with timescale `nudging_tau`; the strength ramps linearly up to `nudging_z_full` (`nudging_z_full` at or below `nudging_z_start` gives a hard switch)

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

- Keep `p_tol` strict enough for long-time runs (`1e-8` works in many LES cases; the built-in default is `1e-10`).
- Increase `p_max_iter` if pressure residuals stagnate before tolerance.
- In parallel runs, monitor both pressure residual trend and `max|divU|`.
- `rhie_chow_dual = .true.` and `cell_projection = 'EXACT'` are the defaults. On meshes with non-orthogonality above
  about 40 degrees (snappyHexMesh or skewed Gmsh meshes), `rhie_chow_dual = .false.` makes the velocity grow by a
  constant factor every step, independent of `dt` (the Rhie-Chow inconsistency); see Theory.
- `cell_projection = 'EXACT'` makes the pressure step energy non-increasing on any mesh and costs about 1.4 to three
  times the step time of `'LS'`. OpenACC builds warn and use `'LS'`. Set `cell_projection = 'LS'` to trade stability
  for speed on well-behaved meshes.
- `QUICK` is not energy neutral: on snappyHexMesh meshes it can still grow slowly next to sliver prism cells at the
  ground (max|U| 13 to 44 over 150 steps in the MichelStadt case). Prefer `KEP` for LES.
- With the default exit test the non-orthogonal corrector loop never runs more than two passes whatever
  `n_nonortho_corr` is, because it compares `max|rhs_p|` (dominated by div u*) between passes. Set
  `nonortho_check = .true.` to stop on the increment of the pressure instead.
- On meshes with faces near 90 degrees of non-orthogonality (snappyHexMesh surfaces), set
  `pressure_min_cos_theta` (default `0` = off, e.g. `0.1`) in `&dopamine_numerics` to floor
  $$\cos\theta$$ in the Poisson coefficient and deferred cross term; the number of limited
  faces is logged at startup. It is not a default because `0.1` made the 45 degree skewed Taylor-Green
  box diverge, so try it only when the pressure solve struggles on such faces.
- With a CUDA-enabled PETSc this linear solve can run on GPU — see [GPU (CUDA) Run](#gpu-cuda-run). Custom `-pc_hypre_boomeramg_interp_type`/`coarsen_type` overrides
  via `&dopamine_PetscParams` should stick to hypre's GPU-supported values
  (`ext+i`/`PMIS`) when GPU mode is active; host-only choices segfault against a
  device-resident matrix.

## Long-Time Stability Checklist

1. Check `n_ranks` in the log header to confirm the intended MPI launch.
2. For wall-model cases, monitor:
   - `max|U|`
   - `E`
   - `max|divU|`
   If only `max|U|` grows while global metrics remain bounded, investigate localized outliers.
3. Keep `p_tol` and `p_max_iter` strict enough for projection quality.
4. Use `solver_stats_frequency` consistently across cases.

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

- Build and library-path problems: [Installation: Troubleshooting](installation.md#troubleshooting).
- Velocity grows by a constant factor every step on a skewed or snappyHexMesh mesh, independent of `dt`: check that `rhie_chow_dual` and `cell_projection` still have their defaults, then see [Pressure Projection and Linear Solve](#pressure-projection-and-linear-solve).
- Step time rose after an upgrade: `cell_projection = 'EXACT'` is now the default; set `'LS'` on orthogonal meshes.
- A run stops at startup with `input namelist parse error`: a parameter name is misspelled or in the wrong group; names are listed in the [Input Reference](input-reference.md).
- Monitor `max|divU|` and the pressure residual; use `solver_stats_frequency` to set how often they print.
