---
title: Input Reference
nav_order: 4
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Input Reference

Every run is configured by one Fortran namelist file passed on the command line:

```sh
mpirun -np 4 ./dopamine my_case.input
```

[`reference.input`](https://github.com/AkshayPatil1994/dopamine/blob/main/reference.input) in the
repository root is an annotated template containing every group. This page is the authoritative list of
**defaults**: the values below are the ones the parser in `src/io/mod_input_parser.f90` applies when a
parameter is omitted. Values written in `reference.input` are example values and may differ.

## Conventions

- Every group is optional. An omitted group, or an omitted parameter inside a group, takes its default.
- Group order in the file does not matter, but each group is written once, as `&name ... /`.
- Comments start with `!`. Strings use single quotes. Logicals are `.true.` / `.false.`.
- Choice strings (`sgs_model`, `advection_scheme`, ...) are case-insensitive. An unknown value stops the
  run with the list of allowed names; an unknown *parameter name* stops the run with a namelist parse error.
- Array parameters are indexed per patch/scalar/object, e.g. `patch_type(2) = 'outlet'`, `patch_U(1,:) = 1.0, 0.0, 0.0`.
- Units are SI unless stated. Lengths are in mesh units.
- Table columns: **Default** is the parser default; `required` means the run is meaningless without it.
- The pressure-velocity coupling defaults (`rhie_chow_dual = .true.`, `cell_projection = 'EXACT'`) favour stability on skewed and non-orthogonal meshes. Inputs written for older versions that need their old results set `rhie_chow_dual = .false.` and `cell_projection = 'LS'`.

## Contents

[&dopamine_domain](#dopamine_domain) ·
[&dopamine_numerics](#dopamine_numerics) ·
[&dopamine_PetscParams](#dopamine_petscparams) ·
[&dopamine_stats](#dopamine_stats) ·
[&dopamine_physics](#dopamine_physics) ·
[&dopamine_abl](#dopamine_abl) ·
[&dopamine_inflow](#dopamine_inflow) ·
[&dopamine_inflow_gen](#dopamine_inflow_gen) ·
[&dopamine_inflow_record](#dopamine_inflow_record) ·
[&dopamine_bc](#dopamine_bc) ·
[&dopamine_initial_condition](#dopamine_initial_condition) ·
[&dopamine_output](#dopamine_output) ·
[&dopamine_checkpoint](#dopamine_checkpoint) ·
[&dopamine_scalars](#dopamine_scalars) ·
[&dopamine_scalar_sources](#dopamine_scalar_sources) ·
[&dopamine_sampling_*](#sampling-groups) ·
[Recipes](#recipes)

---

## &dopamine_domain

Mesh location and reader. See [Pre-processing: Choosing a mesh path](preprocessing.md#choosing-a-mesh-path).

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `mesh_file` | `required` | `gmsh`: path to the `.msh` file. `openfoam`: path to the `polyMesh` directory. `preprocessed`: directory written by `dopamine-decompose`. |
| `mesh_type` | `'gmsh'` | `'gmsh'`, `'openfoam'` (alias `'foam'`), or `'preprocessed'`. `preprocessed` must be run with exactly the rank count it was decomposed for. An unrecognised value prints a warning and aborts at mesh read. |

## &dopamine_numerics

Discretisation, time stepping, pressure solve and SGS model. Theory: [theory.md](theory.md).

### Advection

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `advection_scheme` | `'CENTRAL'` | Momentum advection. `'CENTRAL'`/`'CDS'` lowest dissipation; `'KEP'`/`'SKEW'` kinetic-energy preserving (preferred for LES); `'QUICK'` quadratic upwind, unlimited; `'KOREN_TVD'`, `'VANLEER_TVD'` bounded MUSCL. |
| `temperature_scheme` | `'VANLEER_TVD'` | Temperature advection: `'KOREN_TVD'`, `'VANLEER_TVD'`, `'QUICK'`, `'CDS'`/`'CENTRAL'`. Only active when `boussinesq = .true.`. |
| `use_skew_symmetric_advection` | `.false.` | Deprecated and ignored; set `advection_scheme = 'KEP'` instead. |
| `gradient_scheme` | `'LEAST_SQUARES'` | `'LEAST_SQUARES'` (weighted, accurate on skewed meshes) or `'GREEN_GAUSS'` (cheaper). |
| `gg_iterations` | `2` | Green-Gauss sweeps with skewness correction; only used when `gradient_scheme = 'GREEN_GAUSS'`. `1` = no correction. |

### Time stepping

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `time_integrator` | `'RK3W'` | Wray low-storage RK3 (three stages, three pressure solves). The only integrator. |
| `time_scheme` | `'RK3W'` | Legacy alias of `time_integrator`; do not use in new inputs. |
| `dt` | `1.0e-3` s | Time step. With `adapt_dt = .true.` it is a **ceiling**: the step is `min(dt_cfl, dt)`, so too small a value defeats adaptive stepping. |
| `adapt_dt` | `.false.` | Adaptive step from the CFL and diffusion limits. |
| `cfl_max` | `0.4` | CFL target for adaptive stepping. `0.5` is a typical choice for RK3W. |
| `t_end` | `100.0` s | End time. |
| `nstep_max` | `0` | Maximum number of steps; `0` = unlimited (run to `t_end`). |

### Pressure projection

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `pressure_solver` | `'boomeramg'` | `'boomeramg'` (hypre AMG; needs PETSc built with hypre) or `'gamg'` (PETSc native AMG). Tune further with [`&dopamine_PetscParams`](#dopamine_petscparams). |
| `p_tol` | `1.0e-10` | Poisson relative residual tolerance. |
| `p_max_iter` | `1000` | Maximum Krylov iterations per solve. |
| `p_target_divU` | `1.0e-8` | Parsed but currently has no effect. |
| `use_rotational_projection` | `.true.` | Rotational pressure correction (Brown, Cortez & Minion 2001); lowers the splitting error from O(Δt) to O(Δt²). Recommended for wall-bounded flow. |
| `n_nonortho_corr` | `0` | Non-orthogonal correction passes. `0` = auto (1 on an orthogonal mesh, otherwise 2). |
| `pressure_min_cos_theta` | `0` | Floor on cos(theta) in the pressure Poisson coefficient and cross term. `0` = off, which leaves results unchanged. Meant for meshes with faces near 90 degrees of non-orthogonality (snappyHexMesh surfaces), e.g. `0.1`; the number of limited faces is logged at startup. Not a default because `0.1` made the 45 degree skewed test box diverge. |
| `rhie_chow_dual` | `.true.` | `.true.`: Rhie-Chow restores exactly the cross-free gradient removed by the pre-correction (stable on skewed and non-orthogonal meshes; identical on orthogonal faces). `.false.` restores the pre-default behaviour. |
| `nonortho_check` | `.false.` | `.true.`: the corrector loop stops on the relative increment of the pressure (`1e-3`) or when it stops contracting, and up to `n_nonortho_corr` (max 8) passes can run. Default: stops on `max|rhs_p|`, which includes div u* and ends the loop after two passes. |
| `cell_projection` | `'EXACT'` | `'EXACT'`: the cell velocity is the V-orthogonal projection of the predictor onto fields whose interpolated face flux is divergence free (a second, wide-stencil system solved by preconditioned CG with the pressure solver as preconditioner; about 1.4 to 3 times the step time, also on orthogonal meshes where it brings little, so set `'LS'` there; host builds only: OpenACC builds warn and use `'LS'`). `'LS'`: the cell velocity is corrected by the least-squares reconstruction of the face-normal pressure differences. See Theory. |
| `cell_projection_tol` | `0.1` | `EXACT`: the CG stops when the last velocity-energy step is below this fraction of the accumulated one. |

### Subgrid-scale model and wall-model filtering

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `use_sgs` | `.true.` | `.false.` runs without an SGS model (DNS). |
| `sgs_model` | `'VREMAN'` | `'VREMAN'`, `'SMAGORINSKY'` (Cs fixed at 0.1), `'DYNAMIC_SMAG'` (Germano-Lilly, local), `'AMD'` (anisotropic minimum dissipation). |
| `c_amd` | `1/3` | AMD model constant. Only used by `sgs_model = 'AMD'`. |
| `amd_wall_gradient` | `.true.` | AMD only: use the log-law wall-normal gradient in wall-model cells instead of the resolved one. |
| `wall_model_height` | `0.0` | Wall-model sampling height [mesh units]. `<= 0` samples the wall-adjacent cell. Set above the buffer layer (y+ > 50) to avoid log-layer mismatch. |
| `wall_model_filter_time` | `0.0` s | Time constant for filtering the sampled wall velocity. `<= 0` disables. Try about δ/u_τ. |
| `wall_model_filter_tf_fac` | `0.0` | Adaptive filter `T_f = f·y/(κ·u_τ)`. `> 0` enables and takes precedence over `wall_model_filter_time`. |

### Diagnostics

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `solver_stats_frequency` | `0` | Also accepted here; see [`&dopamine_stats`](#dopamine_stats) (that group wins if both are set). |
| `calc_stats_frequency` | `0` | Same as above. |
| `debug_diagnostics` | `.false.` | Verbose runtime diagnostics: projection-idempotency check, interface KE budget, worst-divergence cell report. |

## &dopamine_PetscParams

Raw PETSc options injected into the PETSc options database before `KSPSetFromOptions`.

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_petsc_options` | `0` | Number of entries used in `petsc_options_tmp` (maximum 32). |
| `petsc_options_tmp(i)` | none | One PETSc option string per entry, exactly as on the command line, e.g. `'-pc_hypre_boomeramg_strong_threshold 0.6'`. Multi-GPU settings: [Usage: Multi-GPU](usage.md#multi-gpu-one-mpi-rank-per-gpu). |

The pressure solve is the largest share of the step on CPU (about 60-65% for a periodic 64^3 box, 1-8 MPI ranks). On that case `'-pc_hypre_boomeramg_agg_nl 1'` with `'-pc_hypre_boomeramg_strong_threshold 0.25'` cut the pressure time by about 26-40% at the same `p_tol` (more CG iterations, cheaper V-cycles). These are not defaults: with `cell_projection = 'EXACT'` they change results at the 1e-5 level, because the cell projection iterates against a loose pressure solve. Keep the smoother symmetric (the default) with CG; `l1-gauss-seidel` raised the iteration count about sixfold.

## &dopamine_stats

Cadence of the on-screen/CSV monitor and default cadence of the sampling objects. A value of `0` disables.

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `solver_stats_frequency` | `0` | Monitor row every N steps (screen and CSV). |
| `calc_stats_frequency` | `0` | Fallback used when `solver_stats_frequency` is `0`. |
| `slice_write_frequency` | `0` | Default write cadence of [sampling planes](#sampling-groups). |
| `probe_write_frequency` | `0` | Default cadence of sampling surfaces, and the fallback for lines and curves when their own cadence is `0`. |
| `line_write_frequency` | `0` | Default cadence of sampling lines. |
| `curve_write_frequency` | `0` | Default cadence of sampling curves. |
| `intersection_decimate_spacing` | `0.0` | Point spacing used when thinning surface-plane intersection curves; `0` keeps all points. |
| `sort_intersection_points` | `.true.` | Order intersection points along the curve. |

## &dopamine_physics

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `nu` | `1.0e-5` m²/s | Kinematic viscosity. |
| `alpha_T` | `1.4e-5` m²/s | Thermal diffusivity. |
| `beta` | `3.4e-3` 1/K | Thermal expansion coefficient (Boussinesq). |
| `T_ref` | `300.0` K | Reference temperature. |
| `Pr_t` | `0.85` | Turbulent Prandtl number for the SGS heat flux. |
| `gx, gy, gz` | `0.0` m/s² | Gravity vector, used by the Boussinesq term. Typical vertical case: `gz = -9.80665`. |
| `fx, fy, fz` | `0.0` m/s² | Constant body force (pressure-gradient drive). Must be zero when `constant_mass_flux_forcing = .true.`. |
| `constant_mass_flux_forcing` | `.false.` | Recompute a uniform body force each step so the volume-averaged velocity tracks `Ux/Uy/Uz_constant`. |
| `Ux_constant, Uy_constant, Uz_constant` | `0.0` m/s | Target bulk velocity for mass-flux forcing. |
| `use_rotation` | `.false.` | Coriolis force `-2 Ω × U`. |
| `omega_x, omega_y, omega_z` | `0.0` 1/s | Rotation vector Ω. |
| `use_harmonic_forcing` | `.false.` | Add `A·cos(ω t + φ)·dir` to the momentum equation. |
| `forcing_amp`, `forcing_omega`, `forcing_phase` | `0.0` | Amplitude A [m/s²], angular frequency ω [rad/s], phase φ [rad]. |
| `forcing_dir_x, forcing_dir_y, forcing_dir_z` | `1.0, 0.0, 0.0` | Direction of the harmonic forcing. |
| `boussinesq` | `.false.` | Buoyancy coupling. Also gates temperature transport: with `.false.` the run is isothermal (T stays at `T_ref`) even in an ABL run. |

## &dopamine_abl

Wall treatment, synthetic-eddy inflow and large-scale nudging. The three are independent.
Workflow: [Usage: ABL module](usage.md#atmospheric-boundary-layer-abl-module).

### Wall model and surface

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `wall_model` | `''` (`'loglaw'`) | Treatment on `wallModel` patches: `'loglaw'` smooth-wall (first cell in the log layer, y+ of about 30 or more), `'most'` Monin-Obukhov rough/stratified (stops if the first-cell height ≤ `z0_m`), `'none'` no wall stress. When empty, `use_abl = .true.` selects `'most'`, otherwise `'loglaw'`. An unknown value prints a warning and uses `'loglaw'`. |
| `use_abl` | `.false.` | Legacy alias for `wall_model = 'most'` when `wall_model` is unset. Does not gate SEM or nudging. |
| `u_ref` | `8.0` m/s | Reference wind speed at `z_ref`. |
| `z_ref` | `100.0` m | Reference height of `u_ref`. |
| `z0_m` | `0.1` m | Momentum roughness length; per-patch `patch_z0_m` overrides. |
| `z0_h` | `0.01` m | Thermal roughness length. |
| `T_surface` | `300.0` K | Surface temperature. |
| `T_ref` | `300.0` K | Background temperature for stability corrections. |
| `L_obukhov` | `0.0` m | Obukhov length seeding the first wall-model call (0 = neutral, > 0 stable, < 0 unstable). |
| `compute_L` | `.true.` | Recompute L per boundary face every step from the wall fluxes. |
| `wind_dir_deg` | `270.0` ° | Meteorological direction (270 = westerly). |

### Synthetic Eddy Method (inflow at `ablInlet` patches)

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `use_sem` | `.true.` | Generate inflow turbulence at `ablInlet` patches. |
| `sem_type` | `'dfsem'` | `'dfsem'` divergence-free eddies with full Reynolds tensor; `'jarrin'` original tent-kernel SEM (not divergence-free). |
| `n_eddies` | `500` | Eddies in the recycling box. Too few makes the field intermittent. |
| `sem_sigma` | `0.0` m | Eddy length scale. `0` = `0.1 × min inlet-plane span`. |
| `sem_seed` | `0` | Random seed; `0` uses a built-in seed. Identical on every rank. |
| `sem_mean_law` | `'log'` | Mean profile when no table is given: `'log'` (from `u_ref`, `z_ref`, `z0_m`) or `'power'`. |
| `sem_alpha` | `0.2` | Exponent for `sem_mean_law = 'power'`. |
| `sem_Iu_ref` | `0.0` | Streamwise intensity at `z_ref`. `> 0` selects the intensity laws; `0` uses Panofsky-Dutton neutral stresses. |
| `sem_Iv_ratio`, `sem_Iw_ratio` | `0.75`, `0.5` | `Iv = ratio·Iu`, `Iw = ratio·Iu`. |
| `sem_rho_uw` | `-0.3` | Correlation `u'w' = ρ σ_u σ_w`. |
| `sem_profile_file` | `''` | CSV table of heights and statistics. Empty = use the laws above. |
| `sem_profile_format` | `'stress'` | `'ti'`: `z U Iu Iv Iw [rho_uw]`. `'stress'`: 6 columns `z U uu vv ww uw` (wind-aligned) or 10 columns `z U V W uu vv ww uv uw vw` (global axes, `wind_dir_deg` not applied). Heights need not be sorted; values are clamped outside the table. |

### Davies nudging

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `use_nudging` | `.false.` | Relax toward time-varying reference profiles. |
| `nudging_file` | `''` | CSV (or netCDF) with rows `time, z, u, v, w, T`: all levels of one time, then the next. |
| `nudging_tau` | `1800.0` s | Relaxation time scale; clamped internally for stability at the current `dt`. |
| `nudging_z_start` | `0.0` m | Nudging acts above this height. |
| `nudging_z_full` | `0.0` m | Height of full strength; the coefficient ramps from `z_start`. At or below `z_start` gives a hard switch. |
| `nudging_horizontal_only` | `.true.` | Leave the vertical velocity un-nudged. |

## &dopamine_inflow

Source of inlet fluctuations at `ablInlet` patches. Details: [Usage: Inflow-plane databases](usage.md#inflow-plane-databases-precursor-and-successor).

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `inflow_source` | `'sem'` | `'sem'` online eddies (see `&dopamine_abl`) or `'planes'` replay of a pre-computed database (needs a spanwise-periodic domain). |
| `inflow_file` | `''` | Database file for `'planes'`. |
| `inflow_seed` | `1` | Seeds segment order and spanwise shifts; restarts reproduce the replay. |
| `inflow_t_offset` | `0.0` s | Start offset into the replay sequence. |
| `inflow_crossfade_time` | `0.0` s | Cross-fade between segments; `0` = 10 % of a segment. |
| `inflow_export_grid` | `''` | Write the inlet plane grid for `dopamine-inflow-gen` to this file and stop. |

## &dopamine_inflow_gen

Read only by `dopamine-inflow-gen`; the eddy description comes from `&dopamine_abl`.

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `gen_file` | `'inflow_planes.nc'` | Output database. |
| `gen_grid_file` | `''` | Plane grid exported by the solver (`inflow_export_grid`). If empty, define the grid below. |
| `gen_W` | `0.0` m | Spanwise period. |
| `gen_y0` | `0.0` m | Spanwise coordinate of the first node. |
| `gen_ny`, `gen_nz` | `0` | Spanwise nodes (uniform over the period) and heights. |
| `gen_z_min`, `gen_z_max` | `0.0` m | Height range of the plane. |
| `gen_dt` | `0.0` s | Frame spacing; `0` = six frames across the shortest streamwise eddy passage. |
| `gen_n_turnovers` | `2.0` | Segment length in eddy turnover times. |
| `gen_T_turnover` | `0.0` s | Override for the turnover time; `0` = `Lx(z_ref)/σ_u(z_ref)`. |
| `gen_n_segments` | `4` | Independent realisations replayed in random order; use at least 4. |
| `gen_n_sigma` | `3.0` | Eddies per kernel volume (≥ 1; 3 is smooth). |
| `gen_seed` | `1` | Random seed. |
| `gen_decorrelate` | `.true.` | Remove the segment-mean spanwise waveform (no spectral line at the replay period). |
| `gen_verify` | `.true.` | Re-read the file and compare statistics with the targets. |

## &dopamine_inflow_record

Precursor run: record a plane of cell centres into an inflow database.

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `record_enabled` | `.false.` | Enable recording. |
| `record_file` | `'inflow_record.nc'` | Output database. |
| `record_point` | `0 0 0` | A point on the recording plane. |
| `record_normal` | `1 0 0` | Flow direction through the plane (horizontal). The plane must form a regular lattice across a spanwise-periodic domain. |
| `record_dt` | `0.0` s | Frame spacing. |
| `record_t_start` | `0.0` s | Skip the spin-up before this time. |
| `record_segment_time` | `0.0` s | Length of one segment. |
| `record_n_segments` | `1` | Number of segments; use at least 4. |

## &dopamine_bc

One entry per boundary patch; names must match the patch names in the mesh. Arrays are indexed by patch (maximum 32).

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_patches` | `0` | Number of entries used below. |
| `patch_name(i)` | `''` | Boundary patch name from the mesh (gmsh physical group or OpenFOAM patch). |
| `patch_type(i)` | `'zeroGradient'` | Boundary type, see the table below. |
| `patch_U(i,1:3)` | `0.0` m/s | Velocity for `velocity`/`inlet`/`movingWall` patches. |
| `patch_T(i)` | `T_ref` K | Temperature for `fixedTemperature` patches. |
| `patch_p(i)` | `0.0` | Pressure for `outlet`/`fixedPressure` patches. |
| `patch_outlet_mode(i)` | `'extrapolate'` | Outlet velocity face value: `'extrapolate'` (linear extrapolation of the cell velocity) or `'convective'` (advective condition for the face velocity, an RK3W state, see [Theory](theory.md#convective-outlet-and-backflow-traction)). Use `'convective'` with turbulent inflow. |
| `patch_backflow_beta(i)` | `0.0`, or `1.0` if `patch_outlet_mode(i) = 'convective'` | Backflow traction coefficient β of an outlet: the face pressure is lowered by `β·½|u|²` where the flow re-enters. `β ≥ 1` bounds the kinetic energy entering through the outlet; `0` switches it off. A negative value selects the default. |
| `patch_z0_m(i)`, `patch_z0_h(i)` | `0.0` m | Per-patch roughness for `wallModel` patches under MOST. `0` inherits `&dopamine_abl`. A negative `patch_z0_h` derives it from `z0_m` with ln(z0_m/z0_h) = 2. |
| `patch_scalar_mode(i,s)` | `'zeroGradient'` | Passive-scalar BC: `'dirichlet'`, `'neumann'`, or `'zeroGradient'`. |
| `patch_scalar_value(i,s)` | `0.0` | Value for `'dirichlet'`. |
| `patch_scalar_flux(i,s)` | `0.0` | Diffusive flux for `'neumann'`. |
| `patch_scalar_flux_omega(i,s)`, `patch_scalar_flux_phase(i,s)` | `0.0` | Harmonic modulation `flux(t) = flux·cos(ω t + φ)` of a Neumann flux. |

`patch_type` names are case-insensitive. An unrecognised name stops the run. Accepted aliases:

| Aliases | Meaning |
|--------------|------------------------------------------------|
| `velocity`, `inlet`, `fixedValue` | Velocity inlet (`patch_U`). |
| `outlet`, `fixedPressure`, `outflow` | Pressure outlet (`patch_p`). |
| `zeroGradient` | Zero gradient. |
| `noSlip` | No-slip wall. |
| `movingWall`, `wallVelocity` | Wall with prescribed velocity. |
| `wallModel`, `wall` | Wall-modelled wall (see `wall_model`). |
| `slip` | Free-slip wall. |
| `symmetry` | Symmetry plane. |
| `periodic`, `cyclic` | Periodic pair. |
| `fixedTemperature` | Fixed temperature (`patch_T`). |
| `ablInlet` | Turbulent ABL inflow (SEM or plane replay). |

## &dopamine_initial_condition

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `ic_type` | `'uniform'` | Initial field: `'uniform'`, `'taylor_green'`, `'linear'`/`'channel_linear'`, `'parabolic'`/`'channel_parabolic'`/`'poiseuille'`, `'loglaw'`/`'channel_loglaw'`, `'boundary_layer'`/`'channel_boundary_layer'`, `'scalar_gauss'`/`'scalar_pulse'`. |
| `ic_uref` | `1.0` m/s | Reference velocity for initialisations that need one. |
| `ic_lref` | `0.0` m | Reference length for initialisations that need one. |
| `noise_level` | `0.0` | Amplitude of random perturbation added where the initialisation supports it. |
| `noise_seed` | `0` | `> 0` makes the noise reproducible for a fixed rank count; `0` seeds from the clock. |
| `wall_normal_direction` | `3` | Wall-normal axis for wall-bounded initialisations (1 = x, 2 = y, 3 = z). |

## &dopamine_output

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `prefix` | `'dopamine'` | Prefix of the output files (XMF, NetCDF) and output directory. |
| `write_interval` | `100` | Field output every N steps. Ignored when negative (see `write_time`). |
| `write_time` | `0.0` s | Physical-time output cadence, active when `write_interval < 0` and `write_time > 0`. The step is snapped down (never up) to land on multiples. |
| `hdf5_use_file_locking` | `.false.` | `.false.` lets live readers open files during the run on most systems. |
| `warn_log_file` | `'solver_warnings.txt'` | File runtime warnings are appended to (eddy-viscosity caps, wall-model matching-point cap, ...). |
| `warn_flush_interval` | `100` | Steps between flushes of the warning log. |

## &dopamine_checkpoint

Per-rank binary checkpoints; restart and repartitioning workflow: [Usage: Checkpoint and restart](usage.md#checkpoint-and-restart).

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `checkpoint_interval` | `0` | Write a checkpoint every N steps; `0` disables. |
| `restart_enabled` | `.false.` | Start from a checkpoint. |
| `restart_dir` | `''` | Directory holding `processorR/checkpoint_XXXXXXXX.bin`. |
| `restart_step` | `-1` | `-1` reads `checkpoint_latest`; `N >= 0` loads step N. |

## &dopamine_scalars

Passive scalars (maximum 16).

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_scalars` | `0` | Number of scalars to solve. |
| `scalar_name(s)` | `''` | Name used for output variables. |
| `scalar_diffusivity(s)` | `0.0` m²/s | Molecular diffusivity. |
| `scalar_pr_t(s)` | `1.0` | Turbulent Prandtl/Schmidt number for the SGS flux. |
| `scalar_tvd_scheme` | `'VANLEER_TVD'` | `'VANLEER_TVD'`, `'KOREN_TVD'`, `'QUICK'`, `'CENTRAL'`/`'CDS'`. Applies to passive scalars only; temperature uses `temperature_scheme`. |

## &dopamine_scalar_sources

Injection sources for passive scalars (maximum 128). Each source is indexed by `i`.

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_scalar_sources` | `0` | Number of sources. |
| `source_scalar_id(i)` | `1` | Scalar the source feeds. |
| `source_name(i)` | `''` | Label. |
| `source_kind(i)` | `'volume'` | `'volume'`, `'point'`, `'line'`, `'surface'`/`'plane'`, or `'obj'`/`'stl'`/`'points'`. |
| `source_shape(i)` | `'sphere'` | Only for `kind = 'volume'`: `'sphere'`, `'box'`/`'aabb'`, `'plane'`/`'surface'`. |
| `source_magnitude(i)` | `0.0` | Source strength (scalar units per second per unit volume). |
| `source_time_mode(i)` | `'constant'` | `'constant'` or `'harmonic'`. |
| `source_omega(i)`, `source_phase(i)` | `0.0` | Angular frequency [rad/s] and phase [rad] of a harmonic source. |
| `source_center(i,1:3)` | `0.0` | Centre of a sphere/point source, or a point on a plane source. |
| `source_radius(i)` | `0.0` | Radius of sphere/point/line (capsule) sources, or the optional disk clip of a plane source. A sphere with radius `0` matches no cells and prints a warning. |
| `source_box_min(i,1:3)`, `source_box_max(i,1:3)` | `0.0` | Corners of a box source. |
| `source_p0(i,1:3)`, `source_p1(i,1:3)` | `0.0` | End points of a line source. |
| `source_normal(i,1:3)`, `source_half_thickness(i)` | `0.0` | Normal and half thickness of a plane source. |
| `source_file(i)`, `source_tag_radius(i)` | `''`, `0.0` | OBJ/STL/xyz file for `kind = 'obj'` and the capture radius around each file point. |

## Sampling groups

Time-resolved sampling of fields along analytic or file-defined geometry. Up to 16 objects per group. A
`write_frequency` of `0` falls back to the matching cadence in [`&dopamine_stats`](#dopamine_stats).

### &dopamine_sampling_planes

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_planes` | `0` | Number of planes. |
| `plane_name(i)` | `''` | Label. |
| `plane_origin(i,1:3)` | `0.0` | A point on the plane. |
| `plane_normal(i,1:3)` | `0.0` | Plane normal. |
| `plane_tangent_hint(i,1:3)` | `0.0` | Hint fixing the in-plane first axis. |
| `plane_half_width(i)`, `plane_half_height(i)` | `0.5` | Half extents of the sampled rectangle. |
| `plane_nx(i)`, `plane_ny(i)` | `64` | Sample points along each in-plane axis. |
| `plane_write_frequency(i)` | `slice_write_frequency` | Steps between writes. |

### &dopamine_sampling_lines

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_lines` | `0` | Number of lines. |
| `line_name(i)` | `''` | Label. |
| `line_p0(i,1:3)`, `line_p1(i,1:3)` | `0.0` | End points. |
| `line_npts(i)` | `128` | Sample points along the line. |
| `line_write_frequency(i)` | `line_write_frequency` | Steps between writes; when `0`, the larger of `line_write_frequency` and `probe_write_frequency`. |

### &dopamine_sampling_curves

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_curves` | `0` | Number of curves. |
| `curve_name(i)` | `''` | Label. |
| `curve_file(i)` | `''` | ASCII xyz (one point per line), OBJ vertices, or ASCII STL vertices. |
| `curve_write_frequency(i)` | `curve_write_frequency` | Steps between writes; when `0`, the larger of `curve_write_frequency` and `probe_write_frequency`. |

### &dopamine_sampling_surfaces

| Parameter | Default | Meaning / how to use |
|--------------------------|----------|------------------------------------------|
| `n_surfaces` | `0` | Number of surfaces. |
| `surface_name(i)` | `surface_<i>` | Label. |
| `surface_file(i)` | `''` | OBJ/STL/xyz point cloud sampled on the fly. |
| `surface_write_frequency(i)` | `probe_write_frequency` | Steps between writes. |
| `surface_intersect_with_planes(i)` | `.true.` | Also write the intersection curves with the sampling planes. |
| `surface_intersection_decimate_spacing(i)` | `intersection_decimate_spacing` | Thinning spacing of intersection curves. |
| `surface_sort_intersection_points(i)` | `sort_intersection_points` | Order intersection points along the curve. |

---

## Recipes

| Case | Parameters to set |
|--------------|------------------------------------------------|
| Periodic DNS (Taylor-Green) | `use_sgs = .false.`, `advection_scheme = 'KEP'`, `ic_type = 'taylor_green'`, all patches `periodic` |
| Wall-modelled channel LES | `sgs_model = 'AMD'` or `'VREMAN'`, `wall_model = 'loglaw'`, `fx` (u_tau drive) or `constant_mass_flux_forcing = .true.` with `Ux_constant`, `wallModel` walls, `periodic` x and y |
| Neutral ABL over roughness | `wall_model = 'most'`, `z0_m`, `u_ref`, `z_ref`, `ablInlet` patch with `use_sem = .true.`, `outlet` patch |
| Stratified ABL | as above plus `boussinesq = .true.`, `gz = -9.80665`, `T_surface`, `compute_L = .true.` |
| Large-scale forced ABL | add `use_nudging = .true.`, `nudging_file`, `nudging_tau`, `nudging_z_start` |
| Precursor/successor inflow | precursor: `&dopamine_inflow_record`; generator: `dopamine-inflow-gen`; successor: `inflow_source = 'planes'`, `inflow_file` |
| Rotating-frame or tidal flow | `use_rotation` with `omega_*`, or `use_harmonic_forcing` with `forcing_*` |
| Tracer release | `n_scalars`, `scalar_diffusivity`, one `&dopamine_scalar_sources` entry |
| Restartable long run | `checkpoint_interval`, then `restart_enabled = .true.`, `restart_dir` |
| Very large mesh | `dopamine-decompose`, then `mesh_type = 'preprocessed'` |
