---
title: Examples
nav_order: 7
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Examples

Each case lives in `examples/<case>/` with a gmsh geometry (`*.geo`) and an input file. Parameter meanings:
[Input Reference](input-reference.md).

## Running a case

```sh
cd examples/taylor_green
gmsh taylor_green.geo -3 -o taylor_green.msh
mpirun -np 1 ../../build/dopamine-decompose decompose.input --nranks 8 --output-dir mesh
mpirun -np 8 ../../build/dopamine taylor_green.input
```

`decompose.input` points `mesh_file` at the `.msh` file with `mesh_type = 'gmsh'`; the run input points at the
decomposed `mesh` directory with `mesh_type = 'preprocessed'`. A single-rank run can use the gmsh form
directly. Output is `<prefix>.nc` plus `<prefix>.xmf` for ParaView.

## Overview

| Case | Physics | Key settings | Reference / check |
|-------|-----------|-------------------|-------------|
| `taylor_green` | Periodic Taylor-Green vortex, Re = 1600 | `nu = 6.25e-4`, `advection_scheme = 'QUICK'`, `ic_type = 'taylor_green'`, `t_end = 20` | Kinetic energy and dissipation against DNS (figure in the repository `figures/tgv_validation.png`) |
| `lid_driven_cavity` | Lid-driven cavity, Re = 1000 | `nu = 1e-3`, `advection_scheme = 'KEP'`, `sgs_model = 'VREMAN'`, `t_end = 100` | Centreline profiles; `compare_grids.py` compares resolutions |
| `channel_flow_retau395` | Wall-modelled LES channel, Re_tau = 395 | `nu = 2.5316e-3`, `fx = 1`, `advection_scheme = 'KEP'`, `ic_type = 'channel_loglaw'`, `t_end = 20` | Mean profile and stresses against channel DNS |
| `RayleighBernard` | Buoyancy-driven convection, Pr = 1, Ra = 2.5e6 | `boussinesq = .true.`, `nu = alpha_T = 6.3246e-4`, `advection_scheme = 'KEP'`, `t_end = 150` | Qualitative plume structure; `RBAnalysis.py` |
| `channelAMDvalidation` | Wall-modelled half channel, Re_tau = 2000 | `sgs_model = 'AMD'`, `wall_model_filter_tf_fac`, 64^3 uniform grid | Mean-velocity profile against Hoyas & Jimenez DNS (data supplied by the user); see its `README.md` |
| `convergence_test` | Spatial and temporal order on Taylor-Green | `CENTRAL` advection, `RK3W`, no SGS | Observed order of accuracy; results in [Theory](theory.md) |

## Regression cases

`ctest -L regression` runs 12 short cases (TGV, channel, Rayleigh-Benard, duct, dynamic Smagorinsky, graded
mesh, hydrostatic and pressure-drop checks, SEM inflow, AMD) at one and two ranks and compares per-step
statistics with stored references. They protect results during refactoring and are not validation cases. See
`tests/regression/README.md`.
