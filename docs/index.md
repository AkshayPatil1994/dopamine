---
title: Home
nav_order: 1
description: "DOPAMINE documentation home"
permalink: /
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# DOPAMINE

**DOPAMINE** (Discrete Operator-based Parallel Algorithm for Modeling
Incompressible Navier-Stokes Equations) is an unstructured, cell-centred
Large Eddy Simulation (LES) code for incompressible flows on general polyhedral
meshes. It is written in modern Fortran with MPI, a PETSc-backed pressure solver,
optional OpenACC/CUDA offload, and NetCDF/XDMF output.

## Where to start

| I want to... | Go to |
|---|---|
| Build the solver (CPU, GPU, offline cluster) | [Installation](installation) |
| Run my first case | [Examples](examples) and [Usage](usage#run) |
| Prepare a mesh or decompose a large one | [Pre-processing](preprocessing) |
| Look up an input parameter, default or unit | [Input Reference](input-reference) |
| Set up an atmospheric boundary layer or urban case | [Usage: ABL module](usage#atmospheric-boundary-layer-abl-module) |
| Understand the numerics and models | [Theory](theory) |
| Use the helper executables and scripts | [Tools](tools) |
| Change or extend the code | [Development](development) |
| Find the literature | [References](references) |

## At a glance

- **Meshes:** gmsh v4 and OpenFOAM `polyMesh`; hex, tet, prism, pyramid and general polyhedra; offline or in-memory partitioning.
- **Numerics:** cell-centred finite volume, Wray RK3, fractional-step projection with rotational correction, least-squares or Green-Gauss gradients, non-orthogonal correction.
- **Advection:** central, kinetic-energy preserving, QUICK, Koren and van Leer TVD.
- **Turbulence:** Vreman, Smagorinsky, dynamic Smagorinsky and AMD subgrid models; log-law and Monin-Obukhov wall models.
- **Atmospheric flows:** synthetic eddy inflow, precursor/successor plane databases, Davies nudging, Boussinesq buoyancy.
- **Extras:** passive scalars with sources, rotation and harmonic forcing, sampling planes/lines/surfaces, checkpoint/restart with repartitioning.
- **Parallel:** MPI domain decomposition with halo exchange, optional OpenMP, OpenACC and multi-GPU.

## Quick start

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cd examples/taylor_green
gmsh taylor_green.geo -3 -o taylor_green.msh
# use mesh_type = 'gmsh' (see decompose.input) for a single rank
mpirun -np 1 ../../build/dopamine decompose.input
```

Details, including PETSc and GPU options, are in [Installation](installation).

## Source, license, citation

Source: [github.com/AkshayPatil1994/dopamine](https://github.com/AkshayPatil1994/dopamine).
DOPAMINE is licensed under [AGPL-3.0-only](https://github.com/AkshayPatil1994/dopamine/blob/main/LICENSE).

## Standalone manual

These pages also build as a single PDF/HTML manual with Pandoc (`bash docs/build.sh`); see
[docs/README.md](https://github.com/AkshayPatil1994/dopamine/tree/main/docs).
