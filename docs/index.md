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
Incompressible Navier–Stokes Equations) is an unstructured, cell-centered
Large Eddy Simulation (LES) code for incompressible flows, written in modern
Fortran with a PETSc-backed pressure solver and NetCDF/XMF output.

This site has three guides:

- **[Theory](theory)** — the governing equations, spatial and temporal
  discretization, advection/gradient/pressure operators, SGS and wall
  modeling, the atmospheric boundary layer (ABL) module, and the
  implementation choices (and their rationale) that affect robustness and
  accuracy on unstructured meshes.
- **[Usage](usage)** — build and run instructions, input namelists and
  parameters, mesh partitioning, output/visualization, the ABL workflow, and
  a long-time stability checklist, with worked examples.
- **[Development](development)** — code layout, build system internals,
  testing, coding style, and how to add a scheme, boundary condition, or ABL
  feature; how to report bugs and issues.

## Getting Started

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/dopamine examples/taylor_green/taylor_green.input
```

See [Usage](usage#build) for the full build matrix (MPI, OpenMP, GPU/CUDA,
OpenACC) and [Development](development#build-system) for build-system
internals.

## Source and License

Source: [github.com/AkshayPatil1994/dopamine](https://github.com/AkshayPatil1994/dopamine).
DOPAMINE is licensed under [AGPL-3.0-only](https://github.com/AkshayPatil1994/dopamine/blob/main/LICENSE).

## Standalone Manual

This documentation also builds as a single PDF/HTML manual via Pandoc,
independent of GitHub Pages:

```bash
bash docs/build.sh
```

See [docs/README.md](https://github.com/AkshayPatil1994/dopamine/tree/main/docs)
in the repository for details.
