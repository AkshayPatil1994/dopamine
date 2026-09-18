---
title: Development
nav_order: 4
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# DOPAMINE Development Guide

This document is for people modifying DOPAMINE itself: adding a numerical
scheme, fixing a bug, extending the ABL module, or touching the build system.
It covers code layout, style, testing, and the contribution workflow. For
governing equations and discretization see [theory.md](theory.md); for
running the solver see [usage.md](usage.md).

## Code Layout

DOPAMINE is modern Fortran (2008+) with a small C shim for the PETSc pressure
backend. Source is organized by responsibility under `src/`:

| Directory | Contents |
|---|---|
| `src/mesh/` | Mesh readers (gmsh, OpenFOAM), partitioner, streaming decomposer, validation, non-conforming interfaces |
| `src/numerics/` | Advection schemes, gradient operators, non-orthogonal correction, viscous fluxes |
| `src/solver/` | Time integration (`mod_fv_solver.f90`), pressure Poisson solve (`mod_pressure_solver.f90`), and the PETSc C backend (`pressure_backend_petsc.c`) |
| `src/bc/` | Boundary condition types and per-patch setup |
| `src/sgs/` | SGS models (Vreman, constant-coefficient and dynamic Smagorinsky), wall models (log-law, MOST) |
| `src/abl/` | Atmospheric boundary layer: MOST profiles/wall coupling, weather nudging, SEM synthetic inflow |
| `src/io/` | Input namelist parsing, NetCDF/XMF output, checkpoint I/O, runtime environment, warning log |
| `src/stats/` | Runtime statistics and advanced sampling (planes/lines/curves/surfaces) |
| `src/init/` | Initial condition assembly |
| `src/tools/` | Standalone executables: `dopamine-decompose`, `dopamine-checkpoint-partition` |

`src/main.f90` is the solver entry point. Everything except `main.f90` and the
`src/tools/*_main.f90` drivers is built into `dopamine_lib`
([CMakeLists.txt:559](../CMakeLists.txt)), which the `dopamine`,
`dopamine-decompose`, and `dopamine-checkpoint-partition` executables link
against.

Module naming follows `mod_<area>.f90` → module `mod_<area>`. Keep a module's
public interface at the top of the file (a `public ::`/`private` block) so a
reader can see the contract before the implementation.

## Build System

CMake drives everything, including fetching and building PETSc, hypre,
METIS/GKlib, and NetCDF when not found on the system
(`-DDOPAMINE_FETCH_PETSC=ON`, etc. — see [usage.md](usage.md#build) for the
full flag list, offline/HPC builds, and the PETSc standalone recipe). Two
points that matter during development:

- `make clean` removes only DOPAMINE's own artifacts; `make cleanall` also
  wipes the fetched external libraries under `DOPAMINE_EXTERNAL_DIR`
  (`external/` by default). Use `cleanall` after changing an external-library
  build option (e.g. `DOPAMINE_PETSC_WITH_CUDA`), not for routine rebuilds —
  rebuilding PETSc from scratch costs minutes.
- The CUDA (`build-gpu/`) and OpenACC (`build-acc26/`) paths use
  separate build/external directories so they don't disturb a CPU build
  sharing the same source tree. Follow [usage.md](usage.md#gpu-cuda-build)'s
  "GPU (CUDA) build" / "OpenACC (nvfortran) build" sections rather than
  reusing `build/`.

A typical development build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DDOPAMINE_ENABLE_PETSC=ON
cmake --build build -j
```

Use `-DCMAKE_BUILD_TYPE=Debug` (bounds checking, no `-O3`) while iterating;
switch to `Release` before checking performance-sensitive changes, since
Fortran array-bounds checking materially changes hot-loop timings.

## Testing

Unit and regression tests live in `tests/` and are registered with CTest via
`tests/CMakeLists.txt`. Each test is a standalone Fortran program
(`test_*.f90`) built against `dopamine_lib`. Run the suite from the build
directory:

```bash
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run a single test by name (names match `add_test(NAME ... COMMAND test_*)` in
[tests/CMakeLists.txt](../tests/CMakeLists.txt), not the source filename):

```bash
ctest --test-dir build -R RhieChowConsistency --output-on-failure
```

### What the existing tests cover

The suite is organized around the invariants theory.md calls out, not just
"does it run":

- **Geometry/gradients**: `MeshGeometry`, `Gradient`, `GradientVector`,
  `GradientRefinement`, `ReconstructionLinearExact` — gradient schemes must
  reproduce a linear field exactly and behave correctly across non-conforming
  refinement interfaces.
- **Advection**: `KEPConservation` (discrete kinetic-energy preservation),
  `QuickBoundedness` (TVD boundedness), `AdvectionKepSoa`/`AdvectionCdsSoa`/
  `AdvectionQuickTvdSoa`/`AdvectionVectorQuickTvdSoa` (SoA/GPU-mirror paths
  match the AoS reference), `AdvectionConvergence`, `Advection2DSnapshot`.
- **Pressure/projection**: `PressureSolver`, `HelmholtzProjection`,
  `RhieChowConsistency` — the last one directly checks the discrete duality
  between the Poisson operator and the velocity-correction stencil described
  in [theory.md](theory.md#velocity-correction-and-discrete-duality); a
  regression here is a silent divergence-control bug, not a crash.
- **SGS/wall models**: `VremanSGS`, `VremanFilter`, `DynamicSmagorinsky`,
  `WallModel`, `WallModelSelect`, `ABLMost`, `SEMProfile`.
- **Parallel consistency**: `PeriodicParallelMapping`,
  `ProcessorFaceConsistency` (MPI-enabled builds only — guarded in
  `tests/CMakeLists.txt`).
- **I/O**: `WriteTimeSnap`.

When you touch a numerics or BC path, look for the test in this list that
already exercises it before writing a new one — several tests were added
specifically because a bug of that shape shipped once (e.g.
`RhieChowConsistency`, `PeriodicParallelMapping`). If none fits, add a new
`test_*.f90` following the pattern of the closest existing one and register it
in `tests/CMakeLists.txt`.

### Case-level regression

`examples/` holds runnable cases used for physical validation beyond unit
tests (Taylor-Green vortex, channel flow at Re_τ=395, lid-driven cavity, ABL
SEM-driven channel, etc. — see [examples/README.md](../examples/README.md)).
There is no CI workflow running these automatically today; validating a
numerics or physics change against the relevant example case before opening a
PR is a manual step. `tools/` has analysis scripts
(`convergence_analysis.py`, `analyze_channel_profiles.py`,
`compare_cavity.py`) for post-processing those runs.

## Style

- **Comments are one line.** No multi-line `!`/`!>` blocks anywhere, Fortran
  or C. A module or subroutine docstring is a single `!>` line. If an inline
  comment needs more than one line to explain, either the code needs a better
  name/structure, or the rationale belongs in a doc under
  `docs/external-lib-guidance/` (GPU porting decisions do this today) with a
  single-line pointer left in the code.
- `!$acc`/`!$omp` directive continuation lines are compiler directives, not
  comments, and are exempt.
- No docstrings on trivial getters/setters.
- Don't repeat a self-explanatory variable name's meaning in a comment.
- `error stop` messages are short (< 80 chars) — put diagnostic detail in a
  preceding `write` to stderr/log if more context is genuinely needed, not in
  the stop message itself.
- Match the surrounding module's conventions for `intent`, `pure`/`elemental`
  use, and SoA vs. AoS field layout — numerics kernels that have an OpenACC
  port (the SGS kernels in `mod_vreman.f90` and
  `mod_dynamic_smagorinsky.f90`, the gradient kernels in `mod_gradient.f90`,
  `viscous_rhs_vector`/`viscous_rhs_scalar` in `mod_viscous.f90`, and the CDS/
  reconstruction/KEP paths in `mod_advection.f90`) carry both an
  AoS and an SoA/mirror code path; changes to the math need to land in both or
  the paths silently diverge (this is exactly what `AdvectionKepSoa` /
  `AdvectionCdsSoa` / `AdvectionQuickTvdSoa` exist to catch).

## Adding a New Advection Scheme

Advection schemes are selected by name (`advection_scheme`,
`temperature_scheme`, `scalar_tvd_scheme`) and implemented in
`src/numerics/mod_advection.f90`. To add one:

1. Add the scheme identifier and dispatch case alongside the existing
   `CENTRAL`/`KEP`/`QUICK`/`VANLEER_TVD`/`KOREN_TVD` cases.
2. If the scheme needs a wider stencil than the immediate face neighbors
   (QUICK and the TVD schemes do), confirm it resolves periodic far-upwind
   neighbors through the periodic-pairing path — see "TVD far-upwind across
   periodic faces" in [theory.md](theory.md); a scheme that silently drops to
   first-order on periodic boundaries is the historical failure mode here.
3. Add a boundedness/conservation test modeled on `QuickBoundedness` or
   `KEPConservation`, and an SoA-consistency test modeled on
   `AdvectionQuickTvdSoa` if you also add a GPU/SoA path.
4. Document the scheme in [theory.md](theory.md) (family + citation) and
   [usage.md](usage.md) (name table + when to prefer it), and add it to the
   advection scheme table in the top-level [README.md](../README.md).

## Adding a Boundary Condition or ABL Feature

BC types live in `src/bc/mod_bc.f90` (types/enumeration) and
`src/bc/mod_bc_setup.f90` (per-patch namelist parsing and setup). ABL-specific
BCs (MOST wall model, SEM inflow) additionally hook into `src/abl/`. Follow
the existing pattern of gating new behavior behind an explicit input flag with
a safe default that preserves old input files' behavior — `use_abl` staying a
legacy alias for `wall_model = 'most'` rather than gating unrelated SEM/
nudging flags is the precedent (see "Enabling ABL Features" in
[usage.md](usage.md)).

## Git Workflow

- Commit messages: short, no long bodies. One
  line stating what changed and, if not obvious, why.
- Prefer focused commits over one commit mixing a numerics change with
  unrelated cleanup — it keeps `git blame` and future audits (like the
  code-review memory entries tracked for this project) attributable to a
  single change.
- The repository carries `main` alongside topic branches such as `PETSc-GPU`;
  check `git log --oneline -5` and the branch name before assuming which
  capabilities are present.

## Reporting Issues and Bugs

Open a bug report against
[.github/ISSUE_TEMPLATE/bug_report.yml](../.github/ISSUE_TEMPLATE/bug_report.yml)
(the form GitHub presents when you click "New issue"). This repository does
not currently run a CI workflow, so a well-specified report is the main
signal a maintainer has to reproduce and triage from. The template asks for:

- The specific branch and commit you observed it on
  (`git rev-parse --abbrev-ref HEAD` / `git rev-parse --short HEAD`), not just
  "main" — GPU/OpenACC behavior in particular depends on which build flags and
  branch you are on.
- The exact input file (or the minimal namelist reproducing it), rank count,
  and build flags (`DOPAMINE_ENABLE_PETSC`, `DOPAMINE_PETSC_WITH_CUDA`,
  compiler). Numerics bugs on this codebase have repeatedly turned out to be
  interpolant/coefficient inconsistencies that only show up at specific mesh
  skewness, rank counts, or scheme combinations (see theory.md's periodic
  pressure coupling and TVD far-upwind sections) — the reproduction details
  are usually the fastest path to root cause.
- For a suspected divergence-control or conservation regression, the relevant
  runtime diagnostics (`max|divU|`, pressure residual, `E`) alongside the
  symptom; these are cheaper to compare across commits than re-deriving the
  numerics from scratch.

Prefer a minimal case built in the style of `examples/` over a large
production mesh — the existing regression tests (`RhieChowConsistency`,
`PeriodicParallelMapping`, etc.) were all born from exactly this kind of
minimal reproduction.

## Documentation

This `docs/` directory is both standalone-buildable (`bash docs/build.sh`,
producing a single PDF/HTML manual via Pandoc — see
[docs/README.md](README.md)) and published as a GitHub Pages site directly
from the `docs/` folder (Jekyll + [just-the-docs](https://just-the-docs.com),
configured in `docs/_config.yml`). When you add a new `##`-level section to
theory.md, usage.md, or this file, no extra step is needed for either output:
Pandoc's `--toc` and just-the-docs' page nav both derive structure from
heading levels and each page's `nav_order` front matter.

If you add a new top-level doc page, give it Jekyll front matter (`title`,
`nav_order`) so it appears in the Pages sidebar, and add it to the `cat`
sequence in [docs/build.sh](build.sh) if it should also appear in the PDF/HTML
manual.
