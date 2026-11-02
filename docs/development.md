---
title: Development
nav_order: 9
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
| `src/abl/` | Atmospheric boundary layer: MOST profiles/wall coupling, weather nudging, divergence-free SEM inflow, inflow-plane database (`mod_inflow_db`), generator (`mod_inflow_gen`, `dopamine-inflow-gen`), replay (`mod_inflow_planes`) and precursor recorder (`mod_inflow_record`) |
| `src/io/` | Input namelist parsing, NetCDF/XMF output, checkpoint I/O, runtime environment, warning log |
| `src/stats/` | Runtime statistics and advanced sampling (planes/lines/curves/surfaces) |
| `src/init/` | Initial condition assembly |
| `src/tools/` | Standalone executables: `dopamine-decompose`, `dopamine-checkpoint-partition`, `dopamine-inflow-gen` |

`src/main.f90` is the solver entry point. Everything except `main.f90` and the
`src/tools/*_main.f90` drivers is built into `dopamine_lib`
(the `dopamine_lib` target in [CMakeLists.txt](../CMakeLists.txt)), which the `dopamine`,
`dopamine-decompose`, `dopamine-checkpoint-partition`, and
`dopamine-inflow-gen` executables link against.

Module naming follows `mod_<area>.f90` → module `mod_<area>`. Keep a module's
public interface at the top of the file (a `public ::`/`private` block) so a
reader can see the contract before the implementation.

## Build System

CMake drives everything, including fetching and building PETSc, hypre,
METIS/GKlib, and NetCDF when not found on the system
(`-DDOPAMINE_FETCH_PETSC=ON`, etc. — see [installation.md](installation.md) for the
full flag list, offline/HPC builds, and the PETSc standalone recipe). Two
points that matter during development:

- `make clean` removes only DOPAMINE's own artifacts; `make cleanall` also
  wipes the fetched external libraries under `DOPAMINE_EXTERNAL_DIR`
  (`external/` by default). Use `cleanall` after changing an external-library
  build option (e.g. `DOPAMINE_PETSC_WITH_CUDA`), not for routine rebuilds —
  rebuilding PETSc from scratch costs minutes.
- The CUDA (`build-gpu/`) and OpenACC (`build-acc26/`) paths use
  separate build/external directories so they don't disturb a CPU build
  sharing the same source tree. Follow [installation.md](installation.md#gpu-cuda-build)'s
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

Tests are registered with CTest in `tests/CMakeLists.txt`. Unit tests are standalone Fortran programs
(`tests/unit/test_*.f90`) built against `dopamine_lib`; regression cases are short solver runs compared
against stored statistics. Unit tests are built on demand (`EXCLUDE_FROM_ALL`), and regression cases need
`gmsh` on `PATH`:

```bash
cmake --build build -j
cmake --build build --target test_amd test_inflow_db test_inflow_gen test_inflow_planes test_sem_dfsem \
  test_mesh_geometry test_nonortho_limiter test_cell_projection
ctest --test-dir build --output-on-failure
```

Run a subset by label or by name (names match `add_test(NAME ...)`, not the source filename):

```bash
ctest --test-dir build -L unit --output-on-failure
ctest --test-dir build -L regression
ctest --test-dir build -R regression_tgv_np2 --output-on-failure
```

### What the existing tests cover

- **Unit** (`unit_*`): `unit_amd` (AMD eddy viscosity), `unit_sem_dfsem` (divergence-free eddy field),
  `unit_inflow_db`, `unit_inflow_gen`, `unit_inflow_planes` (inflow database, generator and replay),
  `unit_mesh_geometry`, `unit_gmsh_open_cell` (mesh geometry and open-cell reading),
  `unit_nonortho_limiter` (`pressure_min_cos_theta`) and `unit_cell_projection` (`cell_projection = 'EXACT'`).
- **Regression** (`regression_<case>_np<N>`): `tgv`, `channel`, `rb`, `duct`, `dynsmag`, `graded`, `pdrop`,
  `hydro`, `hydrowall`, `inflow`, `inflow_jarrin`, `tgv_skew`, `tgv_skew_chk`, `tgv_skew_exact`, `hydro_exact`,
  `duct_exact` and `default_skew` at one and two ranks, and `amd` at one rank. Every case pins `rhie_chow_dual`,
  `cell_projection` and `pressure_min_cos_theta` in its input except `default_skew`, which runs on the defaults. Each run is
  20 steps whose per-step statistics must match `tests/regression/reference/*.csv`; the np2 runs check
  that the result does not depend on the decomposition. The cases are described in
  `tests/regression/README.md`.
- **Inflow chain** (`inflow_chain`): generate planes with `dopamine-inflow-gen`, replay them in a
  successor run.

When you touch a numerics or BC path, find the case above that already exercises it before writing a
new one. If none fits, add a `tests/unit/test_*.f90` program or a regression case following the closest
existing one and register it in `tests/CMakeLists.txt`.

### Case-level regression

`examples/` holds runnable cases used for physical validation beyond unit
tests (Taylor-Green vortex, channel flow at Re_τ=395, lid-driven cavity, ABL
SEM-driven channel, etc. — see [examples/README.md](https://github.com/AkshayPatil1994/dopamine/blob/main/examples/README.md)).
There is no CI workflow running these automatically today; validating a
numerics or physics change against the relevant example case before opening a
PR is a manual step. Each example directory carries its own analysis scripts
(for example `examples/channelAMDvalidation/postprocess.py` and
`examples/convergence_test/plot_convergence.py`).

## Style

- **Keep comments short.** No comment blocks of three or more lines, in Fortran or C. A module or
  subroutine docstring is a short `!>` header. If an inline comment needs a long explanation, the code
  needs a better name or structure, or the rationale belongs in this documentation with a one-line
  pointer left in the code.
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
  the paths silently diverge; the np1 and np2 regression cases exercise both.

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
3. Add a boundedness or conservation test for the scheme, and an SoA-consistency
   check if you also add a GPU/SoA path; model both on the closest test in
   `tests/unit/` and register them in `tests/CMakeLists.txt`.
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
  unrelated cleanup — it keeps `git blame` and bisection attributable to
  a single change.
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
production mesh — the regression cases in `tests/regression/cases/` (for example
`pdrop` and `hydro`) are exactly this kind of minimal reproduction.

## Documentation

This `docs/` directory is both standalone-buildable (`bash docs/build.sh`,
producing a single PDF/HTML manual via Pandoc — see
[docs/README.md](README.md)) and published as a GitHub Pages site directly
from the `docs/` folder (Jekyll + [just-the-docs](https://just-the-docs.com),
configured in `docs/_config.yml`). When you add a new `##`-level section to
any page, no extra step is needed for either output:
Pandoc's `--toc` and just-the-docs' page nav both derive structure from
heading levels and each page's `nav_order` front matter.

If you add a new top-level doc page, give it Jekyll front matter (`title`,
`nav_order`) so it appears in the Pages sidebar, and add it to the `cat`
sequence in [docs/build.sh](build.sh) if it should also appear in the PDF/HTML
manual.
