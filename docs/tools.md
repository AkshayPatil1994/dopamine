---
title: Tools
nav_order: 9
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Tools

Standalone executables are built next to `dopamine` (`build/`) and share the solver library. Python
helpers live in the example directories and under `tests/`.

## dopamine-decompose

Partition a mesh offline so every rank reads only its own piece (`mesh_type = 'preprocessed'`).

```sh
mpirun -np 1 ./build/dopamine-decompose input.input --nranks 128 --output-dir ./mesh_128
```

The input file supplies `mesh_file` and `mesh_type` from `&dopamine_domain`. Flags, memory behaviour and
the load-balance report: [Pre-processing](preprocessing.md#decomposer-options).

## dopamine-checkpoint-partition

Convert checkpoints written by an `N`-rank run into checkpoints for an `M`-rank run.

| Flag | Meaning |
|---|---|
| `--old-dir DIR` | output directory of the original run (holds `processorR/checkpoint_*.bin`) |
| `--old-nranks N` | rank count of the original run |
| `--new-dir DIR` | directory of the mesh decomposed for the new rank count |
| `--new-nranks M` | new rank count |
| `--step S` | checkpoint step to convert |

Full workflow: [Usage: Checkpoint and restart](usage.md#checkpoint-and-restart).

## dopamine-inflow-gen

Generate a database of spanwise-periodic inflow planes from the synthetic eddy description.

```sh
./build/dopamine-inflow-gen generate.input
```

It reads `&dopamine_abl` (the eddy description) and `&dopamine_inflow_gen`
([Input Reference](input-reference.md#dopamine_inflow_gen)) and writes a NetCDF database that a successor run
replays with `inflow_source = 'planes'`. Set `inflow_export_grid` in the solver input first to export the inlet
lattice for the generator. Workflow: [Usage: Inflow-plane databases](usage.md#inflow-plane-databases-precursor-and-successor).

## Example and test scripts

| Script | Purpose |
|---|---|
| `examples/channelAMDvalidation/postprocess.py` | Compare AMD/Vreman channel runs against reference profiles |
| `examples/*/plot_validation.py`, `examples/common/animate_slices.py` | Validation figures and animations for [Validation](validation.md); shared style in `examples/common/plotstyle.py` |
| `examples/convergence_test/run_convergence.py`, `plot_convergence.py` | Run and plot the spatial convergence study |
| `examples/convergence_test/run_advection.py`, `plot_advection.py` | Advection-scheme accuracy study |
| `examples/lid_driven_cavity/compare_grids.py`, `animate_cavity.py` | Cavity grid comparison and animation |
| `examples/RayleighBernard/RBAnalysis.py` | Analysis and animation of the Rayleigh-Benard run |
| `tests/regression/run_case.sh`, `check_case.sh`, `compare_stats.py` | Run a regression case and compare its per-step statistics with the reference CSV |
| `tests/regression/mk_box_mesh.py` | Generate the skewed box mesh used without gmsh |
| `tests/sem_pressure/mk_box.py`, `run_study.py` | SEM inflow pressure study (not part of CTest) |

Regression cases are described in `tests/regression/README.md`.
