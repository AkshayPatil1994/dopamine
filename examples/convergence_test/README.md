# Convergence test suite

Measures the spatial and temporal order of accuracy on the 3D Taylor-Green vortex
(periodic cube, `CENTRAL` advection, `RK3W`, no SGS). Results and discussion:
[docs/theory.md](../../docs/theory.md), section "Verification: Spatial and Temporal Convergence".

## Run

```bash
python3 run_convergence.py all --bin ../../build/dopamine --workdir results
python3 plot_convergence.py results/results.json --out ../../docs/assets/convergence
```

Needs `gmsh`, `mpirun`, Python 3 with `numpy` and `matplotlib`. Options: `--np N` (MPI ranks),
`--force` (redo finished runs); `spatial` / `temporal` run one study. Runs that already
finished are reused. The full study takes about 10 minutes on one GPU.

## Files

- `run_convergence.py`: generates meshes, runs the cases, writes `results/results.json`.
- `plot_convergence.py`: makes the figures and prints the observed orders.
- `warp_mesh.py`: turns a hex mesh into a smoothly warped periodic one (non-orthogonal cells).
- `projection_model.py`: numpy model of one RK3W step used to trace the O(dt h^2) energy loss.
- `inputs/tgv_base.input`: the solver input template.
- `run_advection.py`, `plot_advection.py`, `inputs/advection_base.input`, `meshes/advection_1d.geo`: scalar
  advection scheme comparison on a periodic bar (needs `netCDF4`; the profile is read from the line sampler).
- `meshes/hex.geo`: periodic structured hex cube; `N` is substituted by the driver.

`meshes/tet.geo`, `prism.geo`, `skewed_hex.geo` and `spatial_convergence/` are older and
unused: the tet and prism meshes are rejected by the periodic matching, and the Taylor-Green
initial condition is not periodic on the sheared cell.

## Pitfalls handled by the driver

- `gmsh -setnumber N` does not override `hex.geo`; the driver writes `N` into a copy.
- `prefix` is a `character(32)`; runs use the relative prefix `tgv` inside their own directory.
- The solver clips the step at `cfl_max` even with `adapt_dt = .false.`; the template sets
  `cfl_max = 100` so the requested `dt` is the one used.
- `use_sgs` defaults to `.true.`; the template turns it off.
