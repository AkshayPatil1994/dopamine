# Examples

Each directory is a self-contained case: a gmsh geometry (`*.geo`), an input file, and where relevant
analysis scripts. Generate the mesh, decompose it (or use `mesh_type = 'gmsh'`), then run the solver.
Parameter meanings are in [docs/input-reference.md](../docs/input-reference.md).

| Case | What it demonstrates | Details |
|---|---|---|
| `taylor_green/` | Periodic Taylor-Green vortex at Re = 1600; kinetic-energy decay and dissipation | `tutorial.md` |
| `lid_driven_cavity/` | Wall-bounded laminar flow, grid comparison | `tutorial.md`, `compare_grids.py`, `animate_cavity.py` |
| `channel_flow_retau395/` | Wall-resolved LES channel at Re_tau = 395 | `tutorial.md` |
| `RayleighBernard/` | Buoyancy-driven convection with the Boussinesq term | `tutorial.md`, `RBAnalysis.py` |
| `channelAMDvalidation/` | AMD subgrid model and adaptive wall-model filter in a wall-modelled channel | `README.md`, `postprocess.py` |
| `convergence_test/` | Spatial and temporal order of accuracy on Taylor-Green | `README.md` |

## Common workflow

```sh
gmsh -3 -format msh4 CASE.geo                         # build the mesh
mpirun -np 1 ./dopamine-decompose decompose.input --nranks 4 --output-dir mesh
mpirun -np 4 ./dopamine CASE.input
```

Cases that read `mesh_type = 'gmsh'` skip the decompose step. The regression suite
(`ctest -L regression`) exercises the same physics at small size; see
[docs/development.md](../docs/development.md#testing).
