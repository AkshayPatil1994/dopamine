# Regression cases

Short (20 fixed steps) runs whose per-step integral stats are compared against the
reference CSVs in `reference/`. They exist to catch changes in results while the solver is
being restructured (e.g. the GPU port); they are not validation cases.

| case | mesh | exercises |
|---|---|---|
| `tgv` | 32^3 periodic box | QUICK advection, Green-Gauss gradients, periodic-only BCs |
| `channel` | 32x16x16 (2pi x 2pi x 1, walls in z, Taylor-Green initial field) | KEP, Vreman SGS, log-law wall model, body force, passive scalar with a Dirichlet wall |
| `rb` | 24x24x16 | Boussinesq buoyancy, fixed-temperature walls, van Leer TVD, least-squares gradients, rotational projection |
| `duct` | 32x16x16 (4 x 1 x 1) | velocity inlet, fixed-pressure outlet, no-slip walls, CDS; the Dirichlet-pressure terms of the projection |
| `dynsmag` | 32x16x16 (as `channel`) | Germano-Lilly dynamic Smagorinsky: test filter and the halo exchanges of its packed fields |
| `graded` | 32x16x16, wall bump 0.01 | as `channel`, but the grading creates refinement-interface faces (interface diagnostics on); KEP advection must stay on the device path |
| `pdrop` | skewed 12^3 box (`mk_box_mesh.py`, no gmsh) | outlets at p = 1 and 0 with no body force: Umax must grow as t; Dirichlet outlets with non-orthogonal correctors |
| `hydro` | same mesh | outlets at p = 1 and 0 balancing `fz = -1`: the fluid must stay at rest |
| `hydrowall` | same mesh, all six faces no-slip walls | closed box under `fz = -1`: the wall pressure gradient `dp/dn = f.n`; Umax stays O(dt g) |
| `inflow` | same mesh | divergence-free SEM (`sem_type = 'dfsem'`, the default) `ablInlet` at x0, outlet at xL, periodic y, slip z, CDS: eddies wrap across the periodic seam; the result must not depend on the rank count (np 1, 2, 3 agree to 1e-14) |
| `inflow_jarrin` | same mesh | as `inflow` with the original tent-kernel SEM (`sem_type = 'jarrin'`) |
| `amd` | as `channel` | AMD SGS (`sgs_model = 'AMD'`) with the log-law wall-cell gradient, and the adaptive wall-model filter (`wall_model_filter_tf_fac = 1`); np 1 only |

Run everything through CTest from a build tree (meshes are generated with gmsh on first use):

```sh
ctest --test-dir build -L regression --output-on-failure
```

or one case by hand:

```sh
DOPAMINE_BIN=build/dopamine tests/regression/check_case.sh channel 2 /tmp/regression
```

Each case has one reference per rank count (`<case>_np1.csv`, `<case>_np2.csv`): results
depend on the partition at the 1e-6..1e-3 level (e.g. `nu_sgs_mean` in `channel`), far above
build-to-build noise, so a run is only ever compared with the same rank count. Integral
quantities (`E`, `Umax`, `Ub`, `tauw_mean`, `nu_sgs_mean`) must match to `REGRESSION_RTOL`
(default 1e-6) and the face-flux divergence `max_abs_divU` must stay below 1e-9. CPU and GPU
builds are not bitwise identical (reduction order, FMA, AMG smoother choice), hence the
tolerance.

The stored references come from the CPU build (gfortran, `OMP_NUM_THREADS=1`) of commit
c6f4e6b plus the direct-read `-np > 1` deadlock fix.

The `duct` references were regenerated when the outlet pressure treatment changed (the increment is 0 on
the outlet and the mean-p anchor is gone); `pdrop`, `hydro` and `hydrowall` are new, and `rb` was regenerated when the wall pressure gradient started carrying the body force and buoyancy (the old reference grew spuriously, E 5e-5 -> 1.2e-4 in 20 steps, from a flow that should stay at rest). The np2 references were regenerated after the least-squares stencil started carrying
vertex-adjacent cells across partition cuts: step 1 of every case now agrees with the np1
reference to machine precision (before, `eps` in `rb` differed by 2.5% and `nu_sgs_mean` in
`channel` by 2e-3). The `dynsmag` references were regenerated again after the test filter
started including periodic partners on single-rank meshes.

The `pdrop` and `hydro` references were regenerated after the periodic seam stopped using the owner-only gradient in the deferred non-orthogonal flux (the np1 and np2 results now coincide); `inflow` is new. The `inflow` and `inflow_jarrin` references were regenerated when the analytic stress target changed from $-0.1u_*^2$ on $\overline{u'v'}$ to the momentum flux $\overline{u'w'}=-u_*^2$ (shared by both kernels).

To regenerate a reference after an intentional change in results:

```sh
DOPAMINE_BIN=build/dopamine tests/regression/run_case.sh channel /tmp/ref_channel 1
cp /tmp/ref_channel/case_stats.csv tests/regression/reference/channel_np1.csv
```
