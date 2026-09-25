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

To regenerate a reference after an intentional change in results:

```sh
DOPAMINE_BIN=build/dopamine tests/regression/run_case.sh channel /tmp/ref_channel 1
cp /tmp/ref_channel/case_stats.csv tests/regression/reference/channel_np1.csv
```
