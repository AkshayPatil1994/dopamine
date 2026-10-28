# channelAMDvalidation

Validation of the AMD subgrid model (`sgs_model = 'AMD'`) and the adaptive wall-model input filter
(`wall_model_filter_tf_fac`) in a wall-modelled half channel at Re_tau = 2000, following the finite-difference
DOPAMINE code (`dopamine-fdm`, branch `dev`, Numerics section on the AMD channel validation).

| Item | Value |
|---|---|
| Domain | 2 pi x pi x 1 (x streamwise, y spanwise, z wall-normal), wall at z = 0, free-slip centreline at z = 1 |
| Grid | 64 x 64 x 64, uniform; dx+ = 196, dy+ = 98, first cell centre at y+ ~ 16 |
| Physics | nu = 1/2000, fx = 1 (so u_tau = 1, delta = 1) |
| Wall model | Reichardt log law, sampled in the wall-adjacent cell |
| Averaging | t = 14 .. 24 (delta/u_tau) |
| Reference | Hoyas & Jimenez Re_tau = 2000 DNS (supply it, see `data/README.md`) |

| Input | SGS | Wall-model input |
|---|---|---|
| `amd_filtered.input` | AMD | filtered, T_f = 1 y/(kappa u_tau) |
| `amd_raw.input` | AMD | raw (instantaneous) |
| `vreman_filtered.input` | Vreman | filtered |

## Run

```sh
gmsh channel_amd.geo -3 -o channel_amd.msh
mkdir run_amd_filtered && cd run_amd_filtered
OMP_NUM_THREADS=8 ../../../build/dopamine ../amd_filtered.input        # or decompose + mpirun, see
cd ..                                                                  # ../channel_flow_retau395/tutorial.md
# repeat for amd_raw.input and vreman_filtered.input
python3 postprocess.py amd_filtered=run_amd_filtered/channel.nc amd_raw=run_amd_raw/channel.nc \
    vreman_filtered=run_vreman_filtered/channel.nc --t0 14 --t1 24 \
    --dns data/hoyas_jimenez_retau2000.csv --plot profiles.png
```

## What to look for

The finite-difference code reports `U - U_DNS` at y/delta = 0.1, 0.5, 0.99 of
+2.99 / +2.68 / +2.26 (AMD, raw) against +0.27 / +0.02 / -0.39 (AMD, filtered). The finite-volume
discretisation will not reproduce those numbers, so the acceptance criteria are about the trend:

1. the filtered run is clearly closer to DNS than the raw run (the filter removes the log-layer mismatch);
2. the filtered AMD bulk velocity is within a few percent of the DNS value;
3. `nu_sgs` is positive and bounded, including in the wall-adjacent cell (`amd_wall_gradient = .true.`).

## Differences from the finite-difference reference

* The FD case uses a rough-wall log law with z0 = 6.77e-5 delta. This code only has the smooth Reichardt law.
  At Re_tau = 2000, z0 = 0.135 nu/u_tau, so the two should be close but are not identical.
* The FD log-layer-mismatch strength metric is not reproduced; only `U - U_DNS` and the bulk error are reported.
* C_amd is not tuned for coarse grids (1/3, the Poincare constant of second-order central differences).
