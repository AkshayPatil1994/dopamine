---
title: Validation
nav_order: 7
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# Validation

This page collects the benchmark results used to validate DOPAMINE: order-of-accuracy verification against exact
solutions, and comparisons with published DNS and reference data. Every case lists its setup, the reference, the
figure, the numbers behind it and how to reproduce it. The case inputs live in `examples/<case>/`
([Examples](examples.md)); the regression suite (`ctest -L regression`) protects these results during refactoring
but is not a validation.

Unless noted, all runs use the default numerics of the version that produced them: RK3W time integration,
adaptive time step, `EXACT` cell projection, and the advection scheme named in the table. Quantities are
nondimensional.

## Summary

| Case | Physics | Reference | Grid | Result |
|---|---|---|---|---|
| [Advection order](#advection-order-of-accuracy) | 1D periodic advection | Exact solution | N = 20 to 320 | Observed order 2.0 (CDS, KEP, QUICK), 1.6 to 1.7 (TVD limiters); RK1/2 orders 1.0/2.0 |
| [Taylor-Green vortex](#taylor-green-vortex-re--1600) | Decaying vortex, Re = 1600 | 512^3 DNS (Dairay et al. 2017) | 64^3 to 384^3 | Peak -dE/dt error -1.3% at 384^3; E(t) relative L2 error 1.2% |
| [Lid-driven cavity](#lid-driven-cavity-re--1000) | Re = 1000 | Centreline reference profiles | 40^3, 60^3, 100^3 | RMS error 0.017 (u) and 0.011 (w) at 100^3 |
| [Channel, Re_tau = 395](#wall-modelled-channel-re_tau--395) | Wall-modelled LES | MKM DNS (1999) | 48 x 48 x 32 | Mean velocity within 3% for y+ > 30; stresses within 7% on average |
| [Half channel, Re_tau = 2000](#wall-modelled-half-channel-re_tau--2000) | Wall-modelled LES, AMD | Lee & Moser DNS (2015) | 64^3 | U+ error within 0.9 (4%) everywhere, about 0.7 for y/delta > 0.05 |
| [Michelstadt urban ABL](#michelstadt-urban-abl-in-progress) | Neutral urban LES | OpenFOAM, probes | 1.6 M cells | In progress, no comparison yet |

## About the figures

All figures share one manuscript style (serif text, Computer Modern mathematics, ticks on all sides pointing
inwards, a colour-blind-safe palette, open markers, lowercase panel labels), defined in
`examples/common/plotstyle.py`. Each figure is produced by a script in the matching example directory, reads only
the solver output (NetCDF or the statistics CSV) and the reference file, and also writes the reduced curves and
error metrics to a small JSON file next to the image in `docs/assets/validation/`. The numbers quoted below are
taken from those files. Each script also writes a PDF of the figure for use in manuscripts (PDFs are git-ignored and not committed). Animations are made
by `examples/common/animate_slices.py` and are illustrations of the flow, not evidence of accuracy.

## Advection order of accuracy

Spatial and temporal convergence of the advection schemes and Runge-Kutta integrators on a 1D periodic
problem with a known solution, run on uniform and smoothly stretched meshes (`tests/test_advection_convergence.f90`).
Error norms are taken against the exact solution; the temporal study uses QUICK on N = 400 and error relative
to a same-mesh RK4 reference.

![Advection convergence]({{ site.baseurl }}/assets/validation/advection/advection_validation.png)

*Figure 1. (a) Spatial convergence on a uniform and (b) a stretched mesh; solid lines are the L2 error, dashed
lines the L-infinity error, dotted lines show orders p = 1, 2, 3. (c) Temporal convergence of Forward Euler, Heun
RK2 and classical RK4 against the CFL number.*

| Scheme | Observed order, uniform | Observed order, stretched |
|---|---|---|
| CDS | 2.00 | 2.00 |
| KEP | 2.00 | 2.00 |
| QUICK | 2.01 | 2.01 |
| Koren (TVD) | 1.63 | 1.65 |
| van Leer (TVD) | 1.66 | 1.70 |

Orders are the L2 slope between the two finest meshes. The TVD limiters clip extrema, so their L-infinity error
does not converge (dashed curves in (a), (b)) while the L2 error converges at about 1.6. Temporal orders from the two largest time steps are 1.0 (Euler), 2.0 (Heun) and 3.5 (RK4; at the smallest
steps the RK4 error reaches the spatial error floor of the mesh).

Reproduce:

```sh
python3 examples/convergence_test/plot_validation.py --data <dir with spatial_convergence.csv> \
    --out docs/assets/validation/advection
```

## Taylor-Green vortex, Re = 1600

Decay of the Taylor-Green vortex in a periodic $$(2\pi)^3$$ box, $$\nu = 6.25\times10^{-4}$$, with `KEP` advection
and the default (Vreman) SGS model, on four grids. The reference is the 512^3 spectral-quality DNS of Dairay et al. (2017),
`TGV_Re1600.dat`, whose dissipation $$\epsilon=-\mathrm{d}E_k/\mathrm{d}t$$ is column 3 of the file.

The dissipation rate is shown two ways. The solver's `eps` statistic is gradient based,
$$\epsilon = 2(\nu+\nu_\mathrm{sgs})\langle S_{ij}S_{ij}\rangle$$, evaluated from the reconstructed cell gradients
(`src/stats/mod_stats.f90`). It measures the explicitly modelled viscous and SGS dissipation only. The total
rate of energy loss, $$-\mathrm{d}E_k/\mathrm{d}t$$, taken from the energy history, also contains the
dissipation of the numerical scheme, and is the quantity the DNS reports. A comparison of a coarse-grid solver
with the DNS should therefore use $$-\mathrm{d}E_k/\mathrm{d}t$$.

![Taylor-Green validation]({{ site.baseurl }}/assets/validation/taylor_green/tgv_validation.png)

*Figure 2. (a) Volume-averaged kinetic energy, (b) dissipation rate from the energy history
$$-\mathrm{d}E_k/\mathrm{d}t$$ (central differences of the per-step energy, averaged over 0.2 time units) and
(c) the gradient-based dissipation $$\epsilon$$ reported by the solver, for 64^3, 128^3, 256^3 and 384^3 cells
(symbols) and the 512^3 DNS (line).*

| Grid | Peak $$-\mathrm{d}E_k/\mathrm{d}t$$ ($$\times10^{3}$$) | Error vs DNS | Time of peak | Peak gradient-based $$\epsilon$$ ($$\times10^{3}$$) | Error vs DNS | E(t) relative L2 error |
|---|---|---|---|---|---|---|
| DNS 512^3 | 12.86 | | 8.98 | 12.85 | | |
| 384^3 | 12.69 | -1.3% | 9.02 | 12.05 | -6.2% | 1.2% |
| 256^3 | 11.77 | -8.5% | 9.01 | 10.72 | -16.6% | 3.2% |
| 128^3 | 13.29 | +3.4% | 8.55 | 10.70 | -16.7% | 6.0% |
| 64^3 | 9.93 | -22.7% | 8.54 | 4.99 | -61.1% | 19.6% |

The kinetic energy converges monotonically and tracks the DNS well from 256^3 on. The dissipation rate from the
energy history follows the DNS closely at 384^3 (peak -1.3%, correct peak time and decay), while the
gradient-based $$\epsilon$$ is 6% low at the same grid and 17% low at 128^3 and 256^3. The under-prediction in
panel (c) is therefore largely a property of how $$\epsilon$$ is computed rather than of the energy decay itself.
The gradient-based estimate leaves out the numerical dissipation of the scheme and is evaluated from gradients that
are under-resolved on coarse grids; this reading is consistent with the figure, but the split between the
contributions was not measured separately. The peak of $$-\mathrm{d}E_k/\mathrm{d}t$$ is not monotone in the grid
size: 128^3 overshoots the DNS peak by 3% and peaks about 0.4 time units early, and 256^3 is 8% low. Peak values
come from the smoothed energy-history derivative of each run and depend somewhat on the averaging window.

<img src="{{ site.baseurl }}/assets/validation/taylor_green/tgv_vorticity.gif" alt="Vorticity on two sections of the Taylor-Green vortex" style="max-width:100%">

*Animation. Out-of-plane vorticity on (a) the $$z=L/2$$ plane and (b) the $$y=0.375L$$ plane, 128^3 run, t = 1 to 20,
$$L=2\pi$$.*

Reproduce:

```sh
python3 examples/taylor_green/plot_validation.py --data <dir with 64/ 128/ 256/ 384/> \
    --dns TGV_Re1600.dat --out docs/assets/validation/taylor_green
```

## Lid-driven cavity, Re = 1000

Cubic cavity, lid velocity $$U_\mathrm{lid}=1$$, $$\nu=10^{-3}$$, `KEP` advection with the Vreman SGS model.
Profiles are time averages over the second half of each run (55, 97 and 166 snapshots for 40^3, 60^3 and 100^3),
taken on the mid-plane in the spanwise direction. The reference profiles are the tabulated centreline data
`u.dat`, `v.dat` shipped with the case.

![Cavity centreline profiles]({{ site.baseurl }}/assets/validation/cavity/cavity_validation.png)

*Figure 3. (a) Horizontal velocity along the vertical centreline $$x=0.5$$ and (b) vertical velocity along the
horizontal centreline $$z=0.5$$ for three grids (lines with markers) and the reference (filled circles).*

| Grid | RMS error u | RMS error w | max abs. error u | max abs. error w | min u | min w | max w |
|---|---|---|---|---|---|---|---|
| Reference | | | | | -0.285 | -0.441 | 0.246 |
| 40^3 | 0.038 | 0.020 | 0.153 | 0.041 | -0.258 | -0.408 | 0.227 |
| 60^3 | 0.019 | 0.011 | 0.074 | 0.028 | -0.273 | -0.427 | 0.240 |
| 100^3 | 0.017 | 0.011 | 0.048 | 0.024 | -0.275 | -0.436 | 0.235 |

Errors are the profile interpolated to the reference positions minus the reference value. Errors fall by a factor
of two from 40^3 to 60^3 and then flatten: the remaining difference at 100^3 is a small shift of the shear-layer
profile, around $$z\approx0.6$$, rather than an extremum error, and is not removed by refinement. The cause has
not been isolated.

<img src="{{ site.baseurl }}/assets/validation/cavity/cavity_vorticity.gif" alt="Vorticity on two sections of the lid-driven cavity" style="max-width:100%">

*Animation. Out-of-plane vorticity on (a) the $$y=L/2$$ plane and (b) the $$z=0.85L$$ plane during spin-up, 60^3 run.*

Reproduce:

```sh
python3 examples/lid_driven_cavity/plot_validation.py --nc lid_cavity_20.nc --nc lid_cavity_30.nc \
    --nc lid_cavity_50.nc --u-ref u.dat --v-ref v.dat --out docs/assets/validation/cavity
```

(the NetCDF files are named by the cells per edge divided by two in the case directory; the plot labels use the
true cell count.)

## Wall-modelled channel, Re_tau = 395

Turbulent channel flow driven by a streamwise body force ($$f_x=1$$, so $$u_\tau=1$$ by force balance),
$$h=1$$, $$\nu=1/395$$ in a $$6\times2\times2$$ periodic box with $$48\times48\times32$$ cells, `KEP` advection and
the wall model on both walls. The first cell centre is at $$y^+\approx8$$. Statistics average 500 snapshots over
$$t=121.6$$ to $$300$$ and the two channel halves. The reference is the Moser, Kim and Mansour (1999) DNS at
$$Re_\tau=392.24$$, which is tabulated in wall units and is not rescaled.

![Channel Re_tau 395 profiles]({{ site.baseurl }}/assets/validation/channel395/channel395_validation.png)

*Figure 4. (a) Mean velocity, (b to d) streamwise, wall-normal and spanwise RMS velocity, (e) shear stress budget
(resolved, SGS, viscous and total, against the exact total $$1-y/h$$) and (f) turbulent kinetic energy. Symbols
are DOPAMINE (16 wall-normal levels), lines are DNS.*

| Quantity | Mean abs. relative error, $$y^+>30$$ |
|---|---|
| Mean velocity | 1.7% (max 3.1%) |
| $$u'_\mathrm{rms}$$ | 5.0% |
| $$v'_\mathrm{rms}$$ | 6.7% |
| $$w'_\mathrm{rms}$$ | 4.8% |
| $$-\langle u'v'\rangle$$ | 4.1% |

The mean profile follows the log layer and the DNS to within 3%. The wall-normal stress is under-predicted over the
whole channel (about 10% below the DNS maximum), which is typical of coarse wall-modelled LES. The
first wall-adjacent point lies inside the wall-modelled region and is not meant to match the DNS; it is excluded
from the error measures. The total stress follows $$1-y/h$$, consistent with stationary statistics and
$$u_\tau=1$$.

<img src="{{ site.baseurl }}/assets/validation/channel395/channel395_streaks.gif" alt="Streamwise velocity fluctuations at y+ = 25" style="max-width:100%">

*Animation. Streamwise velocity fluctuation $$u'/u_\tau$$ on the wall-parallel plane of the second cell layer
($$y^+\approx25$$), t = 200 to 300; the streak spacing is visible.*

Reproduce:

```sh
python3 examples/channel_flow_retau395/plot_validation.py --les les_validation_mkm.csv \
    --means chan395.means --stress chan395.reystress --out docs/assets/validation/channel395
```

## Wall-modelled half channel, Re_tau = 2000

Half channel of size $$2\pi\times\pi\times1$$ (wall at $$z=0$$, free-slip centreline at $$z=\delta$$), $$64^3$$
uniform cells ($$\Delta x^+=196$$, first cell at $$y^+\approx16$$), $$\nu=1/2000$$, `AMD` subgrid model with the
adaptive wall-model input filter (`wall_model_filter_tf_fac = 1`) and the Reichardt log-law wall model. The run
covers $$t\,u_\tau/\delta = 0$$ to 116 and the statistics are averaged over $$t=80$$ to 116. The reference is the
Lee & Moser (2015) DNS at $$Re_\tau=1994.756$$.

![Channel Re_tau 2000 mean profile]({{ site.baseurl }}/assets/validation/channel2000/channel2000_profiles.png)

*Figure 5. (a) Mean velocity against $$y/\delta$$ with the DNS and the log law ($$\kappa=0.41$$); (b) difference
from the DNS profile. Markers show a subset of the 64 wall-normal cells.*

| $$y/\delta$$ | $$\langle U\rangle^+ - U^+_\mathrm{DNS}$$ |
|---|---|
| 0.1 | +0.68 |
| 0.5 | -0.20 |
| 0.99 | -0.68 |

The difference stays within about $$\pm0.7$$ wall units (3%) for $$y/\delta>0.05$$, including the centreline,
and reaches $$-0.9$$ at $$y/\delta\approx0.02$$, in the second cell. The profile crosses the DNS near
$$y/\delta\approx0.04$$ and 0.4, the usual log-layer mismatch of a first cell inside the buffer layer.

![Statistical convergence of the Re_tau 2000 run]({{ site.baseurl }}/assets/validation/channel2000/channel2000_convergence.png)

*Figure 6. Convergence evidence. (a) Bulk velocity (thin) and its running mean (thick) against the DNS value; (b) velocity at
$$y/\delta=0.1,\,0.5,\,0.99$$ against time; (c) change of the running-mean profile; (d) probe means against the
averaging start $$t_0$$.*

Panel (d) shows the probe means are insensitive to the averaging start for $$t_0\gtrsim50$$, so the averaging
window $$t=80$$ to 116 is adequate for the mean profile; the run is not long enough to converge higher-order
statistics, which are not claimed here.

Reproduce (the NetCDF output is 4.5 GB; the figures are made with):

```sh
python3 examples/channelAMDvalidation/postprocess.py AMD=run_amd_filtered/channel.nc --t0 80 --t1 116 \
    --dns data/hoyas_jimenez_retau2000.csv --plot profiles.pdf
python3 examples/channelAMDvalidation/plot_convergence.py AMD=run_amd_filtered/channel.nc \
    --dns data/hoyas_jimenez_retau2000.csv --plot convergence.pdf
```

The DNS file `data/hoyas_jimenez_retau2000.csv` currently holds the Lee & Moser profile (see its header).

## Michelstadt urban ABL (in progress)

*This validation is under way. The figures below are snapshots of a run that has not finished, shown to document
the setup and the solver behaviour on a realistic urban geometry. No error measures are claimed and the numbers
will change.*

Neutral atmospheric boundary layer over the Michelstadt urban geometry (buildings meshed with snappyHexMesh, 1.6 M cells, domain $$2400\times1080\times144$$ m). Wind is westerly with
$$U_\mathrm{ref}=5$$ m/s at $$z_\mathrm{ref}=50$$ m, rough-wall Monin-Obukhov wall model ($$z_0=0.01$$ m),
`AMD` subgrid model and `KEP` advection. The inflow is a precomputed set of fluctuation planes
(`dopamine-inflow-gen`; the online SEM is not run) and the upstream fetch is lined with an array of
roughness cubes (small grey squares at $$x<0$$). The case is compared with an OpenFOAM simulation on the same
mesh and with line probes P1 to P8 across the canopy, to be added.

![Michelstadt pedestrian-level velocity magnitude]({{ site.baseurl }}/assets/validation/michelstadt/michelstadt_snapshots.png)

*Figure 7. Velocity magnitude at pedestrian level ($$z\approx1.5$$ m) at (a) t = 400 s and (b) t = 1560 s. Grey
regions contain no cells (solid buildings and roughness cubes). Flow is from left to right.*

<img src="{{ site.baseurl }}/assets/validation/michelstadt/michelstadt_plan.gif" alt="Michelstadt plan view animation" style="max-width:100%">

*Animation. Pedestrian-level velocity magnitude, t = 20 to 1560 s.*

Still to do: finish the 2400 s run, add the OpenFOAM and line-probe comparison (mean and RMS profiles), and report
the time-averaged statistics.

## Known limitations

- The gradient-based dissipation statistic `eps` under-predicts the DNS dissipation on all grids tested because it
  excludes numerical dissipation; use $$-\mathrm{d}E_k/\mathrm{d}t$$ for comparisons with DNS. Refinement beyond
  384^3 was not tested, and the split between numerical and SGS dissipation was not measured.
- The cavity reference file is a tabulated profile without a recorded citation; it is used as a visual and RMS
  benchmark only.
- Channel statistics at Re_tau = 395 are averaged over a finite window and one realisation; the wall-normal stress
  deficit is a property of the coarse wall-modelled grid, not a defect of the averaging.
- Rayleigh-Benard convection (examples/RayleighBernard) and the SEM inflow are exercised by the regression suite
  but have no completed validation yet.
