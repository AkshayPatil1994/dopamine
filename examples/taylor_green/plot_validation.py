#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Taylor-Green vortex Re = 1600: kinetic energy and dissipation against the 512^3 spectral-quality DNS.

Reads <data>/<N>/tgv_stats.csv for each grid (columns time, E, eps) and the Dairay et al. (2017) DNS file
TGV_Re1600.dat. The dissipation is shown two ways: -dE/dt from the energy history (total energy loss, smoothed
over --window time units) and the solver's gradient-based eps = 2 (nu + nu_sgs) <S_ij S_ij>. Writes tgv_validation.{png,pdf} and tgv_validation.json (reduced curves plus error metrics).

  plot_validation.py --data /path/to/taylor_green --dns /path/to/TGV_Re1600.dat --out docs/assets/validation/taylor_green
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402

GRIDS = (64, 128, 256, 384)
STYLE = {64: ("#E69F00", "d"), 128: ("#D55E00", "o"), 256: ("#0072B2", "^"), 384: ("#009E73", "s")}
EPS_SCALE = 1.0e3


def load_run(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["time"], d["E"], d["eps"]


def minus_dEdt(t, E, window):
    """-dE/dt by central differences, averaged over a moving window of `window` time units."""
    g = -np.gradient(E, t)
    k = max(1, int(round(window / np.median(np.diff(t)))))
    pad = np.pad(g, (k // 2, k - 1 - k // 2), mode="edge")
    return np.convolve(pad, np.ones(k) / k, mode="valid")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="directory holding 64/, 128/, 256/, 384/")
    ap.add_argument("--dns", required=True, help="TGV_Re1600.dat")
    ap.add_argument("--out", default=".")
    ap.add_argument("--sample", type=float, default=0.5, help="marker spacing in time units")
    ap.add_argument("--window", type=float, default=0.2, help="smoothing window for -dE/dt in time units")
    args = ap.parse_args()

    dns = np.loadtxt(args.dns, comments="#")
    td, Ed, epsd = dns[:, 0], dns[:, 1], dns[:, 3]
    ip = epsd.argmax()

    plotstyle.apply()
    import matplotlib.pyplot as plt

    fig, (a, b, c) = plt.subplots(3, 1, figsize=(5.2, 8.6), sharex=True, constrained_layout=True)
    a.plot(td, Ed, "k-", lw=2.0, label=r"DNS $512^3$")
    b.plot(td, EPS_SCALE * dns[:, 2], "k-", lw=2.0)
    c.plot(td, EPS_SCALE * epsd, "k-", lw=2.0)

    out = {"dns_peak": {"t": float(td[ip]), "eps": float(epsd[ip])}, "runs": {}}
    for n in GRIDS:
        t, E, eps = load_run(os.path.join(args.data, str(n), "tgv_stats.csv"))
        col, m = STYLE[n]
        ts = np.arange(0.0, t[-1] + 1e-9, args.sample)
        Es, epss = np.interp(ts, t, E), np.interp(ts, t, eps)
        kw = dict(color=col, marker=m, ls="none", markersize=4.5)
        a.plot(ts, Es, label=rf"${n}^3$", **kw)
        c.plot(ts, EPS_SCALE * epss, **kw)
        dE = minus_dEdt(t, E, args.window)
        dEs = np.interp(ts, t, dE)
        b.plot(ts, EPS_SCALE * dEs, **kw)
        ref_E = np.interp(t, td, Ed)
        j = eps.argmax()
        jd = dE.argmax()
        out["runs"][str(n)] = {
            "t": ts.tolist(), "E": Es.tolist(), "eps": epss.tolist(), "minus_dEdt": dEs.tolist(),
            "E_rel_L2": float(np.sqrt(np.mean(((E - ref_E) / ref_E) ** 2))),
            "eps_peak": float(eps[j]), "t_peak": float(t[j]),
            "eps_peak_err_pct": float(100 * (eps[j] / epsd[ip] - 1)),
            "t_peak_shift": float(t[j] - td[ip]),
            "dEdt_peak": float(dE[jd]), "dEdt_t_peak": float(t[jd]),
            "dEdt_peak_err_pct": float(100 * (dE[jd] / dns[ip, 2] - 1)),
            "dEdt_t_peak_shift": float(t[jd] - td[ip]),
            "dEdt_rel_L2": float(np.sqrt(np.mean(((dE - np.interp(t, td, dns[:, 2])) / dns[:, 2].max()) ** 2))),
        }
    out["dns"] = {"t": td[::5].tolist(), "E": Ed[::5].tolist(), "eps": epsd[::5].tolist()}

    a.set_ylabel(r"$E_k$")
    b.set_ylabel(r"$-\mathrm{d}E_k/\mathrm{d}t\times10^{3}$")
    c.set_ylabel(r"$\epsilon\times10^{3}$")
    c.set_xlabel(r"$t\,U_0/L$")
    c.set_xlim(0, 20)
    a.legend(loc="upper right", ncol=1)
    plotstyle.panel_label(a, "(a)")
    plotstyle.panel_label(b, "(b)")
    plotstyle.panel_label(c, "(c)")
    os.makedirs(args.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(args.out, "tgv_validation.pdf"))
    with open(os.path.join(args.out, "tgv_validation.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))
    for n, r in out["runs"].items():
        print(f"{n}^3  E_relL2 = {r['E_rel_L2']:.3e}  eps_peak err = {r['eps_peak_err_pct']:+.1f}%  "
              f"t_peak shift = {r['t_peak_shift']:+.2f}  | -dE/dt peak err = {r['dEdt_peak_err_pct']:+.1f}%  "
              f"t shift = {r['dEdt_t_peak_shift']:+.2f}  L2/max = {r['dEdt_rel_L2']:.3f}")


if __name__ == "__main__":
    main()
