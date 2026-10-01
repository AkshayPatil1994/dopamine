#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
plot_convergence.py
===================
Plot the spatial and temporal convergence figures from run_convergence.py's results.json
and print the observed orders.

    plot_convergence.py results/results.json --out ../../docs/assets/convergence

Writes spatial_convergence.png and temporal_convergence.png into --out.

Metrics
-------
spatial   self-convergence of E(t_end): |E_N - E_2N| against h = L/N (no exact solution is
          needed), and the true error of the initial dissipation rate against the analytic
          eps(0) = 0.75 nu.
temporal  self-convergence |E(dt) - E(dt/2)| against dt on three meshes, and the inviscid
          energy loss (the exact answer is E = 1/8 for nu -> 0).
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito picks: distinguishable under the common colour-vision deficiencies. Series
# are also told apart by marker and line style, never by colour alone.
COL = {"hex": "#0072B2", "warped": "#D55E00", "16": "#0072B2", "32": "#009E73", "64": "#CC79A7"}
MARK = {"hex": "o", "warped": "s", "16": "o", "32": "s", "64": "^"}
INK, GRID = "#222222", "#d9d9d9"


def style(ax, xlabel, ylabel, title):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def slope_line(ax, x, y, p, label, x1, x2, ls="--"):
    """Reference line of slope p passing through (x, y), drawn from x1 to x2."""
    xs = np.array([x1, x2])
    ax.plot(xs, y * (xs / x) ** p, ls, color="#666666", lw=1.1, label=label)


def orders(vals, ratio=2.0):
    return [math.log(vals[i] / vals[i + 1]) / math.log(ratio) for i in range(len(vals) - 1)]


def fmt(o):
    return ", ".join(f"{x:.2f}" for x in o)


def spatial(res, out):
    sp = res["spatial"]
    meta = sp["meta"]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)

    ax = axs[0]
    style(ax, "mesh spacing  h = L/N", r"$|E_N - E_{2N}|$",
          "(a) kinetic energy at t = %g, successive meshes" % meta["t_end"])
    print("SPATIAL  E(t_end) self-convergence, observed order per refinement:")
    for fam in ("hex", "warped"):
        rows = sp[fam]
        h = [r["h"] for r in rows]
        d = [abs(rows[i]["E"] - rows[i + 1]["E"]) for i in range(len(rows) - 1)]
        ax.plot(h[:-1], d, "-", marker=MARK[fam], color=COL[fam], ms=6,
                label=f"{fam}  (orders {fmt(orders(d))})")
        print(f"  {fam:7s} diffs {['%.3e' % x for x in d]} orders {fmt(orders(d))}")
    h0 = sp["hex"][0]["h"]
    d0 = abs(sp["hex"][0]["E"] - sp["hex"][1]["E"])
    slope_line(ax, h0, d0 * 1.6, 2, "slope 2", sp["hex"][2]["h"] * 0.8, h0 * 1.15)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    ax = axs[1]
    style(ax, "mesh spacing  h = L/N", r"$|\varepsilon_0^{num} - 0.75\,\nu|$",
          r"(b) initial dissipation rate vs. analytic value")
    exact = meta["eps0_exact"]
    print("SPATIAL  eps(0) true error vs 0.75*nu, observed order per refinement:")
    for fam in ("hex", "warped"):
        rows = sp[fam]
        h = [r["h"] for r in rows]
        e = [abs(r["eps0"] - exact) for r in rows]
        ax.plot(h, e, "-", marker=MARK[fam], color=COL[fam], ms=6,
                label=f"{fam}  (orders {fmt(orders(e))})")
        print(f"  {fam:7s} errs {['%.3e' % x for x in e]} orders {fmt(orders(e))}")
    e0 = abs(sp["hex"][0]["eps0"] - exact)
    slope_line(ax, sp["hex"][0]["h"], e0 * 1.6, 2, "slope 2",
               sp["hex"][-1]["h"] * 0.8, sp["hex"][0]["h"] * 1.15)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    fig.suptitle(r"Spatial convergence, 3D Taylor-Green vortex, $\nu$ = %g, central advection + RK3W"
                 % meta["nu"], x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.savefig(out / "spatial_convergence.png", dpi=170)
    plt.close(fig)


def temporal(res, out):
    tp = res["temporal"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)

    ax = axs[0]
    style(ax, r"time step  $\Delta t$", r"$|E(\Delta t) - E(\Delta t/2)|$",
          "(a) kinetic energy at t = %g, successive dt" % tp["t_end"])
    print("TEMPORAL  E(t_end) self-convergence, observed order per halving of dt:")
    first = None
    for n, rows in tp["meshes"].items():
        dt = [r["dt"] for r in rows]
        d = [abs(rows[i]["E"] - rows[i + 1]["E"]) for i in range(len(rows) - 1)]
        ax.plot(dt[:-1], d, "-", marker=MARK[n], color=COL[n], ms=6,
                label=rf"${n}^3$ cells")
        print(f"  N={n:>3s} dts {dt[:-1]}\n         diffs {['%.2e' % x for x in d]}\n"
              f"         orders {fmt(orders(d))}")
        if n == "64":
            first = (dt, d)
    dt, d = first
    slope_line(ax, dt[1], d[1] * 0.35, 3, "slope 3 (RK3W)", dt[1] * 0.45, dt[1] * 3)
    slope_line(ax, dt[-2], d[-2] * 0.35, 1, r"slope 1  ($\propto \Delta t\,h^2$)",
               dt[-2] * 0.3, dt[-2] * 6, ls=":")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    inv = tp["inviscid"]
    ax = axs[1]
    style(ax, r"time step  $\Delta t$", r"energy loss  $1/8 - E(t=%g)$" % inv["t_end"],
          r"(b) inviscid run ($\nu \to 0$): spurious energy loss")
    coef = {}
    for n in tp["meshes"]:
        rr = [r for r in inv["runs"] if str(r["N"]) == n]
        dt = [r["dt"] for r in rr]
        loss = [0.125 - r["E"] for r in rr]
        ax.plot(dt, loss, "-", marker=MARK[n], color=COL[n], ms=6, label=rf"${n}^3$ cells")
        coef[n] = loss[-1] / dt[-1]
        print(f"  inviscid N={n}: loss {['%.3e' % x for x in loss]}  orders {fmt(orders(loss))}")
    ref = [r for r in inv["runs"] if r["N"] == 64]
    slope_line(ax, ref[1]["dt"], (0.125 - ref[1]["E"]) * 3, 1, "slope 1",
               ref[-1]["dt"] * 0.8, ref[0]["dt"] * 1.3)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    ax = axs[2]
    style(ax, "mesh spacing  h = L/N", r"loss$/\Delta t$  (coefficient of the $O(\Delta t)$ term)",
          r"(c) the $O(\Delta t)$ coefficient shrinks with $h$")
    hs = [2 * math.pi / int(n) for n in coef]
    cs = list(coef.values())
    ax.plot(hs, cs, "-", marker="o", color="#222222", ms=6, label="measured")
    slope_line(ax, hs[1], cs[1], 2, "slope 2", hs[-1] * 0.8, hs[0] * 1.2)
    print(f"  O(dt) coefficient vs h: {['%.3e' % c for c in cs]}  orders in h {fmt(orders(cs))}")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    fig.suptitle(r"Temporal convergence, 3D Taylor-Green vortex, $\nu$ = %g, central advection + RK3W"
                 % tp["nu"], x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.savefig(out / "temporal_convergence.png", dpi=170)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "figures"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res = json.loads(Path(args.results).read_text())
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": "#888888",
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    if "spatial" in res:
        spatial(res, out)
    if "temporal" in res:
        temporal(res, out)
    print(f"figures in {out}")


if __name__ == "__main__":
    main()
