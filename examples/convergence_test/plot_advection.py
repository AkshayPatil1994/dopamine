#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
plot_advection.py
=================
Figures for the scalar advection scheme comparison from run_advection.py's advection.json:

  advection_schemes.png      final profile of every scheme against the exact one
  advection_convergence.png  error of the smooth Gaussian against the mesh spacing

    plot_advection.py results/advection.json --out ../../docs/assets/convergence
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, GRID, EXACT = "#222222", "#d9d9d9", "#888888"
COL = {"CDS": "#0072B2", "QUICK": "#D55E00", "KOREN_TVD": "#009E73", "VANLEER_TVD": "#CC79A7"}
MARK = {"CDS": "o", "QUICK": "s", "KOREN_TVD": "^", "VANLEER_TVD": "D"}
LABEL = {"CDS": "CDS (central)", "QUICK": "QUICK", "KOREN_TVD": "Koren TVD",
         "VANLEER_TVD": "van Leer TVD"}


def despine(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def profiles(res, out):
    pr = res["profiles"]
    x = np.array(pr["x"])
    ex = np.array(pr["exact"])
    fig, axs = plt.subplots(2, 2, figsize=(10, 6.2), sharex=True, sharey=True,
                            constrained_layout=True)
    print("PROFILES  N =", pr["N"], "(min, max of the final profile; exact is 0 .. 2)")
    for ax, (name, prof) in zip(axs.flat, pr["schemes"].items()):
        prof = np.array(prof)
        ax.step(x, ex, where="mid", color=EXACT, lw=2.4, label="exact", zorder=1)
        ax.plot(x, prof, "-", marker=MARK[name], color=COL[name], ms=3.2, lw=1.0,
                label=LABEL[name], zorder=2)
        ax.set_title(LABEL[name], loc="left", fontsize=10.5, color=INK)
        ax.grid(True, color=GRID, lw=0.6)
        despine(ax)
        ax.text(0.99, 0.97, f"min {prof.min():+.3f}   max {prof.max():.3f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, color=INK)
        print(f"  {name:12s} min {prof.min():+.4f} max {prof.max():.4f}")
    for ax in axs[1]:
        ax.set_xlabel("x")
    for ax in axs[:, 0]:
        ax.set_ylabel("scalar")
    handles = [plt.Line2D([], [], color=EXACT, lw=2.4, label="exact (initial profile)")]
    axs[0, 0].legend(handles=handles, frameon=False, fontsize=8.5, loc="upper left")
    fig.suptitle("Scalar advection schemes: one revolution of a Gaussian and a square pulse, "
                 "%d cells, CFL 0.2" % pr["N"], x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.savefig(out / "advection_schemes.png", dpi=170)
    plt.close(fig)


def convergence(res, out):
    cv = res["convergence"]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    print("CONVERGENCE  observed orders per refinement (N doubles):")
    for ax, key, title in zip(axs, ("L2", "Linf"), (r"(a) $L_2$ error", r"(b) $L_\infty$ error")):
        for name, rows in cv.items():
            h = [r["h"] for r in rows]
            e = [r[key] for r in rows]
            o = [math.log(e[i] / e[i + 1]) / math.log(h[i] / h[i + 1]) for i in range(len(e) - 1)]
            ax.loglog(h, e, "-", marker=MARK[name], color=COL[name], ms=6,
                      label=f"{LABEL[name]}  (last order {o[-1]:.2f})")
            print(f"  {key:4s} {name:12s} orders " + ", ".join(f"{v:.2f}" for v in o))
        ref = cv["CDS"][2][key]
        h2 = cv["CDS"][2]["h"]
        hs = np.array([cv["CDS"][-1]["h"] * 0.85, cv["CDS"][0]["h"] * 1.15])
        for p, ls in ((1, ":"), (2, "--"), (3, "-.")):
            ax.loglog(hs, ref * 0.12 * (hs / h2) ** p, ls, color="#777777", lw=1.0,
                      label=f"slope {p}")
        ax.set_xlabel("mesh spacing  h")
        ax.set_ylabel("error at t = 1")
        ax.set_title(title, loc="left", fontsize=10.5, color=INK)
        ax.grid(True, which="major", color=GRID, lw=0.6)
        despine(ax)
    axs[0].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Scalar advection of a smooth Gaussian, one revolution, CFL 0.2, RK3W",
                 x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.savefig(out / "advection_convergence.png", dpi=170)
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
    profiles(res, out)
    convergence(res, out)
    print(f"figures in {out}")


if __name__ == "__main__":
    main()
