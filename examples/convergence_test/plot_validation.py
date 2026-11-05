#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Advection verification in the shared manuscript style: spatial and temporal order of accuracy.

Reads spatial_convergence.csv (mesh,scheme,N,L2,Linf) and temporal_convergence.csv (integrator,cfl,L2) written
by tests/test_advection_convergence.f90, writes advection_validation.{png,pdf} and advection_validation.json
(observed spatial order from the two finest meshes, temporal order from the two smallest CFL).

  plot_validation.py --data /path/to/advectionTest --out docs/assets/validation/advection
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402

SCHEMES = ["CDS", "KEP", "QUICK", "KOREN_TVD", "VANLEER_TVD"]
STYLE = {"CDS": ("#0072B2", "o"), "KEP": ("#56B4E9", "v"), "QUICK": ("#D55E00", "^"),
         "KOREN_TVD": ("#CC79A7", "D"), "VANLEER_TVD": ("#009E73", "s")}
LABEL = {"KOREN_TVD": "Koren", "VANLEER_TVD": "van Leer"}
INTEGRATORS = [("Forward Euler", "#E69F00", "o", 1), ("Heun RK2", "#0072B2", "s", 2),
               ("Classical RK4", "#000000", "D", 4)]


def order(x, y):
    return float(-np.polyfit(np.log(x[-2:]), np.log(y[-2:]), 1)[0])


def read(path):
    with open(path) as fh:
        return [{k: v.strip() for k, v in r.items()} for r in csv.DictReader(fh)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    sp = defaultdict(lambda: defaultdict(list))
    for r in read(os.path.join(args.data, "spatial_convergence.csv")):
        sp[r["mesh"]][r["scheme"]].append((float(r["N"]), float(r["L2"]), float(r["Linf"])))
    tp = defaultdict(list)
    for r in read(os.path.join(args.data, "temporal_convergence.csv")):
        tp[r["integrator"]].append((float(r["cfl"]), float(r["L2"])))

    plotstyle.apply()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), constrained_layout=True)
    res = {"spatial": {}, "temporal": {}}
    for ax, mesh, tag in zip(axes[:2], ("uniform mesh", "stretched mesh"), "ab"):
        res["spatial"][mesh] = {}
        for s in SCHEMES:
            d = np.array(sorted(sp[mesh][s]))
            c, m = STYLE[s]
            ax.loglog(d[:, 0], d[:, 1], color=c, marker=m, label=LABEL.get(s, s))
            ax.loglog(d[:, 0], d[:, 2], color=c, marker=m, ls="--", lw=0.9, alpha=0.6, markersize=3.5)
            res["spatial"][mesh][s] = {"N": d[:, 0].tolist(), "L2": d[:, 1].tolist(), "Linf": d[:, 2].tolist(),
                                       "order_L2": order(d[:, 0], d[:, 1])}
        n0 = np.array([80.0, 320.0])
        for p in (1, 2, 3):
            ax.loglog(n0, 1.5e-5 * (n0 / 320.0) ** (-p), color=plotstyle.GREY, ls=":", lw=1.0)
            ax.text(n0[0] * 0.93, 1.5e-5 * 0.25 ** (-p), rf"$p={p}$", fontsize=7.5, color=plotstyle.GREY,
                    va="center", ha="right")
        ax.set_xlabel(r"$N$")
        ax.set_ylabel(r"error (solid $L_2$, dashed $L_\infty$)")
        ax.set_xlim(None, 400)
        plotstyle.panel_label(ax, f"({tag})")
    axes[0].legend(loc="upper right", ncol=1)

    ax = axes[2]
    for name, c, m, p in INTEGRATORS:
        d = np.array(sorted(tp[name]))
        ax.loglog(d[:, 0], d[:, 1], color=c, marker=m, label=rf"{name} ($p={p}$)")
        xs = np.array([d[0, 0], d[-1, 0]])
        ax.loglog(xs, 0.3 * d[-1, 1] * (xs / d[-1, 0]) ** p, color=c, ls=":", lw=0.9)
        res["temporal"][name] = {"cfl": d[:, 0].tolist(), "L2": d[:, 1].tolist(), "order": float(np.polyfit(np.log(d[:2, 0]), np.log(d[:2, 1]), 1)[0])}
    ax.set_xlabel(r"CFL$\,\propto\Delta t$")
    ax.set_ylabel(r"$L_2$ error")
    ax.legend(loc="lower right")
    plotstyle.panel_label(ax, "(c)")

    os.makedirs(args.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(args.out, "advection_validation.pdf"))
    with open(os.path.join(args.out, "advection_validation.json"), "w") as fh:
        json.dump(res, fh, separators=(",", ":"))
    for mesh, d in res["spatial"].items():
        print(mesh, {s: round(v["order_L2"], 2) for s, v in d.items()})
    print("temporal", {k: round(v["order"], 2) for k, v in res["temporal"].items()})


if __name__ == "__main__":
    main()
