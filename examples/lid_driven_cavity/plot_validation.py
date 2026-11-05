#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Lid-driven cavity centreline profiles on several grids against reference data, in the manuscript style.

Reads one NetCDF output per grid (coords, offsets, connectivity, time, U), time-averages the last --tavg
fraction of the snapshots, extracts u(z) at x = 0.5 and w(x) at z = 0.5 from the two cell planes straddling the
centre (restricted to the two y-planes straddling the span centre), and plots them against the reference
u.dat (columns u, y) and v.dat (columns x, v). Writes cavity_validation.{png,pdf} and cavity_validation.json.

  plot_validation.py --nc 20=lid_cavity_20.nc --nc 30=lid_cavity_30.nc --u-ref u.dat --v-ref v.dat \\
      --out docs/assets/validation/cavity
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402

COLORS = ["#E69F00", "#D55E00", "#0072B2", "#009E73"]
MARKERS = ["d", "o", "^", "s"]


def planes(c, target=0.5):
    u = np.unique(np.round(c, 6))
    i = int(np.argmin(np.abs(u - target)))
    j = min(max(i + 1 if u[i] < target else i - 1, 0), len(u) - 1)
    return u[[i, j]]


def read_case(path, tavg, chunk=10):
    import netCDF4

    nc = netCDF4.Dataset(path)
    coords, off, conn = nc["coords"][:], nc["offsets"][:], nc["connectivity"][:]
    n = np.diff(off)
    start = off[:-1]
    c = np.stack([np.add.reduceat(coords[conn, k], start) / n for k in range(3)], axis=1)
    nt = nc["U"].shape[0]
    i0 = int(nt * (1 - tavg))
    acc, cnt = 0.0, 0
    for i in range(i0, nt, chunk):
        blk = nc["U"][i:min(i + chunk, nt)]
        acc = acc + blk.sum(axis=0)
        cnt += blk.shape[0]
    U = acc / cnt
    nc.close()
    ny = np.isin(np.round(c[:, 1], 6), planes(c[:, 1]))
    xs, zs = planes(c[:, 0]), planes(c[:, 2])
    mu = ny & np.isin(np.round(c[:, 0], 6), xs)
    mw = ny & np.isin(np.round(c[:, 2], 6), zs)

    def line(mask, ax, comp):
        pos = np.round(c[mask, ax], 6)
        up = np.unique(pos)
        return up, np.array([U[mask, comp][pos == p].mean() for p in up])

    z, u = line(mu, 2, 0)
    x, w = line(mw, 0, 2)
    nx = len(np.unique(np.round(c[:, 0], 6)))
    return dict(nx=nx, snapshots=cnt, z=z, u=u, x=x, w=w)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nc", action="append", required=True, help="path to a NetCDF file (repeatable)")
    ap.add_argument("--u-ref", required=True)
    ap.add_argument("--v-ref", required=True)
    ap.add_argument("--tavg", type=float, default=0.5)
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    ur = np.loadtxt(args.u_ref, delimiter=",", skiprows=1)
    vr = np.loadtxt(args.v_ref, delimiter=",", skiprows=1)
    cases = sorted((read_case(p, args.tavg) for p in args.nc), key=lambda d: d["nx"])

    plotstyle.apply()
    import matplotlib.pyplot as plt

    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True)
    a.plot(ur[:, 0], ur[:, 1], "ko", ms=3.5, mfc="k", label="Reference")
    b.plot(vr[:, 0], vr[:, 1], "ko", ms=3.5, mfc="k", label="Reference")
    out = {"ref": {"u": ur.tolist(), "w": vr.tolist()}, "cases": {}}
    for k, d in enumerate(cases):
        me = max(1, len(d["z"]) // 14)
        kw = dict(color=COLORS[k % 4], marker=MARKERS[k % 4], ls="-", lw=1.2, markersize=4.5, markevery=me)
        lab = rf"${d['nx']}^3$"
        a.plot(d["u"], d["z"], label=lab, **kw)
        b.plot(d["x"], d["w"], label=lab, **kw)
        eu = np.interp(ur[:, 1], d["z"], d["u"]) - ur[:, 0]
        ew = np.interp(vr[:, 0], d["x"], d["w"]) - vr[:, 1]
        out["cases"][str(d["nx"])] = {
            "snapshots_averaged": d["snapshots"], "z": d["z"].tolist(), "u": d["u"].tolist(),
            "x": d["x"].tolist(), "w": d["w"].tolist(),
            "rms_u": float(np.sqrt(np.mean(eu**2))), "rms_w": float(np.sqrt(np.mean(ew**2))),
            "max_u": float(np.abs(eu).max()), "max_w": float(np.abs(ew).max()),
            "u_min": float(d["u"].min()), "w_min": float(d["w"].min()), "w_max": float(d["w"].max())}
    a.set_xlabel(r"$u/U_{\mathrm{lid}}$")
    a.set_ylabel(r"$z/L$")
    a.set_xlim(-0.4, 1.05)
    a.set_ylim(0, 1)
    b.set_xlabel(r"$x/L$")
    b.set_ylabel(r"$w/U_{\mathrm{lid}}$")
    b.set_xlim(0, 1)
    a.legend(loc="upper left")
    plotstyle.panel_label(a, "(a)")
    plotstyle.panel_label(b, "(b)")
    os.makedirs(args.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(args.out, "cavity_validation.pdf"))
    with open(os.path.join(args.out, "cavity_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    for n, r in out["cases"].items():
        print(f"{n}^3: rms u={r['rms_u']:.4f} w={r['rms_w']:.4f} max u={r['max_u']:.4f} w={r['max_w']:.4f}")


if __name__ == "__main__":
    main()
