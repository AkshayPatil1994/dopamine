#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Michelstadt urban LES (work in progress): pedestrian-level velocity-magnitude snapshots and an animation.

Reads the cell-centred U of urban.nc, takes the cells in the slab z0 <= z <= z1 (first cell layer), rasterises
|U| to a regular grid (cells farther than --mask from a cell centre, i.e. buildings, stay blank) and writes
michelstadt_snapshots.{png,pdf}; with --gif it also renders every --stride-th snapshot into
michelstadt_plan.gif.

  plot_snapshots.py --nc urban.nc --out docs/assets/validation/michelstadt --times 400 1560 --gif
"""
import argparse
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402
from gifutil import frames_to_gif  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nc", required=True)
    ap.add_argument("--out", default=".")
    ap.add_argument("--times", type=float, nargs=2, default=(400.0, 1560.0), help="snapshot times for the figure")
    ap.add_argument("--z0", type=float, default=1.0)
    ap.add_argument("--z1", type=float, default=3.0)
    ap.add_argument("--h", type=float, default=3.0, help="raster spacing in m")
    ap.add_argument("--mask", type=float, default=6.0, help="blank raster points farther than this from a cell")
    ap.add_argument("--vmax", type=float, default=6.0)
    ap.add_argument("--gif", action="store_true")
    ap.add_argument("--stride", type=int, default=3)
    args = ap.parse_args()

    import netCDF4
    from scipy.spatial import cKDTree

    d = netCDF4.Dataset(args.nc)
    co, off, cn = d["coords"][:], d["offsets"][:], d["connectivity"][:]
    n = np.diff(off)
    c = np.stack([np.add.reduceat(co[cn, k], off[:-1]) / n for k in range(3)], axis=1)
    slab = np.where((c[:, 2] >= args.z0) & (c[:, 2] <= args.z1))[0]
    xy = c[slab, :2]
    t = np.asarray(d["time"][:])
    gx = np.arange(xy[:, 0].min(), xy[:, 0].max(), args.h)
    gy = np.arange(xy[:, 1].min(), xy[:, 1].max(), args.h)
    X, Y = np.meshgrid(gx, gy)
    dist, idx = cKDTree(xy).query(np.column_stack([X.ravel(), Y.ravel()]))
    blank = (dist > args.mask).reshape(X.shape)

    def field(i):
        u = np.linalg.norm(d["U"][i][slab], axis=1)
        f = u[idx].reshape(X.shape)
        return np.ma.array(f, mask=blank)

    plotstyle.apply()
    import matplotlib.pyplot as plt

    ext = (gx[0] / 1e3, gx[-1] / 1e3, gy[0] / 1e3, gy[-1] / 1e3)

    def draw(ax, f, label, cbar=False, fig=None):
        im = ax.imshow(f, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=args.vmax, aspect="equal",
                       interpolation="nearest")
        ax.set_facecolor("#bdbdbd")
        ax.set_xlabel(r"$x$ (km)")
        ax.set_ylabel(r"$y$ (km)")
        ax.text(0.01, 0.97, label, transform=ax.transAxes, ha="left", va="top",
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
        return im

    ia = [int(np.abs(t - tt).argmin()) for tt in args.times]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), constrained_layout=True, sharex=True)
    for k, (ax, i) in enumerate(zip(axes, ia)):
        im = draw(ax, field(i), rf"$t={t[i]:.0f}$ s")
        plotstyle.panel_label(ax, f"({'ab'[k]})")
    axes[0].set_xlabel("")
    fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02, label=r"$|\mathbf{U}|$ at $z\approx1.5$ m (m s$^{-1}$)")
    os.makedirs(args.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(args.out, "michelstadt_snapshots.pdf"))
    plt.close(fig)

    if args.gif:
        with tempfile.TemporaryDirectory() as tmp:
            sel = list(range(0, len(t), args.stride))
            for k, i in enumerate(sel):
                f2, ax = plt.subplots(figsize=(7.2, 3.4), constrained_layout=True)
                im = draw(ax, field(i), rf"$t={t[i]:.0f}$ s")
                f2.colorbar(im, ax=ax, pad=0.02, label=r"$|\mathbf{U}|$ (m s$^{-1}$)")
                f2.savefig(os.path.join(tmp, f"f_{k:04d}.png"), dpi=130)
                plt.close(f2)
            frames_to_gif(os.path.join(tmp, "f_%04d.png"), os.path.join(args.out, "michelstadt_plan.gif"), fps=8,
                          width=720)


if __name__ == "__main__":
    main()
