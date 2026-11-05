#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Validation animations in the shared manuscript style (structured-grid cases, NetCDF input).

  --case tgv      out-of-plane vorticity on two sections of the Taylor-Green vortex (z = L/2 and y = 3L/8)
  --case cavity   out-of-plane vorticity on two sections of the lid-driven cavity (y = 0.5 and z = 0.85)
  --case channel  streamwise fluctuation u'/u_tau on a wall-parallel plane (channel, Re_tau = 395)

  animate_slices.py --case tgv --nc tgv.nc --out docs/assets/validation/taylor_green/tgv_vorticity.gif
"""
import argparse
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle  # noqa: E402
from gifutil import frames_to_gif  # noqa: E402


def grid_index(nc):
    """Return (ix, iy, iz, shape, axes): the structured index of every cell from its centroid."""
    co, off, cn = nc["coords"][:], nc["offsets"][:], nc["connectivity"][:]
    n = np.diff(off)
    c = np.stack([np.add.reduceat(co[cn, k], off[:-1]) / n for k in range(3)], axis=1)
    out, shape, axes = [], [], []
    for k in range(3):
        u, inv = np.unique(np.round(c[:, k], 6), return_inverse=True)
        out.append(inv)
        shape.append(len(u))
        axes.append(u)
    return out, tuple(shape), axes


def to_grid(f, idx, shape):
    g = np.empty(shape, dtype=np.float64)
    g[idx[0], idx[1], idx[2]] = f
    return g


def nearest(a, v):
    return int(np.abs(a - v).argmin())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, choices=("tgv", "cavity", "channel"))
    ap.add_argument("--nc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--t-min", type=float, default=0.0)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--width", type=int, default=720)
    ap.add_argument("--layer", type=int, default=1, help="channel: wall-normal cell layer to show")
    ap.add_argument("--pct", type=float, default=99.0, help="percentile of |omega| that sets the colour limit")
    args = ap.parse_args()

    import netCDF4

    nc = netCDF4.Dataset(args.nc)
    idx, shape, ax = grid_index(nc)
    t = np.asarray(nc["time"][:])
    frames = [i for i in range(len(t)) if t[i] >= args.t_min][:: args.stride]
    plotstyle.apply()
    import matplotlib.pyplot as plt

    L = 2 * np.pi
    if args.case == "tgv":
        dx = [a[1] - a[0] for a in ax]
        kz, jy = nearest(ax[2], 0.5 * L), nearest(ax[1], 0.375 * L)
        cmap, sym = plotstyle.vorticity_cmap(), True
        label = r"$\omega\,L/U_0$"
        panels = [dict(ext=(0, 1, 0, 1), xl=r"$x/L$", yl=r"$y/L$", txt=r"$\omega_z$, $z=L/2$"),
                  dict(ext=(0, 1, 0, 1), xl=r"$x/L$", yl=r"$z/L$", txt=rf"$\omega_y$, $y={ax[1][jy] / L:.3f}L$")]
        scale = L

        def frame(i):
            u = [to_grid(nc["U"][i][:, k], idx, shape) for k in range(3)]
            d = lambda f, a: (np.roll(f, -1, a) - np.roll(f, 1, a)) / (2 * dx[a])
            wz = (d(u[1], 0) - d(u[0], 1))[:, :, kz].T
            wy = (d(u[0], 2) - d(u[2], 0))[:, jy, :].T
            return [wz, wy]
    elif args.case == "cavity":
        jy, kz = nearest(ax[1], 0.5), nearest(ax[2], 0.85)
        cmap, sym = plotstyle.vorticity_cmap(), True
        label = r"$\omega\,L/U_{\mathrm{lid}}$"
        panels = [dict(ext=(0, 1, 0, 1), xl=r"$x/L$", yl=r"$z/L$", txt=r"$\omega_y$, $y=L/2$"),
                  dict(ext=(0, 1, 0, 1), xl=r"$x/L$", yl=r"$y/L$", txt=rf"$\omega_z$, $z={ax[2][kz]:.2f}L$")]
        scale = 1.0

        def frame(i):
            u = [to_grid(nc["U"][i][:, k], idx, shape) for k in range(3)]
            g = lambda f, a: np.gradient(f, ax[a], axis=a)
            wy = (g(u[0], 2) - g(u[2], 0))[:, jy, :].T
            wz = (g(u[1], 0) - g(u[0], 1))[:, :, kz].T
            return [wy, wz]
    else:
        cmap, sym = "RdBu_r", True
        label = r"$u'/u_\tau$"
        panels = [dict(ext=(ax[0][0], ax[0][-1], ax[1][0], ax[1][-1]), xl=r"$x/h$", yl=r"spanwise coordinate$/h$",
                       txt=None)]

        def frame(i):
            u = to_grid(nc["U"][i][:, 0], idx, shape)
            return [(u[:, :, args.layer] - u.mean(axis=(0, 1))[args.layer]).T]

    fields = [frame(i) for i in frames]
    if args.case == "channel":
        vmax = 5.0
    else:
        vmax = float(np.percentile(np.abs(np.stack([np.stack(f) for f in fields])), args.pct))
    npan = len(panels)
    wide = args.case == "channel"
    with tempfile.TemporaryDirectory() as tmp:
        for k, (i, fl) in enumerate(zip(frames, fields)):
            if wide:
                fig, axs = plt.subplots(1, 1, figsize=(7.0, 2.8), constrained_layout=True)
                axs = [axs]
            else:
                fig, axs = plt.subplots(1, npan, figsize=(7.4, 3.9), constrained_layout=True)
            for a, f, p in zip(axs, fl, panels):
                im = a.imshow(f, origin="lower", extent=p["ext"], cmap=cmap, vmin=-vmax, vmax=vmax, aspect="equal",
                              interpolation="bilinear")
                a.set_xlabel(p["xl"])
                a.set_ylabel(p["yl"])
                if p["txt"]:
                    a.text(0.03, 0.97, p["txt"], transform=a.transAxes, ha="left", va="top",
                           bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
            if not wide:
                for a, tag in zip(axs, "ab"):
                    plotstyle.panel_label(a, f"({tag})")
                axs[1].set_ylabel("")
                fig.colorbar(im, ax=axs, orientation="horizontal", pad=0.04, shrink=0.7, aspect=40, label=label)
                fig.text(0.99, 0.985, rf"$t={t[i]:.1f}$", ha="right", va="top")
            else:
                axs[0].text(0.02, 0.96, rf"$t={t[i]:.1f}$", transform=axs[0].transAxes, ha="left", va="top",
                            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
                fig.colorbar(im, ax=axs[0], pad=0.02, label=label)
            fig.savefig(os.path.join(tmp, f"f_{k:04d}.png"), dpi=130)
            plt.close(fig)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        frames_to_gif(os.path.join(tmp, "f_%04d.png"), args.out, fps=args.fps, width=args.width)
    print(f"{args.out}: {len(frames)} frames, vmax={vmax:.3g}, {os.path.getsize(args.out) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
