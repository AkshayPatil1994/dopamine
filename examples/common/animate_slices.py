#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Validation animations in the shared manuscript style (structured-grid cases, NetCDF input).

  --case tgv      vorticity magnitude |omega| on the z = pi plane (Taylor-Green, Re = 1600)
  --case cavity   velocity magnitude |U|/U_lid on the y = 0.5 mid-plane (lid-driven cavity)
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
    """Return (ix, iy, iz, shape): the structured index of every cell from its centroid."""
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
    args = ap.parse_args()

    import netCDF4

    nc = netCDF4.Dataset(args.nc)
    idx, shape, ax = grid_index(nc)
    t = np.asarray(nc["time"][:])
    frames = [i for i in range(len(t)) if t[i] >= args.t_min][:: args.stride]
    plotstyle.apply()
    import matplotlib.pyplot as plt

    if args.case == "tgv":
        dx = [a[1] - a[0] for a in ax]
        kmid = int(np.abs(ax[2] - np.pi).argmin())
        cmap, label, vmax, vmin, extent = "viridis", r"$|\omega|\,L/U_0$", None, 0.0, None
        xl, yl, xs = r"$x/L$", r"$y/L$", 1 / (2 * np.pi)
        ext = (0, 2 * np.pi * xs, 0, 2 * np.pi * xs)

        def frame(i):
            u = [to_grid(nc["U"][i][:, k], idx, shape) for k in range(3)]
            d = lambda f, a: (np.roll(f, -1, a) - np.roll(f, 1, a)) / (2 * dx[a])
            wx = d(u[2], 1) - d(u[1], 2)
            wy = d(u[0], 2) - d(u[2], 0)
            wz = d(u[1], 0) - d(u[0], 1)
            return np.sqrt(wx**2 + wy**2 + wz**2)[:, :, kmid].T
    elif args.case == "cavity":
        jmid = int(np.abs(ax[1] - 0.5).argmin())
        cmap, label, vmax, vmin = "viridis", r"$|\mathbf{U}|/U_{\mathrm{lid}}$", 1.0, 0.0
        xl, yl = r"$x/L$", r"$z/L$"
        ext = (0, 1, 0, 1)

        def frame(i):
            u = [to_grid(nc["U"][i][:, k], idx, shape) for k in (0, 2)]
            return np.sqrt(u[0] ** 2 + u[1] ** 2)[:, jmid, :].T
    else:
        cmap, label, vmax, vmin = "RdBu_r", r"$u'/u_\tau$", 5.0, -5.0
        xl, yl = r"$x/h$", r"spanwise coordinate$/h$"
        ext = (ax[0][0], ax[0][-1], ax[1][0], ax[1][-1])

        def frame(i):
            u = to_grid(nc["U"][i][:, 0], idx, shape)
            return (u[:, :, args.layer] - u.mean(axis=(0, 1))[args.layer]).T

    fields = [frame(i) for i in frames]
    if vmax is None:
        vmax = float(np.percentile(np.stack(fields), 99.5))
    aspect = (ext[1] - ext[0]) / (ext[3] - ext[2])
    h = 3.4 if aspect < 1.5 else 2.6
    with tempfile.TemporaryDirectory() as tmp:
        for k, (i, f) in enumerate(zip(frames, fields)):
            fig, a = plt.subplots(figsize=(h * min(aspect, 3.0) + 1.3, h), constrained_layout=True)
            im = a.imshow(f, origin="lower", extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal",
                          interpolation="bilinear")
            a.set_xlabel(xl)
            a.set_ylabel(yl)
            a.text(0.02, 0.96, rf"$t={t[i]:.1f}$", transform=a.transAxes, ha="left", va="top",
                   bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
            fig.colorbar(im, ax=a, pad=0.02, label=label)
            fig.savefig(os.path.join(tmp, f"f_{k:04d}.png"), dpi=130)
            plt.close(fig)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        frames_to_gif(os.path.join(tmp, "f_%04d.png"), args.out, fps=args.fps, width=args.width)
    print(f"{args.out}: {len(frames)} frames, {os.path.getsize(args.out) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
