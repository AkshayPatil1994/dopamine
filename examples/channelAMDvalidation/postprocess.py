#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Mean streamwise-velocity profile of a wall-modelled half-channel run, compared with DNS.

Reads the NetCDF field output (cell-centred U, nu_sgs), averages Ux over x, y and the snapshots in
[t0, t1] for every wall-normal cell layer, writes profiles_<label>.csv, and, given a DNS profile,
prints U - U_DNS at y/delta = 0.1, 0.5, 0.99 and the bulk-velocity error.

  postprocess.py amd_filtered=run_amd_filtered/channel.nc amd_raw=run_amd_raw/channel.nc \\
      --t0 14 --t1 24 --dns data/hoyas_jimenez_retau2000.csv

The wall-normal direction is z (see channel_amd.geo). The runs are driven with fx = 1, so u_tau = 1
and U is already U+; y/delta = z because delta = 1.
"""
import argparse
import csv
import sys

import numpy as np
import netCDF4 as nc

PROBES = (0.1, 0.5, 0.99)
_trapz = getattr(np, "trapezoid", None) or np.trapz


def layer_profile(path, t0, t1):
    """Return z_c, dz, <Ux>, <nu_sgs> per layer plus the number of snapshots averaged."""
    d = nc.Dataset(path)
    t = np.asarray(d["time"][:])
    sel = np.where((t >= t0) & (t <= t1))[0]
    if sel.size == 0:
        sys.exit(f"{path}: no snapshots in [{t0}, {t1}] (run covers {t.min():.3g}..{t.max():.3g})")
    xyz = np.asarray(d["coords"][:])
    conn = np.asarray(d["connectivity"][:])
    off = np.asarray(d["offsets"][:])
    n_cells = off.size - 1
    # Cell z-centre and thickness from the vertex z of each cell (hexes: 8 vertices).
    zc = np.empty(n_cells)
    dz = np.empty(n_cells)
    for c in range(n_cells):
        zv = xyz[conn[off[c]:off[c + 1]], 2]
        zc[c] = zv.mean()
        dz[c] = zv.max() - zv.min()
    key = np.round(zc / dz.min() * 1e3).astype(np.int64)  # layer id; layers are z-aligned
    layers, inv = np.unique(key, return_inverse=True)
    n_l = layers.size
    cnt = np.bincount(inv, minlength=n_l)
    u_sum = np.zeros(n_l)
    nu_sum = np.zeros(n_l)
    for s in sel:
        u_sum += np.bincount(inv, weights=np.asarray(d["U"][s, :, 0]), minlength=n_l)
        nu_sum += np.bincount(inv, weights=np.asarray(d["nu_sgs"][s, :]), minlength=n_l)
    u = u_sum / (cnt * sel.size)
    nus = nu_sum / (cnt * sel.size)
    zl = np.bincount(inv, weights=zc, minlength=n_l) / cnt
    dzl = np.bincount(inv, weights=dz, minlength=n_l) / cnt
    order = np.argsort(zl)
    return zl[order], dzl[order], u[order], nus[order], sel.size


def read_dns(path):
    """Two-column CSV/whitespace file: y/delta, U+ ('#' comments allowed)."""
    arr = np.loadtxt(path, delimiter=None if "," not in open(path).read() else ",", comments="#")
    return arr[:, 0], arr[:, 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=path/to/channel.nc")
    ap.add_argument("--t0", type=float, default=14.0)
    ap.add_argument("--t1", type=float, default=24.0)
    ap.add_argument("--dns", help="DNS profile: y/delta, U+ (columns)")
    ap.add_argument("--plot", metavar="FILE", help="save a U+ vs y/delta figure (needs matplotlib)")
    args = ap.parse_args()

    dns = read_dns(args.dns) if args.dns else None
    results = {}
    for item in args.runs:
        label, _, path = item.partition("=")
        if not path:
            sys.exit(f"expected label=path, got '{item}'")
        z, dz, u, nus, n_snap = layer_profile(path, args.t0, args.t1)
        with open(f"profiles_{label}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["y_over_delta", "dz", "U_plus", "nu_sgs"])
            w.writerows(zip(z, dz, u, nus))
        ub = float(np.sum(u * dz) / np.sum(dz))
        results[label] = (z, u, ub)
        print(f"{label}: {n_snap} snapshots in [{args.t0}, {args.t1}], Ub+ = {ub:.3f}, "
              f"first-cell U+ = {u[0]:.3f} at y+ = {z[0] * 2000:.1f}")

    if dns is None:
        print("No --dns given: profiles written, no comparison made.")
        return
    yd, ud = dns
    ub_dns = float(_trapz(ud, yd) / (yd[-1] - yd[0]))
    print(f"\nDNS bulk U+ (trapezoid over the supplied profile) = {ub_dns:.3f}")
    print(f"{'run':<18}" + "".join(f"U-U_DNS@{p:<6}" for p in PROBES) + "  Ub error [%]")
    for label, (z, u, ub) in results.items():
        diffs = [float(np.interp(p, z, u) - np.interp(p, yd, ud)) for p in PROBES]
        print(f"{label:<18}" + "".join(f"{x:+13.2f}" for x in diffs) + f"  {100 * (ub - ub_dns) / ub_dns:+8.2f}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.semilogx(yd * 2000, ud, "k-", label="DNS")
        for label, (z, u, _) in results.items():
            ax.semilogx(z * 2000, u, "o-", ms=3, label=label)
        ax.set_xlabel("$y^+$")
        ax.set_ylabel("$U^+$")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)


if __name__ == "__main__":
    main()
