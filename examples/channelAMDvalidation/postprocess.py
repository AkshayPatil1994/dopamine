#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Mean streamwise-velocity profile of a wall-modelled half-channel run, compared with DNS.

Reads the NetCDF field output (cell-centred U, nu_sgs), averages Ux over x, y and the snapshots in
[t0, t1] for every wall-normal cell layer, writes profiles_<label>.csv, and, given a DNS profile,
prints U - U_DNS at y/delta = 0.1, 0.5, 0.99 and the bulk-velocity error. With --plot it writes a
two-panel manuscript figure: (a) <U>/u_tau vs y/delta, (b) the difference from DNS.

  postprocess.py amd_filtered=run_amd_filtered/channel.nc amd_raw=run_amd_raw/channel.nc \\
      --t0 14 --t1 24 --dns data/hoyas_jimenez_retau2000.csv --plot profiles.pdf

A .pdf --plot also writes a .png next to it. Labels amd_raw, amd_filtered, vreman_raw and
vreman_filtered get the paper's colours and line styles; other labels get a default style.
"""
import argparse
import sys
import os
import csv

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402
from ncprofiles import RE_TAU, mean_profile, read_dns, split_run

PROBES = (0.1, 0.5, 0.99)
_trapz = getattr(np, "trapezoid", None) or np.trapz


def marker_indices(z, n=9):
    """Indices of about n points spread evenly in log(z), for sparse markers on a log axis."""
    targets = np.logspace(np.log10(z[0]), np.log10(z[-1]), n)
    return np.unique(np.abs(np.log(z)[None, :] - np.log(targets)[:, None]).argmin(axis=1))


def make_figure(results, dns, args):
    plotstyle.apply()
    import matplotlib.pyplot as plt

    yd, ud = dns
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    xmin = min(z[0] for z, _ in results.values())

    ys = np.logspace(np.log10(xmin), 0, 200)
    a.semilogx(ys, np.log(ys * RE_TAU) / args.kappa + args.B, color=plotstyle.GREY, ls=":", lw=1.3,
               label=rf"Log law, $\kappa = {args.kappa:g}$")
    a.semilogx(yd, ud, "k-", lw=2.4, label="DNS (Lee & Moser)")
    b.axhline(0, color="k", lw=0.8)

    for k, (label, (z, u)) in enumerate(results.items()):
        st = plotstyle.run_style(label, k)
        mk = marker_indices(z) if st["marker"] else None
        kw = dict(color=st["color"], ls=st["ls"], marker=st["marker"], markevery=mk, label=st["label"])
        a.semilogx(z, u, **kw)
        b.semilogx(z, u - np.interp(z, yd, ud), **kw)

    for ax in (a, b):
        ax.set_xlim(xmin, 1)
        ax.set_xlabel(r"$y/\delta$")
    a.set_ylabel(r"$\langle U\rangle/u_\tau$")
    b.set_ylabel(r"$\langle U\rangle/u_\tau - U^+_{\mathrm{DNS}}$")
    a.legend(loc="lower right")
    plotstyle.panel_label(a, "(a)")
    plotstyle.panel_label(b, "(b)")
    if args.title:
        fig.suptitle(args.title)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=path/to/channel.nc")
    ap.add_argument("--t0", type=float, default=0.0, help="start of averaging window (default: 0)")
    ap.add_argument("--t1", type=float, default=np.inf, help="end of averaging window (default: latest snapshot)")
    ap.add_argument("--dns", help="DNS profile: y/delta, U+ (columns)")
    ap.add_argument("--plot", metavar="FILE", help="save the figure (.pdf also writes a .png; needs --dns)")
    ap.add_argument("--kappa", type=float, default=0.41, help="log-law slope in the figure")
    ap.add_argument("--B", type=float, default=4.8, help="log-law intercept in the figure")
    ap.add_argument("--title", help="figure title")
    args = ap.parse_args()

    dns = read_dns(args.dns) if args.dns else None
    results = {}
    bulk = {}
    for item in args.runs:
        label, path = split_run(item)
        z, dz, u, nus, ta, tb = mean_profile(path, args.t0, args.t1)
        with open(f"profiles_{label}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["y_over_delta", "dz", "U_plus", "nu_sgs"])
            w.writerows(zip(z, dz, u, nus))
        ub = float(np.sum(u * dz) / np.sum(dz))
        results[label] = (z, u)
        bulk[label] = ub
        print(f"{label}: averaged over t = {ta:.3g}..{tb:.3g}, Ub+ = {ub:.3f}, "
              f"first-cell U+ = {u[0]:.3f} at y+ = {z[0] * RE_TAU:.1f}")

    if dns is None:
        print("No --dns given: profiles written, no comparison made.")
        return
    yd, ud = dns
    ub_dns = float(_trapz(ud, yd) / (yd[-1] - yd[0]))
    print(f"\nDNS bulk U+ (trapezoid over the supplied profile) = {ub_dns:.3f}")
    print(f"{'run':<18}" + "".join(f"U-U_DNS@{p:<6}" for p in PROBES) + "  Ub error [%]")
    for label, (z, u) in results.items():
        diffs = [float(np.interp(p, z, u) - np.interp(p, yd, ud)) for p in PROBES]
        print(f"{label:<18}" + "".join(f"{x:+13.2f}" for x in diffs)
              + f"  {100 * (bulk[label] - ub_dns) / ub_dns:+8.2f}")

    if args.plot:
        for f in plotstyle.save(make_figure(results, dns, args), args.plot):
            print(f"wrote {f}")


if __name__ == "__main__":
    main()
