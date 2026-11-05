#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Time evolution of the x-y averaged mean velocity, to judge whether the statistics have settled.

For every snapshot in the NetCDF output it computes the plane-averaged Ux per wall-normal layer
(as postprocess.py does) and writes a four-panel manuscript figure:

  (a) the bulk velocity Ub+(t) and its running mean from t_start, against the DNS bulk velocity,
  (b) Ux+(t) at y/delta = 0.1, 0.5, 0.99 together with the DNS value (black),
  (c) the change of the running-mean profile, max_z |<U>_[t_start,t] - <U>_[t_start,t_end]|, versus t
      (last snapshot omitted, it is zero by construction),
  (d) the mean over [t0, t_end] at the probes versus the averaging start t0.

  plot_convergence.py amd_filtered=run_amd_filtered/channel.nc --dns data/hoyas_jimenez_retau2000.csv \\
      --plot convergence.pdf

A .pdf --plot also writes a .png next to it. Line colour identifies the run, line style the probe.
"""
import argparse
import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402
from ncprofiles import read_dns, snapshot_profiles, split_run

PROBES = (0.1, 0.5, 0.99)
PROBE_LS = ("-", "--", ":")
_trapz = getattr(np, "trapezoid", None) or np.trapz


def running_mean(x, t):
    """Time-weighted (trapezoid) running mean of x[t, ...] from t[0]; first entry is x[0]."""
    out = np.empty_like(x)
    out[0] = x[0]
    acc = np.zeros_like(x[0])
    for i in range(1, t.size):
        acc += 0.5 * (x[i] + x[i - 1]) * (t[i] - t[i - 1])
        out[i] = acc / (t[i] - t[0])
    return out


def window_mean(x, t, t0):
    """Trapezoid mean of x over [t0, t_end]; needs at least two snapshots in the window."""
    m = t >= t0
    if m.sum() < 2:
        return np.full(x.shape[1:], np.nan)
    return _trapz(x[m], t[m], axis=0) / (t[m][-1] - t[m][0])


def probe_values(z, u, probes=PROBES):
    """Interpolate u[..., layer] to the probe heights."""
    u = np.atleast_2d(u)
    return np.array([[np.interp(p, z, row) for p in probes] for row in u])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=path/to/channel.nc")
    ap.add_argument("--dns", help="DNS profile: y/delta, U+ (columns); adds reference lines")
    ap.add_argument("--plot", metavar="FILE", default="convergence.pdf")
    args = ap.parse_args()

    plotstyle.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    dns = None
    if args.dns:
        yd, ud = read_dns(args.dns)
        dns = [float(np.interp(p, yd, ud)) for p in PROBES]
        ub_dns = float(_trapz(ud, yd) / (yd[-1] - yd[0]))

    fig, ax = plt.subplots(2, 2, figsize=(9, 6.6), constrained_layout=True)
    (a, b), (c, d) = ax
    run_handles = []

    for k, item in enumerate(args.runs):
        label, path = split_run(item)
        t, z, dz, u = snapshot_profiles(path)
        if t.size < 3:
            raise SystemExit(f"{path}: only {t.size} snapshots, nothing to check")
        st = plotstyle.run_style(label, k)
        col = st["color"]
        run_handles.append(Line2D([], [], color=col, label=st["label"]))

        ub = (u * dz).sum(axis=1) / dz.sum()
        a.plot(t, ub, color=col, lw=0.8, alpha=0.45)
        a.plot(t, running_mean(ub, t), color=col)

        for j, pr in enumerate(probe_values(z, u).T):
            b.plot(t, pr, color=col, ls=PROBE_LS[j], lw=1.2)

        run_u = running_mean(u, t)
        dev = np.abs(run_u - run_u[-1]).max(axis=1)
        c.semilogy(t[:-1], dev[:-1],  # last point is zero by construction
                    color=col, marker="o", ms=3, mfc=col)

        t0s = t[:-2]
        wm = np.array([[np.interp(p, z, window_mean(u, t, t0)) for p in PROBES] for t0 in t0s])
        for j in range(len(PROBES)):
            d.plot(t0s, wm[:, j], color=col, ls=PROBE_LS[j], lw=1.2)

        print(f"{label}: {t.size} snapshots, t = {t[0]:.3g}..{t[-1]:.3g}")
        print(f"  {'t0':>7}" + "".join(f"  <U>@{p:<5}" for p in PROBES) + "      Ub")
        for t0 in np.linspace(t[0], t[-3], 6):
            m = window_mean(u, t, t0)
            print(f"  {t0:7.2f}" + "".join(f"  {np.interp(p, z, m):10.3f}" for p in PROBES)
                  + f"  {np.sum(m * dz) / dz.sum():7.3f}")

    probe_handles = [Line2D([], [], color="k", ls=ls, lw=1.2, label=rf"$y/\delta = {p:g}$")
                     for p, ls in zip(PROBES, PROBE_LS)]
    if dns:
        a.axhline(ub_dns, color="k", lw=0.8)
        run_handles.append(Line2D([], [], color="k", lw=0.8, label=rf"DNS $U_b^+ = {ub_dns:.2f}$"))
        for j, v in enumerate(dns):
            for ax_ in (b, d):
                ax_.axhline(v, color="k", ls=PROBE_LS[j], lw=0.8)
        probe_handles.append(Line2D([], [], color="k", lw=0.8, label="DNS"))

    a.set_ylabel(r"$U_b/u_\tau$")
    b.set_ylabel(r"$\langle U\rangle_{xy}/u_\tau$")
    c.set_ylabel(r"$\max|\Delta\langle U\rangle|/u_\tau$")
    d.set_ylabel(r"$\langle U\rangle_{[t_0,t_e]}/u_\tau$")
    for ax_ in (a, b, c):
        ax_.set_xlabel(r"$t\,u_\tau/\delta$")
    d.set_xlabel(r"averaging start $t_0\,u_\tau/\delta$")
    a.legend(handles=run_handles, loc="lower right")
    b.legend(handles=probe_handles, loc="best", ncol=2)
    for ax_, s in zip(ax.flat, "abcd"):
        plotstyle.panel_label(ax_, f"({s})")

    for f in plotstyle.save(fig, args.plot):
        print(f"wrote {f}")


if __name__ == "__main__":
    main()
