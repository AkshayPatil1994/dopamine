#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Wall-modelled LES channel at Re_tau = 395 against the Moser-Kim-Mansour (1999) DNS.

Reads the wall-normal LES statistics (les_validation_mkm.csv, columns z, y_plus, U_plus, u/v/w_rms_plus,
minus_uv_plus, tau_*_plus) and the MKM tables chan395.means / chan395.reystress, writes
channel395_validation.{png,pdf} and channel395_validation.json (reduced profiles plus error metrics).

  plot_validation.py --les les_validation_mkm.csv --means MKM_MEANS.dat --stress MKM_REYSTRESS.dat \\
      --out docs/assets/validation/channel395
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import plotstyle  # noqa: E402

KAPPA, B = 0.41, 5.2
BLUE, ORANGE, GREEN = plotstyle.BLUE, plotstyle.ORANGE, "#009E73"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--les", required=True)
    ap.add_argument("--means", required=True)
    ap.add_argument("--stress", required=True)
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    names = open(args.les).readline().strip().split(",")
    les = {k: v for k, v in zip(names, np.loadtxt(args.les, delimiter=",", skiprows=1).T)}
    mm = np.loadtxt(args.means, comments="#")
    ms = np.loadtxt(args.stress, comments="#")
    yd, ypd, Ud = mm[:, 0], mm[:, 1], mm[:, 2]
    rms = {"u": np.sqrt(ms[:, 2]), "v": np.sqrt(ms[:, 3]), "w": np.sqrt(ms[:, 4])}
    uvd = -ms[:, 5]
    kd = 0.5 * (ms[:, 2] + ms[:, 3] + ms[:, 4])
    yp, z = les["y_plus"], les["z"]
    k_les = 0.5 * (les["u_rms_plus"] ** 2 + les["v_rms_plus_wall_normal"] ** 2 + les["w_rms_plus_spanwise"] ** 2)

    plotstyle.apply()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.6), constrained_layout=True)
    a, b, c, d, e, f = ax.ravel()
    kw = dict(color=BLUE, ls="none", marker="o")
    dns = dict(color="k", lw=2.0)

    ys = np.logspace(0, np.log10(ypd[-1]), 200)
    a.semilogx(ypd[1:], Ud[1:], label="DNS (MKM)", **dns)
    a.semilogx(ys[ys < 6], ys[ys < 6], color=plotstyle.GREY, ls=":", lw=1.3, label=r"$U^+=y^+$")
    ysl = ys[ys > 25]
    a.semilogx(ysl, np.log(ysl) / KAPPA + B, color=plotstyle.GREY, ls="--", lw=1.3,
               label=rf"$\kappa^{{-1}}\ln y^+ + B$")
    a.semilogx(yp, les["U_plus"], label="DOPAMINE", **kw)
    a.set_xlim(1, 500)
    a.set_ylim(0, 22)
    a.set_ylabel(r"$\langle U\rangle^+$")
    a.legend(loc="upper left")

    for axis, key, lab, ref in ((b, "u_rms_plus", r"$u'^+_{\mathrm{rms}}$", "u"),
                                (c, "v_rms_plus_wall_normal", r"$v'^+_{\mathrm{rms}}$", "v"),
                                (d, "w_rms_plus_spanwise", r"$w'^+_{\mathrm{rms}}$", "w")):
        axis.semilogx(ypd[1:], rms[ref][1:], **dns)
        axis.semilogx(yp, les[key], **kw)
        axis.set_xlim(1, 500)
        axis.set_ylim(0, None)
        axis.set_ylabel(lab)

    e.plot(yd[1:], uvd[1:], label=r"$-\langle u'v'\rangle^+$ DNS", **dns)
    e.plot(z, les["minus_uv_plus"], label=r"$-\langle u'v'\rangle^+$ resolved", **kw)
    e.plot(z, les["tau_total_plus"], color=ORANGE, ls="none", marker="s", label=r"$\tau^+_{\mathrm{tot}}$")
    e.plot(z, les["tau_sgs_plus"], color=GREEN, ls="none", marker="^", label=r"$\tau^+_{\mathrm{sgs}}$")
    e.plot(z, les["tau_visc_plus"], color=plotstyle.GREY, ls="none", marker="v", label=r"$\tau^+_{\nu}$")
    e.plot(yd, 1 - yd, color=plotstyle.GREY, ls="--", lw=1.3, label=r"$1-y/h$")
    e.set_xlim(0, 1)
    e.set_ylim(0, 1.2)
    e.set_ylabel(r"stress$^+$")
    e.legend(loc="upper right", fontsize=7.5, ncol=1)

    f.semilogx(ypd[1:], kd[1:], **dns)
    f.semilogx(yp, k_les, **kw)
    f.set_xlim(1, 500)
    f.set_ylim(0, None)
    f.set_ylabel(r"$k^+$")

    for axis in (a, b, c, d, f):
        axis.set_xlabel(r"$y^+$")
    e.set_xlabel(r"$y/h$")
    for axis, tag in zip(ax.ravel(), "abcdef"):
        plotstyle.panel_label(axis, f"({tag})")

    os.makedirs(args.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(args.out, "channel395_validation.pdf"))

    sel = yp > 30
    interp = lambda ref: np.interp(yp, ypd, ref)
    metrics = {
        "n_points": int(len(yp)),
        "U_rel_err_pct_max_y30": float(100 * np.max(np.abs(les["U_plus"][sel] / interp(Ud)[sel] - 1))),
        "U_rel_err_pct_mean_y30": float(100 * np.mean(np.abs(les["U_plus"][sel] / interp(Ud)[sel] - 1))),
    }
    for axis_key, key, ref in (("u", "u_rms_plus", "u"), ("v", "v_rms_plus_wall_normal", "v"),
                               ("w", "w_rms_plus_spanwise", "w")):
        metrics[f"{axis_key}_rms_rel_err_pct_mean_y30"] = float(
            100 * np.mean(np.abs(les[key][sel] / interp(rms[ref])[sel] - 1)))
    metrics["uv_rel_err_pct_mean_y30"] = float(
        100 * np.mean(np.abs(les["minus_uv_plus"][sel] / interp(uvd)[sel] - 1)))
    out = {"les": {k: v.tolist() for k, v in les.items()}, "metrics": metrics}
    with open(os.path.join(args.out, "channel395_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()
