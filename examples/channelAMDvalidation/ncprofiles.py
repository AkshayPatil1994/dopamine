# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Read plane-averaged (x-y) wall-normal profiles from the DOPAMINE NetCDF field output.

The wall-normal direction is z; the runs are driven with fx = 1, so u_tau = 1, delta = 1 and U = U+.
"""
import sys

import netCDF4 as nc
import numpy as np

RE_TAU = 2000.0


def split_run(item):
    """Parse 'label=path' from the command line."""
    label, _, path = item.partition("=")
    if not path:
        sys.exit(f"expected label=path, got '{item}'")
    return label, path


def read_dns(path):
    """CSV with columns y/delta, U+ ('#' comments allowed)."""
    arr = np.loadtxt(path, delimiter=",", comments="#")
    return arr[:, 0], arr[:, 1]


def _layers(d):
    """Group cells into z-aligned layers. Returns inv (layer of each cell), cnt, z_l, dz_l, order (sorted by z)."""
    xyz = np.asarray(d["coords"][:])
    conn = np.asarray(d["connectivity"][:])
    off = np.asarray(d["offsets"][:])
    zv = xyz[conn, 2]
    start = off[:-1]
    zc = np.add.reduceat(zv, start) / np.diff(off)
    dz = np.maximum.reduceat(zv, start) - np.minimum.reduceat(zv, start)
    key = np.round(zc / dz.min() * 1e3).astype(np.int64)
    _, inv = np.unique(key, return_inverse=True)
    n_l = inv.max() + 1
    cnt = np.bincount(inv, minlength=n_l)
    zl = np.bincount(inv, weights=zc, minlength=n_l) / cnt
    dzl = np.bincount(inv, weights=dz, minlength=n_l) / cnt
    order = np.argsort(zl)
    return inv, cnt, zl[order], dzl[order], order


def layer_average(d, var, inv, cnt, order, i, comp=None):
    """Plane average of snapshot i of `var` per layer, sorted by z."""
    f = d[var][i, :, comp] if comp is not None else d[var][i, :]
    return (np.bincount(inv, weights=np.asarray(f), minlength=cnt.size) / cnt)[order]


def snapshot_profiles(path):
    """Return t, z, dz, U[t, layer] for every snapshot, sorted in time."""
    d = nc.Dataset(path)
    t = np.asarray(d["time"][:])
    inv, cnt, z, dz, order = _layers(d)
    u = np.stack([layer_average(d, "U", inv, cnt, order, i, 0) for i in range(t.size)])
    idx = np.argsort(t)
    return t[idx], z, dz, u[idx]


def mean_profile(path, t0, t1):
    """Return z, dz, <U>, <nu_sgs> averaged over snapshots in [t0, t1], and the first and last snapshot times used."""
    d = nc.Dataset(path)
    t = np.asarray(d["time"][:])
    sel = np.where((t >= t0) & (t <= t1))[0]
    if sel.size == 0:
        sys.exit(f"{path}: no snapshots in [{t0}, {t1}] (run covers {t.min():.3g}..{t.max():.3g})")
    inv, cnt, z, dz, order = _layers(d)
    u = sum(layer_average(d, "U", inv, cnt, order, s, 0) for s in sel) / sel.size
    nus = sum(layer_average(d, "nu_sgs", inv, cnt, order, s) for s in sel) / sel.size
    return z, dz, u, nus, t[sel[0]], t[sel[-1]]
