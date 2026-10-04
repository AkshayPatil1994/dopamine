#!/usr/bin/env python3
"""Cross-check a recorded database against the solver's own NetCDF output (frame differences,
so the removed mean cancels). Usage: check_recording.py case.nc rec.nc x_layer"""
import sys
import numpy as np
import netCDF4 as nc

case, rec, xl = sys.argv[1], sys.argv[2], float(sys.argv[3])
d, db = nc.Dataset(case), nc.Dataset(rec)
X = np.array(d['coords'][:]); conn = np.array(d['connectivity'][:]); off = np.array(d['offsets'][:])
cent = np.array([X[conn[off[i]:off[i + 1]]].mean(axis=0) for i in range(len(off) - 1)])
t = np.array(d['time'][:]); U = d['U']
y = np.array(db['y'][:]); z = np.array(db['z'][:]); u = db['u']      # u(seg, time, comp, z, y)
dt = float(db.getncattr('dt_p'))
layer = np.where(np.abs(cent[:, 0] - xl) < 1e-6)[0]
assert len(layer) == len(y) * len(z), (len(layer), len(y), len(z))

def nc_plane(ti):
    P = np.zeros((len(z), len(y), 3))
    for c in layer:
        P[np.argmin(abs(z - cent[c, 2])), np.argmin(abs(y - cent[c, 1])), :] = U[ti, c, :]
    return P

def db_plane(n):
    return np.transpose(np.array(u[0, n, :, :, :]), (1, 2, 0))

worst, scale = 0.0, 0.0
for ia, ib in [(0, 1), (0, 2), (1, 3)]:
    na, nb = int(round(t[ia] / dt)), int(round(t[ib] / dt))
    dnc = nc_plane(ia) - nc_plane(ib); ddb = db_plane(na) - db_plane(nb)
    worst = max(worst, np.abs(dnc - ddb).max()); scale = max(scale, np.abs(dnc).max())
mean = np.abs(np.array(u[0]).mean(axis=(0, 3))).max()
print(f'recording check: max difference {worst:.2e} (field differences up to {scale:.3f}); residual mean {mean:.1e}')
sys.exit(0 if worst < 1e-5 * max(scale, 1.0) and mean < 1e-6 else 1)
