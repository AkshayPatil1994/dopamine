"""
Animate the lid-driven cavity centre plane (y=0.5) over all snapshots.

Reads every time record in the dopamine NetCDF output, extracts the x-z plane
closest to y=0.5, and writes an animated GIF of velocity magnitude with
streamlines overlaid.

Usage:
    python3 animate_cavity.py [--nc PATH] [--out GIF] [--fps FPS] [--stride N]
                              [--vmax VMAX] [--dpi DPI] [--no-streamlines]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.interpolate import griddata


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--nc', default=None, help='Path to the .nc output (default: auto-detect)')
    p.add_argument('--out', default='cavity_plane.gif', help='Output GIF (default: cavity_plane.gif)')
    p.add_argument('--fps', type=float, default=4.0, help='Frames per second (default: 4)')
    p.add_argument('--stride', type=int, default=1, help='Use every Nth snapshot (default: 1)')
    p.add_argument('--vmax', type=float, default=None,
                   help='Colour-scale maximum (default: max over all snapshots)')
    p.add_argument('--dpi', type=int, default=110, help='GIF resolution (default: 110)')
    p.add_argument('--no-streamlines', action='store_true', help='Skip the streamline overlay')
    return p.parse_args()


def find_nc(path_hint):
    if path_hint:
        if not os.path.isfile(path_hint):
            sys.exit(f'ERROR: {path_hint} not found')
        return path_hint
    here = os.path.dirname(os.path.abspath(__file__))
    cands = sorted(f for f in os.listdir(here) if f.endswith('.nc'))
    if not cands:
        sys.exit('ERROR: no .nc file found — pass --nc /path/to/file')
    if len(cands) > 1:
        print(f'  Multiple .nc files found, using {cands[0]}')
    return os.path.join(here, cands[0])


def compute_centroids(connectivity, offsets, coords):
    """Vectorised cell-centroid computation from CSR vertex list."""
    starts = offsets[:-1]
    n_per = offsets[1:] - starts
    cx = np.add.reduceat(coords[connectivity, 0], starts) / n_per
    cy = np.add.reduceat(coords[connectivity, 1], starts) / n_per
    cz = np.add.reduceat(coords[connectivity, 2], starts) / n_per
    return cx, cy, cz


def read_nu(nc_path):
    """Parse nu from cavity.input next to the .nc file (for the Re label)."""
    input_path = os.path.join(os.path.dirname(os.path.abspath(nc_path)), 'cavity.input')
    if os.path.isfile(input_path):
        for line in open(input_path):
            s = line.strip()
            if s.startswith('!') or '=' not in s:
                continue
            if s.split('=')[0].strip() == 'nu':
                try:
                    return float(s.split('=')[1].split()[0].rstrip(','))
                except ValueError:
                    pass
    return None


def build_plane_binning(cx, cy, cz):
    """Select the cells nearest y=0.5 and map them onto a structured x-z grid.

    Returns (cell_idx, row, col, x_grid, z_grid) where row/col give the grid
    slot of each selected cell, so per-snapshot fields are binned with a single
    np.add.at instead of a per-grid-point search.
    """
    y_unique = np.array(sorted(set(np.round(cy, 4))))
    k = int(np.argmin(np.abs(y_unique - 0.5)))
    y1 = y_unique[max(0, k - 1)]
    y2 = y_unique[min(len(y_unique) - 1, k + 1)]

    mask = (np.abs(cy - y1) < 0.01) | (np.abs(cy - y2) < 0.01)
    cell_idx = np.flatnonzero(mask)
    if cell_idx.size == 0:
        sys.exit('ERROR: no cells found near the y=0.5 plane')
    print(f'  Plane cells: {cell_idx.size} between y={y1:.4f} and y={y2:.4f}')

    x_sel = np.round(cx[cell_idx], 4)
    z_sel = np.round(cz[cell_idx], 4)
    x_grid, col = np.unique(x_sel, return_inverse=True)
    z_grid, row = np.unique(z_sel, return_inverse=True)
    print(f'  Structured grid: {len(x_grid)} x {len(z_grid)} (x by z)')
    return cell_idx, row, col, x_grid, z_grid


def bin_snapshot(U_snap, row, col, shape):
    """Average cell values onto the (nz, nx) grid; empty slots become NaN."""
    counts = np.zeros(shape)
    np.add.at(counts, (row, col), 1.0)
    out = []
    for comp in range(U_snap.shape[1]):
        acc = np.zeros(shape)
        np.add.at(acc, (row, col), U_snap[:, comp])
        with np.errstate(invalid='ignore', divide='ignore'):
            out.append(np.where(counts > 0, acc / np.maximum(counts, 1), np.nan))
    return out


def main():
    args = parse_args()
    nc_path = find_nc(args.nc)
    print(f'Loading: {nc_path}')

    try:
        import netCDF4
    except ImportError:
        sys.exit('ERROR: netCDF4 not installed. Run: pip install netCDF4')

    nc = netCDF4.Dataset(nc_path)
    coords = nc.variables['coords'][:]
    offsets = nc.variables['offsets'][:]
    conn = nc.variables['connectivity'][:]
    times = np.asarray(nc.variables['time'][:])
    U_var = nc.variables['U']

    n_t = U_var.shape[0]
    print(f'  n_cells={offsets.shape[0] - 1}, n_snapshots={n_t}, '
          f't=[{times[0]:.3g}, {times[-1]:.3g}]')

    nu = read_nu(nc_path)
    re_label = f'Re={1.0 / nu:.0f}' if nu else 'Re=?'

    cx, cy, cz = compute_centroids(conn, offsets, coords)
    cell_idx, row, col, x_grid, z_grid = build_plane_binning(cx, cy, cz)
    shape = (len(z_grid), len(x_grid))
    X_grid, Z_grid = np.meshgrid(x_grid, z_grid)

    frames = list(range(0, n_t, max(1, args.stride)))
    print(f'Reading {len(frames)} snapshots of the centre plane...')

    plane_frames = []
    for t in frames:
        U_snap = np.asarray(U_var[t, :, :])[cell_idx, :]
        u, v, w = bin_snapshot(U_snap, row, col, shape)
        plane_frames.append((u, v, w, np.sqrt(u ** 2 + v ** 2 + w ** 2)))
    nc.close()

    vmax = args.vmax
    if vmax is None:
        vmax = float(np.nanmax([np.nanmax(f[3]) for f in plane_frames]))
        vmax = max(vmax, 1e-6)
    print(f'  Colour scale: 0 to {vmax:.3f}')

    # Uniform grid for the streamline overlay (mesh is graded)
    x_uniform = np.linspace(x_grid.min(), x_grid.max(), 40)
    z_uniform = np.linspace(z_grid.min(), z_grid.max(), 40)
    X_uniform, Z_uniform = np.meshgrid(x_uniform, z_uniform)
    points = np.column_stack([X_grid.ravel(), Z_grid.ravel()])

    levels = np.linspace(0, vmax, 21)
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.colorbar(ScalarMappable(norm=Normalize(0, vmax), cmap='magma'),
                        ax=ax, label='Velocity magnitude')

    writer = PillowWriter(fps=args.fps)
    print(f'Rendering {len(frames)} frames -> {args.out}')
    with writer.saving(fig, args.out, dpi=args.dpi):
        for n, (t, (u, v, w, vmag)) in enumerate(zip(frames, plane_frames)):
            ax.clear()
            ax.contourf(X_grid, Z_grid, vmag, levels=levels, cmap='magma', extend='max')
            ax.contour(X_grid, Z_grid, vmag, levels=levels,
                       colors='black', linewidths=0.3, alpha=0.3)

            if not args.no_streamlines:
                U_uniform = griddata(points, u.ravel(), (X_uniform, Z_uniform), method='cubic')
                W_uniform = griddata(points, w.ravel(), (X_uniform, Z_uniform), method='cubic')
                U_uniform = np.nan_to_num(U_uniform)
                W_uniform = np.nan_to_num(W_uniform)
                ax.streamplot(x_uniform, z_uniform, U_uniform, W_uniform,
                              color='white', linewidth=0.8, density=1.5, arrowsize=1.2)

            ax.set_xlabel('$x/L$', fontsize=12)
            ax.set_ylabel('$z/L$', fontsize=12)
            ax.set_title(f'Centre plane (y=0.5), {re_label}\n'
                         f'snapshot {t + 1}/{n_t}   t = {times[t]:.2f}',
                         fontsize=12, fontweight='bold')
            ax.set_aspect('equal')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            writer.grab_frame()
            print(f'  frame {n + 1}/{len(frames)} (t={times[t]:.2f})')

    plt.close(fig)
    print(f'Animation saved: {args.out}')


if __name__ == '__main__':
    main()
