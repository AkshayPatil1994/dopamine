"""
Lid-driven cavity grid comparison: overlay all available .nc cases on one figure.

Discovers every *.nc file in a directory (default: script directory), extracts the
centreline profiles exactly as compare_cavity.py does, and plots them together with
one open (no facecolor) marker style per case, against the u.dat/v.dat reference.

Usage:
    python compare_grids.py [--dir DIR] [--nc PATH ...] [--u-ref PATH] [--v-ref PATH]
                            [--tavg-frac FRAC] [--out PNG]

Arguments:
    --dir       Directory to scan for .nc files (default: script directory)
    --nc        Explicit .nc paths (repeatable); overrides --dir scanning
    --u-ref     Path to u.dat reference profile (u vs y on vertical centerline)
    --v-ref     Path to v.dat reference profile (v vs x on horizontal centerline)
    --tavg-frac Fraction of time record to include in time average (default: 0.5)
    --out       Output figure file (default: cavity_grid_comparison.png)

Cavity orientation assumed (same as compare_cavity.py):
    Lid at z=1 moves in +x  →  ref y ↔ dopamine z, ref v ↔ dopamine w (U[:,2])
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from compare_cavity import (
    compute_centroids,
    extract_profile_z,
    extract_profile_x,
    load_reference_data,
    rms_error,
    time_average,
)

# Open-marker cycle: (marker, colour). Facecolors are always 'none'.
MARKERS = ['o', 's', '^', 'v', 'D', '<', '>', 'p', 'h', '*']
COLORS = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e',
          '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#7f7f7f']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    p.add_argument('--dir', default=here, help='Directory to scan for .nc files')
    p.add_argument('--nc', action='append', default=None,
                   help='Explicit .nc path (repeatable); overrides --dir scanning')
    p.add_argument('--u-ref', default='u.dat', help='Reference u profile (vertical centreline)')
    p.add_argument('--v-ref', default='v.dat', help='Reference v profile (horizontal centreline)')
    p.add_argument('--tavg-frac', type=float, default=0.5,
                   help='Fraction of time record to time-average (default 0.5)')
    p.add_argument('--out', default='cavity_grid_comparison.png', help='Output figure path')
    return p.parse_args()


def find_nc_files(args):
    """Return the list of .nc files to process, sorted by trailing integer if present."""
    if args.nc:
        files = [f for f in args.nc if os.path.isfile(f)]
        missing = [f for f in args.nc if not os.path.isfile(f)]
        for m in missing:
            print(f'  WARNING: skipping missing file {m}')
    else:
        files = sorted(glob.glob(os.path.join(args.dir, '*.nc')))

    if not files:
        sys.exit(f'ERROR: no .nc files found (dir={args.dir})')

    def sort_key(path):
        stem = os.path.splitext(os.path.basename(path))[0]
        tail = stem.rsplit('_', 1)[-1]
        return (0, int(tail)) if tail.isdigit() else (1, 0)

    return sorted(files, key=sort_key)


def read_nu(nc_path):
    """Parse nu from cavity.input next to the .nc file; None if unavailable."""
    input_path = os.path.join(os.path.dirname(os.path.abspath(nc_path)), 'cavity.input')
    if not os.path.isfile(input_path):
        return None
    for line in open(input_path):
        if 'nu' in line and '=' in line and not line.strip().startswith('!'):
            try:
                return float(line.split('=')[1].split()[0].rstrip(','))
            except Exception:
                pass
    return None


def straddling_planes(coord, target=0.5):
    """Return the two unique cell planes lying either side of target."""
    uniq = np.unique(np.round(coord, 6))
    i = int(np.argmin(np.abs(uniq - target)))
    j = i + 1 if uniq[i] < target else i - 1
    j = min(max(j, 0), len(uniq) - 1)
    return tuple(sorted((uniq[i], uniq[j])))


def midplane_mask(cy, target=0.5):
    """Mask selecting the two cell planes straddling y=target."""
    y1, y2 = straddling_planes(cy, target)
    return (np.abs(cy - y1) < 1e-6) | (np.abs(cy - y2) < 1e-6)


def process_case(nc_path, tavg_frac):
    """
    Extract the centreline profiles for one case.

    Returns a dict with sim_z/sim_u (vertical centreline), sim_x/sim_w (horizontal
    centreline), the grid dimensions and Re, or None if the file cannot be read.
    """
    import netCDF4

    print(f'\nLoading: {nc_path}')
    try:
        nc = netCDF4.Dataset(nc_path)
    except Exception as exc:
        print(f'  WARNING: cannot open ({exc}) — skipping')
        return None

    try:
        coords = nc.variables['coords'][:]
        offsets = nc.variables['offsets'][:]
        conn = nc.variables['connectivity'][:]
        times = nc.variables['time'][:]
    except KeyError as exc:
        print(f'  WARNING: missing variable {exc} — skipping')
        nc.close()
        return None

    print(f'  n_cells={offsets.shape[0]-1}, n_snapshots={len(times)}, '
          f't=[{times[0]:.3g}, {times[-1]:.3g}]')

    cx, cy, cz = compute_centroids(conn, offsets, coords)

    print('  Time-averaging velocity field...')
    U_mean = time_average(nc.variables['U'], tavg_frac)
    nc.close()

    # The y=0 / y=1 walls are no-slip, so span-averaging mixes the end-wall
    # boundary layers into the profile; restrict to the y=0.5 midplane instead.
    ymask = midplane_mask(cy)
    print(f'  Restricting to y=0.5 midplane: {ymask.sum()} of {len(cy)} cells')
    U_mean, cx, cy, cz = U_mean[ymask], cx[ymask], cy[ymask], cz[ymask]

    # Cells are not necessarily centred on 0.5 (graded meshes); take the two cell
    # planes that straddle the centre.
    x1, x2 = straddling_planes(cx)
    z1, z2 = straddling_planes(cz)

    tol = 1.0e-5  # x1/x2/z1/z2 are exact cell-plane coordinates
    z1_arr, u1 = extract_profile_z(U_mean, cx, cy, cz, x_target=x1, avg_y=True, x_tol=tol)
    z2_arr, u2 = extract_profile_z(U_mean, cx, cy, cz, x_target=x2, avg_y=True, x_tol=tol)
    x1_arr, w1 = extract_profile_x(U_mean, cx, cy, cz, z_target=z1, avg_y=True, z_tol=tol)
    x2_arr, w2 = extract_profile_x(U_mean, cx, cy, cz, z_target=z2, avg_y=True, z_tol=tol)

    if len(z1_arr) > 0 and len(z2_arr) > 0:
        sim_z = z1_arr
        sim_u = (u1 + np.interp(sim_z, z2_arr, u2, left=np.nan, right=np.nan)) / 2
    else:
        sim_z, sim_u = (z1_arr, u1) if len(z1_arr) > 0 else (z2_arr, u2)

    if len(x1_arr) > 0 and len(x2_arr) > 0:
        sim_x = x1_arr
        sim_w = (w1 + np.interp(sim_x, x2_arr, w2, left=np.nan, right=np.nan)) / 2
    else:
        sim_x, sim_w = (x1_arr, w1) if len(x1_arr) > 0 else (x2_arr, w2)

    valid_u = ~np.isnan(sim_u)
    sim_z, sim_u = sim_z[valid_u], sim_u[valid_u]
    valid_w = ~np.isnan(sim_w)
    sim_x, sim_w = sim_x[valid_w], sim_w[valid_w]

    if len(sim_z) == 0 or len(sim_x) == 0:
        print('  WARNING: empty centreline profile — skipping')
        return None

    lid_vel = 1.0
    nu = read_nu(nc_path)
    if nu is None:
        nu = 0.001
        print('  WARNING: could not parse nu from cavity.input; defaulting to 0.001')
    sim_u = sim_u / lid_vel
    sim_w = sim_w / lid_vel

    nx = len(np.unique(np.round(cx, 8)))
    ny = len(np.unique(np.round(cy, 8)))
    nz = len(np.unique(np.round(cz, 8)))
    print(f'  Grid: {nx} x {ny} x {nz} cells, Re={lid_vel/nu:.0f}')
    print(f'  u-profile: {len(sim_z)} points, w-profile: {len(sim_x)} points')

    return dict(
        name=os.path.splitext(os.path.basename(nc_path))[0],
        sim_z=sim_z, sim_u=sim_u, sim_x=sim_x, sim_w=sim_w,
        re=lid_vel / nu, nx=nx, ny=ny, nz=nz,
    )


def main():
    args = parse_args()

    try:
        import netCDF4  # noqa: F401
    except ImportError:
        sys.exit('ERROR: netCDF4 not installed. Run: pip install netCDF4')

    files = find_nc_files(args)
    print(f'Found {len(files)} .nc file(s):')
    for f in files:
        print(f'  {f}')

    # Reference data: resolve relative to cwd first, then to the first case's directory.
    base_dir = os.path.dirname(os.path.abspath(files[0]))
    u_ref_path = args.u_ref if os.path.isfile(args.u_ref) else os.path.join(base_dir, args.u_ref)
    v_ref_path = args.v_ref if os.path.isfile(args.v_ref) else os.path.join(base_dir, args.v_ref)
    ref = load_reference_data(u_ref_path, v_ref_path)

    cases = [c for c in (process_case(f, args.tavg_frac) for f in files) if c is not None]
    if not cases:
        sys.exit('ERROR: no case produced a usable profile')

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    re_vals = sorted({round(c['re']) for c in cases})
    re_txt = f'Re={re_vals[0]:.0f}' if len(re_vals) == 1 else 'Re=' + ', '.join(f'{r:.0f}' for r in re_vals)
    fig.suptitle(f'Lid-driven cavity: grid comparison  |  {re_txt}',
                 fontsize=13, fontweight='bold')

    ax_u, ax_w = axes
    ax_u.plot(ref['u_u'], ref['u_y'], 'ko', ms = 3, lw=1.2, label='Reference', zorder=1)
    ax_w.plot(ref['v_x'], ref['v_v'], 'ko', ms = 3, lw=1.2, label='Reference', zorder=1)

    print('\n--- RMS errors ---')
    for i, case in enumerate(cases):
        marker = MARKERS[i % len(MARKERS)]
        color = COLORS[i % len(COLORS)]
        rms_u = rms_error(ref['u_y'], ref['u_u'], case['sim_z'], case['sim_u'])
        rms_v = rms_error(ref['v_x'], ref['v_v'], case['sim_x'], case['sim_w'])
        case['rms_u'], case['rms_v'] = rms_u, rms_v
        label = f"{case['nx']}$^3$"
        style = dict(marker=marker, mfc='none', mec=color, mew=1.0,
                     ms=4, ls='none', zorder=3 + i)
        ax_u.plot(case['sim_u'], case['sim_z'], label=label, **style)
        ax_w.plot(case['sim_x'], case['sim_w'], label=label, **style)
        print(f"  {case['name']:<24} grid {case['nx']}x{case['ny']}x{case['nz']}   "
              f'RMS u={rms_u:.4f}   RMS w={rms_v:.4f}')

    ax_u.axhline(0.5, color='gray', lw=0.6, ls='--')
    ax_u.axvline(0.0, color='gray', lw=0.6, ls='--')
    ax_u.set_xlabel('$u/U_{lid}$', fontsize=12)
    ax_u.set_ylabel('$z/L$  (vertical position)', fontsize=12)
    ax_u.set_title('$u$-velocity along $x=0.5$', fontsize=11)
    ax_u.legend(fontsize=8)
    ax_u.set_xlim(-0.7, 1.1)
    ax_u.set_ylim(0, 1)
    ax_u.grid(True, alpha=0.3)

    ax_w.axhline(0.0, color='gray', lw=0.6, ls='--')
    ax_w.axvline(0.5, color='gray', lw=0.6, ls='--')
    ax_w.set_xlabel('$x/L$  (horizontal position)', fontsize=12)
    ax_w.set_ylabel('$w/U_{lid}$  (vertical velocity)', fontsize=12)
    ax_w.set_title('$w$-velocity along $z=0.5$', fontsize=11)
    ax_w.legend(fontsize=8)
    ax_w.set_xlim(0, 1)
    ax_w.set_ylim(-0.6, 0.6)
    ax_w.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.out, dpi=150, bbox_inches='tight')
    print(f'\nFigure saved: {args.out}')


if __name__ == '__main__':
    main()
