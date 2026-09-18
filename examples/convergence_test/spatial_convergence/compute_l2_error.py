#!/usr/bin/env python3
"""
Compute L2 error of DOPAMINE Taylor-Green vortex solution against exact solution.

Exact solution for the 2.5D Taylor-Green vortex (w=0) on [0,2pi]^3:
  u(x,y,z,t) = u0 * sin(x) * cos(y) * cos(z) * exp(-2*nu*t)
  v(x,y,z,t) = -u0 * cos(x) * sin(y) * cos(z) * exp(-2*nu*t)
  w(x,y,z,t) = 0

Computes: L2 = sqrt( sum((U_numerical - U_exact)^2) / n_cells )
"""

import sys
import numpy as np
from pathlib import Path

def exact_taylor_green(x, y, z, t, nu, u0=1.0):
    """Compute exact Taylor-Green solution at given coordinates and time."""
    decay = np.exp(-2.0 * nu * t)

    u = u0 * np.sin(x) * np.cos(y) * np.cos(z) * decay
    v = -u0 * np.cos(x) * np.sin(y) * np.cos(z) * decay
    w = np.zeros_like(x)

    return u, v, w


def compute_l2_error(output_file, nu):
    """
    Compute L2 error from DOPAMINE NetCDF output.

    Args:
        output_file: Path to output NetCDF file (e.g., "convergence_tg_16x16.nc")
        nu: Kinematic viscosity

    Returns:
        L2 error (float)
    """

    try:
        import netCDF4 as nc
    except ImportError:
        print("ERROR: netCDF4 not available", file=sys.stderr)
        print("Install with: pip install netCDF4", file=sys.stderr)
        return None

    # Check if file exists
    if not Path(output_file).exists():
        print(f"ERROR: Output file not found: {output_file}", file=sys.stderr)
        return None

    try:
        # Open NetCDF file
        ds = nc.Dataset(output_file, 'r')

        # Get dimensions
        n_cells = ds.dimensions['n_cells'].size
        n_time = ds.dimensions['time'].size

        # Get final time step
        time_var = ds.variables['time']
        final_time = float(time_var[-1])  # Last time step

        # Read final velocity field
        # Note: NetCDF shape is (n_time, n_cells, 3) for U
        U_var = ds.variables['U']
        U_final = U_var[-1, :, :]  # Shape: (n_cells, 3)

        u_num = U_final[:, 0]  # u-component
        v_num = U_final[:, 1]  # v-component
        w_num = U_final[:, 2]  # w-component

        # Read cell coordinates (centroids)
        # Mesh is unstructured, so we need to reconstruct cell centroids
        # from vertex coordinates and connectivity
        coords = ds.variables['coords'][:]  # Shape: (n_verts, 3)
        connectivity = ds.variables['connectivity'][:]  # 0-based CSR indices
        offsets = ds.variables['offsets'][:]  # CSR row offsets

        # Compute cell centroids
        cell_centroids = np.zeros((n_cells, 3))
        for c in range(n_cells):
            idx_start = offsets[c]
            idx_end = offsets[c + 1]
            cell_verts = connectivity[idx_start:idx_end].astype(int)
            cell_coords = coords[cell_verts, :]
            cell_centroids[c, :] = np.mean(cell_coords, axis=0)

        x = cell_centroids[:, 0]
        y = cell_centroids[:, 1]
        z = cell_centroids[:, 2]

        # Compute exact solution at final time
        u_exact, v_exact, w_exact = exact_taylor_green(x, y, z, final_time, nu, u0=1.0)

        # Compute L2 error
        du = u_num - u_exact
        dv = v_num - v_exact
        dw = w_num - w_exact

        error_sq = du**2 + dv**2 + dw**2
        l2_error = np.sqrt(np.mean(error_sq))

        ds.close()

        return l2_error

    except Exception as e:
        print(f"ERROR: Failed to compute L2 error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return None


def main():
    if len(sys.argv) < 3:
        print("Usage: python compute_l2_error.py <output_file> <nu>")
        print("  output_file: DOPAMINE output NetCDF file (e.g., convergence_tg_16x16.nc)")
        print("  nu: kinematic viscosity")
        sys.exit(1)

    output_file = sys.argv[1]
    try:
        nu = float(sys.argv[2])
    except ValueError:
        print(f"ERROR: Invalid viscosity: {sys.argv[2]}")
        sys.exit(1)

    l2_error = compute_l2_error(output_file, nu)

    if l2_error is not None:
        print(f"{l2_error:.6e}")
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
