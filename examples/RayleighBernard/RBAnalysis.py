#!/usr/bin/env python3
"""
analyze_snapshots.py
---------------------
Analyze RB.nc 3D snapshots from a DOPAMINE run on the structured
128x128x64 hex mesh (Lx=Ly=2, Lz=1, z wall-clustered).

Subcommands
-----------
profile   Average T over x, y and time (t >= --start-time) to build the
          mean vertical profile T(z), and compare it against
          reference_dns.dat.

              python analyze_snapshots.py profile --start-time 50

animate   Animate x/y/z-normal slices of T or |U| from --start-time to the
          end of the run and save as a .gif. Each of --x/--y/--z takes a
          LOCATION and a PARAM (T or U) and adds one subplot, e.g.:

              python analyze_snapshots.py animate --start-time 100 \\
                  --x 1.0 U --y 1.0 T --save slices.gif

          --x 1.0 U --y 1.0 T plots the x=1.0-normal slice of |U| next to
          the y=1.0-normal slice of T.
"""

import argparse
import sys

import matplotlib
import numpy as np
import netCDF4 as nc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


# --------------------------------------------------------------------------
# Structured-mesh helpers
# --------------------------------------------------------------------------

def build_structured_grid(ds: nc.Dataset, decimals: int = 4):
    """Recover the (Nx, Ny, Nz) structured indexing of the hex cells from
    their vertex connectivity, since the file stores an unstructured CSR
    mesh even though the underlying grid is a structured box.

    Returns
    -------
    cell_idx : (Nx, Ny, Nz) int array mapping structured indices -> flat cell id
    ux, uy, uz : 1D arrays of the unique cell-center coordinates along each axis
    ix, iy, iz : (n_cells,) int arrays, the structured index of each cell along
                 each axis (useful for bincount-style reductions)
    """
    coords = ds.variables["coords"][:]                       # (n_verts, 3)
    connectivity = ds.variables["connectivity"][:]
    offsets = ds.variables["offsets"][:]
    n_cells = offsets.shape[0] - 1
    stride = int(np.diff(offsets[:2])[0])
    conn = connectivity.reshape(n_cells, stride)

    centers = coords[conn].mean(axis=1)                      # (n_cells, 3)
    centers = np.round(centers, decimals)

    ux, ix = np.unique(centers[:, 0], return_inverse=True)
    uy, iy = np.unique(centers[:, 1], return_inverse=True)
    uz, iz = np.unique(centers[:, 2], return_inverse=True)

    cell_idx = np.empty((len(ux), len(uy), len(uz)), dtype=np.int64)
    cell_idx[ix, iy, iz] = np.arange(n_cells)

    if len(ux) * len(uy) * len(uz) != n_cells:
        raise RuntimeError(
            f"Recovered grid {len(ux)}x{len(uy)}x{len(uz)} does not match "
            f"n_cells={n_cells}; increase rounding tolerance in "
            f"build_structured_grid().")

    return cell_idx, ux, uy, uz, ix, iy, iz


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


# --------------------------------------------------------------------------
# 1. + 3. Mean temperature profile T(z), compared to reference_dns.dat
# --------------------------------------------------------------------------

def compute_temperature_profile(ds: nc.Dataset, start_time: float):
    cell_idx, ux, uy, uz, ix, iy, iz = build_structured_grid(ds)
    times = ds.variables["time"][:]

    t0 = int(np.searchsorted(times, start_time))
    if t0 >= len(times):
        sys.exit(f"Error: --start-time {start_time} is beyond the last "
                  f"snapshot time {times[-1]}.")

    nz = len(uz)
    counts = np.bincount(iz, minlength=nz).astype(float)     # = Nx*Ny for every z
    tsum = np.zeros(nz)

    n_used = len(times) - t0
    for t in range(t0, len(times)):
        T_t = ds.variables["T"][t, :]
        tsum += np.bincount(iz, weights=T_t, minlength=nz)

    profile = tsum / n_used / counts
    return uz, profile, times[t0], times[-1], n_used


def load_reference(path: str):
    data = np.loadtxt(path, delimiter=",", comments="#")
    z_ref, T_ref = data[:, 0], data[:, 1]
    return z_ref, T_ref


def cmd_profile(args):
    ds = nc.Dataset(args.nc)
    uz, profile, t_start, t_end, n_used = compute_temperature_profile(ds, args.start_time)
    ds.close()

    print(f"Averaged {n_used} snapshots over t in [{t_start:.3f}, {t_end:.3f}]")

    fig, ax = plt.subplots(figsize=(5, 6))
    ax.plot(profile, uz, "k-", lw=1.5, label=f"DOPAMINE")

    if args.reference:
        try:
            z_ref, T_ref = load_reference(args.reference)
            ax.plot(z_ref, T_ref, "o", ms=4, mfc="none", mec="C3",
                    label="DNS")
        except OSError:
            print(f"Warning: could not read reference file '{args.reference}', "
                  f"skipping comparison.")

    ax.set_xlabel(r"$T$")
    ax.set_ylabel(r"$z$")
    ax.set_title("Horizontally- and time-averaged temperature profile")
    ax.set_ylim(uz.min(), uz.max())
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(frameon=False, loc="lower left", fontsize=15)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=200, bbox_inches="tight")
        print(f"Figure saved to: {args.save}")
    else:
        plt.show()

    if args.csv:
        np.savetxt(args.csv, np.column_stack([uz, profile]),
                   header="z,T", delimiter=",", comments="")
        print(f"Profile data saved to: {args.csv}")


# --------------------------------------------------------------------------
# 2. Slice animation
# --------------------------------------------------------------------------

AXES = ("x", "y", "z")
PARAM_LABELS = {"T": r"$T$", "U": r"$|\mathbf{u}|$"}


def parse_slice_specs(args):
    """Turn --x/--y/--z LOC PARAM options into an ordered list of slice specs."""
    specs = []
    for axis in AXES:
        val = getattr(args, axis)
        if val is None:
            continue
        loc_str, param = val
        param = param.upper()
        if param not in ("T", "U"):
            sys.exit(f"Error: --{axis} parameter must be 'T' or 'U', got '{param}'.")
        specs.append({"axis": axis, "loc": float(loc_str), "param": param})
    if not specs:
        sys.exit("Error: at least one of --x, --y, --z must be given for 'animate'.")
    return specs


def extract_slice(field: np.ndarray, cell_idx: np.ndarray, spec: dict) -> np.ndarray:
    axis, idx = spec["axis"], spec["grid_idx"]
    if axis == "x":
        return field[cell_idx[idx, :, :]]        # (Ny, Nz)
    if axis == "y":
        return field[cell_idx[:, idx, :]]         # (Nx, Nz)
    return field[cell_idx[:, :, idx]]             # (Nx, Ny)


def cmd_animate(args):
    specs = parse_slice_specs(args)

    ds = nc.Dataset(args.nc)
    cell_idx, ux, uy, uz, ix, iy, iz = build_structured_grid(ds)
    axis_coords = {"x": ux, "y": uy, "z": uz}
    plane_coords = {"x": (uy, uz), "y": (ux, uz), "z": (ux, uy)}
    plane_labels = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}

    for spec in specs:
        coords = axis_coords[spec["axis"]]
        spec["grid_idx"] = nearest_index(coords, spec["loc"])
        spec["actual_loc"] = coords[spec["grid_idx"]]

    times = ds.variables["time"][:]
    t0 = int(np.searchsorted(times, args.start_time))
    if t0 >= len(times):
        sys.exit(f"Error: --start-time {args.start_time} is beyond the last "
                  f"snapshot time {times[-1]}.")
    frame_idx = np.arange(t0, len(times))

    need_T = any(s["param"] == "T" for s in specs)
    need_U = any(s["param"] == "U" for s in specs)

    frames = [[] for _ in specs]
    vmin = {"T": np.inf, "U": np.inf}
    vmax = {"T": -np.inf, "U": -np.inf}

    print(f"Reading {len(frame_idx)} snapshots (t in [{times[t0]:.3f}, "
          f"{times[-1]:.3f}]) for {len(specs)} slice(s)...")

    for t in frame_idx:
        T_full = ds.variables["T"][t, :] if need_T else None
        if need_U:
            U_full = ds.variables["U"][t, :, :]
            Umag_full = np.sqrt(np.sum(U_full ** 2, axis=1))
        for i, spec in enumerate(specs):
            field = T_full if spec["param"] == "T" else Umag_full
            sl = extract_slice(field, cell_idx, spec)
            frames[i].append(sl)
            p = spec["param"]
            vmin[p] = min(vmin[p], sl.min())
            vmax[p] = max(vmax[p], sl.max())
    ds.close()

    n_frames = len(frame_idx)
    n_sub = len(specs)
    fig, axes = plt.subplots(1, n_sub, figsize=(5.2 * n_sub, 4.6), squeeze=False)
    axes = axes[0]

    images = []
    for i, spec in enumerate(specs):
        ax = axes[i]
        c0, c1 = plane_coords[spec["axis"]]
        l0, l1 = plane_labels[spec["axis"]]
        p = spec["param"]
        data0 = frames[i][0]
        im = ax.imshow(data0.T, origin="lower",
                       extent=[c0.min(), c0.max(), c1.min(), c1.max()],
                       aspect="equal", cmap="inferno_r" if p == "T" else "viridis",
                       vmin=vmin[p], vmax=vmax[p])
        ax.set_xlabel(f"${l0}$")
        ax.set_ylabel(f"${l1}$")
        ax.set_title(f"{spec['axis']} = {spec['actual_loc']:.3f}, "
                     f"{PARAM_LABELS[p]}")
        fig.colorbar(im, ax=ax, shrink=0.85, label=PARAM_LABELS[p])
        images.append(im)

    suptitle = fig.suptitle(f"t = {times[frame_idx[0]]:.2f}")
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    def update(frame_i):
        for i in range(n_sub):
            images[i].set_data(frames[i][frame_i].T)
        suptitle.set_text(f"t = {times[frame_idx[frame_i]]:.2f}")
        return images + [suptitle]

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    anim.save(args.save, writer=PillowWriter(fps=args.fps))
    print(f"Animation saved to: {args.save}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def add_common_args(p):
    p.add_argument("--nc", default="RB.nc",
                   help="Path to the 3D snapshot NetCDF file (default: RB.nc)")
    p.add_argument("--start-time", type=float, default=0.0,
                   help="Skip snapshots with t < this value (default: 0.0)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_profile = sub.add_parser("profile", help="Mean T(z) profile vs. reference_dns.dat")
    add_common_args(p_profile)
    p_profile.add_argument("--reference", default="reference_dns.dat",
                           help="Reference DNS profile file (default: reference_dns.dat)")
    p_profile.add_argument("--save", default=None, help="Save figure to this path")
    p_profile.add_argument("--csv", default=None,
                           help="Also save the computed profile (z,T) to this CSV path")
    p_profile.set_defaults(func=cmd_profile)

    p_anim = sub.add_parser("animate", help="Animate T or |U| slices as a .gif")
    add_common_args(p_anim)
    p_anim.add_argument("--x", nargs=2, metavar=("LOC", "PARAM"),
                        help="x-normal slice at x=LOC of PARAM (T or U)")
    p_anim.add_argument("--y", nargs=2, metavar=("LOC", "PARAM"),
                        help="y-normal slice at y=LOC of PARAM (T or U)")
    p_anim.add_argument("--z", nargs=2, metavar=("LOC", "PARAM"),
                        help="z-normal slice at z=LOC of PARAM (T or U)")
    p_anim.add_argument("--save", default="animation.gif", help="Output .gif path")
    p_anim.add_argument("--fps", type=int, default=10, help="Animation frame rate")
    p_anim.set_defaults(func=cmd_animate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
