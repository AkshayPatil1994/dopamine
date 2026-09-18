#!/usr/bin/env python3
"""
compute_convergence.py
======================
Build a convergence CSV (compatible with tools/convergence_analysis.py) from
dopamine run directories using statistics written to solver.log.

Error metric
------------
The volume-averaged kinetic energy ``E_k`` from the final statistics row in
``solver.log`` is used as the convergence observable:

    error(N)  = | E_k(N, t_end)  -  E_k(N_ref, t_end) |   [spatial mode]
    error(dt) = | E_k(dt, t_end) -  E_k(dt_ref, t_end) |  [temporal mode]

This avoids mesh-interpolation entirely and works for any cell type (hex, tet,
prism, skewed).  Kinetic energy is a smooth global scalar, so it converges at
the scheme's design order (2nd for CDS advection, 3rd for RK3W time
integration) provided the flow remains smooth over [0, t_end].

Output CSV columns (matching convergence_analysis.py)
------------------------------------------------------
    resolution, nx, ny, nz, dx, dt, l2_error, reference

Usage
-----
Spatial convergence (called by run_spatial_convergence.sh):

    python3 compute_convergence.py \\
        --mode spatial \\
        --mesh-type hex \\
        --n-values 8 16 32 64 \\
        --n-ref 64 \\
        --results-dir results/spatial \\
        --domain-length 6.283185307179586 \\
        --output results/spatial/hex_spatial.csv

Temporal convergence (called by run_temporal_convergence.sh):

    python3 compute_convergence.py \\
        --mode temporal \\
        --scheme rk3w \\
        --dt-values 0.04 0.02 0.01 0.005 0.0025 \\
        --dt-ref 0.0025 \\
        --results-dir results/temporal \\
        --output results/temporal/rk3w_temporal.csv
"""

import argparse
import csv
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_final_ke_from_solver_log(solver_log: Path) -> float:
    """Return kinetic energy E from the last stats line in solver.log.

    Stats line format from mod_stats.f90 (rank 0 screen output):
      step dt time E eps max|U| max|divU| <|divU|>V [optional output tag]

    We parse lines that start with an integer step and extract scientific-
    notation floats. Depending on whether dt is printed as '-' or numeric,
    the float layout is either:
      [time, E, eps, umax, max_div, abs_div]           (6 floats)
      [dt, time, E, eps, umax, max_div, abs_div]       (7 floats)
    """
    if not solver_log.exists():
        raise FileNotFoundError(f"solver.log not found: {solver_log}")

    float_pat = re.compile(r"[-+]?\d+(?:\.\d*)?[Ee][-+]?\d+")
    step_pat = re.compile(r"^\s*\d+")

    last_e = None
    with open(solver_log, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not step_pat.match(line):
                continue

            vals = [float(x) for x in float_pat.findall(line)]
            if len(vals) >= 7:
                # [dt, time, E, eps, umax, max_div, abs_div]
                last_e = vals[2]
            elif len(vals) >= 6:
                # [time, E, eps, umax, max_div, abs_div]
                last_e = vals[1]

    if last_e is None:
        raise ValueError(
            f"No parseable stats lines found in {solver_log}. "
            "Ensure solver_stats_frequency > 0 and rank-0 output is captured."
        )

    return last_e


def find_solver_log(run_dir: Path) -> Path:
    """Locate solver.log produced by dopamine in run_dir."""
    candidates = list(run_dir.glob("solver.log"))
    if not candidates:
        candidates = list(run_dir.glob("**/solver.log"))
    if not candidates:
        raise FileNotFoundError(
            f"No solver.log found in {run_dir}. "
            "Check that the solver ran successfully and logging was redirected."
        )
    if len(candidates) > 1:
        # Pick the most recently modified one
        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


# ---------------------------------------------------------------------------
# Spatial convergence
# ---------------------------------------------------------------------------

def run_spatial(args) -> None:
    results_dir = Path(args.results_dir)
    L = args.domain_length
    n_ref = int(args.n_ref)
    n_values = [int(n) for n in args.n_values]

    # Resolve reference kinetic energy
    ref_dir = results_dir / f"{args.mesh_type}_N{n_ref}"
    ref_log = find_solver_log(ref_dir)
    e_ref = read_final_ke_from_solver_log(ref_log)
    print(f"  Reference: {args.mesh_type} N={n_ref}  E_ref = {e_ref:.8e}")

    rows = []
    for n in n_values:
        run_dir = results_dir / f"{args.mesh_type}_N{n}"
        try:
            solver_log = find_solver_log(run_dir)
            e = read_final_ke_from_solver_log(solver_log)
        except (FileNotFoundError, ValueError) as exc:
            print(f"  WARNING: skipping N={n}: {exc}", file=sys.stderr)
            continue

        error = abs(e - e_ref)
        dx = L / n
        # nx = ny = nz = n for isotropic mesh families
        rows.append({
            "resolution": n,
            "nx": n,
            "ny": n,
            "nz": n,
            "dx": f"{dx:.6e}",
            "dt": "0.0",
            "l2_error": f"{error:.6e}",
            "reference": "reference" if n == n_ref else "no",
        })
        print(f"  N={n:3d}  dx={dx:.4f}  E={e:.8e}  error={error:.4e}")

    _write_csv(args.output, rows)
    print(f"  Written: {args.output}")


# ---------------------------------------------------------------------------
# Temporal convergence
# ---------------------------------------------------------------------------

def run_temporal(args) -> None:
    results_dir = Path(args.results_dir)
    dt_values = [float(dt) for dt in args.dt_values]
    dt_ref = float(args.dt_ref)

    def dt_to_tag(dt: float) -> str:
        return str(dt).replace(".", "_")

    # Resolve reference kinetic energy
    ref_run_id = f"{args.scheme}_dt{dt_to_tag(dt_ref)}"
    ref_dir = results_dir / ref_run_id
    ref_log = find_solver_log(ref_dir)
    e_ref = read_final_ke_from_solver_log(ref_log)
    print(f"  Reference: {args.scheme} dt={dt_ref}  E_ref = {e_ref:.8e}")

    rows = []
    for dt in dt_values:
        run_id = f"{args.scheme}_dt{dt_to_tag(dt)}"
        run_dir = results_dir / run_id
        try:
            solver_log = find_solver_log(run_dir)
            e = read_final_ke_from_solver_log(solver_log)
        except (FileNotFoundError, ValueError) as exc:
            print(f"  WARNING: skipping dt={dt}: {exc}", file=sys.stderr)
            continue

        error = abs(e - e_ref)
        # convergence_analysis.py expects an integer resolution column.
        # Use inverse timestep as a temporal resolution proxy.
        temporal_resolution = int(round(1.0 / dt))
        rows.append({
            "resolution": temporal_resolution,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "dx": "0.0",
            "dt": f"{dt:.6e}",
            "l2_error": f"{error:.6e}",
            "reference": "reference" if dt == dt_ref else "no",
        })
        print(f"  dt={dt:.5f}  E={e:.8e}  error={error:.4e}")

    _write_csv(args.output, rows)
    print(f"  Written: {args.output}")


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def _write_csv(output_path: str, rows: list) -> None:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["resolution", "nx", "ny", "nz", "dx", "dt",
                  "l2_error", "reference"]
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build convergence CSVs from dopamine solver.log files."
    )

    # Compatibility shim: some wrapper scripts may accidentally forward
    # convergence_analysis.py arguments (e.g., positional results file,
    # --expected-order). Accept and ignore them to avoid hard failures.
    parser.add_argument("results_file", nargs="?", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--expected-order", default=None,
                        help=argparse.SUPPRESS)

    # ---- common arguments ---
    parser.add_argument("--mode", choices=["spatial", "temporal"],
                        required=True)
    parser.add_argument("--results-dir", default=None,
                        help="Base directory containing per-run subdirectories. "
                            "Default: <script_dir>/results/{spatial|temporal}")
    parser.add_argument("--output", default=None,
                        help="Output CSV path. "
                            "Default: <results_dir>/{mesh_type}_spatial.csv or "
                            "<results_dir>/{scheme}_temporal.csv")

    # ---- spatial-specific ---
    parser.add_argument("--mesh-type",
                        help="Mesh family name (hex, tet, prism, skewed_hex).")
    parser.add_argument("--n-values", nargs="+", type=str,
                        help="Resolution values (e.g. 8 16 32 64).")
    parser.add_argument("--n-ref", type=str, default="64",
                        help="Reference resolution (finest mesh).")
    parser.add_argument("--domain-length", type=float,
                        default=6.283185307179586,
                        help="Domain side length L (default: 2*pi).")

    # ---- temporal-specific ---
    parser.add_argument("--scheme",
                        help="Time scheme name (rk3w).")
    parser.add_argument("--dt-values", nargs="+", type=str,
                        help="dt values (e.g. 0.04 0.02 0.01 0.005 0.0025).")
    parser.add_argument("--dt-ref", type=str,
                        help="Reference dt (finest time step).")

    args, unknown = parser.parse_known_args()
    if unknown:
        print("  WARNING: ignoring unknown arguments:", " ".join(unknown),
              file=sys.stderr)

    # Resolve defaults lazily so mode-specific naming can be used.
    script_dir = Path(__file__).resolve().parent
    if args.results_dir is None:
        default_subdir = "spatial" if args.mode == "spatial" else "temporal"
        args.results_dir = str(script_dir / "results" / default_subdir)

    if args.mode == "spatial":
        if not args.mesh_type or not args.n_values:
            parser.error("--mode spatial requires --mesh-type and --n-values")
        if args.output is None:
            args.output = str(Path(args.results_dir) / f"{args.mesh_type}_spatial.csv")
        run_spatial(args)
    elif args.mode == "temporal":
        if not args.scheme or not args.dt_values or not args.dt_ref:
            parser.error(
                "--mode temporal requires --scheme, --dt-values, and --dt-ref")
        if args.output is None:
            args.output = str(Path(args.results_dir) / f"{args.scheme}_temporal.csv")
        run_temporal(args)


if __name__ == "__main__":
    main()
