#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
run_convergence.py
==================
Spatial and temporal convergence study on the 3D Taylor-Green vortex.

Spatial: hex and warped-hex families, N = 8, 16, 32, 64, one fixed dt.
Temporal: hex N = 16, 32, 64, dt halved from 0.2 to 0.003125 at fixed t_end, plus a
inviscid energy-drift check.

Every run is a normal dopamine run in <workdir>/<study>/<run>/ and the observables are
read from its <prefix>_stats.csv: E (volume-mean kinetic energy) and eps (viscous
dissipation rate).  Results go to <workdir>/results.json, which plot_convergence.py reads.

    run_convergence.py all --bin ../../build/dopamine --workdir results
    run_convergence.py spatial --np 2
    run_convergence.py temporal --force        # redo runs that already finished

Needs: gmsh on PATH (or --gmsh), mpirun on PATH (or --mpirun), numpy not required.
"""
import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2.0 * math.pi

NU = 0.01                       # Re = 100 for U = L = 1: coarse meshes are in the asymptotic range
SPATIAL_N = [8, 16, 32, 64]
SPATIAL_DT = 4.0e-3             # same dt for every mesh: its error is mesh-independent and cancels
SPATIAL_TEND = 1.0
FAMILIES = ["hex", "warped"]
INITIAL_DT = 1.0e-4             # one tiny step: eps(t ~ 0) against the analytic 0.75*nu

TEMPORAL_N = [16, 32, 64]
TEMPORAL_DT = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
CFL_MAX_TEST = 0.9              # largest U dt / h attempted (U ~ 1)
TEMPORAL_TEND = 1.0


def sh(cmd, cwd=None, log=None):
    """Run cmd; on failure show the tail of the log and exit."""
    with open(log, "w") if log else open(os.devnull, "w") as fh:
        rc = subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        tail = ""
        if log:
            tail = "".join(open(log, errors="replace").readlines()[-6:])
        raise RuntimeError(f"command failed ({rc}): {' '.join(map(str, cmd))}\n{tail}")


def make_mesh(args, family, n, mesh_dir):
    """Return the path of the (cached) mesh for a family and resolution."""
    mesh_dir.mkdir(parents=True, exist_ok=True)
    hexmesh = mesh_dir / f"hex_N{n}.msh"
    if not hexmesh.exists():
        # `gmsh -setnumber N` does not override DefineNumber in hex.geo, so write N into a copy
        geo = (HERE / "meshes" / "hex.geo").read_text()
        geo, hits = re.subn(r"^N\s*=.*$", f"N = {n};", geo, flags=re.M)
        if hits != 1:
            sys.exit("hex.geo: expected exactly one 'N = ...' line")
        geo_n = mesh_dir / f"hex_N{n}.geo"
        geo_n.write_text(geo)
        sh([args.gmsh, "-3", "-format", "msh4", str(geo_n), "-o", str(hexmesh)],
           log=mesh_dir / f"hex_N{n}.gmsh.log")
    if family == "hex":
        return hexmesh
    warped = mesh_dir / f"warped_N{n}.msh"
    if not warped.exists():
        sh([sys.executable, str(HERE / "warp_mesh.py"), str(hexmesh), str(warped)])
    return warped


def read_row(run_dir, t_target, dt):
    """Stats row whose time is closest to t_target -> (time, E, eps, nu_sgs_over_nu)."""
    path = run_dir / "tgv_stats.csv"
    best = None
    with open(path) as fh:
        for r in csv.DictReader(fh, skipinitialspace=True):
            t = float(r["time"])
            if best is None or abs(t - t_target) < abs(best[0] - t_target):
                best = (t, float(r["E"]), float(r["eps"]))
    if best is None or abs(best[0] - t_target) > 1e-8 * max(1.0, t_target):
        raise RuntimeError(f"{path}: no stats row at t={t_target} (closest {best})")
    return best


def run_case(args, mesh, run_dir, dt, t_end, nu=NU):
    """Run dopamine once (skipped if the run finished earlier); return its stats row."""
    run_dir.mkdir(parents=True, exist_ok=True)
    stats = run_dir / "tgv_stats.csv"
    if not args.force and stats.exists():
        try:
            return read_row(run_dir, t_end, dt)
        except RuntimeError:
            pass                                        # incomplete run: redo it
    text = (HERE / "inputs" / "tgv_base.input").read_text()
    for key, val in (("__MESH_FILE__", str(mesh)), ("__DT__", f"{dt:.10g}"),
                     ("__T_END__", f"{t_end - 0.1 * dt:.12g}"), ("__NU__", f"{nu:.10g}")):
        text = text.replace(key, val)
    (run_dir / "run.input").write_text(text)
    stats.unlink(missing_ok=True)
    # t_end - dt/10 above makes the last step land on t_end whatever the rounding of t += dt
    cmd = [args.mpirun, "-np", str(args.np)] + args.mpi_flags.split() + \
          [args.bin, "run.input"]
    t0 = time.time()
    sh(cmd, cwd=run_dir, log=run_dir / "solver.log")
    row = read_row(run_dir, t_end, dt)
    print(f"    {run_dir.name:<22s} t={row[0]:.5f}  E={row[1]:.10e}  "
          f"eps={row[2]:.8e}  ({time.time() - t0:.0f} s)", flush=True)
    return row


def cell_count(run_dir):
    for line in open(run_dir / "solver.log", errors="replace"):
        m = re.search(r"Read mesh:\s+(\d+) cells", line)
        if m:
            return int(m.group(1))
    return None


def study_spatial(args, res):
    print("== spatial ==", flush=True)
    out = res.setdefault("spatial", {})
    for fam in FAMILIES:
        rows = []
        for n in SPATIAL_N:
            mesh = make_mesh(args, fam, n, args.workdir / "meshes")
            d = args.workdir / "spatial" / f"{fam}_N{n}"
            _, e, eps = run_case(args, mesh, d, SPATIAL_DT, SPATIAL_TEND)
            d0 = args.workdir / "spatial" / f"{fam}_N{n}_t0"
            # t_end = dt: a single step, so eps is the initial dissipation rate
            _, _, eps0 = run_case(args, mesh, d0, INITIAL_DT, INITIAL_DT)
            rows.append({"N": n, "h": TWO_PI / n, "cells": cell_count(d),
                         "E": e, "eps": eps, "eps0": eps0})
        out[fam] = rows
    out["meta"] = {"nu": NU, "dt": SPATIAL_DT, "t_end": SPATIAL_TEND,
                   "initial_dt": INITIAL_DT, "eps0_exact": 0.75 * NU,
                   "E0_exact": 0.125}


def study_temporal(args, res):
    print("== temporal ==", flush=True)
    out = {"nu": NU, "t_end": TEMPORAL_TEND, "meshes": {}}
    for n in TEMPORAL_N:
        mesh = make_mesh(args, "hex", n, args.workdir / "meshes")
        rows = []
        for dt in TEMPORAL_DT:
            if dt * 1.0 / (TWO_PI / n) > CFL_MAX_TEST:       # beyond RK3W + central stability
                continue
            d = args.workdir / "temporal" / f"hex_N{n}_dt{dt:g}"
            _, e, eps = run_case(args, mesh, d, dt, TEMPORAL_TEND)
            rows.append({"dt": dt, "E": e, "eps": eps})
        out["meshes"][str(n)] = rows

    # Inviscid energy drift: E should stay at 1/8 when nu -> 0; what it loses is numerical.
    print("  inviscid drift", flush=True)
    drift = []
    for n in TEMPORAL_N:
        mesh = make_mesh(args, "hex", n, args.workdir / "meshes")
        for dt in (0.02, 0.01, 0.005):
            d = args.workdir / "temporal" / f"inviscid_N{n}_dt{dt:g}"
            _, e, _ = run_case(args, mesh, d, dt, 0.2, nu=1.0e-8)
            drift.append({"N": n, "dt": dt, "E": e})
    out["inviscid"] = {"t_end": 0.2, "nu": 1.0e-8, "runs": drift}
    res["temporal"] = out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study", choices=["spatial", "temporal", "all"])
    ap.add_argument("--bin", default=os.environ.get(
        "DOPAMINE_BIN", str(HERE.parent.parent / "build" / "dopamine")))
    ap.add_argument("--workdir", default=str(HERE / "results"))
    ap.add_argument("--np", type=int, default=int(os.environ.get("NP", "1")))
    ap.add_argument("--gmsh", default=os.environ.get("GMSH", "gmsh"))
    ap.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun"))
    ap.add_argument("--mpi-flags", default=os.environ.get("MPIRUN_FLAGS", "--bind-to none"))
    ap.add_argument("--force", action="store_true", help="rerun finished cases")
    args = ap.parse_args()

    args.bin = str(Path(args.bin).resolve())
    args.workdir = Path(args.workdir).resolve()
    if not os.access(args.bin, os.X_OK):
        sys.exit(f"dopamine binary not found: {args.bin} (set --bin or DOPAMINE_BIN)")
    if shutil.which(args.gmsh) is None:
        sys.exit(f"gmsh not found: {args.gmsh}")
    args.workdir.mkdir(parents=True, exist_ok=True)

    path = args.workdir / "results.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    if args.study in ("spatial", "all"):
        study_spatial(args, res)
        path.write_text(json.dumps(res, indent=1))
    if args.study in ("temporal", "all"):
        study_temporal(args, res)
        path.write_text(json.dumps(res, indent=1))
    print(f"results: {path}")


if __name__ == "__main__":
    main()
