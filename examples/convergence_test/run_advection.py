#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
run_advection.py
================
Comparison of the scalar advection schemes on a periodic 1D bar (unit length, uniform
velocity U = 1, one revolution, so the exact solution at t = 1 is the initial profile).

profiles     N = 100 cells, Gaussian + square pulse (ic_type = 'scalar_pulse'), every scheme.
convergence  smooth Gaussian (ic_type = 'scalar_gauss'), N = 25 ... 400, every scheme,
             fixed CFL = U dt / h = 0.2.

Each run is a normal dopamine run; the final profile is read from the line sampler
(<prefix>_sample_line_x.nc, one point per cell centre).  Results go to
<workdir>/advection.json, which plot_advection.py reads.  Needs gmsh, mpirun and netCDF4.

    run_advection.py --bin ../../build/dopamine --workdir results
"""
import argparse
import json
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np

from run_convergence import HERE, sh

SCHEMES = ["CDS", "QUICK", "KOREN_TVD", "VANLEER_TVD"]
PROFILE_N = 100
CONV_N = [25, 50, 100, 200, 400]
CFL = 0.2


def exact(x, kind):
    """Initial profile (= exact solution after one revolution)."""
    d = x - 0.25
    d -= np.round(d)
    f = np.exp(-0.5 * (d / 0.05) ** 2)
    if kind == "scalar_pulse":
        f = f + ((x >= 0.55) & (x < 0.75))
    return f


def make_mesh(args, n):
    mesh = args.workdir / "meshes" / f"adv_N{n}.msh"
    if mesh.exists():
        return mesh
    mesh.parent.mkdir(parents=True, exist_ok=True)
    geo = (HERE / "meshes" / "advection_1d.geo").read_text()
    geo, hits = re.subn(r"^N\s*=.*$", f"N = {n};", geo, flags=re.M)
    if hits != 1:
        sys.exit("advection_1d.geo: expected exactly one 'N = ...' line")
    geo_n = mesh.with_suffix(".geo")
    geo_n.write_text(geo)
    sh([args.gmsh, "-3", "-format", "msh4", str(geo_n), "-o", str(mesh)],
       log=mesh.with_suffix(".gmsh.log"))
    return mesh


def run_case(args, scheme, kind, n):
    """One revolution; returns the final profile at the cell centres."""
    from netCDF4 import Dataset
    h = 1.0 / n
    nsteps = int(round(n / CFL))                     # dt = CFL h
    dt = CFL * h
    run_dir = args.workdir / "advection" / f"{kind}_{scheme}_N{n}"
    run_dir.mkdir(parents=True, exist_ok=True)
    out = run_dir / "adv_sample_line_x.nc"
    if not args.force and out.exists():
        try:
            with Dataset(out) as d:
                if abs(d["time"][-1] - 1.0) < 1e-6:
                    return np.array(d["T"][-1])
        except OSError:
            pass
    out.unlink(missing_ok=True)
    text = (HERE / "inputs" / "advection_base.input").read_text()
    yz = 1.5 * h                                       # a cell centre in y and z
    for key, val in (("__MESH_FILE__", str(make_mesh(args, n))), ("__SCHEME__", scheme),
                     ("__DT__", f"{dt:.12g}"), ("__T_END__", f"{(nsteps - 0.1) * dt:.12g}"),
                     ("__IC__", kind), ("__N__", str(n)), ("__NSTEPS__", str(nsteps)),
                     ("__P0__", f"{h / 2:.12g}, {yz:.12g}, {yz:.12g}"),
                     ("__P1__", f"{1 - h / 2:.12g}, {yz:.12g}, {yz:.12g}")):
        text = text.replace(key, val)
    (run_dir / "run.input").write_text(text)
    t0 = time.time()
    cmd = [args.mpirun, "-np", str(args.np)] + args.mpi_flags.split() + [args.bin, "run.input"]
    sh(cmd, cwd=run_dir, log=run_dir / "solver.log")
    with Dataset(out) as d:
        t_final = float(d["time"][-1])
        prof = np.array(d["T"][-1])
    if abs(t_final - 1.0) > 1e-6:
        sys.exit(f"{run_dir}: final sample at t = {t_final}, expected 1")
    print(f"    {run_dir.name:<34s} ({time.time() - t0:.0f} s)", flush=True)
    return prof


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bin", default=os.environ.get(
        "DOPAMINE_BIN", str(HERE.parent.parent / "build" / "dopamine")))
    ap.add_argument("--workdir", default=str(HERE / "results"))
    ap.add_argument("--np", type=int, default=int(os.environ.get("NP", "1")))
    ap.add_argument("--gmsh", default=os.environ.get("GMSH", "gmsh"))
    ap.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun"))
    ap.add_argument("--mpi-flags", default=os.environ.get("MPIRUN_FLAGS", "--bind-to none"))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    args.bin = str(Path(args.bin).resolve())
    args.workdir = Path(args.workdir).resolve()
    if not os.access(args.bin, os.X_OK):
        sys.exit(f"dopamine binary not found: {args.bin}")
    if shutil.which(args.gmsh) is None:
        sys.exit(f"gmsh not found: {args.gmsh}")

    res = {"profiles": {}, "convergence": {}}
    x = (np.arange(PROFILE_N) + 0.5) / PROFILE_N
    print("== profiles ==", flush=True)
    res["profiles"] = {"N": PROFILE_N, "x": x.tolist(), "exact": exact(x, "scalar_pulse").tolist(),
                       "schemes": {s: run_case(args, s, "scalar_pulse", PROFILE_N).tolist()
                                   for s in SCHEMES}}
    print("== convergence ==", flush=True)
    for s in SCHEMES:
        rows = []
        for n in CONV_N:
            xn = (np.arange(n) + 0.5) / n
            err = run_case(args, s, "scalar_gauss", n) - exact(xn, "scalar_gauss")
            rows.append({"N": n, "h": 1.0 / n, "L2": float(np.sqrt(np.mean(err ** 2))),
                         "Linf": float(np.abs(err).max())})
        res["convergence"][s] = rows
    args.workdir.mkdir(parents=True, exist_ok=True)
    (args.workdir / "advection.json").write_text(json.dumps(res, indent=1))
    print(f"results: {args.workdir / 'advection.json'}")


if __name__ == "__main__":
    main()
