#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Compare two dopamine *_stats.csv files step by step.

Bitwise agreement is not expected between CPU/GPU builds or rank counts (reduction order,
AMG hierarchy), so integral quantities are compared with a relative tolerance, and the
divergence columns are only required to stay below an absolute ceiling in the test run.

Usage: compare_stats.py reference.csv test.csv [--rtol 1e-6] [--div-max 1e-9]
Exit status 0 on pass, 1 on mismatch.
"""
import argparse
import csv
import sys

# Integral quantities that must track the reference.
COMPARED = ["time", "E", "Umax", "Ub", "tauw_mean", "nu_sgs_mean"]
# Face-flux divergence: checked against an absolute ceiling, not against the reference. (The
# *_cell columns are the cell-centred divergence, O(h^2) by construction, so not a check.)
DIVERGENCE = ["max_abs_divU"]


def load(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {int(r["step"]): {k: float(v) for k, v in r.items() if k != "step"} for r in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("test")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="absolute floor for quantities that are ~0 (e.g. nu_sgs without SGS)")
    ap.add_argument("--div-max", type=float, default=1e-9)
    args = ap.parse_args()

    ref, tst = load(args.reference), load(args.test)
    missing = sorted(set(ref) - set(tst))
    if missing:
        print(f"FAIL: test run is missing steps {missing[:5]}{'...' if len(missing) > 5 else ''}")
        return 1

    failures, worst = [], {}
    for step in sorted(ref):
        for col in COMPARED:
            if col not in ref[step]:
                continue
            a, b = ref[step][col], tst[step][col]
            err = abs(a - b) / max(abs(a), args.atol / args.rtol)
            if err > worst.get(col, (0.0, 0))[0]:
                worst[col] = (err, step)
            if err > args.rtol:
                failures.append(f"step {step} {col}: ref={a:.12e} test={b:.12e} rel={err:.2e}")
        for col in DIVERGENCE:
            if col in tst[step] and tst[step][col] > args.div_max:
                failures.append(f"step {step} {col}={tst[step][col]:.3e} exceeds {args.div_max:.1e}")

    for col, (err, step) in sorted(worst.items()):
        print(f"  {col:12s} max rel diff {err:.2e} (step {step})")
    if failures:
        print(f"FAIL: {len(failures)} mismatches (rtol={args.rtol:g})")
        for line in failures[:20]:
            print("  " + line)
        return 1
    print(f"PASS: {len(ref)} steps within rtol={args.rtol:g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
