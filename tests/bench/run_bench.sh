#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
#
# Periodic-box timing run: generate a binary OpenFOAM hex box, run, report s/step.
#
#   run_bench.sh <N> <workdir> [nstep] [freq]
#
# Mesh is N^3 on [0,2pi]^3; dt scales as 0.32/N so CFL stays fixed. The monitor prints
# wall time per <freq> steps; the first row (warm-up, first-touch) is dropped.
# Environment: DOPAMINE_BIN (default build/dopamine), MPIRUN (default mpirun), NP (ranks, default 1), BIND, MESH_CACHE.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
n="${1:?usage: run_bench.sh <N> <workdir> [nstep] [freq]}"
workdir="${2:?usage: run_bench.sh <N> <workdir> [nstep] [freq]}"
nstep="${3:-60}"
freq="${4:-10}"
bin="${DOPAMINE_BIN:-$root/build/dopamine}"

mkdir -p "$workdir"
workdir="$(cd "$workdir" && pwd)"
cache="${MESH_CACHE:-$workdir/../meshes}"
mkdir -p "$cache"
mesh="$(cd "$cache" && pwd)/box$n/constant/polyMesh"
L=6.283185307179586
[[ -f "$mesh/boundary" ]] || python3 "$here/mk_foam_box.py" "$n" "$n" "$n" "$mesh" $L $L $L

dt="$(python3 -c "print(0.32*$L/$n/ (6.283185307179586*1.0) )")"
sed -e "s#@MESH@#$mesh#" -e "s#@DT@#$dt#" -e "s#@NSTEP@#$nstep#" -e "s#@FREQ@#$freq#" \
  "$here/bench_tgv.input" > "$workdir/bench.input"

cd "$workdir"
"${MPIRUN:-mpirun}" -np "${NP:-1}" --bind-to "${BIND:-none}" "$bin" bench.input > run.log 2>&1 || {
  echo "run failed; see $workdir/run.log" >&2; tail -20 run.log >&2; exit 1; }

grep "Setup time" run.log
awk -v n="$n" -v f="$freq" '
  /^ +[0-9]+ +[0-9.E+-]+ +[0-9.E+-]+ +/ && NF >= 8 { w[++k] = $NF }
  END {
    if (k < 2) { print "too few monitor rows"; exit 1 }
    s = 0; for (i = 2; i <= k; i++) s += w[i]
    spp = s / ((k - 1) * f)
    printf "N=%d cells=%d s/step=%.5f Mcell-steps/s=%.3f (rows %d, first dropped)\n", n, n^3, spp, n^3 / spp / 1e6, k
  }' run.log
