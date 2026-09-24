#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
#
# Run one regression case and leave its per-step stats CSV at <workdir>/case_stats.csv.
#
#   run_case.sh <case> <workdir> [np]
#
# Environment:
#   DOPAMINE_BIN   solver executable            (default: build/dopamine)
#   MPIRUN         MPI launcher                 (default: mpirun)
#   MPIRUN_FLAGS   extra launcher flags         (default: --bind-to none)
#   GMSH           gmsh executable              (default: gmsh on PATH)
#   MESH_CACHE     directory for generated meshes (default: <workdir>/../meshes)
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"

case_name="${1:?usage: run_case.sh <case> <workdir> [np]}"
workdir="${2:?usage: run_case.sh <case> <workdir> [np]}"
np="${3:-1}"

case_dir="$here/cases/$case_name"
[[ -f "$case_dir/case.input" ]] || { echo "unknown case: $case_name" >&2; exit 2; }

bin="${DOPAMINE_BIN:-$root/build/dopamine}"
mpirun="${MPIRUN:-mpirun}"
mpirun_flags="${MPIRUN_FLAGS:---bind-to none}"
gmsh="${GMSH:-gmsh}"

mkdir -p "$workdir"
workdir="$(cd "$workdir" && pwd)"
mesh_cache="${MESH_CACHE:-$workdir/../meshes}"
mkdir -p "$mesh_cache"
mesh_cache="$(cd "$mesh_cache" && pwd)"

# shellcheck source=/dev/null
source "$case_dir/case.env"
mesh="$mesh_cache/$case_name.msh"
if [[ ! -f "$mesh" ]]; then
  sed -e "$GEO_SED" "$root/$GEO_SRC" > "$mesh_cache/$case_name.geo"
  "$gmsh" -3 -format msh4 -o "$mesh" "$mesh_cache/$case_name.geo" > "$mesh_cache/$case_name.gmsh.log" 2>&1
fi

cp "$case_dir/case.input" "$workdir/case.input"
ln -sf "$mesh" "$workdir/mesh.msh"
rm -f "$workdir/case_stats.csv"

cd "$workdir"
# shellcheck disable=SC2086
if ! "$mpirun" -np "$np" $mpirun_flags "$bin" case.input > run.log 2>&1; then
  echo "dopamine failed for $case_name (np=$np); see $workdir/run.log" >&2
  tail -20 run.log >&2
  exit 1
fi
[[ -s case_stats.csv ]] || { echo "no case_stats.csv produced; see $workdir/run.log" >&2; exit 1; }
