#!/usr/bin/env bash
# run_temporal_convergence.sh
# ==============================================================================
# Temporal convergence study for the Taylor-Green vortex.
#
# Uses a fixed N=32 hex mesh and sweeps dt ∈ {0.04, 0.02, 0.01, 0.005, 0.0025}
# for the RK3W time integration scheme:
#   RK3W    — Wray low-storage 3rd-order Runge-Kutta (Williamson 1980)
#
# Error metric: |E_k(t_end, dt) - E_k(t_end, dt_ref)| where E_k is the
# volume-averaged kinetic energy read from <prefix>_stats.csv (final row).
# The finest dt (0.0025) is used as the reference.
#
# Expected convergence rate: ~3rd order in dt.
#
# Prerequisites: same as run_spatial_convergence.sh
#
# Environment variables:
#   DOPAMINE_BIN  — path to dopamine executable  [default: ../../build/dopamine]
#   NP            — number of MPI ranks           [default: 1]
#   MESH_N        — hex mesh resolution           [default: 32]
#
# Usage:
#   chmod +x run_temporal_convergence.sh
#   ./run_temporal_convergence.sh
#   ./run_temporal_convergence.sh 2>&1 | tee temporal_run.log
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOPAMINE_BIN="${DOPAMINE_BIN:-${SCRIPT_DIR}/../../build/dopamine}"
NP="${NP:-1}"
MESH_N="${MESH_N:-32}"

# dt sweep — ratio 2:1 between consecutive entries
DT_VALUES=(0.04 0.02 0.01 0.005 0.0025)
DT_REF="0.0025"    # finest dt used as reference

SCHEMES=(rk3w)

MESH_DIR="${SCRIPT_DIR}/meshes"
RESULTS_DIR="${SCRIPT_DIR}/results/temporal"

# --------------------------------------------------------------------------- #
# Sanity checks
# --------------------------------------------------------------------------- #
if [[ ! -x "${DOPAMINE_BIN}" ]]; then
  echo "ERROR: dopamine binary not found at '${DOPAMINE_BIN}'"
  echo "       Set DOPAMINE_BIN or build the project first."
  exit 1
fi

if ! command -v gmsh &>/dev/null; then
  echo "ERROR: gmsh not found in PATH. Install from https://gmsh.info"
  exit 1
fi

mkdir -p "${RESULTS_DIR}"

# --------------------------------------------------------------------------- #
# Generate the single fixed mesh (hex N=MESH_N) once
# --------------------------------------------------------------------------- #
MESH_FILE="${MESH_DIR}/hex_N${MESH_N}.msh"
if [[ ! -f "${MESH_FILE}" ]]; then
  echo "[mesh] Generating hex N=${MESH_N} for temporal sweep ..."
  gmsh -3 -format msh4 -setnumber N "${MESH_N}" \
       "${MESH_DIR}/hex.geo" \
       -o "${MESH_FILE}" \
       > "${MESH_DIR}/hex_N${MESH_N}_gmsh.log" 2>&1
  echo "       -> ${MESH_FILE}"
else
  echo "[mesh] Reusing existing ${MESH_FILE}"
fi

# --------------------------------------------------------------------------- #
# Solver runs
# --------------------------------------------------------------------------- #
for SCHEME in "${SCHEMES[@]}"; do
  echo ""
  echo "======================================================================"
  echo " Scheme: ${SCHEME^^}"
  echo "======================================================================"

  INPUT_TEMPLATE="${SCRIPT_DIR}/inputs/temporal_${SCHEME}_base.input"

  for DT in "${DT_VALUES[@]}"; do
    # Replace '.' with 'p' in DT for safe directory/file names
    DT_TAG="${DT//./_}"
    RUN_ID="${SCHEME}_dt${DT_TAG}"
    RUN_DIR="${RESULTS_DIR}/${RUN_ID}"

    mkdir -p "${RUN_DIR}"

    sed \
      -e "s|__MESH_FILE__|${MESH_FILE}|g" \
      -e "s|__DT__|${DT}|g" \
      -e "s|__PREFIX__|${RUN_DIR}/tgv|g" \
      "${INPUT_TEMPLATE}" > "${RUN_DIR}/run.input"

    echo "  [run]  ${RUN_ID} (dt=${DT}) ..."
    mpirun -np "${NP}" "${DOPAMINE_BIN}" "${RUN_DIR}/run.input" \
      > "${RUN_DIR}/solver.log" 2>&1 \
      && echo "         -> done" \
      || { echo "  ERROR: solver failed — see ${RUN_DIR}/solver.log"; continue; }
  done

  # ---- compute convergence metrics for this scheme ------------------------
  echo "  [conv] Computing convergence for ${SCHEME} ..."
  python3 "${SCRIPT_DIR}/compute_convergence.py" \
    --mode temporal \
    --scheme "${SCHEME}" \
    --dt-values "${DT_VALUES[@]}" \
    --dt-ref "${DT_REF}" \
    --results-dir "${RESULTS_DIR}" \
    --output "${RESULTS_DIR}/${SCHEME}_temporal.csv"
done

# --------------------------------------------------------------------------- #
# Print convergence summary
# --------------------------------------------------------------------------- #
echo ""
echo "======================================================================"
echo " Temporal Convergence Summary"
echo "======================================================================"
for SCHEME in "${SCHEMES[@]}"; do
  CSV="${RESULTS_DIR}/${SCHEME}_temporal.csv"
  if [[ -f "${CSV}" ]]; then
    echo ""
    echo "--- ${SCHEME^^} ---"
    python3 "${SCRIPT_DIR}/../../tools/convergence_analysis.py" \
      "${CSV}" --mode temporal --expected-order 3
  fi
done
