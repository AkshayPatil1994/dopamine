#!/usr/bin/env bash
# run_spatial_convergence.sh
# ==============================================================================
# Spatial convergence study for the Taylor-Green vortex.
#
# Sweeps mesh resolution N ∈ {8, 16, 32, 64} for four mesh types:
#   hex         — structured orthogonal hexahedra
#   tet         — unstructured Delaunay tetrahedra
#   prism       — 2D triangles extruded in z (triangular prisms)
#   skewed_hex  — structured non-orthogonal hexahedra (sheared parallelepiped)
#
# Error metric: |E_k(t_end, N) - E_k(t_end, N_ref)| where E_k is the
# volume-averaged kinetic energy read from <prefix>_stats.csv (final row).
# The finest mesh in each family (N=64 by default) is the reference.
#
# Prerequisites:
#   gmsh          — mesh generation (https://gmsh.info)
#   dopamine      — solver binary (build/dopamine or set DOPAMINE_BIN)
#   python3       — for compute_convergence.py and convergence_analysis.py
#   numpy, scipy  — Python packages
#
# Environment variables:
#   DOPAMINE_BIN  — path to dopamine executable  [default: ../../build/dopamine]
#   NP            — number of MPI ranks           [default: 1]
#   N_REF         — reference resolution          [default: 64]
#
# Usage:
#   chmod +x run_spatial_convergence.sh
#   ./run_spatial_convergence.sh
#   ./run_spatial_convergence.sh 2>&1 | tee spatial_run.log
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOPAMINE_BIN="${DOPAMINE_BIN:-${SCRIPT_DIR}/../../build/dopamine}"
NP="${NP:-1}"
N_REF="${N_REF:-64}"

MESH_TYPES=(hex tet prism skewed_hex)
N_VALUES=(8 16 32 64)

MESH_DIR="${SCRIPT_DIR}/meshes"
INPUT_TEMPLATE="${SCRIPT_DIR}/inputs/spatial_base.input"
RESULTS_DIR="${SCRIPT_DIR}/results/spatial"

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
# Mesh generation + solver runs
# --------------------------------------------------------------------------- #
for MESH_TYPE in "${MESH_TYPES[@]}"; do
  echo ""
  echo "======================================================================"
  echo " Mesh type: ${MESH_TYPE}"
  echo "======================================================================"

  for N in "${N_VALUES[@]}"; do
    RUN_ID="${MESH_TYPE}_N${N}"
    RUN_DIR="${RESULTS_DIR}/${RUN_ID}"
    MESH_FILE="${MESH_DIR}/${MESH_TYPE}_N${N}.msh"

    # ---- mesh generation ---------------------------------------------------
    if [[ ! -f "${MESH_FILE}" ]]; then
      echo "  [mesh] Generating ${MESH_TYPE} N=${N} ..."
      gmsh -3 -format msh4 -setnumber N "${N}" \
           "${MESH_DIR}/${MESH_TYPE}.geo" \
           -o "${MESH_FILE}" \
           > "${MESH_DIR}/${RUN_ID}_gmsh.log" 2>&1
      echo "         -> ${MESH_FILE}"
    else
      echo "  [mesh] Reusing existing ${MESH_FILE}"
    fi

    # ---- input file --------------------------------------------------------
    mkdir -p "${RUN_DIR}"
    sed \
      -e "s|__MESH_FILE__|${MESH_FILE}|g" \
      -e "s|__PREFIX__|${RUN_DIR}/tgv|g" \
      "${INPUT_TEMPLATE}" > "${RUN_DIR}/run.input"

    # ---- solver run --------------------------------------------------------
    echo "  [run]  ${RUN_ID} ..."
    mpirun -np "${NP}" "${DOPAMINE_BIN}" "${RUN_DIR}/run.input" \
      > "${RUN_DIR}/solver.log" 2>&1 \
      && echo "         -> done" \
      || { echo "  ERROR: solver failed — see ${RUN_DIR}/solver.log"; continue; }
  done

  # ---- compute convergence metrics for this mesh type ---------------------
  echo "  [conv] Computing convergence for ${MESH_TYPE} ..."
  python3 "${SCRIPT_DIR}/compute_convergence.py" \
    --mode spatial \
    --mesh-type "${MESH_TYPE}" \
    --n-values "${N_VALUES[@]}" \
    --n-ref "${N_REF}" \
    --results-dir "${RESULTS_DIR}" \
    --domain-length 6.283185307179586 \
    --output "${RESULTS_DIR}/${MESH_TYPE}_spatial.csv"
done

# --------------------------------------------------------------------------- #
# Print convergence summary
# --------------------------------------------------------------------------- #
echo ""
echo "======================================================================"
echo " Spatial Convergence Summary"
echo "======================================================================"
for MESH_TYPE in "${MESH_TYPES[@]}"; do
  CSV="${RESULTS_DIR}/${MESH_TYPE}_spatial.csv"
  if [[ -f "${CSV}" ]]; then
    echo ""
    echo "--- ${MESH_TYPE} ---"
    python3 "${SCRIPT_DIR}/../../tools/convergence_analysis.py" \
      "${CSV}" --mode spatial --expected-order 2
  fi
done
