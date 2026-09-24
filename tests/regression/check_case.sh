#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
#
# Run a regression case and compare it against the stored reference for the same rank count
# (results depend on the partition at the ~1e-5 level, so np=1 and np=2 have separate references).
#
#   check_case.sh <case> <np> <workroot> [compare_stats.py flags...]
#
# Same environment variables as run_case.sh. REGRESSION_RTOL overrides the default tolerance.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case_name="$1"; np="$2"; workroot="$3"; shift 3
workdir="$workroot/${case_name}_np${np}"
"$here/run_case.sh" "$case_name" "$workdir" "$np"
exec python3 "$here/compare_stats.py" "$here/reference/${case_name}_np${np}.csv" "$workdir/case_stats.csv" \
  --rtol "${REGRESSION_RTOL:-1e-6}" "$@"
