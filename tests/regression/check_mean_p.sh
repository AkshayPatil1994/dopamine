#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
#
# Run a regression case and require a bounded mean pressure in its last output frame.
#
#   check_mean_p.sh <case> <np> <workroot> <bound>
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case_name="$1"; np="$2"; workroot="$3"; bound="$4"
workdir="$workroot/${case_name}_np${np}"
"$here/run_case.sh" "$case_name" "$workdir" "$np"
exec python3 "$here/check_mean_p.py" "$workdir/case.nc" "$bound"
