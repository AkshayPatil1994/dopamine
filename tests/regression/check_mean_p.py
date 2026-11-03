#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Fail if the mean pressure in the last frame of <case.nc> exceeds a bound (exit 77 if netCDF4 is missing)."""
import sys

try:
    import netCDF4
    import numpy as np
except ImportError:
    print("SKIP: python netCDF4/numpy not available")
    sys.exit(77)

path, bound = sys.argv[1], float(sys.argv[2])
with netCDF4.Dataset(path) as nc:
    p = np.array(nc["p"][-1])
    t = float(nc["time"][-1])
mean_p = float(p.mean())
print(f"t = {t:.3f}  mean p = {mean_p:+.4e}  (bound {bound:g})")
sys.exit(0 if abs(mean_p) <= bound else 1)
