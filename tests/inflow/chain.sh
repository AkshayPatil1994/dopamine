#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
#
# End-to-end check of the inflow-plane chain on a 32x16x16 duct:
#   export the inlet grid -> generate planes on it -> replay them in the solver (np 1 and 2),
#   stop checks, restart reproducibility, and the precursor recorder against the solver's own output.
#
#   chain.sh <workdir>      (DOPAMINE_BIN, DOPAMINE_GEN, MPIRUN, GMSH as in tests/regression)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
work="${1:?usage: chain.sh <workdir>}"
bin="${DOPAMINE_BIN:-$root/build/dopamine}"
gen="${DOPAMINE_GEN:-$root/build/dopamine-inflow-gen}"
mpirun="${MPIRUN:-mpirun}"; flags="${MPIRUN_FLAGS:---bind-to none}"; gmsh="${GMSH:-gmsh}"
mkdir -p "$work"; work="$(cd "$work" && pwd)"; cd "$work"
fail=0
ok()   { echo "  ok    $*"; }
bad()  { echo "  FAIL  $*"; fail=1; }
run()  { local np=$1; shift; $mpirun -np "$np" $flags "$bin" "$@" > "${1%.input}.log" 2>&1; }
mk()   { # mk <template> <out> KEY=VAL...
  local t=$1 o=$2; shift 2; python3 - "$t" "$o" "$@" <<'PY'
import sys
s = open(sys.argv[1]).read()
kv = dict(a.split('=', 1) for a in sys.argv[3:])
for k in ('NSTEP','YTYPE','WIND','INFLOW','EXTRA','WRITE','RECFILE','NSEG'):
    s = s.replace('@%s@' % k, kv.get(k, {'YTYPE':'periodic','WIND':'270.0','INFLOW':'','EXTRA':'','NSTEP':'1','WRITE':'100000','RECFILE':'rec.nc','NSEG':'1'}[k]).replace('\\n', '\n'))
open(sys.argv[2], 'w').write(s)
PY
}
umax() { python3 - "$1" "$2" <<'PY'
import csv, sys
with open(sys.argv[1]) as f:
    r = csv.reader(f); h = [x.strip() for x in next(r)]
    rows = {int(float(x[0])): x for x in r}
print(float(rows[int(sys.argv[2])][h.index('Umax')]))
PY
}
rel() { python3 -c "import sys; a,b=float(sys.argv[1]),float(sys.argv[2]); print(abs(a-b)/abs(a))" "$1" "$2"; }

sed 's/^Nx = .*/Nx = 32;/;s/^Ny = .*/Ny = 16;/;s/^Nz = .*/Nz = 16;/;s/^Lx = .*/Lx = 4.0;/;s/^Ly = .*/Ly = 1.0;/' "$root/examples/RayleighBernard/RB.geo" > m.geo
"$gmsh" -3 -format msh4 -o mesh.msh m.geo > gmsh.log 2>&1

echo "inflow chain: export grid, generate, replay"
mk "$here/duct_abl.input.in" export.input "INFLOW=&dopamine_inflow\n  inflow_export_grid = 'inlet_grid.txt'\n/"
run 2 export.input && grep -q "^ny 16" inlet_grid.txt && ok "inlet grid exported on 2 ranks (16 spanwise nodes, graded heights)" || bad "grid export"
mk "$here/gen.input.in" gen.input
"$gen" gen.input > gen.log 2>&1 && grep -q "worst deviation" gen.log && ok "planes generated on the exported grid, statistics verified" || bad "generator"
cp "$work/planes.nc" planes_a.nc

pl() { mk "$here/duct_abl.input.in" "$1" "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes.nc'\n  inflow_seed = $2\n/" NSTEP=$3 ${4:-}; }
pl s5_np1.input 5 1; run 1 s5_np1.input; cp case_stats.csv s5_np1.csv
pl s5_np2.input 5 1; run 2 s5_np2.input; cp case_stats.csv s5_np2.csv
pl s6_np1.input 6 1; run 1 s6_np1.input; cp case_stats.csv s6_np1.csv
a=$(umax s5_np1.csv 1); b=$(umax s5_np2.csv 1); c=$(umax s6_np1.csv 1)
[[ $(python3 -c "print(1 if $(rel $a $b) < 1e-4 else 0)") == 1 ]] && ok "inflow independent of the rank count (Umax $a vs $b)" || bad "np1 vs np2 differ: $a $b"
[[ $(python3 -c "print(1 if $(rel $a $c) > 1e-5 else 0)") == 1 ]] && ok "a different replay seed changes the inflow ($a vs $c)" || bad "seed has no effect"
grep -q "coinciding with database nodes - spanwise  100.0%, vertical  100.0%" s5_np1.log && ok "mesh-matched grid: no spatial interpolation" || bad "grid not recognised as matched"

echo "database on a different lattice than the mesh (interpolation path)"
cat > gen2.input <<'GEN'
&dopamine_abl
  wall_model = 'none'
  u_ref = 1.0
  z_ref = 1.0
  z0_m = 0.001
  wind_dir_deg = 270.0
  sem_sigma = 0.15
  sem_mean_law = 'power'
  sem_alpha = 0.15
  sem_Iu_ref = 0.1
/
&dopamine_inflow_gen
  gen_file = 'planes_off.nc'
  gen_W = 1.0
  gen_y0 = 0.01
  gen_ny = 24
  gen_nz = 20
  gen_z_min = 0.005
  gen_z_max = 0.995
  gen_n_segments = 4
  gen_n_turnovers = 1.0
/
GEN
"$gen" gen2.input > gen2.log 2>&1 && ok "planes generated on an unrelated 24 x 20 lattice" || bad "generator (off-lattice)"
mk "$here/duct_abl.input.in" off.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes_off.nc'\n/" NSTEP=100
run 1 off.input && grep -q "spanwise    0.0%, vertical    0.0%" off.log && ok "unmatched lattice: interpolated (0% node matches), run completes" || bad "off-lattice replay"
umax_off=$(umax case_stats.csv 100); [[ $(python3 -c "print(1 if 0.9 < $umax_off < 1.6 else 0)") == 1 ]] && ok "interpolated inflow is well behaved (Umax $umax_off)" || bad "off-lattice Umax $umax_off"

echo "stop checks"
mk "$here/duct_abl.input.in" neg1.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes.nc'\n/" YTYPE=slip
run 1 neg1.input && bad "non-periodic span was accepted" || { grep -q "spanwise periodicity required" neg1.log && ok "non-periodic spanwise domain stops the run" || bad "wrong message"; }
mk "$here/duct_abl.input.in" neg2.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes.nc'\n/" WIND=250.0
run 1 neg2.input && bad "wind mismatch accepted" || { grep -q "wind direction differs" neg2.log && ok "wind direction mismatch stops the run" || bad "wrong message"; }

echo "restart reproducibility"
mk "$here/duct_abl.input.in" full.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes.nc'\n  inflow_seed = 5\n/" NSTEP=200 "EXTRA=&dopamine_checkpoint\n  checkpoint_interval = 100\n/"
rm -rf case; run 1 full.input; cp case_stats.csv full.csv; cp case.nc full_case.nc
mk "$here/duct_abl.input.in" rst.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'planes.nc'\n  inflow_seed = 5\n/" NSTEP=200 "EXTRA=&dopamine_checkpoint\n  restart_enabled = .true.\n  restart_dir = 'case'\n  restart_step = 100\n/"
run 1 rst.input
a=$(umax full.csv 200); b=$(umax case_stats.csv 200)
[[ $(python3 -c "print(1 if $(rel $a $b) < 1e-4 else 0)") == 1 ]] && ok "restart reproduces the continuous run (Umax $a vs $b)" || bad "restart differs: $a $b"

echo "precursor recorder"
if python3 -c "import netCDF4" 2>/dev/null; then
  rm -rf prec && mkdir prec && cd prec && cp ../mesh.msh .
  mk "$here/precursor.input.in" p1.input NSTEP=300 WRITE=50 RECFILE=rec1.nc NSEG=1
  $mpirun -np 1 $flags "$bin" p1.input > p1.log 2>&1
  python3 "$here/check_recording.py" case.nc rec1.nc 2.0625 && ok "recorded frames equal the solver's own cell data (float32), mean removed" || bad "recording differs from the solver output"
  mk "$here/precursor.input.in" p4.input NSTEP=900 WRITE=100000 RECFILE=rec4.nc NSEG=4
  $mpirun -np 2 $flags "$bin" p4.input > p4.log 2>&1
  grep -q "removed the segment-mean spanwise-mean waveform" p4.log && ok "4 segments recorded on 2 ranks and decorrelated" || bad "4-segment recording"
  cd ..
  cp prec/rec4.nc rec4.nc
  mk "$here/duct_abl.input.in" succ.input "INFLOW=&dopamine_inflow\n  inflow_source = 'planes'\n  inflow_file = 'rec4.nc'\n/" NSTEP=100
  run 2 succ.input && grep -q "coinciding with database nodes - spanwise  100.0%" succ.log && ok "successor driven by the precursor recording (matched lattice, 2 ranks)" || bad "successor from recording"
else
  echo "  skip  python netCDF4 not available: recorder cross-check skipped"
fi
exit $fail
