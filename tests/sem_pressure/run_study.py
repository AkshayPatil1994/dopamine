#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Inlet pressure-fluctuation study for the two synthetic-eddy kernels ('jarrin', 'dfsem').

A box (periodic y, slip z, SEM ablInlet at x0, outlet at xL) is driven by each kernel for several
flow-through times; pressure and velocity are sampled on lines along x at cell centres. For every
station the script reports the RMS of the pressure fluctuation, the RMS velocity, and their ratio

    Cp = p_rms / (u_rms^2),   u_rms^2 = <u'_i u'_i> / 3,

(for developed turbulence Cp is O(0.5)). The spike in the first cells at the inlet is quantified as
Cp(first 3 cells) / Cp(reference band), and the mean-flow acceleration (sum of U along the line) is
reported to confirm that the inflow does not distort the mean.

    run_study.py --workdir W [--np 4] [--cases iso,shear,abl] [--kernels jarrin,dfsem] [--nx 128 ...]

Needs the dopamine binary (DOPAMINE_BIN, default build/dopamine), mpirun, python3 + numpy + netCDF4.
"""
import argparse, os, subprocess, sys, textwrap, glob
import numpy as np
import netCDF4

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, '..', '..'))

CASES = {
    # uniform mean, 10 % isotropic turbulence: isolates the eddy kernel from mean shear and anisotropy
    'iso':   dict(profile=[(0.0, 1.0, 0.01, 0.01, 0.01, 0.0), (1.0, 1.0, 0.01, 0.01, 0.01, 0.0)]),
    # uniform mean with anisotropic stresses and a shear stress: needs the principal-axis rotation
    'shear': dict(profile=[(0.0, 1.0, 0.016, 0.009, 0.004, -0.004), (1.0, 1.0, 0.016, 0.009, 0.004, -0.004)]),
    # analytic neutral log law and Panofsky-Dutton stresses (the default when no profile is given)
    'abl':   dict(profile=None),
}


def write_input(path, a, kernel, case, np_):
    prof = CASES[case]['profile']
    if prof is not None:  # stress columns scaled by --stress-scale (turbulence intensity squared); mean U = --uref
        prof = [(z, a.uref, *(a.stress_scale * v for v in r)) for (z, U, *r) in prof]
    prof_line = ''
    if prof is not None:
        with open(os.path.join(os.path.dirname(path), 'profile.csv'), 'w') as f:
            f.write('# z U uu vv ww uw\n')
            for r in prof:
                f.write(' '.join(f'{v:.6g}' for v in r) + '\n')
        prof_line = "  sem_profile_file = 'profile.csv'\n"
    a.ztype = 'periodic' if a.zperiodic else 'slip'
    a.numerics = a.numerics.replace(';', '\n  ')
    dx, dy, dz = a.lx / a.nx, a.ly / a.ny, a.lz / a.nz
    nl = 0
    lines = ''
    for (fy, fz) in [(0.2, 0.3), (0.5, 0.5), (0.8, 0.7), (0.35, 0.8), (0.65, 0.2), (0.1, 0.6)]:
        nl += 1
        y = (int(fy * a.ny) + 0.5) * dy
        z = (int(fz * a.nz) + 0.5) * dz
        lines += (f"  line_name({nl}) = 'ln{nl}'\n  line_p0({nl},1) = {0.5*dx:.9g}\n  line_p0({nl},2) = {y:.9g}\n"
                  f"  line_p0({nl},3) = {z:.9g}\n  line_p1({nl},1) = {a.lx-0.5*dx:.9g}\n  line_p1({nl},2) = {y:.9g}\n"
                  f"  line_p1({nl},3) = {z:.9g}\n  line_npts({nl}) = {a.nx}\n  line_write_frequency({nl}) = {a.every}\n")
    txt = f"""\
&dopamine_domain
  mesh_file = 'mesh.msh'
  mesh_type = 'gmsh'
/
&dopamine_numerics
  advection_scheme = '{a.scheme}'
  time_integrator  = 'RK3W'
  use_sgs          = .true.
  sgs_model        = 'VREMAN'
  adapt_dt         = .true.
  dt               = {a.dt}
  cfl_max          = {a.cfl}
  t_end            = {a.t_end}
  p_tol            = 1.0e-9
  p_max_iter       = 500
{a.numerics}
/
&dopamine_physics
  nu = {a.nu}
/
&dopamine_abl
  wall_model = 'none'
  u_ref = {a.uref}
  z_ref = 0.5
  z0_m  = {a.z0}
  wind_dir_deg = 270.0
  use_sem = .true.
  sem_type = '{kernel}'
  n_eddies = {a.neddies}
  sem_sigma = {a.sigma}
{prof_line}/
&dopamine_bc
  n_patches = 6
  patch_name(1) = 'x0'
  patch_type(1) = 'ablInlet'
  patch_name(2) = 'xL'
  patch_type(2) = 'outlet'
  patch_p(2)    = 0.0
  patch_name(3) = 'y0'
  patch_type(3) = 'periodic'
  patch_name(4) = 'yL'
  patch_type(4) = 'periodic'
  patch_name(5) = 'z0'
  patch_type(5) = '{a.ztype}'
  patch_name(6) = 'zL'
  patch_type(6) = '{a.ztype}'
/
&dopamine_initial_condition
  ic_type = 'uniform'
  ic_uref = {a.uref}
  noise_level = 0.0
/
&dopamine_output
  prefix         = 'case'
  write_interval = 100000000
/
&dopamine_stats
  solver_stats_frequency = 50
/
&dopamine_sampling_lines
  n_lines = {nl}
{lines}/
"""
    open(path, 'w').write(txt)


def run_case(a, kernel, case):
    wd = os.path.join(a.workdir, f'{case}_{kernel}_np{a.np}{a.tag}')
    os.makedirs(wd, exist_ok=True)
    mesh = os.path.join(a.workdir, f'box_{a.nx}x{a.ny}x{a.nz}.msh')
    if not os.path.exists(mesh):
        subprocess.check_call([sys.executable, os.path.join(here, 'mk_box.py'), str(a.nx), str(a.ny), str(a.nz),
                               str(a.lx), str(a.ly), str(a.lz), mesh])
    link = os.path.join(wd, 'mesh.msh')
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(mesh, link)
    write_input(os.path.join(wd, 'case.input'), a, kernel, case, a.np)
    if not a.reuse or not glob.glob(os.path.join(wd, 'case_sample_line_ln1*.nc')):
        for f in glob.glob(os.path.join(wd, 'case_sample_line_*.nc')) + glob.glob(os.path.join(wd, 'case_stats.csv')):
            os.remove(f)
        binary = os.environ.get('DOPAMINE_BIN', os.path.join(root, 'build', 'dopamine'))
        with open(os.path.join(wd, 'run.log'), 'w') as log:
            env = dict(os.environ, OMP_NUM_THREADS='1')  # one thread per rank: ranks already fill the cores
            rc = subprocess.call([os.environ.get('MPIRUN', 'mpirun'), '-np', str(a.np), '--bind-to', 'none',
                                  binary, 'case.input'], cwd=wd, stdout=log, stderr=subprocess.STDOUT, env=env)
        if rc != 0:
            print(open(os.path.join(wd, 'run.log')).read()[-2000:])
            raise SystemExit(f'run failed: {wd}')
    return wd


def analyse(a, wd):
    """Return x, p_rms(x), u_rms(x), mean U(x) averaged over the sampled lines after the spin-up."""
    pv, uv, um = [], [], []
    for f in sorted(glob.glob(os.path.join(wd, 'case_sample_line_ln*.nc'))):
        ds = netCDF4.Dataset(f)
        t = np.array(ds['time'][:])
        x = np.array(ds['coords'][:])
        p = np.array(ds['p'][:])   # (time, point) or (point, time) depending on writer
        U = np.array(ds['U'][:])
        if p.shape[0] != len(t):
            p = p.T
        if U.shape[0] != len(t):
            U = np.moveaxis(U, -1, 0)
        if U.shape[-1] != 3 and U.shape[1] == 3:
            U = np.moveaxis(U, 1, 2)
        sel = t >= a.t_spin
        p, U = p[sel], U[sel]
        xs = x[0] if x.shape[0] == 3 else x[:, 0]
        pv.append(p.var(axis=0)); um.append(U.mean(axis=0))
        uv.append(U.var(axis=0).sum(axis=-1) / 3.0)
        ds.close()
    return xs, np.sqrt(np.mean(pv, axis=0)), np.sqrt(np.mean(uv, axis=0)), np.mean(um, axis=0)[:, 0], int(sel.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--np', type=int, default=4)
    ap.add_argument('--cases', default='iso,shear,abl')
    ap.add_argument('--kernels', default='jarrin,dfsem')
    ap.add_argument('--nx', type=int, default=80); ap.add_argument('--ny', type=int, default=32)
    ap.add_argument('--nz', type=int, default=32)
    ap.add_argument('--lx', type=float, default=2.5); ap.add_argument('--ly', type=float, default=1.0)
    ap.add_argument('--lz', type=float, default=1.0)
    ap.add_argument('--sigma', type=float, default=0.15); ap.add_argument('--neddies', type=int, default=400)
    ap.add_argument('--scheme', default='KEP'); ap.add_argument('--nu', type=float, default=1e-4)
    ap.add_argument('--z0', type=float, default=0.005)
    ap.add_argument('--dt', type=float, default=5e-3); ap.add_argument('--cfl', type=float, default=0.4)
    ap.add_argument('--t-end', dest='t_end', type=float, default=9.0)
    ap.add_argument('--t-spin', dest='t_spin', type=float, default=2.5)
    ap.add_argument('--every', type=int, default=4, help='line sampling interval [steps]')
    ap.add_argument('--stress-scale', dest='stress_scale', type=float, default=1.0,
                    help='multiply the prescribed Reynolds stresses (profile cases) by this factor')
    ap.add_argument('--uref', type=float, default=1.0, help='mean speed and eddy convection speed')
    ap.add_argument('--numerics', default='', help='extra &dopamine_numerics lines, separated by ;')
    ap.add_argument('--zperiodic', action='store_true', help='periodic z (no slip walls, no cut eddies)')
    ap.add_argument('--tag', default='', help='suffix for the run directories')
    ap.add_argument('--reuse', action='store_true', help='analyse existing runs instead of rerunning')
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    dx = a.lx / a.nx
    print(f'# box {a.lx}x{a.ly}x{a.lz}, {a.nx}x{a.ny}x{a.nz}, sigma = {a.sigma} ({a.sigma/dx:.1f} cells), '
          f'{a.neddies} eddies, scheme {a.scheme}, t in [{a.t_spin}, {a.t_end}]')
    res = {}
    for case in a.cases.split(','):
        for kernel in a.kernels.split(','):
            wd = run_case(a, kernel, case)
            res[(case, kernel)] = analyse(a, wd)
    for case in a.cases.split(','):
        print(f'\n== {case}: Cp = p_rms / u_rms^2 (u_rms^2 = <u\'u\'>/3) ==')
        print('  x/sigma   ' + '  '.join(f'{k:>22s}' for k in a.kernels.split(',')))
        print('            ' + '  '.join(f'{"p_rms":>8s} {"u_rms":>6s} {"Cp":>5s}' for k in a.kernels.split(',')))
        ref = {}
        xs = None
        for k in a.kernels.split(','):
            xs, pr, ur, um_, ns = res[(case, k)]
            ref[k] = (pr, ur)
        stations = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 26.5, 36.5, 50.5, 70.5, 100.5]
        for st in stations:
            i = int(round(st - 0.5)) if st < a.nx else a.nx - 1
            if i >= a.nx:
                continue
            row = f'  {xs[i]/a.sigma:7.3f}   '
            for k in a.kernels.split(','):
                pr, ur = ref[k]
                row += f'{pr[i]:8.4f} {ur[i]:6.4f} {pr[i]/max(ur[i]**2,1e-30):5.2f}   '
            print(row)
        for k in a.kernels.split(','):
            pr, ur = ref[k]
            cp = pr / np.maximum(ur ** 2, 1e-30)
            band = (xs > 2.0 * a.sigma) & (xs < 2.5)
            print(f'  {k:>7s}: Cp[first 3 cells] = {cp[:3].mean():.2f}, Cp[max, first 0.5 sigma] = '
                  f'{cp[xs < 0.5*a.sigma].max():.2f}, Cp[2 sigma..2.5] = {cp[band].mean():.2f}, '
                  f'p_rms ratio (first 3 cells / band) = {pr[:3].mean()/pr[band].mean():.2f}, '
                  f'u_rms first/band = {ur[:3].mean()/ur[band].mean():.2f}, samples = {res[(case,k)][4]}')


if __name__ == '__main__':
    main()
