#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
projection_model.py
===================
Minimal model of advance_rk3w on a uniform periodic grid, used to trace the O(dt h^2)
energy loss seen in the temporal study.  numpy only.

It keeps what matters: cell velocities and a separate face flux, central advection with the
stored flux, the RK3W weights, pre-correction with the wide (cell) gradient of p^{k-1/2},
the Rhie-Chow face flux, a compact-stencil Poisson solve (FFT), flux correction with the
compact gradient and cell-velocity correction with the wide gradient.  For the inviscid
Taylor-Green vortex it reproduces the solver's energy loss to 0.03%
(N=16, dt=0.01, t=0.2: 5.024e-6 here, 5.025e-6 in dopamine).

Variants (options of run()):
  rc=False                  drop the Rhie-Chow term from the face flux
  precorr=False             no pressure pre-correction (non-incremental)
  convflux="u"              advect with interp(u) instead of the stored flux
  cellgrad="adjoint"        cell velocity projected with the operator adjoint to the wide
                            gradient (an exact orthogonal projection of the cell field)

    python3 projection_model.py
"""
import math

import numpy as np

AL = (8 / 15, 5 / 12, 3 / 4)
BE = (0, -17 / 60, -5 / 12)
GA = (8 / 15, 2 / 15, 1 / 3)


def run(N, dt, tend, opts=None):
    """Inviscid Taylor-Green vortex; returns (kinetic energy at tend, max |div q|)."""
    o = dict(cellgrad="wide", rc=True, convflux="q", precorr=True)
    o.update(opts or {})
    L = 2 * np.pi
    h = L / N
    x = (np.arange(N) + .5) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = np.array([np.sin(X) * np.cos(Y) * np.cos(Z), -np.cos(X) * np.sin(Y) * np.cos(Z), 0 * X])

    def sh(a, d, s):                                    # a[i+s] along axis d
        return np.roll(a, -s, axis=d)

    def wide(p):                                        # cell gradient (p[i+1]-p[i-1])/2h
        return np.array([(sh(p, d, 1) - sh(p, d, -1)) / (2 * h) for d in range(3)])

    def comp(p):                                        # compact face gradient
        return np.array([(sh(p, d, 1) - p) / h for d in range(3)])

    def interp(v):                                      # cell -> face average
        return np.array([0.5 * (v[d] + sh(v[d], d, 1)) for d in range(3)])

    def div(q):
        return sum((q[d] - sh(q[d], d, -1)) / h for d in range(3))

    k = 2 * np.pi * np.fft.fftfreq(N, d=h)
    kk = np.meshgrid(k, k, k, indexing="ij")
    lam = sum((2 * np.cos(a * h) - 2) / h ** 2 for a in kk)
    lam[0, 0, 0] = 1.0
    lw = sum(-(np.sin(a * h) / h) ** 2 for a in kk)
    bad = np.abs(lw) < 1e-12
    lw[bad] = 1.0

    def conv(u, q):
        H = np.zeros_like(u)
        for c in range(3):
            for d in range(3):
                uf = 0.5 * (u[c] + sh(u[c], d, 1))
                H[c] -= (q[d] * uf - sh(q[d] * uf, d, -1)) / h
        return H

    p = np.zeros((N, N, N))
    q = interp(u)
    H0 = np.zeros_like(u)
    for _ in range(int(round(tend / dt))):
        for s in range(3):
            H = conv(u, q if o["convflux"] == "q" else interp(u))
            u = u + dt * (AL[s] * H + BE[s] * H0)
            gk = GA[s] * dt
            pold = p.copy() if o["precorr"] else 0 * p
            gp = wide(pold)
            us = u - gk * gp
            qs = interp(us)
            if o["rc"]:
                qs = qs + gk * (interp(gp) - comp(pold))
            r = np.fft.fftn(div(qs) / gk)
            r[0, 0, 0] = 0
            phi = np.real(np.fft.ifftn(r / lam))
            q = qs - gk * comp(phi)
            if o["cellgrad"] == "adjoint":
                r = np.fft.fftn(sum((sh(us[d], d, 1) - sh(us[d], d, -1)) / (2 * h) for d in range(3)) / gk) / lw
                r[bad] = 0
                u = us - gk * wide(np.real(np.fft.ifftn(r)))
            else:
                u = us - gk * wide(phi)
            p = pold + phi if o["precorr"] else phi
            H0 = H
    return 0.5 * np.mean(np.sum(u ** 2, axis=0)), np.abs(div(q)).max()


def order(a, b):
    return math.log2(abs(a / b))


def main():
    print("Inviscid Taylor-Green, t = 0.2: energy loss 1/8 - E (exact answer: 0)")
    cases = [("baseline (as in dopamine)", {}),
             ("no Rhie-Chow term", dict(rc=False)),
             ("no pre-correction", dict(precorr=False)),
             ("advect with interp(u)", dict(convflux="u")),
             ("cell projection made orthogonal", dict(cellgrad="adjoint"))]
    for name, opt in cases:
        print(f"\n{name}")
        for N in (8, 16, 32):
            a = 0.125 - run(N, 0.01, 0.2, opt)[0]
            b = 0.125 - run(N, 0.005, 0.2, opt)[0]
            print(f"  N={N:2d}  dt=0.01: {a:+.3e}   dt=0.005: {b:+.3e}   order in dt: {order(a, b):.2f}")


if __name__ == "__main__":
    main()
