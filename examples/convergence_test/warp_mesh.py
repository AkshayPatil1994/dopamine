#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""
warp_mesh.py
============
Turn a structured periodic hex mesh (gmsh MSH 4.1, ASCII) into a smoothly warped,
non-orthogonal one while keeping it periodic.

Every node is displaced by a smooth field that is periodic with period L in each
direction:

    x' = x + A sin(k x) sin(k (y + z))
    y' = y + A sin(k y) sin(k (z + x))
    z' = z + A sin(k z) sin(k (x + y)),        k = 2 pi / L,  A = amp * L

The factor sin(k x) in the x-displacement (and likewise for y and z) makes the box faces
stay flat: nodes never leave the planes x = 0, L etc., so the bounding box is exactly
[0, L]^3 (the solver derives the Taylor-Green wavelength from it) and the opposite faces
move identically, so the periodic face pairs still match.
The exact solution does not depend on the mesh, only the cell shapes do: the warped
family exercises the non-orthogonal correction and the least-squares gradients while
the expected order stays second.

Usage:  warp_mesh.py in.msh out.msh [--amp 0.01] [--length 6.283185307179586]
"""
import argparse
import math
import sys


def warp(src, dst, amp, length):
    k = 2.0 * math.pi / length
    a = amp * length
    with open(src) as f:
        lines = f.read().split("\n")

    out = []
    i = 0
    n_nodes = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        i += 1
        if line.strip() != "$Nodes":
            continue
        n_blocks, n_nodes, _, _ = (int(t) for t in lines[i].split())
        out.append(lines[i]); i += 1
        for _ in range(n_blocks):
            hdr = lines[i].split()
            out.append(lines[i]); i += 1
            if int(hdr[2]) != 0:
                sys.exit("parametric nodes are not supported")
            nb = int(hdr[3])
            out.extend(lines[i:i + nb])        # node tags
            i += nb
            for _ in range(nb):                # coordinates
                x, y, z = (float(t) for t in lines[i].split()[:3])
                i += 1
                xw = x + a * math.sin(k * x) * math.sin(k * (y + z))
                yw = y + a * math.sin(k * y) * math.sin(k * (z + x))
                zw = z + a * math.sin(k * z) * math.sin(k * (x + y))
                out.append(f"{xw:.16e} {yw:.16e} {zw:.16e}")
    with open(dst, "w") as f:
        f.write("\n".join(out))
    return n_nodes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--amp", type=float, default=0.01,
                    help="displacement amplitude as a fraction of L (default 0.01, ~10 deg max non-orthogonality)")
    ap.add_argument("--length", type=float, default=2.0 * math.pi)
    args = ap.parse_args()
    n = warp(args.src, args.dst, args.amp, args.length)
    print(f"warped {n} nodes: {args.src} -> {args.dst}")


if __name__ == "__main__":
    main()
