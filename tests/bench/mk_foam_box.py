#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Uniform hex box as a binary OpenFOAM polyMesh, for single-GPU benchmarks (no gmsh needed).

    mk_foam_box.py NX NY NZ OUTDIR [LX LY LZ]

Writes OUTDIR/{points,faces,owner,neighbour,boundary} (binary, label=32, scalar=64, compact faces).
Patches 1..6 are x0, xL, y0, yL, z0, zL (type 'patch'); pair them via &dopamine_bc patch_type.
"""
import os
import sys

import numpy as np

HEADER = """FoamFile
{{
    version     2.0;
    format      binary;
    arch        "LSB;label=32;scalar=64";
    class       {cls};
    location    "constant/polyMesh";
    object      {obj};
}}

"""


def write_list(path, cls, obj, arrays):
    with open(path, "wb") as f:
        f.write(HEADER.format(cls=cls, obj=obj).encode())
        for a in arrays:
            f.write(f"{a.shape[0]}\n(".encode())
            f.write(np.ascontiguousarray(a).tobytes())
            f.write(b")\n")


def main():
    nx, ny, nz = (int(v) for v in sys.argv[1:4])
    out = sys.argv[4]
    lx, ly, lz = (float(v) for v in sys.argv[5:8]) if len(sys.argv) >= 8 else (1.0, 1.0, 1.0)
    os.makedirs(out, exist_ok=True)

    def pid(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    # cell index c = i + nx*(j + ny*k): flatten with k slowest, i fastest
    i, j, k = (a.transpose(2, 1, 0).ravel() for a in (i, j, k))
    cell = i + nx * (j + ny * k)

    px = (np.arange(nx + 1) * lx / nx)
    py = (np.arange(ny + 1) * ly / ny)
    pz = (np.arange(nz + 1) * lz / nz)
    pts = np.empty(((nx + 1) * (ny + 1) * (nz + 1), 3))
    pts[:, 0] = np.tile(px, (ny + 1) * (nz + 1))
    pts[:, 1] = np.tile(np.repeat(py, nx + 1), nz + 1)
    pts[:, 2] = np.repeat(pz, (nx + 1) * (ny + 1))

    def quad(a, b, c, d):
        return np.stack([a, b, c, d], axis=1).astype(np.int32)

    # +x, +y, +z faces of every cell, outward normals (right-hand rule)
    fx = quad(pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i + 1, j + 1, k + 1), pid(i + 1, j, k + 1))
    fy = quad(pid(i, j + 1, k), pid(i, j + 1, k + 1), pid(i + 1, j + 1, k + 1), pid(i + 1, j + 1, k))
    fz = quad(pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1))
    # -x, -y, -z faces: reversed point order
    gx = quad(pid(i, j, k), pid(i, j, k + 1), pid(i, j + 1, k + 1), pid(i, j + 1, k))
    gy = quad(pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j, k + 1), pid(i, j, k + 1))
    gz = quad(pid(i, j, k), pid(i, j + 1, k), pid(i + 1, j + 1, k), pid(i + 1, j, k))

    # internal faces: per owner cell in order +x, +y, +z (neighbour index increasing)
    nc = nx * ny * nz
    stack = np.stack([fx, fy, fz], axis=1)  # (nc, 3, 4)
    nbr = np.stack([cell + 1, cell + nx, cell + nx * ny], axis=1)
    mask = np.stack([i < nx - 1, j < ny - 1, k < nz - 1], axis=1)
    int_faces = stack[mask]
    int_owner = np.repeat(cell, 3).reshape(nc, 3)[mask]
    int_nbr = nbr[mask]

    # boundary patches, each ordered by owner cell
    patches = [
        ("x0", i == 0, gx), ("xL", i == nx - 1, fx),
        ("y0", j == 0, gy), ("yL", j == ny - 1, fy),
        ("z0", k == 0, gz), ("zL", k == nz - 1, fz),
    ]
    bfaces, bowner, info = [], [], []
    start = len(int_faces)
    for name, m, f in patches:
        bfaces.append(f[m])
        bowner.append(cell[m])
        info.append((name, int(m.sum()), start))
        start += int(m.sum())

    faces = np.concatenate([int_faces] + bfaces).astype(np.int32)
    owner = np.concatenate([int_owner] + bowner).astype(np.int32)
    nf = faces.shape[0]
    offsets = (np.arange(nf + 1, dtype=np.int32) * 4)

    write_list(f"{out}/points", "vectorField", "points", [pts])
    write_list(f"{out}/faces", "faceCompactList", "faces", [offsets, faces.ravel()])
    write_list(f"{out}/owner", "labelList", "owner", [owner])
    write_list(f"{out}/neighbour", "labelList", "neighbour", [int_nbr.astype(np.int32)])

    with open(f"{out}/boundary", "w") as f:
        f.write(HEADER.format(cls="polyBoundaryMesh", obj="boundary").replace("binary", "ascii"))
        f.write(f"{len(info)}\n(\n")
        for name, n, s in info:
            f.write(f"    {name}\n    {{\n        type            patch;\n        nFaces          {n};\n"
                    f"        startFace       {s};\n    }}\n")
        f.write(")\n")
    print(f"{nc} cells, {nf} faces ({len(int_faces)} internal) -> {out}")


if __name__ == "__main__":
    main()
