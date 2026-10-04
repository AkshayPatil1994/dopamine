#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Orthogonal uniform hex box (Gmsh v4 ASCII) with named faces x0, xL, y0, yL, z0, zL.

    mk_box.py NX NY NZ LX LY LZ OUT.msh

Used by the SEM pressure-fluctuation study; gmsh is not needed.
"""
import sys
NX,NY,NZ=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
LX,LY,LZ=float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]); out=sys.argv[7]
def nid(i,j,k): return 1+i+(NX+1)*(j+(NY+1)*k)
nodes=[(LX*i/NX,LY*j/NY,LZ*k/NZ) for k in range(NZ+1) for j in range(NY+1) for i in range(NX+1)]
hexes=[[nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)]
       for k in range(NZ) for j in range(NY) for i in range(NX)]
quads={s:[] for s in range(1,7)}
for a in range(NY):
  for b in range(NZ):
    quads[1].append([nid(0,a,b),nid(0,a+1,b),nid(0,a+1,b+1),nid(0,a,b+1)])
    quads[2].append([nid(NX,a,b),nid(NX,a+1,b),nid(NX,a+1,b+1),nid(NX,a,b+1)])
for a in range(NX):
  for b in range(NZ):
    quads[3].append([nid(a,0,b),nid(a+1,0,b),nid(a+1,0,b+1),nid(a,0,b+1)])
    quads[4].append([nid(a,NY,b),nid(a+1,NY,b),nid(a+1,NY,b+1),nid(a,NY,b+1)])
for a in range(NX):
  for b in range(NY):
    quads[5].append([nid(a,b,0),nid(a+1,b,0),nid(a+1,b+1,0),nid(a,b+1,0)])
    quads[6].append([nid(a,b,NZ),nid(a+1,b,NZ),nid(a+1,b+1,NZ),nid(a,b+1,NZ)])
names=['x0','xL','y0','yL','z0','zL']
f=open(out,'w')
f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
f.write("$PhysicalNames\n6\n")
for s in range(6): f.write(f'2 {s+1} "{names[s]}"\n')
f.write("$EndPhysicalNames\n")
f.write("$Entities\n0 0 6 1\n")
for s in range(1,7): f.write(f"{s} 0 0 0 1 1 1 1 {s} 0\n")
f.write("1 0 0 0 1 1 1 0 0\n$EndEntities\n")
nn=len(nodes)
f.write(f"$Nodes\n1 {nn} 1 {nn}\n3 1 0 {nn}\n")
for n in range(1,nn+1): f.write(f"{n}\n")
for (x,y,z) in nodes: f.write(f"{x:.16g} {y:.16g} {z:.16g}\n")
f.write("$EndNodes\n")
ne=len(hexes)+sum(len(v) for v in quads.values())
f.write(f"$Elements\n7 {ne} 1 {ne}\n")
t=1
for s in range(1,7):
  f.write(f"2 {s} 3 {len(quads[s])}\n")
  for q in quads[s]:
    f.write(f"{t} {q[0]} {q[1]} {q[2]} {q[3]}\n"); t+=1
f.write(f"3 1 5 {len(hexes)}\n")
for h in hexes:
  f.write(f"{t} "+" ".join(map(str,h))+"\n"); t+=1
f.write("$EndElements\n")
