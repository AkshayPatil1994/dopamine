#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Curvilinear unit-box hex mesh (Gmsh v4 ASCII) for the skewed-mesh regression cases.

    mk_box_mesh.py N AMP OUT.msh

Nodes of an N^3 grid are displaced by AMP*sin(pi X)sin(2 pi Y)sin(pi Z) (cyclically in x, y, z), so the six
faces stay planar and the unit-cube domain is unchanged while interior faces are skewed
(AMP = 0.04 at N = 12 gives about 30 degrees of non-orthogonality). Physical surfaces 1..6 are
x0, xL, y0, yL, z0, zL; gmsh is not needed.
"""
import sys, math
N=int(sys.argv[1]); amp=float(sys.argv[2]); out=sys.argv[3]
def nid(i,j,k): return 1+i+(N+1)*(j+(N+1)*k)
nodes=[]
for k in range(N+1):
  for j in range(N+1):
    for i in range(N+1):
      X,Y,Z=i/N,j/N,k/N
      x=X+amp*math.sin(math.pi*X)*math.sin(2*math.pi*Y)*math.sin(math.pi*Z)
      y=Y+amp*math.sin(math.pi*Y)*math.sin(2*math.pi*Z)*math.sin(math.pi*X)
      z=Z+amp*math.sin(math.pi*Z)*math.sin(2*math.pi*X)*math.sin(math.pi*Y)
      nodes.append((x,y,z))
hexes=[]
for k in range(N):
  for j in range(N):
    for i in range(N):
      hexes.append([nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)])
# boundary quads by surface: 1 x0,2 xL,3 y0,4 yL,5 z0,6 zL
quads={s:[] for s in range(1,7)}
for a in range(N):
  for b in range(N):
    quads[1].append([nid(0,a,b),nid(0,a+1,b),nid(0,a+1,b+1),nid(0,a,b+1)])
    quads[2].append([nid(N,a,b),nid(N,a+1,b),nid(N,a+1,b+1),nid(N,a,b+1)])
    quads[3].append([nid(a,0,b),nid(a+1,0,b),nid(a+1,0,b+1),nid(a,0,b+1)])
    quads[4].append([nid(a,N,b),nid(a+1,N,b),nid(a+1,N,b+1),nid(a,N,b+1)])
    quads[5].append([nid(a,b,0),nid(a+1,b,0),nid(a+1,b+1,0),nid(a,b+1,0)])
    quads[6].append([nid(a,b,N),nid(a+1,b,N),nid(a+1,b+1,N),nid(a,b+1,N)])
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
