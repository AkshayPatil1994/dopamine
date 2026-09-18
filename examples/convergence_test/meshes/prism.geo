// Prismatic mesh for Taylor-Green convergence testing
// 2D unstructured triangles extruded uniformly in z -> triangular prisms
// Periodic cube [0, 2*pi]^3
//
// Usage:
//   gmsh -3 -format msh4 -setnumber N 16 -o prism_N16.msh prism.geo
//   gmsh -3 -format msh4 -setnumber N 32 -o prism_N32.msh prism.geo
//
// Physical surface order from Extrude{}{Surface{1}; Layers{N};}:
//   out[0] = Volume,  out[1] = top face (zL),
//   out[2] = y=0 wall (from Line 1),  out[3] = x=L wall (from Line 2),
//   out[4] = y=L wall (from Line 3),  out[5] = x=0 wall (from Line 4)

L = 6.283185307179586;

// Resolution parameter (override from command line with -setnumber N <value>)
N = DefineNumber[ 16, Name "N" ];

// Bottom face (z=0)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};

Line(1) = {1, 2};  // y=0 bottom edge
Line(2) = {2, 3};  // x=L bottom edge
Line(3) = {3, 4};  // y=L bottom edge
Line(4) = {4, 1};  // x=0 bottom edge

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Transfinite edges control face-to-face node count for proper periodic pairing
Transfinite Line{1, 2, 3, 4} = N + 1;

// Extrude in z-direction with N uniform layers; no Recombine -> prism cells
out[] = Extrude{0, 0, L} {
  Surface{1};
  Layers{N};
};

Physical Surface("periodic_z0", 1) = {1};        // bottom face (original)
Physical Surface("periodic_zL", 2) = {out[1]};   // top face (extruded replica)
Physical Surface("periodic_y0", 3) = {out[2]};   // y=0 wall  (from Line 1)
Physical Surface("periodic_xL", 4) = {out[3]};   // x=L wall  (from Line 2)
Physical Surface("periodic_yL", 5) = {out[4]};   // y=L wall  (from Line 3)
Physical Surface("periodic_x0", 6) = {out[5]};   // x=0 wall  (from Line 4)
Physical Volume("fluid", 7) = {out[0]};

Mesh.ElementOrder = 1;
Mesh.Algorithm   = 5;   // 2D: Delaunay (triangles)
Mesh.Optimize = 1;
