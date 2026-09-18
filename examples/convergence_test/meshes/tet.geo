// Unstructured tetrahedral mesh for Taylor-Green convergence testing
// Periodic cube [0, 2*pi]^3 with tetrahedral cells
// Target cell spacing h ~= L/N (Delaunay refinement)
//
// Usage:
//   gmsh -3 -format msh4 -setnumber N 16 -o tet_N16.msh tet.geo
//   gmsh -3 -format msh4 -setnumber N 32 -o tet_N32.msh tet.geo

L = 6.283185307179586;

// Resolution parameter (override from command line with -setnumber N <value>)
N = DefineNumber[ 16, Name "N" ];

Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};
Point(5) = {0, 0, L, 1.0};
Point(6) = {L, 0, L, 1.0};
Point(7) = {L, L, L, 1.0};
Point(8) = {0, L, L, 1.0};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 5};
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};      Plane Surface(1) = {1}; // z=0
Line Loop(2) = {5, 6, 7, 8};      Plane Surface(2) = {2}; // z=L
Line Loop(3) = {-1, 9, 5, -10};   Plane Surface(3) = {3}; // y=0
Line Loop(4) = {-3, 11, 7, -12};  Plane Surface(4) = {4}; // y=L
Line Loop(5) = {-4, 12, 8, -9};   Plane Surface(5) = {5}; // x=0
Line Loop(6) = {-2, 10, 6, -11};  Plane Surface(6) = {6}; // x=L

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Control edge density (no Transfinite Surface/Volume - leaves meshing free)
Transfinite Line{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = N + 1;

// No Recombine: pure tetrahedral decomposition
Physical Surface("periodic_z0", 1) = {1};
Physical Surface("periodic_zL", 2) = {2};
Physical Surface("periodic_y0", 3) = {3};
Physical Surface("periodic_yL", 4) = {4};
Physical Surface("periodic_x0", 5) = {5};
Physical Surface("periodic_xL", 6) = {6};
Physical Volume("fluid", 7) = {1};

Mesh.ElementOrder = 1;
Mesh.Algorithm   = 5;   // 2D: Delaunay
Mesh.Algorithm3D = 4;   // 3D: Frontal-Delaunay
Mesh.Optimize = 1;
