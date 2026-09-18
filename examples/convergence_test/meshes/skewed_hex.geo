// Skewed hexahedral mesh for Taylor-Green convergence testing
// Sheared parallelepiped: the z-edges are tilted in x by `shear = 0.3*L`
// producing non-orthogonal hex cells while remaining fully periodic.
//
// Domain:   z=0 face is the standard square [0,L]x[0,L]
//           z=L face is shifted +shear in x -> [shear, L+shear] x [0,L]
// All six boundary patches remain planar and form valid periodic pairs:
//   x0 <-> xL  (shift = (L, 0, 0))
//   y0 <-> yL  (shift = (0, L, 0))
//   z0 <-> zL  (shift = (shear, 0, L))
// DOPAMINE's geometric periodic pairing resolves all three automatically.
//
// Usage:
//   gmsh -3 -format msh4 -setnumber N 16 -o skewed_hex_N16.msh skewed_hex.geo
//   gmsh -3 -format msh4 -setnumber N 32 -o skewed_hex_N32.msh skewed_hex.geo

L     = 6.283185307179586;
shear = 0.3 * L;   // ~16.7 degree skew angle; adjust to taste

// Resolution parameter (override from command line with -setnumber N <value>)
N = DefineNumber[ 16, Name "N" ];

// Bottom face (z=0): standard square
Point(1) = {0,         0, 0, 1.0};
Point(2) = {L,         0, 0, 1.0};
Point(3) = {L,         L, 0, 1.0};
Point(4) = {0,         L, 0, 1.0};

// Top face (z=L): shifted by +shear in x
Point(5) = {shear,     0, L, 1.0};
Point(6) = {shear + L, 0, L, 1.0};
Point(7) = {shear + L, L, L, 1.0};
Point(8) = {shear,     L, L, 1.0};

// Bottom and top face edges
Line(1)  = {1, 2};  // z=0, y=0
Line(2)  = {2, 3};  // z=0, x=L
Line(3)  = {3, 4};  // z=0, y=L
Line(4)  = {4, 1};  // z=0, x=0
Line(5)  = {5, 6};  // z=L, y=0
Line(6)  = {6, 7};  // z=L, x=L+shear
Line(7)  = {7, 8};  // z=L, y=L
Line(8)  = {8, 5};  // z=L, x=shear

// Skewed vertical (z) edges
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Six planar surfaces (all four vertices of each surface are coplanar)
Line Loop(1) = {1, 2, 3, 4};      Plane Surface(1) = {1};  // z=0
Line Loop(2) = {5, 6, 7, 8};      Plane Surface(2) = {2};  // z=L
Line Loop(3) = {-1, 9, 5, -10};   Plane Surface(3) = {3};  // y=0
Line Loop(4) = {-3, 11, 7, -12};  Plane Surface(4) = {4};  // y=L
Line Loop(5) = {-4, 12, 8, -9};   Plane Surface(5) = {5};  // x=0 (tilted)
Line Loop(6) = {-2, 10, 6, -11};  Plane Surface(6) = {6};  // x=L (tilted)

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line{1, 3, 5, 7} = N + 1;
Transfinite Line{2, 4, 6, 8} = N + 1;
Transfinite Line{9, 10, 11, 12} = N + 1;
Transfinite Surface{1, 2, 3, 4, 5, 6} Alternated;
Transfinite Volume{1};
Recombine Surface{1, 2, 3, 4, 5, 6};
Recombine Volume{1};

Physical Surface("periodic_z0", 1) = {1};
Physical Surface("periodic_zL", 2) = {2};
Physical Surface("periodic_y0", 3) = {3};
Physical Surface("periodic_yL", 4) = {4};
Physical Surface("periodic_x0", 5) = {5};
Physical Surface("periodic_xL", 6) = {6};
Physical Volume("fluid", 7) = {1};

Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
