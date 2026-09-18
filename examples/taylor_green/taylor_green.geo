// Taylor-Green vortex periodic box mesh
// Generates a structured hex mesh over [0, 2*pi]^3
// Usage:
//   gmsh -3 -format msh4 -o taylor_green.msh taylor_green.geo

L = 6.283185307179586;

// Resolution (increase for convergence studies)
n_points = 64;
Nx = n_points;
Ny = n_points;
Nz = n_points;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};
Point(5) = {0, 0, L, 1.0};
Point(6) = {L, 0, L, 1.0};
Point(7) = {L, L, L, 1.0};
Point(8) = {0, L, L, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // z=0
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2}; // z=L
Line Loop(3) = {-1, 9, 5, -10};
Plane Surface(3) = {3}; // y=0
Line Loop(4) = {-3, 11, 7, -12};
Plane Surface(4) = {4}; // y=L
Line Loop(5) = {-4, 12, 8, -9};
Plane Surface(5) = {5}; // x=0
Line Loop(6) = {-2, 10, 6, -11};
Plane Surface(6) = {6}; // x=L

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line{1, 3, 5, 7} = Nx + 1;
Transfinite Line{2, 4, 6, 8} = Ny + 1;
Transfinite Line{9, 10, 11, 12} = Nz + 1;
Transfinite Surface{1, 2, 3, 4, 5, 6} Alternated;
Transfinite Volume{1};
Recombine Surface{1, 2, 3, 4, 5, 6};
Recombine Volume{1};

Physical Surface("periodic_x0", 1) = {5};
Physical Surface("periodic_xL", 2) = {6};
Physical Surface("periodic_y0", 3) = {3};
Physical Surface("periodic_yL", 4) = {4};
Physical Surface("periodic_z0", 5) = {1};
Physical Surface("periodic_zL", 6) = {2};
Physical Volume("fluid", 7) = {1};

Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
