// Structured box mesh for Rayleigh-Benard convection.
// Coordinates: x,y=horizontal (periodic), z=vertical/wall-normal (no-slip, heated/cooled).
// Aspect ratio Lx=Ly=2*Lz approx. matches Kunnen et al. (2009) / Lavezzo et al. (2011) Table 1.

SetFactory("OpenCASCADE");

// Domain size
Lx = 2.0;
Ly = 2.0;
Lz = 1.0;
// Number of grid points
Nx = 128;
Ny = 128;
Nz = 64;
// Stretching ratio (< 0.0 for clustering near walls, > 0.0 for clustering near center)
BumpZ = 0.25;

Point(1) = {0,  0,  0, 1.0};
Point(2) = {Lx, 0,  0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0,  Ly, 0, 1.0};
Point(5) = {0,  0,  Lz, 1.0};
Point(6) = {Lx, 0,  Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0,  Ly, Lz, 1.0};

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

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {4, 9, -8, -12};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Curve Loop(5) = {1, 10, -5, -9};
Plane Surface(5) = {5};
Curve Loop(6) = {3, 12, -7, -11};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1, 3, 5, 7} = Nx + 1;
Transfinite Line {2, 4, 6, 8} = Ny + 1;
Transfinite Line {9, 10, 11, 12} = Nz + 1 Using Bump BumpZ;

Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};
Recombine Volume {1};

Periodic Surface {4} = {3} Translate {Lx, 0, 0};
Periodic Surface {6} = {5} Translate {0, Ly, 0};

Physical Surface("bottom_wall") = {1};
Physical Surface("top_wall") = {2};
Physical Surface("periodic_x0") = {3};
Physical Surface("periodic_xL") = {4};
Physical Surface("periodic_y0") = {5};
Physical Surface("periodic_yL") = {6};
Physical Volume("fluid") = {1};
