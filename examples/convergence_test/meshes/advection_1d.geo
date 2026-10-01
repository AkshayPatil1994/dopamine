// Thin periodic bar for 1D advection tests: [0, 1] x [0, 4/N] x [0, 4/N]
// with N x 4 x 4 cubic hex cells (cell size h = 1/N). run_advection.py writes N
// into a copy of this file (gmsh -setnumber does not override DefineNumber).

Lx = 1.0;

// Resolution parameter (override from command line with -setnumber N <value>)
N = DefineNumber[ 16, Name "N" ];

L = Lx;
Nx = N;
Ny = 4;
Nz = 4;
Ly = 4.0 / N;
Lz = 4.0 / N;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};
Point(5) = {0, 0, Lz, 1.0};
Point(6) = {Lx, 0, Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0, Ly, Lz, 1.0};

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

Transfinite Line{1, 3, 5, 7} = Nx + 1;
Transfinite Line{2, 4, 6, 8} = Ny + 1;
Transfinite Line{9, 10, 11, 12} = Nz + 1;
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
