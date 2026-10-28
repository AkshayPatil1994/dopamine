// Half channel (wall at z=0, free-slip centreline at z=Lz=delta) for the AMD validation at Re_tau = 2000.
// Structured channel mesh for LES wall-model validation.
// Coordinates: x=streamwise, y=spanwise, z=wall-normal (vertical).

SetFactory("OpenCASCADE");

Lx = 2*Pi;  
Ly = Pi;  
Lz = 1.0;                // half-channel height: wall at z=0, free-slip centreline at z=Lz

Nx = 64;   // streamwise cells: dx+ = 196
Ny = 64;   // spanwise cells: dy+ = 98
Nz = 64;   // wall-normal cells: first cell centre at y+ = 16
xPeriodic = 1;  // 1 - periodic in streamwise direction, 0 for inflow/outflow

// BumpZ < 0 - refines at one wall | -1.05 is 5% stretching at bottom wall, -0.95 is 5% stretching at top wall
// BumpZ > 0 - refines at both walls | 2.0 is 5% stretching in the center, 0.5 is 5% stretching at the walls
BumpZ = -1.0;   // uniform wall-normal spacing (Progression 1.0)



// Define all points
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
Plane Surface(1) = {1}; // z=0 bottom wall

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2}; // z=Lz top wall

Curve Loop(3) = {4, 9, -8, -12};
Plane Surface(3) = {3}; // x=0 periodic

Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4}; // x=Lx periodic

Curve Loop(5) = {1, 10, -5, -9};
Plane Surface(5) = {5}; // y=0 periodic

Curve Loop(6) = {3, 12, -7, -11};
Plane Surface(6) = {6}; // y=Ly periodic

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1, 3, 5, 7} = Nx + 1;
Transfinite Line {2, 4, 6, 8} = Ny + 1;

If (BumpZ > 0)
   Transfinite Line {9, 10, 11, 12} = Nz + 1 Using Bump BumpZ;
Else
   temp_BumpZ = Fabs(BumpZ);
   Transfinite Line {9, 10, 11, 12} = Nz + 1 Using Progression temp_BumpZ;
EndIf

Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};
Recombine Volume {1};

// Periodic boundaries along x and y directions
If (xPeriodic == 1)
   Periodic Surface {4} = {3} Translate {Lx, 0, 0};   
Else
   print("X is not Periodic")
EndIf

Periodic Surface {6} = {5} Translate {0, Ly, 0};

Physical Surface("bottom_wall") = {1};
Physical Surface("top_slip") = {2};
Physical Surface("periodic_x0") = {3};
Physical Surface("periodic_xL") = {4};
Physical Surface("periodic_y0") = {5};
Physical Surface("periodic_yL") = {6};
Physical Volume("fluid") = {1};