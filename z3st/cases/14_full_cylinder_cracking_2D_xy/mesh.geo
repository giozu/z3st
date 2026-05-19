// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D circular pellet cross-section (transverse cut),
//  with a cold-contact arc on the right side.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

// Geometry
R = 10.0e-3;                          // radius (m) = 10 mm
half_contact_deg = 30.0;
half_contact = half_contact_deg * Pi / 180.0;

// Pre-crack: radial slit at theta = 15 deg (middle of upper-half contact arc)
crack_theta_deg  = 15.0;
crack_theta      = crack_theta_deg * Pi / 180.0;
crack_length     = 2.5e-4;            // 250 um = 5 * lc;
crack_R_inner    = R - crack_length;  // r-coordinate of the crack tip

// Mesh sizing
lc_outer  = 1.25e-5;                  // 12.5 um (lc/4 with phase-field lc = 50 um)
lc_center = 2.0e-4;                   // 200 um (4 lc; coarsens away from rim)
pin_size  = 5.0e-5;                   // 50 um pin segment at (-R, 0)

// Points
Point(1) = {0, 0, 0, lc_center};
Point(2) = {R, 0, 0, lc_outer};
Point(3) = {R*Cos(crack_theta), R*Sin(crack_theta), 0, lc_outer};
Point(4) = {R*Cos(half_contact), R*Sin(half_contact), 0, lc_outer};
Point(5) = {-R, 0, 0, lc_outer};
Point(6) = {-R + pin_size, 0, 0, lc_outer}; 
Point(7) = {crack_R_inner*Cos(crack_theta), crack_R_inner*Sin(crack_theta), 0, lc_outer};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Line(4)   = {5, 6}; 
Line(5)   = {6, 1};
Line(6)   = {1, 2};
Line(7)   = {3, 7};

// Curve loop
Curve Loop(1) = {1, 7, -7, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Physical groups (tags must match those in geometry.yaml)
Physical Surface("uo2", 10) = {1};
Physical Curve("contact_wall", 20)   = {1, 2};     // both halves of the cold contact arc
Physical Curve("insulated_wall", 21) = {3};        // remaining 5/6 of the perimeter
Physical Curve("pin", 30)            = {4};        // 50-um Clamp_x segment
Physical Curve("symmetry", 40)       = {5, 6};     // y=0 mirror symmetry plane (Clamp_y)
Physical Curve("crack_seed", 50)     = {7};        // pre-crack slit; Dirichlet D = 1 applied here

// Refinement field
Field[1] = Distance;
Field[1].CurvesList = {1, 2, 7};

Field[2] = Threshold;
Field[2].InField   = 1;
Field[2].SizeMin   = lc_outer;
Field[2].SizeMax   = lc_center;
Field[2].DistMin   = 0.0;
Field[2].DistMax   = 5.0e-3;

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm                  = 6;

// Usage:
//   gmsh -2 mesh.geo -format msh2
