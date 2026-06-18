// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a rectangular plate
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

L   = 0.025;        // plate length x (m) = 25 mm
H   = 0.0098;       // plate height y (m) = 9.8 mm
pin = 5.0e-5;       // 50-um Clamp_y segment at the bottom of the right edge

// Phase-field internal length lc = 9.2e-5 m.
lc_fine      = 3.0e-5;   // ~ lc/3 along the quenched edges (refine to ~2e-5 for converged crack counts)
lc_coarse    = 5.0e-4;   // interior
keep_fine    = 1.0e-3;   // hold lc_fine through the crack-penetration depth (~0.93 mm at 10 ms)
grade_to     = 2.5e-3;   // then coarsen out to lc_coarse

Point(1) = {0, 0,   0, lc_coarse};   // bottom-left
Point(2) = {L, 0,   0, lc_coarse};   // bottom-right
Point(3) = {L, pin, 0, lc_coarse};   // top of pin segment
Point(4) = {L, H,   0, lc_coarse};   // top-right
Point(5) = {0, H,   0, lc_coarse};   // top-left

Line(1) = {1, 2};   // bottom   (quench)
Line(2) = {2, 3};   // pin      (Clamp_y)
Line(3) = {3, 4};   // symmetry (Clamp_x)
Line(4) = {4, 5};   // top      (quench)
Line(5) = {5, 1};   // left     (quench)

Curve Loop(1)    = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// Physical groups (tags must match geometry.yaml)
Physical Surface("plate", 10)  = {1};
Physical Curve("quench", 20)   = {1, 4, 5};   // bottom + top + left
Physical Curve("symmetry", 30) = {3};
Physical Curve("pin", 40)      = {2};

// Refinement: fine within refine_depth of the three quenched edges
Field[1] = Distance;
Field[1].CurvesList = {1, 4, 5};
Field[1].Sampling = 2000;   // straight Lines sample only endpoints by default; force dense sampling

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_fine;
Field[2].SizeMax = lc_coarse;
Field[2].DistMin = keep_fine;   // fine out to the crack depth, then grade
Field[2].DistMax = grade_to;

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm                  = 6;

// Usage:
//   gmsh -2 mesh.geo -format msh2
