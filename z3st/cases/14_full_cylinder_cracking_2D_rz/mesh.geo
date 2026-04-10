// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D Axisymmetric cylinder section (Rectangle)
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

R = 10.0e-3;      // Radius (m)
H = 10.0e-3;      // Height (m)

lc_outer  = 1.0e-5;   // Refined mesh near contact surface (10 um)
lc_center = 1.0e-3;   // Coarser near axis (1 mm)

Point(1) = {0, 0, 0, lc_center};
Point(2) = {R, 0, 0, lc_outer};
Point(3) = {R, H, 0, lc_outer};
Point(4) = {0, H, 0, lc_center};

Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Outer wall (Quench)
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Symmetry axis

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("uo2", 10) = {1};
Physical Curve("bottom", 30) = {1};
Physical Curve("contact_wall", 20) = {2};
Physical Curve("top", 40) = {3};
Physical Curve("axis", 50) = {4};

// Mesh settings
Mesh.Algorithm = 6;

// Refine near the outer boundary (where cracks start)
Field[1] = Distance;
Field[1].CurvesList = {2};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_outer;
Field[2].SizeMax = lc_center;
Field[2].DistMin = 0.0;
Field[2].DistMax = 4.0e-3; // Match 3D refinement distance (4mm)

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
