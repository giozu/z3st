// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2D elliptical cavity in a rectangular plate (cross-section)
//  Coherent with a given semi-dihedral angle
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// --- Nouveaux paramètres pour le cas Jiang ---
scale = 1e-6; // (meters)
Lx = 60 * scale;
Ly = 60 * scale;
ax = 7.9 * scale; // Semi-axis X









ay = 7.9 * scale; // Semi-axis Y











h_plate = 4.0 * scale;    // coarse mesh size
h_cavity = 0.25 * scale;  // fine mesh size (article de Jiang: h = 0.25 µm)

// --- Anciens paramètres (mis en commentaire) ---
// scale = 1; // (micrometers)
// Lx = 2.0 * scale;
// Ly = 2.0 * scale;
// ax = 0.06 * scale; // Semi-axis X
// theta = 60 * Pi / 180; // Semi-dihedral angle in radiants 
// ay = ax * (1 - Cos(theta)) / Sin(theta); // ax * tan(theta/2)
// rho = ay*ay/ax;
// h_plate = scale * 0.0400; // coarse mesh size
// h_cavity = rho / 3.0;     // fine mesh size, adaptive to curvature

// Rectangular plate (Quart de symétrie : coin inférieur gauche en 0,0)
Rectangle(1) = {0, 0, 0, Lx, Ly};
// // Ancien rectangle centré (mis en commentaire)
// Rectangle(1) = {-Lx/2, -Ly/2, 0, Lx, Ly};

// Elliptical hole
Ellipse(10) = {0, 0, 0, ax, ay};

Curve Loop(10) = {10};
Plane Surface(10) = {10};

// Subtract the ellipse from the rectangle
BooleanDifference{ Surface{1}; Delete; }{ Surface{10}; Delete; }

// --- Groupes Physiques ---
Physical Surface("uo2") = {1};
// // Anciens groupes physiques statiques (mis en commentaire)
// Physical Curve("cavity") = {10};  // Internal elliptical boundary
// Physical Curve("ymin") = {11};    // Bottom
// Physical Curve("xmin") = {12};    // Left
// Physical Curve("xmax") = {13};    // Right
// Physical Curve("ymax") = {14};    // Top

// Sélection dynamique des frontières par boîte englobante pour éviter les erreurs d'IDs post-opération booléenne
eps = 2e-6;
c_ymin() = Curve In BoundingBox{ax-eps, -eps, -eps, Lx+eps, eps, eps};
c_xmin() = Curve In BoundingBox{-eps, ay-eps, -eps, eps, Ly+eps, eps};
c_xmax() = Curve In BoundingBox{Lx-eps, -eps, -eps, Lx+eps, Ly+eps, eps};
c_ymax() = Curve In BoundingBox{-eps, Ly-eps, -eps, Lx+eps, Ly+eps, eps};
c_cavity() = Curve In BoundingBox{-eps, -eps, -eps, ax+eps, ay+eps, eps};

Physical Curve("ymin") = {c_ymin()};
Physical Curve("xmin") = {c_xmin()};
Physical Curve("xmax") = {c_xmax()};
Physical Curve("ymax") = {c_ymax()};
Physical Curve("cavity") = {c_cavity()};

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {c_cavity(), c_ymin()};
// Field[1].CurvesList = {10}; // Ancienne liste de courbes statique
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 0.5 * scale;   // Taille minimale jusqu'à 0.5 µm des courbes
Field[2].DistMax = 10.0 * scale;  // Transition progressive vers la taille grossière sur 10 µm

Background Field = 2;

// Mesh generation options
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Algorithm = 6; // Frontal-Delaunay for better quality in 2D
// Mesh 2;
// Save "mesh.msh";
