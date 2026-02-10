// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2D elliptical cavity in a rectangular plate (cross-section)
//  Coherent with a given semi-dihedral angle
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

scale = 1; // (micrometers)

// Parameters
Lx = 2.0 * scale;
Ly = 2.0 * scale;
ax = 0.06 * scale; // Semi-axis X

theta = 60 * Pi / 180; // Semi-dihedral angle in radiants 
ay = ax * (1 - Cos(theta)) / Sin(theta); // ax * tan(theta/2)

Printf("The value of ax (major semi-axis) is: %g", ax);
Printf("The value of ay (minor semi-axis) is: %g", ay);

rho = ay*ay/ax;
Printf("Computed curvature radius rho: %g", rho);

h_plate = scale * 0.0400; // coarse mesh size
h_cavity = rho / 3.0;     // fine mesh size, adaptive to curvature

// Rectangular plate
Rectangle(1) = {-Lx/2, -Ly/2, 0, Lx, Ly};

// Elliptical hole
Ellipse(10) = {0, 0, 0, ax, ay};

Curve Loop(10) = {10};
Plane Surface(10) = {10};

// Subtract the ellipse from the rectangle
BooleanDifference{ Surface{1}; Delete; }{ Surface{10}; Delete; }

// Physical groups
Physical Surface("solid") = {1};
Physical Curve("cavity") = {10};  // Internal elliptical boundary
Physical Curve("ymin") = {11};    // Bottom
Physical Curve("xmin") = {12};    // Left
Physical Curve("xmax") = {13};    // Right
Physical Curve("ymax") = {14};    // Top

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {10};
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 0.02;
Field[2].DistMax = Lx/2;

Background Field = 2;

// Mesh generation options
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Algorithm = 6; // Frontal-Delaunay for better quality in 2D
// Mesh 2;
// Save "mesh.msh";
