// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2 2D elliptical cavities in a rectangular plate (cross-section)
//  Coherent with a semi-dihedral angle of 50 degrees
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

scale = 1; // (micrometers)

// Parameters
ax = 0.100;          // Semi-axis X (100 nm)
Fc_target = 0.40;    // Target coverage (40%)

Lx = 2 * ax * Sqrt(Pi / Fc_target);
Ly = Lx;
x_pos = Lx / 4; 

Printf("Target Fc: %g", Fc_target);
Printf("Lx: %g", Lx);

theta = 50 * Pi / 180; // Semi-dihedral angle in radiants 
ay = ax * (1 - Cos(theta)) / Sin(theta); // Semi-axis Y coherent

Printf("The value of ax (major semi-axis) is: %g", ax);
Printf("The value of ay (minor semi-axis) is: %g", ay);

rho = ay*ay/ax;
Printf("Computed curvature radius rho: %g", rho);

h_plate = scale * 0.00400; // coarse mesh size
h_cavity = rho / 3.0;      // fine mesh size, adaptive to curvature

// Rectangular plate
Rectangle(1) = {-Lx/2, -Ly/2, 0, Lx, Ly};

// Elliptical hole
Ellipse(10) = {-Lx/4, 0, 0, ax, ay};
Curve Loop(10) = {10};
Plane Surface(10) = {10};

Ellipse(11) = {+Lx/4, 0, 0, ax, ay};
Curve Loop(11) = {11};
Plane Surface(11) = {11};

// Subtract the ellipse from the rectangle
BooleanDifference{ Surface{1}; Delete; }{ Surface{10, 11}; Delete; }

// Physical groups
Physical Surface("solid") = {1};
Physical Curve("cavity") = {10, 11};  // Internal elliptical boundary
Physical Curve("ymin") = {12};    // Bottom
Physical Curve("xmin") = {13};    // Left
Physical Curve("xmax") = {14};    // Right
Physical Curve("ymax") = {15};    // Top

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {10, 11};
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
