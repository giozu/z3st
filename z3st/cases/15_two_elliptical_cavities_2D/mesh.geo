// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2 2D elliptical cavities in a rectangular plate (cross-section)
//  Coherent with a semi-dihedral angle of 50 degrees
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

scale = 1; // Work in micrometers for stable OpenCASCADE booleans

// Parameters
ax = 11.2 * scale;  // Semi-axis X in micrometers
Fc_target = 0.40;    // Target coverage (40%)

Lx = 2 * ax * Sqrt(Pi / Fc_target);
Ly = Lx;
x_pos = Lx / 4; 

Printf("Target Fc: %g", Fc_target);
Printf("Lx: %g", Lx);

theta = 71.3 * Pi / 180; // Semi-dihedral angle in radiants 
ay = ax * (1 - Cos(theta)) / Sin(theta); // Semi-axis Y coherent

Printf("The value of ax (major semi-axis) is: %g", ax);
Printf("The value of ay (minor semi-axis) is: %g", ay);

rho = ay*ay/ax;
Printf("Computed curvature radius rho: %g", rho);

h_plate = 1.5 * scale;     // coarse mesh size (1.5 microns)
h_cavity = 0.15 * scale;    // extremely fine mesh size to prevent crack locking

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

// We capture the new curve IDs using Bounding Boxes to avoid OpenCASCADE renumbering issues
eps = 1e-4;
c_ymin[] = Curve In BoundingBox{-Lx/2-eps, -Ly/2-eps, -eps, Lx/2+eps, -Ly/2+eps, eps};
c_ymax[] = Curve In BoundingBox{-Lx/2-eps,  Ly/2-eps, -eps, Lx/2+eps,  Ly/2+eps, eps};
c_xmin[] = Curve In BoundingBox{-Lx/2-eps, -Ly/2-eps, -eps, -Lx/2+eps,  Ly/2+eps, eps};
c_xmax[] = Curve In BoundingBox{ Lx/2-eps, -Ly/2-eps, -eps,  Lx/2+eps,  Ly/2+eps, eps};
c_cavity[] = Curve In BoundingBox{-Lx/2+eps, -Ly/2+eps, -eps, Lx/2-eps, Ly/2-eps, eps};

// Physical groups
Physical Surface("solid") = {1}; // This is a surface, so it's just the tag
Physical Curve("cavity") = {c_cavity[]};
Physical Curve("ymin") = {c_ymin[]};
Physical Curve("xmin") = {c_xmin[]};
Physical Curve("xmax") = {c_xmax[]};
Physical Curve("ymax") = {c_ymax[]};

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {c_cavity[]};
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 6.0 * scale;
Field[2].DistMax = Lx/2;

Background Field = 2;

// Mesh generation options
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Algorithm = 6; // Frontal-Delaunay for better quality in 2D

// Scale the final generated mesh into meters for the SI physics solver
Mesh.ScalingFactor = 1e-6;
// Mesh 2;
// Save "mesh.msh";
