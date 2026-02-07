// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  3D lenticular bubble
//  Coherent with a semi-dihedral angle of 50 degrees
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// Parameters
Lx = 0.5;          // plate length (m)
ax = 0.1072;       // ellipse semi-axis X (half-width on grain boundary)

h_plate = 0.0400;  // coarse mesh size
h_cavity = 0.0015; // fine mesh size at bubble tip

theta = 50 * Pi / 180;
R = ax / Sin(theta);    // Curvature radius of the bubble 
d = R * Cos(theta);     // Distance of the bubble center from the median plane

Sphere(2) = {0, 0, -d, R};
Sphere(3) = {0, 0,  d, R};

// Intersection of the two spheres
BooleanIntersection(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };

// Physical groups
Physical Surface("solid") = {1};
Physical Curve("cavity") = {10};  // Internal elliptical boundary
Physical Curve("xmin") = {11};    // Left
Physical Curve("xmax") = {13};    // Right
Physical Curve("ymin") = {12};    // Bottom
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
Mesh 2;
Save "mesh.msh";
