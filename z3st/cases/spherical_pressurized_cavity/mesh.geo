// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Elliptical cavity inside a rectangular box
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// --- Parameters ---
Lx = 1.0;          // box length  (m)
Ly = 1.0;          // box width   (m)
Lz = 1.0;          // box height  (m)

ax  = 0.04;        // ellipse semi-axis X (m)
ay  = 0.04;        // ellipse semi-axis Y (m)
az  = 0.04;        // ellipse semi-axis Y (m)

h_box = 0.0500;    // target element size (m)
h_sph = 0.0025;    // target element size (m)

Point(1) = { 0, 0, 0 };       // center

Box(1) = { -Lx/2, -Ly/2, -Lz/2,  Lx, Ly, Lz };
Sphere(2) = { 0, 0, 0, 1.0 };

Dilate {{0, 0, 0}, {ax, ay, az}} { Volume{2}; }

BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// --- Physical groups ---
Physical Volume("solid") = {1};
Physical Surface("cavity") = {7};
Physical Surface("xmin") = {8};
Physical Surface("ymin") = {9};
Physical Surface("zmax") = {10};
Physical Surface("ymax") = {11};
Physical Surface("zmin") = {12};
Physical Surface("xmax") = {13};

// --- Global mesh options ---
Mesh.CharacteristicLengthMin = h_sph;   // global minimum element size
Mesh.CharacteristicLengthMax = h_box;   // global maximum element size
Mesh.Optimize = 1;                      // optimize mesh quality
Mesh.ElementOrder = 1;                  // linear elements (order 1)

// --- Local mesh refinement using fields ---
Field[1] = Distance;
Field[1].SurfacesList = {7};            // surface id of the cavity
Field[1].NumPointsPerCurve = 100;       // accuracy of distance computation

// Field[2]: map distance → element size
Field[2] = Threshold;
Field[2].InField = 1;                   // use distance field[1]
Field[2].SizeMin = h_sph;               // fine size at cavity surface
Field[2].SizeMax = h_box;               // coarse size far away
Field[2].DistMin = 0.0;                 // up to this distance → SizeMin
Field[2].DistMax = Lx / 4;                 // beyond this distance → SizeMax

// set background mesh field
Background Field = 2;

Mesh 3;
Save "mesh.msh";
