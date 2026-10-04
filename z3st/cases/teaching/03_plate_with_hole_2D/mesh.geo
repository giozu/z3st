// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a quarter of a plate with a circular hole
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

L = 0.200;        // half-side of the plate (m)
a = 0.010;        // hole radius (m)
h_hole = a/40;    // element size on the hole edge
h_far  = L/15;    // element size far from it

Rectangle(1) = {0, 0, 0, L, L};
Disk(2) = {0, 0, 0, a};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

e = 1e-6;
hole[] = Curve In BoundingBox{-e, -e, -e, a+e, a+e, e};
xmin[] = Curve In BoundingBox{-e, a-e, -e, e, L+e, e};
xmax[] = Curve In BoundingBox{L-e, -e, -e, L+e, L+e, e};
ymin[] = Curve In BoundingBox{a-e, -e, -e, L+e, e, e};
ymax[] = Curve In BoundingBox{-e, L-e, -e, L+e, L+e, e};

// Fine on the hole, where the stress is concentrated, coarse far away.
Field[1] = Distance;
Field[1].CurvesList = {hole[]};
Field[1].Sampling = 200;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_hole;
Field[2].SizeMax = h_far;
Field[2].DistMin = 0.2*a;
Field[2].DistMax = 6*a;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Labels match geometry.yaml.
Physical Surface("plate", 1) = {1};
Physical Curve("hole", 2)    = {hole[]};
Physical Curve("xmin", 3)    = {xmin[]};
Physical Curve("xmax", 4)    = {xmax[]};
Physical Curve("ymin", 5)    = {ymin[]};
Physical Curve("ymax", 6)    = {ymax[]};
