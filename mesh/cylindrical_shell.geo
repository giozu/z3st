// Gmsh project created on Mon Jul 29 17:10:09 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.8, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Extrude {0, 0, 4} {
  Curve{2}; Curve{1}; Surface{1}; // Layers {5};
}
//+
Physical Surface("top", 11) = {6};
//+
Physical Surface("bottom", 12) = {1};
//+
Physical Surface("inner", 13) = {3};
//+
Physical Surface("outer", 14) = {2};
//+
Physical Volume("volume", 15) = {1};
