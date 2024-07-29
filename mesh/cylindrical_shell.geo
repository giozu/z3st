// Gmsh project created on Mon Jul 29 17:10:09 2024
SetFactory("OpenCASCADE");

//+
Circle(1) = {0, 0, 0, 5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 4.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Extrude {0, 0, 4} {
  Surface{1}; Layers {5}; 
}
//+
Physical Surface("top", 7) = {4};
//+
Physical Surface("bottom", 8) = {1};
//+
Physical Surface("external_lateral", 9) = {2};
//+
Physical Surface("internal_lateral", 10) = {3};
//+
Physical Volume("domain", 11) = {1};
