SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.8, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {2};
//+
Curve Loop(5) = {1};
//+
Curve Loop(6) = {2};
//+
Plane Surface(1) = {5, 6};
//+
Extrude {0, 0, 5} {
  Curve{2}; Curve{1}; Surface{1}; Point{2}; Point{3}; Layers {10};
}
//+
Physical Surface("top", 13) = {6};
//+
Physical Surface("bottom", 14) = {1};
//+
Physical Surface("outer", 15) = {2};
//+
Physical Surface("inner", 16) = {3};
//+
Physical Volume("shell", 17) = {1};
