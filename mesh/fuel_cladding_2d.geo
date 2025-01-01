SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1.8, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1.5, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 2, 0, 2*Pi};
//+
Curve Loop(1) = {3};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {2};
//+
Plane Surface(2) = {3};
//+
Physical Surface("cladding", 2) = {1};
//+
Physical Surface("fuel", 1) = {2};
//+
Recombine Surface {2};
//+
Recombine Surface {1};
