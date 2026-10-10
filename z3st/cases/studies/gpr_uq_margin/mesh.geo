SetFactory("OpenCASCADE");

Ri = 0.000;
Ro = 0.0041;
Lz = 0.020;

nx = 31;
ny = 7;

Point(1) = {Ri, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = nx;
Transfinite Line {2, 4} = ny;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("volume", 10) = {1};
Physical Curve("inner", 1) = {4};
Physical Curve("outer", 2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top", 4) = {3};
