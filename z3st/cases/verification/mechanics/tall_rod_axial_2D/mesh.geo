// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) rod
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 0.000;   // inner radius (m)
Ro = 0.005;   // outer radius (m)
Lz = 0.400;   // rod height   (m)

nx = 12;      // radial divisions
ny = 161;     // axial divisions

Point(1) = {Ri, 0,  0};
Point(2) = {Ro, 0,  0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2};   // bottom (z = 0)
Line(2) = {2, 3};   // outer  (r = Ro)
Line(3) = {3, 4};   // top    (z = Lz)
Line(4) = {4, 1};   // axis   (r = 0)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = nx;
Transfinite Line {2, 4} = ny;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("fuel", 10) = {1};
Physical Curve("axis", 1)   = {4};
Physical Curve("outer", 2)  = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top", 4)    = {3};
