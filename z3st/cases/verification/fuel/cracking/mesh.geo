// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric block
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 0.000;   // inner radius (m)
Ro = 0.0045;  // outer radius (m)
Lz = 0.010;   // height (m)

nx = 3;       // radial
ny = 5;       // axial

Point(1) = {Ri, 0,  0};
Point(2) = {Ro, 0,  0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2}; // bottom (z=0)
Line(2) = {2, 3}; // outer  (r=Ro)
Line(3) = {3, 4}; // top    (z=Lz)
Line(4) = {4, 1}; // inner  (r=0, axis)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = nx;
Transfinite Line {2, 4} = ny;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("fuel", 10) = {1};
Physical Curve("inner", 1)  = {4};
Physical Curve("outer", 2)  = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top", 4)    = {3};
