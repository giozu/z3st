// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  LFR reactor vessel, 2D axisymmetric (r = x, z = y):
//  cylindrical shell + hemispherical bottom head, same wall thickness.
//
//  Size: ELSY, Waltar-Todd-Tsvetkov (2012) Table 18.1, H/D = 9.5/12.5 m.
//  Thickness: 50 mm, an assumption (the order of the Super Phenix tank,
//  same book, Sect. 12.3.1); Table 18.1 does not give it.
//
//  Inner bottom pole at z = 0, equator (head/shell junction) at z = Ri,
//  top edge (support ring) at z = H.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 6.25;    // inner radius (m)
t  = 0.05;    // wall thickness (m)
H  = 9.5;     // inner height, pole to top edge (m)
Ro = Ri + t;

n_t   = 9;    // nodes through the wall
n_arc = 601;  // nodes along the head meridian (~16 mm)
n_cyl = 201;  // nodes along the shell (~16 mm)

Point(1) = {0,  Ri,     0};   // centre of the sphere
Point(2) = {0,  0,      0};   // inner pole
Point(3) = {0,  -t,     0};   // outer pole
Point(4) = {Ri, Ri,     0};   // inner equator
Point(5) = {Ro, Ri,     0};   // outer equator
Point(6) = {Ri, H,      0};   // inner top
Point(7) = {Ro, H,      0};   // outer top

Circle(1) = {2, 1, 4};        // inner head arc
Circle(2) = {3, 1, 5};        // outer head arc
Line(3)   = {3, 2};           // wall at the pole (on the axis)
Line(4)   = {4, 5};           // head/shell junction
Line(5)   = {4, 6};           // inner shell
Line(6)   = {5, 7};           // outer shell
Line(7)   = {6, 7};           // top edge

Curve Loop(1) = {3, 1, 4, -2};
Plane Surface(1) = {1};
Curve Loop(2) = {-4, 5, 7, -6};
Plane Surface(2) = {2};

Transfinite Line {3, 4, 7} = n_t;
Transfinite Line {1, 2} = n_arc;
Transfinite Line {5, 6} = n_cyl;
Transfinite Surface {1, 2};
Recombine Surface {1, 2};

Physical Surface("steel", 1) = {1, 2};
Physical Curve("inner", 10) = {1, 5};
Physical Curve("outer", 11) = {2, 6};
Physical Curve("symmetry_axis", 12) = {3};
Physical Curve("top_edge", 13) = {7};
