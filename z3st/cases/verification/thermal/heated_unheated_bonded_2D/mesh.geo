// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: slab of a heated layer (0 < x < a) bonded to an unheated one
//  (a < x < Lx). One conforming structured mesh: the interface nodes
//  belong to both materials.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

a  = 0.02;  // interface (m)
Lx = 0.10;  // slab thickness (m)
Ly = 0.02;  // slab height (m)

// Coarse on purpose: the spurious source an unheated element picks up next
// to the interface scales with the element width.
n_h = 5;    // nodes across the heated layer
n_c = 17;   // nodes across the unheated layer
n_y = 3;    // nodes along y

Point(1) = {0,  0,  0};
Point(2) = {a,  0,  0};
Point(3) = {Lx, 0,  0};
Point(4) = {Lx, Ly, 0};
Point(5) = {a,  Ly, 0};
Point(6) = {0,  Ly, 0};

Line(1) = {1, 2};  // ymin, heated
Line(2) = {2, 3};  // ymin, cold
Line(3) = {3, 4};  // xmax
Line(4) = {4, 5};  // ymax, cold
Line(5) = {5, 6};  // ymax, heated
Line(6) = {6, 1};  // xmin
Line(7) = {2, 5};  // interface x = a

Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Transfinite Line {1, 5} = n_h;
Transfinite Line {2, 4} = n_c;
Transfinite Line {3, 6, 7} = n_y;
Transfinite Surface {1, 2};
Recombine Surface {1, 2};

Physical Surface("heated", 10) = {1};
Physical Surface("cold", 11) = {2};
Physical Curve("xmin", 1) = {6};
Physical Curve("xmax", 2) = {3};
Physical Curve("ymin", 3) = {1, 2};
Physical Curve("ymax", 4) = {4, 5};
