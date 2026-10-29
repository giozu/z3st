// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 1D alloy slab, graded toward the salt-facing end (xmax)
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 150.0e-6;  // alloy depth (m)

nx = 601;       // nodes -> 600 elements, 0.71 um at xmin, 0.05 um at xmax

Point(1) = {0,  0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};

Line(1) = {1, 2};

Transfinite Curve {1} = nx Using Progression 0.99557;

Physical Point("xmin") = {1};
Physical Point("xmax") = {2};
Physical Curve("alloy") = {1};
