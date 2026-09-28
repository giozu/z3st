// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for the cohesive 1D bar (Vicentini et al. 2026, Fig. 2)
//
//  Bar of length 2L clamped at the left end, prescribed displacement at
//  the right end. The element count is odd so that no node sits at the
//  mid-point: a node there would let the localizing element be picked by
//  the mesh rather than by the solution.
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 2.0e-3;    // bar length 2L (m), L = 1 mm

nx = 202;       // nodes -> 201 line elements, h = ell/5

Point(1) = {0,  0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};

Line(1) = {1, 2};

Transfinite Curve {1} = nx Using Progression 1;

Physical Point("xmin") = {1};
Physical Point("xmax") = {2};
Physical Curve("bar") = {1};
