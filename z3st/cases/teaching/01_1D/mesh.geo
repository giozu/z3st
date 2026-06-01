// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 1D bar
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 1.000;     // bar length (m)

nx = 4;         // nodes along x -> (nx - 1) line elements

// Endpoints
Point(1) = {0,  0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};

// Bar (1D element collection)
Line(1) = {1, 2};

// Structured 1D mesh
Transfinite Curve {1} = nx Using Progression 1;

// Physical groups
Physical Point("xmin") = {1};
Physical Point("xmax") = {2};
Physical Curve("steel") = {1};
