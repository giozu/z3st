// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a box Lx Ly with a structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 1.000;
Ly = 0.500;
Dn = 0.25;

nx = 81;
ny = 41;

// Corner points
Point(1) = {0,  0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0,  Ly, 0, 1.0};

Point(5) = {Dn, Ly/2, 0, 1.0};
Point(6) = {0, Ly/2, 0, 1.0};

// Edges
Line(1) = {1, 2}; // y-min
Line(2) = {2, 3}; // x-max
Line(3) = {3, 4}; // y-max
Line(4) = {4, 6}; // x-min
Line(5) = {6, 1}; // x-min
Line(6) = {6, 5};

// Surface
Line Loop(1) = {1, 2, 3, 4, 6, -6, 5};
Plane Surface(1) = {1};

// Transfinite structure for quad
Transfinite Curve {1, 2, 3} = nx Using Progression 1;
Transfinite Curve {4, 5, 6} = nx Using Progression 1;

// Transfinite Surface {1};
// Recombine Surface {1};

// Physical groups
Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3};
Physical Curve("xmin") = {4,5};
Physical Curve("crack") = {6};
Physical Surface("steel") = {1};

Mesh 2;
Save "mesh.msh";
