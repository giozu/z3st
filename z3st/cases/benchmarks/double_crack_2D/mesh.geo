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

lc_damage = 0.002;
h_fine = lc_damage / 5.0;
h_coarse = 0.01;

// Corner points
Point(1) = {0,  0,  0, h_coarse};
Point(2) = {Lx, 0,  0, h_coarse};
Point(3) = {Lx, Ly, 0, h_coarse};
Point(4) = {0,  Ly, 0, h_coarse};
Point(5) = {Dn, Ly/2, 0, h_fine};
Point(6) = {0,  Ly/2, 0, h_fine};
Point(7) = {Lx-Dn, Ly/2, 0, h_fine};
Point(8) = {Lx,  Ly/2, 0, h_fine};

// Edges
Line(1) = {1, 2};       // y-min
Line(2) = {2, 8};       // x-max
Line(3) = {8, 7};       // CRACK
Line(4) = {8, 3};       // x-max
Line(5) = {3, 4};       // y-max
Line(6) = {4, 6};       // x-min
Line(7) = {6, 5};       // CRACK
Line(8) = {6, 1};       // x-min

// Surface
Line Loop(1) = {1, 2, 3, -3, 4, 5, 6, 7, -7, 8};
Plane Surface(1) = {1};

// Physical groups
Physical Curve("ymin")  = {1};
Physical Curve("xmax")  = {2, 4};
Physical Curve("ymax")  = {5};
Physical Curve("xmin")  = {6, 8};
Physical Curve("crack") = {3, 7};
Physical Surface("steel") = {1};

// Mesh 2;
// Save "mesh.msh";
