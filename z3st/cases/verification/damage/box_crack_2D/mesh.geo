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

// Recombine Surface {1};

// Physical groups
Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3};
Physical Curve("xmin") = {4,5};
Physical Curve("crack") = {6};
Physical Surface("steel") = {1};

// Mesh 2;
// Save "mesh.msh";

// crack-path band ahead of the pre-crack: elements below lc so the crack can propagate
h_path = lc_damage / 2.5;
w_path = 2 * lc_damage;
Field[1] = Box;
Field[1].VIn = h_path;
Field[1].VOut = h_coarse;
Field[1].XMin = Dn;
Field[1].XMax = Lx;
Field[1].YMin = Ly/2 - w_path;
Field[1].YMax = Ly/2 + w_path;
Field[1].Thickness = 0.05;
Background Field = 1;
