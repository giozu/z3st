// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for two 2D coaxial cylinders
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.05;
r_2_i = 0.06;
r_2_o = 0.065;

n_f = 81; // circumferential divisions

Point(100) = {0, 0, 0};

// (1)
Circle(1) = {0, 0, 0, r_1_o};
Curve Loop(1) = {1};
Plane Surface(1) = {1};

// (2)
Circle(2) = {0, 0, 0, r_2_i};
Circle(3) = {0, 0, 0, r_2_o};
Curve Loop(2) = {2};
Curve Loop(3) = {3};
Plane Surface(2) = {3, 2}; // in-between surface

Transfinite Curve {1, 2, 3} = n_f;

Physical Surface("cyl_1", 1) = {1};
Physical Surface("cyl_2", 2) = {2};

Physical Curve("interface_1", 10) = {1};
Physical Curve("interface_2", 11) = {2};
Physical Curve("outer_boundary", 12) = {3};

Mesh 2;
Save "mesh.msh";