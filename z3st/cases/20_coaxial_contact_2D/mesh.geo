// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for two 2D coaxial cylinders
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


SetFactory("OpenCASCADE");

r_int = 0.055;
r_ext = 0.065; 

Circle(1) = {0, 0, 0, r_int};
Circle(2) = {0, 0, 0, r_ext};

Curve Loop(1) = {1};
Plane Surface(1) = {1};

Curve Loop(2) = {2};
Curve Loop(3) = {1};
Plane Surface(2) = {2, 3};

BooleanFragments{ Surface{1, 2}; Delete; }{}

Physical Surface("cyl_1", 1) = {1};
Physical Surface("cyl_2", 2) = {2};
Physical Curve("contact_interface", 10) = {1};

Mesh 2;
Save "mesh.msh";