// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Lz = 0.050;
r  = 0.040;
r_hole = 0.004;

n_z = 10;

// 2D quarter annulus
Point(0) = {0, 0, 0};
Point(1) = {r, 0, 0};
Point(2) = {0, r, 0};
Circle(1) = {1, 0, 2};
Line(2) = {0, 1};
Line(3) = {2, 0};
Curve Loop(1) = {2, 1, 3};
Plane Surface(1) = {1};

// inner hole
Disk(2) = {3*r_hole, 3*r_hole, 0, r_hole, r_hole};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

// Extrude
out[] = Extrude {0, 0, Lz} {
  Surface{1};
  Layers{n_z};
};

// Physical groups
Physical Surface("zmin") = {1};
Physical Surface("ymin") = {2};
Physical Surface("xmin") = {3};
Physical Surface("arc") = {4};
Physical Surface("hole") = {5};
Physical Surface("zmax") = {6};
Physical Volume("volume") = {1};

// Mesh
Mesh.CharacteristicLengthMax = 0.003;
Mesh 3;

Save "mesh.msh";
