// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a unit cube
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 1.0;
Ly = 1.0;
Lz = 1.0;

n = 5; // Mesh divisions

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};
Point(5) = {0, 0, Lz, 1.0};
Point(6) = {Lx, 0, Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0, Ly, Lz, 1.0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Surfaces
Line Loop(1) = {1, 2, 3, 4}; Plane Surface(1) = {1}; // zmin
Line Loop(2) = {5, 6, 7, 8}; Plane Surface(2) = {2}; // zmax
Line Loop(3) = {1, 10, -5, -9}; Plane Surface(3) = {3}; // ymin
Line Loop(4) = {2, 11, -6, -10}; Plane Surface(4) = {4}; // xmax
Line Loop(5) = {3, 12, -7, -11}; Plane Surface(5) = {5}; // ymax
Line Loop(6) = {4, 9, -8, -12}; Plane Surface(6) = {6}; // xmin

// Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Transfinite mesh
Transfinite Curve "*" = n;
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";

// Physical Groups
Physical Volume("grain", 1) = {1};
Physical Surface("zmin", 10) = {1};
Physical Surface("zmax", 11) = {2};
Physical Surface("ymin", 12) = {3};
Physical Surface("xmax", 13) = {4};
Physical Surface("ymax", 14) = {5};
Physical Surface("xmin", 15) = {6};

Mesh 3;
Save "mesh.msh";
