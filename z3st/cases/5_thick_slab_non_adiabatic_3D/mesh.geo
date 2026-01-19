// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a box Lx Ly Lz with a structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 0.400; // (m)
Ly = 2.000; // (m)
Lz = 2.000; // (m)

nx = 81;  
ny = 11;
nz = 11;

// Corner points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};
Point(5) = {0, 0, Lz, 1.0};
Point(6) = {Lx, 0, Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0, Ly, Lz, 1.0};

// Edges
Line(1) = {1, 2}; // x-direction
Line(2) = {2, 3}; // y-direction
Line(3) = {3, 4}; // x-direction
Line(4) = {4, 1}; // y-direction
Line(5) = {1, 5}; // z-direction
Line(6) = {2, 6}; // z-direction
Line(7) = {3, 7}; // z-direction
Line(8) = {4, 8}; // z-direction
Line(9) = {5, 6}; // x-direction
Line(10) = {6, 7}; // y-direction
Line(11) = {7, 8}; // x-direction
Line(12) = {8, 5}; // y-direction

// Faces (Loops)
Line Loop(1) = {1, 2, 3, 4};    Plane Surface(1) = {1}; // zmin
Line Loop(2) = {5, 9, -6, -1};  Plane Surface(2) = {2}; // ymin
Line Loop(3) = {6, 10, -7, -2}; Plane Surface(3) = {3}; // xmax
Line Loop(4) = {7, 11, -8, -3}; Plane Surface(4) = {4}; // ymax
Line Loop(5) = {8, 12, -5, -4}; Plane Surface(5) = {5}; // xmin
Line Loop(6) = {9, 10, 11, 12}; Plane Surface(6) = {6}; // zmax

// Cube volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Curve {1, 3, 9, 11} = nx Using Progression 1;
Transfinite Curve {2, 4, 10, 12} = ny Using Progression 1;
Transfinite Curve {5, 6, 7, 8} = nz Using Progression 1;

Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};

Physical Surface("zmin") = {1};
Physical Surface("ymin") = {2};
Physical Surface("xmax") = {3};
Physical Surface("ymax") = {4};
Physical Surface("xmin") = {5};
Physical Surface("zmax") = {6};
Physical Volume("steel") = {1};

Mesh.RecombineAll = 1;
Mesh 3;
Save "mesh.msh";