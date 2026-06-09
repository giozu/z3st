// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a box Lx Ly with a structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 0.100;
Ly = 1.000;

nx = 11;
ny = 41;

// Corner points
Point(1) = {0,  0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0,  Ly, 0, 1.0};

// Edges
Line(1) = {1, 2}; // Bottom (y-min)
Line(2) = {2, 3}; // Right  (x-max)
Line(3) = {3, 4}; // Top    (y-max)
Line(4) = {4, 1}; // Left   (x-min)

// Surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Transfinite structure for quad mesh ---
Transfinite Curve {1, 3} = nx Using Progression 1;
Transfinite Curve {2, 4} = ny Using Progression 1;

Transfinite Surface {1};
Recombine Surface {1};

// --- Physical Groups for Z3ST ---
Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3};
Physical Curve("xmin") = {4};
Physical Surface("steel") = {1};

// Meshing settings
// Mesh 2;

// Save mesh file
// Save "mesh.msh";
