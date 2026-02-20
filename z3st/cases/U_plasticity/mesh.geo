// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a box Lx Ly with a structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 0.100;
Ly = 0.100;

nxy = 41;   // Lx, Ly

// Corner points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};

// Edges
Line(1) = {1, 2}; // ymin
Line(2) = {2, 3}; // xmax
Line(3) = {3, 4}; // ymax
Line(4) = {4, 1}; // xmin

// Surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Transfinite structure: ensure structured meshing
Transfinite Curve {1, 2, 3, 4} = nxy Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

// Define physical groups for boundary conditions
Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3};
Physical Curve("xmin") = {4};
Physical Surface("steel") = {1};

// Generate structured mesh
Mesh.RecombineAll = 1;
// Mesh 2;

// Save mesh file
// Save "mesh.msh";
