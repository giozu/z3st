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

n = 20; // Number of divisions per Lx/y/z

// Corner points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};

// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Faces
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(10) = {10};

// Transfinite structure
Transfinite Curve {1} = n Using Progression 1.2;
Transfinite Curve {-3} = n Using Progression 1.2;
Transfinite Curve {2, 4} = n Using Progression 1;
Transfinite Surface {10};
Recombine Surface {10};

// Extrude
out[] = Extrude {0, 0, Lz} {  Surface{10}; Layers{n}; Recombine; };

// Physical groups for boundary conditions
Physical Surface("zmin") = {10};
Physical Surface("zmax") = {out[0]};
Physical Surface("ymin") = {out[2]};    // y = 0
Physical Surface("xmax") = {out[3]};    // x = Lx
Physical Surface("ymax") = {out[4]};    // y = Ly
Physical Surface("xmin") = {out[5]};    // x = 0
Physical Volume("steel") = {out[1]};

// Generate structured mesh
Mesh.RecombineAll = 1;
Mesh 3;

// Save mesh file
Save "mesh.msh";
