// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a cube LxLxL with a structured mesh
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// Use built-in kernel instead of OpenCASCADE (important for transfinite meshing)
SetFactory("Built-in");

L = 1.0;  // Cube side length
n = 20;   // Number of divisions per edge (adjust as needed)

// Define cube corner points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};
Point(5) = {0, 0, L, 1.0};
Point(6) = {L, 0, L, 1.0};
Point(7) = {L, L, L, 1.0};
Point(8) = {0, L, L, 1.0};

// Define cube edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

// Define cube faces
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 9, -6, -1};
Plane Surface(2) = {2};

Line Loop(3) = {6, 10, -7, -2};
Plane Surface(3) = {3};

Line Loop(4) = {7, 11, -8, -3};
Plane Surface(4) = {4};

Line Loop(5) = {8, 12, -5, -4};
Plane Surface(5) = {5};

Line Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};

// Define cube volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Transfinite structure: ensure structured meshing
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12} = n Using Progression 1;
Transfinite Curve {5, 6, 7, 8} = n Using Progression 1;
Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};

// Define physical groups for boundary conditions
Physical Surface("bottom") = {1};
Physical Surface("face1") = {2};
Physical Surface("face2") = {3};
Physical Surface("face3") = {4};
Physical Surface("face4") = {5};
Physical Surface("top") = {6};
Physical Volume("cube") = {1};

// Generate structured mesh
Mesh 3;

// Save mesh file
Save "cube.msh";
