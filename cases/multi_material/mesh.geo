// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a cube LxLxL with a structured mesh
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// Use built-in kernel instead of OpenCASCADE (important for transfinite meshing)
SetFactory("Built-in");

L = 1.0;  // Side 
H = 1.0;  // Height
n = 20;   // Number of divisions per edge (adjust as needed)

// First solid

// Define corner points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};
Point(5) = {0, 0, H, 1.0};
Point(6) = {L, 0, H, 1.0};
Point(7) = {L, L, H, 1.0};
Point(8) = {0, L, H, 1.0};

// Define edges
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

// Define faces
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

// Define volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Second solid

dx = L+0.01;  // offset in x

// Define corner points
Point(9) = {dx + 0, 0, 0, 1.0};
Point(10) = {dx + L, 0, 0, 1.0};
Point(11) = {dx + L, L, 0, 1.0};
Point(12) = {dx + 0, L, 0, 1.0};
Point(13) = {dx + 0, 0, H, 1.0};
Point(14) = {dx + L, 0, H, 1.0};
Point(15) = {dx + L, L, H, 1.0};
Point(16) = {dx + 0, L, H, 1.0};

// Define edges
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};
Line(17) = {9, 13};
Line(18) = {10, 14};
Line(19) = {11, 15};
Line(20) = {12, 16};
Line(21) = {13, 14};
Line(22) = {14, 15};
Line(23) = {15, 16};
Line(24) = {16, 13};

// Define faces
Line Loop(7) = {13, 14, 15, 16};     // bottom
Plane Surface(7) = {7};

Line Loop(8) = {17, 21, -18, -13};   // face1
Plane Surface(8) = {8};

Line Loop(9) = {18, 22, -19, -14};   // face2
Plane Surface(9) = {9};

Line Loop(10) = {19, 23, -20, -15};  // face3
Plane Surface(10) = {10};

Line Loop(11) = {20, 24, -17, -16};  // face4
Plane Surface(11) = {11};

Line Loop(12) = {21, 22, 23, 24};    // top
Plane Surface(12) = {12};

// Define volume
Surface Loop(2) = {7, 8, 9, 10, 11, 12};
Volume(2) = {2};

// Transfinite mesh for first cube
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12} = n Using Progression 1;
Transfinite Curve {5, 6, 7, 8} = n Using Progression 1;
Transfinite Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};

// Transfinite mesh for second cube
Transfinite Curve {13, 14, 15, 16, 21, 22, 23, 24} = n Using Progression 1;
Transfinite Curve {17, 18, 19, 20} = n Using Progression 1;
Transfinite Surface {7, 8, 9, 10, 11, 12};
Transfinite Volume {2};

// Physical groups for boundary conditions (first cube)
Physical Surface("bottom1") = {1};
Physical Surface("face_yy_0_1") = {2};
Physical Surface("face_xx_1_1") = {3};
Physical Surface("face_yy_1_1") = {4};
Physical Surface("face_xx_0_1") = {5};
Physical Surface("top1")    = {6};
Physical Volume("fuel")    = {1};

// Physical groups for second cube
Physical Surface("bottom2") = {7};
Physical Surface("face_yy_0_2") = {8};
Physical Surface("face_xx_1_2") = {9};
Physical Surface("face_yy_1_2") = {10};
Physical Surface("face_xx_0_2") = {11};
Physical Surface("top2")    = {12};
Physical Volume("cladding")    = {2};

// Generate structured mesh
Mesh.RecombineAll = 1;
Mesh 3;

// Save mesh
Save "mesh.msh";
