// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 3D thick-walled cylinder, radial-structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 0.02;  // Inner radius (m)
Ro = 0.03;  // Outer radius (m)
Lz = 0.01;  // Height

num_points_radial = 10;      // Number of points in the radial direction
num_points_circ_qtr = 15;    // Number of points on each 90-degree arc
num_layers_height = 21;      // Number of element layers along the height

// 2D slice
Point(1) = {0, 0, 0};        // Center point

// Points on inner circle (at 0, 90, 180, 270 degrees)
Point(2) = {Ri, 0, 0};
Point(3) = {0, Ri, 0};
Point(4) = {-Ri, 0, 0};
Point(5) = {0, -Ri, 0};

// Points on outer circle (at 0, 90, 180, 270 degrees)
Point(6) = {Ro, 0, 0};
Point(7) = {0, Ro, 0};
Point(8) = {-Ro, 0, 0};
Point(9) = {0, -Ro, 0};

// Define arcs for inner circle (Tags 1-4)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Define arcs for outer circle (Tags 5-8)
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Radial lines connecting inner and outer circles (Tags 9-12)
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Define curve loops for the four annular sectors
Curve Loop(1) = {1, 10, -5, -9};      // First quadrant
Curve Loop(2) = {2, 11, -6, -10};     // Second quadrant
Curve Loop(3) = {3, 12, -7, -11};     // Third quadrant
Curve Loop(4) = {4, 9, -8, -12};      // Fourth quadrant

// Define the four plane surfaces from the curve loops
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

// 2D Structured Mesh
Transfinite Line {9, 10, 11, 12} = num_points_radial;

// Circumferential arcs (inner and outer):
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8} = num_points_circ_qtr;

// Apply transfinite algorithm to surfaces and recombine into quadrilaterals
Transfinite Surface {1, 2, 3, 4};
Recombine Surface {1, 2, 3, 4};

// Extrude the 2D surfaces to create 3D volumes.
out[] = Extrude {0, 0, Lz} {
    Surface {1, 2, 3, 4};
    Layers{num_layers_height};
    Recombine;
};

Coherence;

Physical Surface("bottom") = {1, 2, 3, 4};
Physical Surface("top") = {out[0], out[6], out[12], out[18]};
Physical Surface("inner_radius") = {out[2], out[8], out[14], out[20]};
Physical Surface("outer_radius") = {out[4], out[10], out[16], out[22]};
Physical Volume("steel")         = {out[1],  out[7],  out[13], out[19]};

// --- Mesh Generation ---
Mesh.ElementOrder = 1;      // Use linear elements
Mesh.Optimize = 1;          // Improve mesh quality

// Generate the 3D mesh
Mesh 3;

// Save the mesh
Save "mesh.msh";
