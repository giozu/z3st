// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for two 2D coaxial cylinders
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

// 1D mesh for cluster dynamics
// x-coordinate represents cluster size n

// Parameters
n_ele = 100;     // number of elements
L = 100.0;       // domain length

// Points
Point(1) = {1, 0, 0};
Point(2) = {L, 0, 0};

// Line
Line(1) = {1, 2};

// Mesh control
Transfinite Line(1) = n_ele + 1;

// Physical groups (explicit tag numbers)
Physical Point("left") = {1};  // left boundary
Physical Point("right") = {2};  // right boundary
Physical Curve("domain") = {1};  // domain
