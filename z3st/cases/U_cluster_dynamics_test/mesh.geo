// 1D mesh for cluster dynamics
// x-coordinate represents cluster size n

// Parameters
n_ele = 10000;  // number of elements
L = 100.0;     // domain length

// Points
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};

// Line
Line(1) = {1, 2};

// Mesh control
Transfinite Line(1) = n_ele + 1;

// Physical groups (explicit tag numbers)
Physical Point(1) = {1};  // boundary_left
Physical Point(2) = {2};  // boundary_right
Physical Curve(3) = {1};  // domain
