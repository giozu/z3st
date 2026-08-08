// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) block
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.0045;     // pellet outer radius (m)
h1    = 0.010;      // axial segment height (m)

n_r1 = 25;          // radial divisions
n_z  = 31;          // axial divisions (shared)

// --. rectangle [0, r_1_o] x [0, h1] --..
Point(1) = {0,     0, 0};
Point(2) = {r_1_o, 0, 0};
Point(3) = {r_1_o, h1, 0};
Point(4) = {0,     h1, 0};
Line(1) = {1, 2};   // bottom_1  (z = 0)
Line(2) = {2, 3};   // lateral_1 (r = r_1_o)
Line(3) = {3, 4};   // top_1     (z = h1)
Line(4) = {4, 1};   // axis_1    (r = 0)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- structured quad mesh ---
Transfinite Line {1, 3} = n_r1;
Transfinite Line {2, 4, 6, 8} = n_z;
Transfinite Surface {1};
Recombine Surface {1, 2};

// --- physical groups (ids must match geometry.yaml) ---
Physical Surface("fuel", 1) = {1};

Physical Curve("bottom_1", 1)  = {1};
Physical Curve("lateral_1", 2) = {2};
Physical Curve("top_1", 3)     = {3};
Physical Curve("axis_1", 4)    = {4};

Mesh.ElementOrder = 1;
