// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) rod
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.0045;     // pellet outer radius (m)
r_2_i = 0.004565;   // clad inner radius (m)
r_2_o = 0.005315;   // clad outer radius (m)
h1    = 0.010;      // axial segment height (m)
h2    = 0.010;      // axial segment height (m)

n_r1 = 25;          // radial divisions, pellet
n_r2 = 7;           // radial divisions, clad
n_z  = 31;          // axial divisions (shared)

// --. fuel : rectangle [0, r_1_o] x [0, h1] --..
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

// --. clad : rectangle [r_2_i, r_2_o] x [0, h2] --..
Point(5) = {r_2_i, 0, 0};
Point(6) = {r_2_o, 0, 0};
Point(7) = {r_2_o, h2, 0};
Point(8) = {r_2_i, h2, 0};
Line(5) = {5, 6};   // bottom_2 (z = 0)
Line(6) = {6, 7};   // outer_2  (r = r_2_o)
Line(7) = {7, 8};   // top_2    (z = 2)
Line(8) = {8, 5};   // inner_2  (r = r_2_i)
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// --- structured quad mesh ---
Transfinite Line {1, 3} = n_r1;
Transfinite Line {5, 7} = n_r2;
Transfinite Line {2, 4, 6, 8} = n_z;
Transfinite Surface {1};
Transfinite Surface {2};
Recombine Surface {1, 2};

// --- physical groups (ids must match geometry.yaml) ---
Physical Surface("fuel", 1) = {1};
Physical Surface("clad", 2) = {2};

Physical Curve("bottom_1", 1)  = {1};
Physical Curve("bottom_2", 2)  = {5};
Physical Curve("lateral_1", 3) = {2};
Physical Curve("outer_2", 4)   = {6};
Physical Curve("inner_2", 5)   = {8};
Physical Curve("top_2", 6)     = {7};
Physical Curve("top_1", 7)     = {3};
Physical Curve("axis_1", 8)    = {4};

Mesh.ElementOrder = 1;
