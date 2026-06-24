// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) segment
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.0041;     // pellet outer radius (m)        4.10 mm
r_2_i = 0.00413;    // clad inner radius (m)          4.13 mm  -> gap = 30 um
r_2_o = 0.00475;    // clad outer radius (m)          4.75 mm
h     = 0.010;      // axial segment height (m)       10 mm

n_r1 = 8;           // radial divisions, pellet
n_r2 = 4;           // radial divisions, clad
n_z  = 11;          // axial divisions (shared)

// --- cyl_1 : fuel pellet, rectangle [0, r_1_o] x [0, h] ---
Point(1) = {0,     0, 0};
Point(2) = {r_1_o, 0, 0};
Point(3) = {r_1_o, h, 0};
Point(4) = {0,     h, 0};
Line(1) = {1, 2};   // bottom_1  (z = 0)
Line(2) = {2, 3};   // lateral_1 (r = r_1_o, gap-facing)
Line(3) = {3, 4};   // top_1     (z = h)
Line(4) = {4, 1};   // axis_1    (r = 0, symmetry axis)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- cyl_2 : cladding, rectangle [r_2_i, r_2_o] x [0, h] ---
Point(5) = {r_2_i, 0, 0};
Point(6) = {r_2_o, 0, 0};
Point(7) = {r_2_o, h, 0};
Point(8) = {r_2_i, h, 0};
Line(5) = {5, 6};   // bottom_2 (z = 0)
Line(6) = {6, 7};   // outer_2  (r = r_2_o, coolant)
Line(7) = {7, 8};   // top_2    (z = h)
Line(8) = {8, 5};   // inner_2  (r = r_2_i, gap-facing)
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
Physical Surface("cyl_1", 1) = {1};
Physical Surface("cyl_2", 2) = {2};

Physical Curve("bottom_1", 1)  = {1};
Physical Curve("bottom_2", 2)  = {5};
Physical Curve("lateral_1", 3) = {2};
Physical Curve("outer_2", 4)   = {6};
Physical Curve("inner_2", 5)   = {8};
Physical Curve("top_2", 6)     = {7};
Physical Curve("top_1", 7)     = {3};
Physical Curve("axis_1", 8)    = {4};

Mesh.ElementOrder = 1;
