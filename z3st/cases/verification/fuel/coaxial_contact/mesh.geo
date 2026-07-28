// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D axisymmetric (r-z) segment
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.0041;     // outer radius (m)
r_2_i = 0.00411;    // inner radius (m)
r_2_o = 0.00475;    // outer radius (m)
h     = 0.010;      // axial segment height (m)

n_r1 = 20;           // radial divisions, region 1
n_r2 = 20;           // radial divisions, region 2 
n_z  = 11;          // axial divisions

// --- cyl_1 : rectangle [0, r_1_o] x [0, h] ---
Point(1) = {0,     0, 0};
Point(2) = {r_1_o, 0, 0};
Point(3) = {r_1_o, h, 0};
Point(4) = {0,     h, 0};

Line(1) = {1, 2};   // ymin_1  (z = 0)
Line(2) = {2, 3};   // xmax_1  (r = r_1_o, gap-facing)
Line(3) = {3, 4};   // ymax_1  (z = h)
Line(4) = {4, 1};   // xmin_1  (r = 0, symmetry axis)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- cyl_2 : rectangle [r_2_i, r_2_o] x [0, h] ---
Point(5) = {r_2_i, 0, 0};
Point(6) = {r_2_o, 0, 0};
Point(7) = {r_2_o, h, 0};
Point(8) = {r_2_i, h, 0};

Line(5) = {5, 6};   // ymin_2  (z = 0)
Line(6) = {6, 7};   // xmax_2  (r = r_2_o)
Line(7) = {7, 8};   // ymax_2  (z = h)
Line(8) = {8, 5};   // xmin_2  (r = r_2_i, gap-facing)

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
Physical Surface("cyl_1", 10) = {1};
Physical Surface("cyl_2", 20) = {2};

// Physical Curve("label", chosen_number) = {gmsh line number};
Physical Curve("ymin_1", 1)  = {1};
Physical Curve("xmin_1", 5)  = {4};
Physical Curve("xmax_1", 7)  = {2};
Physical Curve("ymax_1", 6)  = {3};
Physical Curve("ymin_2", 2)  = {5};
Physical Curve("ymax_2", 3)  = {7};
Physical Curve("xmax_2", 4)  = {6};
Physical Curve("xmin_2", 8)  = {8};

Mesh.ElementOrder = 1;
