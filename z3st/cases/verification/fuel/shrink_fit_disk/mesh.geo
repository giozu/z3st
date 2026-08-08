// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 2D in-plane (r-theta) quarter disk, plane strain
//
//  Quarter symmetry: theta in [0, 90 deg], symmetry planes on y = 0
//  (bottom_*) and x = 0 (axis_*).
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

r_1_o = 0.0041;     // pellet outer radius (m)        4.10 mm
r_2_i = 0.00413;    // clad inner radius (m)          4.13 mm  -> gap = 30 um
r_2_o = 0.00475;    // clad outer radius (m)          4.75 mm

lc_f = 0.20e-3;     // pellet element size (m)
lc_c = 0.10e-3;     // clad element size (m)

// --- cyl_1 : fuel pellet, quarter disk [0, r_1_o] x [0, 90 deg] ---
Point(1) = {0,     0,     0, lc_f};
Point(2) = {r_1_o, 0,     0, lc_f};
Point(3) = {0,     r_1_o, 0, lc_f};
Line(1)   = {1, 2};       // bottom_1  (y = 0, symmetry)
Circle(2) = {2, 1, 3};    // lateral_1 (r = r_1_o, gap-facing arc)
Line(3)   = {3, 1};       // axis_1    (x = 0, symmetry)
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

// --- cyl_2 : cladding, quarter annulus [r_2_i, r_2_o] x [0, 90 deg] ---
Point(4) = {r_2_i, 0,     0, lc_c};
Point(5) = {r_2_o, 0,     0, lc_c};
Point(6) = {0,     r_2_o, 0, lc_c};
Point(7) = {0,     r_2_i, 0, lc_c};
Line(4)   = {4, 5};       // bottom_2 (y = 0, symmetry)
Circle(5) = {5, 1, 6};    // outer_2  (r = r_2_o, coolant arc)
Line(6)   = {6, 7};       // axis_2   (x = 0, symmetry)
Circle(7) = {7, 1, 4};    // inner_2  (r = r_2_i, gap-facing arc)
Curve Loop(2) = {4, 5, 6, 7};
Plane Surface(2) = {2};

// --- physical groups (ids must match geometry.yaml) ---
Physical Surface("cyl_1", 1) = {1};
Physical Surface("cyl_2", 2) = {2};

Physical Curve("bottom_1", 1)  = {1};
Physical Curve("bottom_2", 2)  = {4};
Physical Curve("lateral_1", 3) = {2};
Physical Curve("outer_2", 4)   = {5};
Physical Curve("inner_2", 5)   = {7};
Physical Curve("axis_2", 6)    = {6};
Physical Curve("axis_1", 7)    = {3};

Mesh.ElementOrder = 1;
