// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 3D quarter-disk segment (r-theta extruded along z)
//
//  The 2D r-theta quarter disk of coaxial_contact_disk extruded by a thin
//  height h. Quarter symmetry: theta in [0, 90 deg], symmetry planes on
//  y = 0 (bottom_*) and x = 0 (axis_*); both z faces (z0_*, z1_*) are
//  clamped axially so the state is exactly plane strain.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

r_1_o = 0.0041;     // pellet outer radius (m)        4.10 mm
r_2_i = 0.00413;    // clad inner radius (m)          4.13 mm  -> gap = 30 um
r_2_o = 0.00475;    // clad outer radius (m)          4.75 mm
h     = 0.001;      // extrusion height (m)           1 mm

lc_f = 0.30e-3;     // pellet element size (m)
lc_c = 0.15e-3;     // clad element size (m)

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

// --- extrude along z (unstructured tets) ---
// ext[0] = top surface (z = h), ext[1] = volume, ext[2..] = lateral faces
// in curve-loop order.
ext1[] = Extrude {0, 0, h} { Surface{1}; };
ext2[] = Extrude {0, 0, h} { Surface{2}; };

// --- physical groups (ids must match geometry.yaml) ---
Physical Volume("cyl_1", 1) = {ext1[1]};
Physical Volume("cyl_2", 2) = {ext2[1]};

Physical Surface("bottom_1", 1)  = {ext1[2]};   // y = 0 (from Line 1)
Physical Surface("bottom_2", 2)  = {ext2[2]};   // y = 0 (from Line 4)
Physical Surface("lateral_1", 3) = {ext1[3]};   // arc   (from Circle 2)
Physical Surface("outer_2", 4)   = {ext2[3]};   // arc   (from Circle 5)
Physical Surface("inner_2", 5)   = {ext2[5]};   // arc   (from Circle 7)
Physical Surface("axis_2", 6)    = {ext2[4]};   // x = 0 (from Line 6)
Physical Surface("axis_1", 7)    = {ext1[4]};   // x = 0 (from Line 3)
Physical Surface("z0_1", 8)      = {1};         // pellet z = 0 face
Physical Surface("z1_1", 9)      = {ext1[0]};   // pellet z = h face
Physical Surface("z0_2", 10)     = {2};         // clad z = 0 face
Physical Surface("z1_2", 11)     = {ext2[0]};   // clad z = h face

Mesh.ElementOrder = 1;
