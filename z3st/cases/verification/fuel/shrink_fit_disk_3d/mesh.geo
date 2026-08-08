// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: 3D quarter-disk segment (r-theta extruded along z)
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

r_1_o = 0.0041;     // cylinder 1 outer radius (m)        4.10 mm
r_2_i = 0.00413;    // cylinder 2 inner radius (m)        4.13 mm  -> gap = 30 um
r_2_o = 0.00475;    // cylinder 2 outer radius (m)        4.75 mm
h     = 0.001;      // extrusion height (m)               1 mm

lc_f = 0.30e-3;     // cylinder 1 element size (m)
lc_c = 0.15e-3;     // cylinder 2 element size (m)

// cyl_1 : quarter disk [0, r_1_o] x [0, 90 deg]
Point(1) = {0,     0,     0, lc_f};
Point(2) = {r_1_o, 0,     0, lc_f};
Point(3) = {0,     r_1_o, 0, lc_f};
Line(1)   = {1, 2};
Circle(2) = {2, 1, 3};
Line(3)   = {3, 1};
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1}; // z0_1

// cyl_2 : cladding, quarter annulus [r_2_i, r_2_o] x [0, 90 deg]
Point(4) = {r_2_i, 0,     0, lc_c};
Point(5) = {r_2_o, 0,     0, lc_c};
Point(6) = {0,     r_2_o, 0, lc_c};
Point(7) = {0,     r_2_i, 0, lc_c};
Line(4)   = {4, 5};
Circle(5) = {5, 1, 6};
Line(6)   = {6, 7};
Circle(7) = {7, 1, 4};
Curve Loop(2) = {4, 5, 6, 7};
Plane Surface(2) = {2}; // z0_2

// extrude along z (unstructured tets)
ext1[] = Extrude {0, 0, h} { Surface{1}; };
ext2[] = Extrude {0, 0, h} { Surface{2}; };

// physical groups (ids must match geometry.yaml)
Physical Volume("cyl_1", 1) = {ext1[1]};
Physical Volume("cyl_2", 2) = {ext2[1]};

Physical Surface("ymin_1", 1)    = {ext1[2]};
Physical Surface("lateral_1", 3) = {ext1[3]};
Physical Surface("xmin_1", 7)    = {ext1[4]};
Physical Surface("zmin_1", 8)    = {1};

Physical Surface("ymin_2", 2)    = {ext2[2]};
Physical Surface("outer_2", 4)   = {ext2[3]};
Physical Surface("inner_2", 5)   = {ext2[5]};
Physical Surface("xmin_2", 6)    = {ext2[4]};
Physical Surface("zmin_2", 10)   = {2};

Mesh.ElementOrder = 1;
