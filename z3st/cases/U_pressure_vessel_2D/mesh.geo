// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D section of a pressure vessel with hemi-spherical heads
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// lengths in m
scale = 1e-3;

// Parameters
Ri = 1520 * scale;              // Inner radius
t_head = 26 * scale;            // Head thickness
t_cyl = 52 * scale;             // Cylindrical part thickness
L_cyl = 2000 * scale;           // Length of the cylindrical part
L_taper = 150 * scale;          // Length of the transition zone (tapering)

// Mesh Size
lc = 5.0 * scale;              // general fineness
lc_fine = 1.0 * scale;          // refinement in the transition zone

// Center of curvature of the head
Point(1) = {0, 0, 0, lc};

// Inner profile
Point(2) = {0, Ri, 0, lc};                  // Top inner (on symmetry axis)
Point(3) = {Ri, 0, 0, lc_fine};             // Inner junction head-cylinder
Point(4) = {Ri, -L_cyl, 0, lc};             // Inner bottom

// Outer profile
Point(5) = {0, Ri + t_head, 0, lc};         // Top outer
Point(6) = {Ri + t_head, 0, 0, lc_fine};    // Outer junction head-cylinder (thickness 26)
Point(7) = {Ri + t_cyl, -L_taper, 0, lc_fine}; // End of tapering (thickness 52)
Point(8) = {Ri + t_cyl, -L_cyl, 0, lc};     // Outer bottom

// Lines
Circle(1) = {2, 1, 3};        // Inner head arc
Line(2)   = {3, 4};           // Inner cylinder wall
Circle(3) = {5, 1, 6};        // Outer head arc
Line(4)   = {6, 7};           // Taper (outer)
Line(5)   = {7, 8};           // Outer cylinder wall
Line(6)   = {4, 8};           // Bottom closure
Line(7)   = {2, 5};           // Top closure (symmetry axis)

// Surface
Curve Loop(1) = {7, 3, 4, 5, -6, -2, -1};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("steel", 1) = {1};
Physical Line("inner", 10) = {1, 2};
Physical Line("symmetry_axis", 11) = {7};
Physical Line("bottom_boundary", 12) = {6};

