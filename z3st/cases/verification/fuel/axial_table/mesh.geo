// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D axisymmetric fuel stack (axial-power verification)
//
//  Tall (r-z) column: the radial source is flat, the axial chopped-cosine
//  profile is what the case verifies, so the axial direction carries the
//  resolution.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 0.000;   // Inner radius (m)
Ro = 0.0041;  // Outer radius (m)
Lz = 0.400;   // Active fuel height L (m)

// Divisions
nx = 6;    // radial (source is radially flat)
ny = 161;  // axial (resolves the cosine; odd -> a node sits exactly at z_mid)

Point(1) = {Ri, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2}; // Bottom (z=0)
Line(2) = {2, 3}; // Outer radius (r=Ro)
Line(3) = {3, 4}; // Top (z=Lz)
Line(4) = {4, 1}; // Inner radius (r=Ri)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured
Transfinite Line {1, 3} = nx;
Transfinite Line {2, 4} = ny;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("volume", 10) = {1};
Physical Curve("inner", 1) = {4};
Physical Curve("outer", 2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top", 4) = {3};
