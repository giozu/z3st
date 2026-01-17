// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D section of a thin-walled cylinder, structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ri = 2.000;  // Inner radius (m)
Ro = 2.100;  // Outer radius (m)
Lz = 20.00;  // Height

// Divisions
nx = 41;     // 40 radial elements
ny = 91;     // 90 axial elements

Point(1) = {Ri, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2}; // Bottom (z=0)
Line(2) = {2, 3}; // Outer radius (r=Ro)
Line(3) = {3, 4}; // Top (z=H)
Line(4) = {4, 1}; // Inner radius (r=Ri)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured
Transfinite Line {1, 3} = nx;
Transfinite Line {2, 4} = ny;
Transfinite Surface {1};
Recombine Surface {1};

// Hex
Recombine Surface{1};

Physical Surface("steel", 10) = {1};
Physical Curve("inner_radius", 1) = {4};
Physical Curve("outer_radius", 2) = {2};
Physical Curve("bottom", 3) = {1};
Physical Curve("top", 4) = {3};

// Generate the 2D mesh
Mesh 2;

// Save the mesh
Save "mesh.msh";
