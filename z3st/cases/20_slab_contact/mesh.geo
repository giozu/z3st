// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
// Gmsh GEO for r-z section of two coaxial cylinders (adjacent blocks)
//
// Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

H = 0.1;           
r_min = 0.0;     
r_inter = 0.055;   
r_max = 0.065;     

// --- Mesh Density Parameters ---
n_radial = 15;  // Nodes along each radial segment (thickness)
n_vertical = 50; // Nodes along the height (z-axis)

// Create the two surfaces
Rectangle(1) = {r_min, 0, 0, (r_inter - r_min), H};
Rectangle(2) = {r_inter, 0, 0, (r_max - r_inter), H};

// Merge them to share the interface line
BooleanFragments{ Surface{1, 2}; Delete; }{}

// --- Structured Mesh Commands ---
// Radial lines (bottom and top segments)
Transfinite Curve {1, 5, 3, 7} = n_radial;

// Vertical lines (inner wall, contact interface, outer wall)
Transfinite Curve {4, 2, 6} = n_vertical;

// Force structured grid on surfaces
Transfinite Surface {1};
Transfinite Surface {2};

// Convert triangles to quadrilaterals (better for interface stress)
Recombine Surface {1, 2};

// --- Physical Surfaces ---
Physical Surface("cyl_inner", 1) = {1};
Physical Surface("cyl_outer", 2) = {2};

// --- Physical Curves ---
Physical Curve("bottom", 11) = {1, 5};
Physical Curve("top", 12) = {3, 7};
Physical Curve("contact", 13) = {2};
Physical Curve("outer_wall", 14) = {6};
Physical Curve("inner_wall", 15) = {4};

Mesh.ElementOrder = 1;
Mesh 2;
Save "mesh.msh";