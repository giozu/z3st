// -----------------------------------------------------------------------------
//
//  Gmsh GEO for fuel and cladding 3D slice (Improved Version)
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// Use the OpenCASCADE geometry kernel.
SetFactory("OpenCASCADE");

// --- Parameters ---
h = 0.01;         // Height of the cylinder
r = 0.004;        // Radius of the fuel

n_r = 20;         // Number of divisions in the radial direction
n_theta = 32;     // Number of divisions along the circumference (use a multiple of 4)
n_z = 20;         // Number of divisions along the height

// For N elements, we need N+1 points in a transfinite definition.
n_theta_quarter = n_theta / 4;

// --- 2D Base Geometry (built from 4 quadrants for full control) ---
// Points
p_center = newp; Point(p_center) = {0, 0, 0, 0};
p_E = newp; Point(p_E) = {r, 0, 0};
p_N = newp; Point(p_N) = {0, r, 0};
p_W = newp; Point(p_W) = {-r, 0, 0};
p_S = newp; Point(p_S) = {0, -r, 0};

// Radial Lines
l_E = newl; Line(l_E) = {p_center, p_E};
l_N = newl; Line(l_N) = {p_center, p_N};
l_W = newl; Line(l_W) = {p_center, p_W};
l_S = newl; Line(l_S) = {p_center, p_S};

// Circumferential Arcs
c_NE = newl; Circle(c_NE) = {p_E, p_center, p_N};
c_NW = newl; Circle(c_NW) = {p_N, p_center, p_W};
c_SW = newl; Circle(c_SW) = {p_W, p_center, p_S};
c_SE = newl; Circle(c_SE) = {p_S, p_center, p_E};

// Create the four quadrant surfaces
s_list[0] = news; Curve Loop(s_list[0]) = {l_N, c_NW, -l_W}; Plane Surface(s_list[0]) = {s_list[0]};
s_list[1] = news; Curve Loop(s_list[1]) = {l_W, c_SW, -l_S}; Plane Surface(s_list[1]) = {s_list[1]};
s_list[2] = news; Curve Loop(s_list[2]) = {l_S, c_SE, -l_E}; Plane Surface(s_list[2]) = {s_list[2]};
s_list[3] = news; Curve Loop(s_list[3]) = {l_E, c_NE, -l_N}; Plane Surface(s_list[3]) = {s_list[3]};

// --- 2D Structured Mesh ---
// Define the number of points for radial and circumferential lines
Transfinite Curve {l_E, l_N, l_W, l_S} = n_r + 1;
Transfinite Curve {c_NE, c_NW, c_SW, c_SE} = n_theta_quarter + 1;

// Apply transfinite algorithm to all four surfaces
For i In {0:3}
  Transfinite Surface {s_list[i]};
  Recombine Surface {s_list[i]};
EndFor

// --- 3D Extrusion ---
// Extrude the four surfaces together.
// The output 'geom' array is more complex here, but we can access it predictably.
// For each of the 4 surfaces, it returns [top_surface, volume, 3_lateral_surfaces].
// So each block is 5 elements long.
geom[] = Extrude {0, 0, h} {
    Surface{s_list[]};
    Layers{n_z};
    Recombine;
};

// --- Physical Groups ---
// We need to collect the tags for each part from the extrusion output.
// The original four surfaces form the bottom face.
Physical Surface("bottom") = {s_list[]};

// The top surfaces are the 1st element of each 5-element block.
Physical Surface("top") = {geom[0], geom[5], geom[10], geom[15]};

// The lateral surfaces are created from extruding the outer arcs (c_NW, c_SW, etc).
// In our loops, these are the 2nd curve, so they become the 3rd lateral surface
// in the output list (index 3 in a 0-indexed block of 5).
Physical Surface("lateral") = {geom[3], geom[8], geom[13], geom[18]};

// The volumes are the 2nd element of each 5-element block.
Physical Volume("fuel") = {geom[1], geom[6], geom[11], geom[16]};

// --- Mesh Generation ---
Mesh.ElementOrder = 1;
Mesh 3;
Mesh.SaveAll = 0;

// Save the mesh file (uncomment to use)
Save "full_cylinder.msh";
