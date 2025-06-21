// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a full cylinder
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// Use the OpenCASCADE geometry kernel.
SetFactory("OpenCASCADE");

// --- Parameters ---
h = 0.01;         // Height of the cylinder
r = 0.004;        // Radius of the cylinder

n_theta = 128;     // Number of divisions along the circumference (use a multiple of 4)
n_z = 20;         // Number of divisions along the height

// Calculate divisions for each quarter-circle arc.
n_theta_quarter = n_theta / 4;

// --- 2D Base Geometry (a single circular surface) ---
Point(0) = {0, 0, 0};
Point(1) = {r, 0, 0};
Point(2) = {0, r, 0};
Point(3) = {-r, 0, 0};
Point(4) = {0, -r, 0};

// Arcs forming the circle
Circle(1) = {1, 0, 2}; // NE quadrant
Circle(2) = {2, 0, 3}; // NW quadrant
Circle(3) = {3, 0, 4}; // SW quadrant
Circle(4) = {4, 0, 1}; // SE quadrant

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- 2D Structured Mesh ---
// Apply transfinite meshing to the boundary curves.
Transfinite Curve {1, 2, 3, 4} = n_theta_quarter + 1;

// Make the surface transfinite, providing the 4 boundary corners and the center point.
Transfinite Surface {1} = {1, 2, 3, 4, 0};
Recombine Surface {1};

// --- 3D Extrusion ---
geom[] = Extrude {0, 0, h} {
    Surface{1};
    Layers{n_z};
    Recombine;
};

// --- Physical Groups ---
Physical Surface("bottom") = {1};
Physical Surface("top") = {geom[0]};
Physical Surface("lateral") = {geom[2], geom[3], geom[4], geom[5]};
Physical Volume("fuel") = {geom[1]};

// --- Mesh Generation ---
Mesh.ElementOrder = 1;
Mesh 3;

// Save the mesh file in the same directory as this .geo script.
Save "full_cylinder_structured.msh";
