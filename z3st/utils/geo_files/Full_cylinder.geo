// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a full cylinder
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Lz = 0.100;        // Height of the cylinder (m)
r  = 0.040;        // Radius of the cylinder (m)

n_theta = 50;
n_z = 31;         // Number of divisions along the height

// Calculate divisions for each quarter-circle arc.

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
Transfinite Curve {1, 2, 3, 4} = n_theta + 1;

// Make the surface transfinite, providing the 4 boundary corners and the center point.
Transfinite Surface {1} = {1, 2, 3, 4}; // uncomment for structured meshing
Recombine Surface {1};

// --- 3D Extrusion ---
geom[] = Extrude {0, 0, Lz} {
    Surface{1};
    Layers{n_z};
    Recombine;
};

// --- Physical Groups ---
Physical Surface("bottom") = {1};
Physical Surface("top") = {geom[0]};
Physical Surface("lateral") = {geom[2], geom[3], geom[4], geom[5]};
Physical Volume("oxide") = {geom[1]};

// --- Mesh Generation ---
Mesh.ElementOrder = 1;
Mesh 3;

// Save the mesh file in the same directory as this .geo script.
Save "mesh.msh";
