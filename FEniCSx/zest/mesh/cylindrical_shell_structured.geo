// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a 3D thick-walled cylinder, radial-structured mesh
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// Inner and outer radii
ri = 0.0042;
ro = 0.0047;
n_r = 10;
n_h = 10;

h = 0.01; // Height of the 3D extrusion

// --------
// 2D slice
// --------
Point(1) = {0, 0, 0};        // Center point

// inner circle
Point(2) = {ri, 0, 0};       // Point on inner radius
Point(3) = {0, ri, 0};       // Point on inner radius (90°)
Point(4) = {-ri, 0, 0};      // Point on inner radius (180°)
Point(5) = {0, -ri, 0};      // Point on inner radius (270°)

// outer circle
Point(6) = {ro, 0, 0};       // Point on outer radius
Point(7) = {0, ro, 0};       // Point on outer radius (90°)
Point(8) = {-ro, 0, 0};      // Point on outer radius (180°)
Point(9) = {0, -ro, 0};      // Point on outer radius (270°)

// Define arcs for inner
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Define arcs for inner and outer circles
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Radial lines connecting inner and outer circles
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Curve loops and surface
Curve Loop(1) = {1, 2, 3, 4};         // Inner circle
Curve Loop(2) = {5, 6, 7, 8};         // Outer circle

Curve Loop(3) = {1, 10, -5, -9};      // First radial section
Curve Loop(4) = {2, 11, -6, -10};     // Second radial section
Curve Loop(5) = {3, 12, -7, -11};     // Third radial section
Curve Loop(6) = {4, 9, -8, -12};      // Fourth radial section

Plane Surface(1) = {3}; // Radial sector 1
Plane Surface(2) = {4}; // Radial sector 2
Plane Surface(3) = {5}; // Radial sector 3
Plane Surface(4) = {6}; // Radial sector 4

// Transfinite commands for structured radial mesh
Transfinite Line {1:12} = n_r;       // Uniform divisions for radial and angular lines
Transfinite Surface {1, 2, 3, 4};   // Structured quadrilateral mesh

Recombine Surface {1, 2, 3, 4};     // Hexahedral elements

// -------
// Extrude
// -------
Extrude {0, 0, h} {
    Surface {1, 2, 3, 4}; 
    Layers {n_h};
    Recombine;
}

// Define physical surfaces
Physical Surface("cladding_bottom") = {1, 2, 3, 4};
Physical Surface("cladding_inner_lateral") = {5, 10, 14, 18};
Physical Surface("cladding_outer_lateral") = {7, 12, 16, 19};
Physical Surface("cladding_top") = {9, 13, 17, 20};

// Define physical volume
Physical Volume("cladding") = {1, 2, 3, 4};

// Mesh
Mesh.ElementOrder = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MeshSizeFromPoints = 0.5;
Mesh.MeshSizeExtendFromBoundary = 0.5;
Mesh 3;

// Save the mesh
Save "cylindrical_shell_structured.msh";
