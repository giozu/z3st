// -----------------------------------------------------------------------------
//
//  Gmsh GEO for fuel and cladding 3D slice
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// The OpenCASCADE geometry kernel supports several useful features for solid
// modelling.

SetFactory("OpenCASCADE");

// --------
// 2D slice
// --------
Circle(1) = {0,0,0, 1.0};       Curve Loop(1) = 1;
Circle(2) = {0,0,0, 1.5};       Curve Loop(2) = 2;
Circle(3) = {0,0,0, 2.0};       Curve Loop(3) = 3;

Plane Surface(1) = {1};
Plane Surface(2) = {2, 3};

Physical Surface("fuel_bottom") = {1};
Physical Surface("cladding_bottom") = {2};

// -------
// Extrude
// -------
// Extrusion along the z axis. The 2D mesh is also extruded. This is done 
// with the same `Extrude' command, but by specifying element 'Layers' 
// (10 layers, for a total height of {1} * h):

h = 2;

// Extrude {0, 0, h} {
//    Surface{1,2}; Layers{ {2, 8}, {0.5, 1} };
// }

Extrude {0, 0, h} {
    Surface{1,2}; Layers{ {8}, {1} }; Recombine; // Recombine to create quadrangles instead of triangles
}

Physical Surface("fuel_lateral") = {3};
Physical Surface("fuel_top") = {4};
Physical Surface("cladding_inner_lateral") = {5};
Physical Surface("cladding_outer_lateral") = {6};
Physical Surface("cladding_top") = {7};

Physical Volume ("fuel", 1) = {1};
Physical Volume ("cladding", 2) = {5, 6, 2, 7};

Recombine Surface{1, 2}; // Recombine to create quadrangles instead of triangles

// Mesh
MeshSize { PointsOf{ Volume{1, 2}; } } = 0.2;
Mesh.ElementOrder = 1;
Mesh 3; // 3 = volumes are meshed

Save "fuel_cladding_3d_slice.msh";