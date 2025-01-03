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
Circle(1) = {0,0,0, 1.5};       Curve Loop(1) = 1;
Circle(2) = {0,0,0, 2.0};       Curve Loop(2) = 2;

Plane Surface(1) = {1, 2};

Physical Surface("cladding_bottom") = {1};

// -------
// Extrude
// -------

h = 2;

Extrude {0, 0, h} {
    Surface{1};
}

Physical Surface("cladding_inner_lateral") = {2};
Physical Surface("cladding_outer_lateral") = {3};
Physical Surface("cladding_top") = {4};

Physical Volume ("cladding", 1) = {3, 4, 1, 5};

// Mesh
MeshSize { PointsOf{ Volume{1}; } } = 0.4;
// Mesh.ElementOrder = 1;
Mesh 3; // 3 = volumes are meshed

Save "cladding_3d_slice.msh";