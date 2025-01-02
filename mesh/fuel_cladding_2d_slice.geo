// -----------------------------------------------------------------------------
//
//  Gmsh GEO for fuel and cladding 2D slice
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

// The OpenCASCADE geometry kernel supports several useful features for solid
// modelling.

SetFactory("OpenCASCADE");

Circle(1) = {0,0,0, 1.0};       Curve Loop(1) = 1;
Circle(2) = {0,0,0, 1.5};       Curve Loop(2) = 2;
Circle(3) = {0,0,0, 2.0};       Curve Loop(3) = 3;

Plane Surface(1) = {1};
Plane Surface(2) = {2, 3};

Physical Surface("fuel_bottom") = {1};
Physical Surface("cladding_bottom") = {2};

// Mesh
// mesh_size = 0.25;
// MeshSize{ PointsOf{ Surface{1,2}; } } = mesh_size;

// Recombine Surface{1, 2}; // to create quadrangles instead of triangles

// Mesh.Algorithm = 5; // Delaunay for 2D meshes
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes

// Mesh.ElementOrder = 2; // Quadratic elements

// Mesh.HighOrderOptimize = 2;

// Mesh.MeshSizeExtendFromBoundary = 0;
// Mesh.MeshSizeFromPoints = 0;
// Mesh.MeshSizeFromCurvature = 20; // number of elements per 2 * pi rad

Mesh.MeshSizeMin = 0.001;
Mesh.MeshSizeMax = 1.0;

Mesh 2; // 2 = surfaces are meshed

RecombineMesh;
Mesh.SubdivisionAlgorithm = 1;
RefineMesh;

Save "fuel_cladding_2d_slice.msh";