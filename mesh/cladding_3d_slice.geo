// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a 3D cylindrical shell
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// Inner and outer radii
ri = 0.5;
ro = 1.0;
h = 2;

n_r = 64; // radial divisions
n_h = 40; // axial divisions

// --------
// 2D slice
// --------
Circle(1) = {0, 0, 0, ri};
Circle(2) = {0, 0, 0, ro};

// Curve loops for the annular surface
Curve Loop(1) = {1};
Curve Loop(2) = {2};

// Annular plane surface
Plane Surface(1) = {2, 1};

// Transfinite commands for structured mesh
Transfinite Line {1, 2} = n_r;
Recombine Surface {1};

// -------
// Extrude
// -------
Extrude {0, 0, h} {
    Surface{1}; Layers{n_h}; Recombine;
}

// Define physical surfaces
Physical Surface("cladding_bottom") = {1};
Physical Surface("cladding_inner_lateral") = {2};
Physical Surface("cladding_outer_lateral") = {3};
Physical Surface("cladding_top") = {4};

// Define physical volume
Physical Volume("cladding") = {1};

// Mesh generation
Mesh 3; // Generate a 3D mesh

// Save the mesh
Save "cladding_3d_slice.msh";
