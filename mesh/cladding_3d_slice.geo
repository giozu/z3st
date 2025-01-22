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

ri = 0.5;
ro = 1.0;

// --------
// 2D slice
// --------
Circle(1) = {0,0,0,ri};       Curve Loop(1) = 1;
Circle(2) = {0,0,0,ro};       Curve Loop(2) = 2;

Plane Surface(1) = {1, 2};

// Transfinite commands for structured mesh
Transfinite Line {1, 2} = 8;

Recombine Surface {1};

// -------
// Extrude
// -------
h = 2;
Extrude {0, 0, h} {
    Surface{1}; Layers{10}; Recombine;
}

Physical Surface("cladding_bottom") = {1};
Physical Surface("cladding_inner_lateral") = {2};
Physical Surface("cladding_outer_lateral") = {3};
Physical Surface("cladding_top") = {4};

Physical Volume ("cladding", 1) = {3, 4, 1, 5};

Mesh 3; // 3 = volumes are meshed

Save "cladding_3d_slice.msh";
