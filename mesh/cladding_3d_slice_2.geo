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

ri = 0.5; // inner radius of the cylinder
ro = 1.5; // outer radius of the cylinder
h = 10.; // height of the cylinder

n = 20; // number of radial divisions
n_l = 20; // number of extruded layers
ms = 1.; // mesh element size

// Inner circle
Point(1) = {ri * Sqrt(2) / 2, ri * Sqrt(2) / 2, 0, ms };
Point(2) = {ri * Sqrt(2) / 2, -ri * Sqrt(2) / 2, 0, ms };
Point(3) = {-ri * Sqrt(2) / 2, -ri * Sqrt(2) / 2, 0, ms };
Point(4) = {-ri * Sqrt(2) / 2, ri * Sqrt(2) / 2, 0, ms };

// Outer circle
Point(5) = {ro * Sqrt(2) / 2, ro * Sqrt(2) / 2, 0, ms };
Point(6) = {ro * Sqrt(2) / 2, -ro * Sqrt(2) / 2, 0, ms };
Point(7) = {-ro * Sqrt(2) / 2, -ro * Sqrt(2) / 2, 0, ms };
Point(8) = {-ro * Sqrt(2) / 2, ro * Sqrt(2) / 2, 0, ms };

Point(9) = {0, 0, 0};

Circle(1) = {1, 9, 2};
Circle(2) = {2, 9, 3};
Circle(3) = {3, 9, 4};
Circle(4) = {4, 9, 1};

Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};

// Curve loops
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};

// Transfinite commands for structured mesh
Transfinite Line {1, 2, 3, 4} = n Using Progression 1.;
Transfinite Line {5, 6, 7, 8} = 2*n Using Progression 1.;

Transfinite Surface;

Recombine Surface {1};

// Physical surfaces
Physical Surface("bottom") = {1};
Physical Surface("inner") = {2, 3, 4, 5};
Physical Surface("outer") = {6, 7, 8, 9};
Physical Surface("top") = {10};

// Extrusion
Extrude {0, 0, h} {
  Surface{1}; Layers {n_l}; Recombine;
}

// Physical volume
Physical Volume("volume", 1) = {1};

Mesh 3;

Save "cladding_3d_slice_2.msh";
