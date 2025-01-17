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

r = 1.; // outer radius of the cylinder
l = 0.5; // edge of the inner square
h = 2.; // height of the cylinder

n = 12; // number of elements on the square edge
n_l = 12; // number of extruded layers
ms = 1.; // mesh element size

// Inner square
Point(1) = {l/2, l/2, 0, ms};
Point(2) = {l/2, -l/2, 0, ms};
Point(3) = {-l/2, -l/2, 0, ms};
Point(4) = {-l/2, l/2, 0, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Circle
Point(5) = {r * Sqrt(2) / 2, r * Sqrt(2) / 2, 0, ms };
Point(6) = {r * Sqrt(2) / 2, -r * Sqrt(2) / 2, 0, ms };
Point(7) = {-r * Sqrt(2) / 2, -r * Sqrt(2) / 2, 0, ms };
Point(8) = {-r * Sqrt(2) / 2, r * Sqrt(2) / 2, 0, ms };

Point(9) = {0, 0, 0};

Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};

// Lines connecting inner square and circle
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Curve loops for outer sectors
Curve Loop(2) = {9, 5, -10, 1};
Curve Loop(3) = {10, 6, -11, 2};
Curve Loop(4) = {11, 7, -12, 3};
Curve Loop(5) = {12, 8, -9, 4};

Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// Transfinite commands for structured mesh - inner square
Transfinite Line {1} = n Using Progression 1.;
Transfinite Line {2} = n Using Progression 1.;
Transfinite Line {3} = n Using Progression 1.;
Transfinite Line {4} = n Using Progression 1.;

Transfinite Surface {1} = {1, 2, 3, 4};

// Transfinite commands for structured mesh - outer sectors
Transfinite Line {5, 6, 7, 8, 9, 10, 11, 12} = n;
Transfinite Surface {2, 3, 4, 5} = {9, 5, -10, 1, 10, 6, -11, 2, 11, 7, -12, 3, 12, 8, -9, 4};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};

// Physical surfaces
Physical Surface("bottom") = {1, 2, 3, 4, 5};
Physical Surface("lateral") = {12, 15, 18, 21};
Physical Surface("top") = {10, 14, 17, 20, 22};

// Extrusion
Extrude {0, 0, h} {
  Surface{1, 2, 3, 4, 5}; Layers {n_l}; Recombine;
}

// Physical volume
// Physical Volume("volume1") = {1};
// Physical Volume("volume2") = {2};
// Physical Volume("volume3") = {3};
// Physical Volume("volume4") = {4};
// Physical Volume("volume5") = {5};
Physical Volume("volume") = {1, 2, 3, 4, 5};

Mesh 3;

Save "cylindrical_bar.msh";
