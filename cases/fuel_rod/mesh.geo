// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a fuel rod
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

h = 0.01; // height
n_h = 20;
n_r = 10;
n_f = 30;

Point(100) = {0, 0, 0};       // Center of the bottom face
Point(101) = {0, 0, h};       // Center of the top face

r_fuel_o = 0.004;    // Outer radius of the fuel
r_clad_i = 0.0042;    // Inner radius of the cladding
r_clad_o = 0.0047;    // Outer radius of the cladding

// Fuel boundary points
Point(1) = {r_fuel_o, 0, 0};
Point(2) = {0, r_fuel_o, 0};
Point(3) = {-r_fuel_o, 0, 0};
Point(4) = {0, -r_fuel_o, 0};

// Cladding inner boundary points
Point(5)  = {r_clad_i, 0, 0};
Point(6)  = {0, r_clad_i, 0};
Point(7)  = {-r_clad_i, 0, 0};
Point(8)  = {0, -r_clad_i, 0};

// Cladding outer boundary points
Point(9)  = {r_clad_o, 0, 0};
Point(10) = {0, r_clad_o, 0};
Point(11) = {-r_clad_o, 0, 0};
Point(12) = {0, -r_clad_o, 0};

// Fuel outer boundary
Circle(1) = {1, 100, 2};  
Circle(2) = {2, 100, 3};
Circle(3) = {3, 100, 4};
Circle(4) = {4, 100, 1};

// Cladding inner boundary
Circle(5)  = {5, 100, 6};  
Circle(6)  = {6, 100, 7};
Circle(7)  = {7, 100, 8};
Circle(8)  = {8, 100, 5};

// Cladding outer boundary
Circle(9)   = {9, 100, 10};  
Circle(10)  = {10, 100, 11};
Circle(11)  = {11, 100, 12};
Circle(12)  = {12, 100, 9};

// Radial connections between inner and outer cladding
Line(13) = {5, 9};
Line(14) = {6, 10};
Line(15) = {7, 11};
Line(16) = {8, 12};

Curve Loop(1) = {1, 2, 3, 4};

Curve Loop(2) = {13, 9, -14, -5}; 
Curve Loop(3) = {14, 10, -15, -6}; 
Curve Loop(4) = {15, 11, -16, -7}; 
Curve Loop(5) = {16, 12, -13, -8}; 

Plane Surface(1) = {1};

Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Extrude {0, 0, h} {
    Surface {1, 2, 3, 4, 5};
    Layers {n_h}; // Number of layers in the height direction
    Recombine;
}

Transfinite Line {1:12} = n_f;
Transfinite Surface {1, 2, 3, 4, 5};  
Recombine Surface {1, 2, 3, 4, 5};

Transfinite Line {13:16} = n_r; 

Physical Surface("fuel_bottom") = {1};
Physical Surface("cladding_bottom") = {2, 3, 4, 5};

Physical Surface("fuel_lateral") = {6, 7, 8, 9};
Physical Surface("cladding_outer") = {12, 16, 20, 24};
Physical Surface("cladding_inner") = {14, 18, 22, 25};
Physical Surface("cladding_top") = {15, 19, 23, 26};
Physical Surface("fuel_top") = {10};

Physical Volume ("fuel", 1) = {1};
Physical Volume ("cladding", 2) = {2, 3, 4, 5};

Mesh.ElementOrder = 1;
Mesh 3;

Save "mesh.msh";
