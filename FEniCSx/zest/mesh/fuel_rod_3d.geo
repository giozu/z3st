// -----------------------------------------------------------------------------
//
//  Gmsh GEO for a fuel rod
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE");

h = 2; // Altezza del cilindro

// ------------------------
// 1️⃣ DEFINIZIONE DEI PUNTI CENTRALI
// ------------------------
Point(100) = {0, 0, 0};       // Centro della base
Point(101) = {0, 0, h};       // Centro della cima

// ------------------------
// 2️⃣ DEFINIZIONE DEI PUNTI SUL CONTORNO
// ------------------------
r_fuel = 1.0;      // Raggio Fuel
r_clad_inner = 1.5;  // Raggio interno Cladding
r_clad_outer = 2.0;  // Raggio esterno Cladding

// Punti per il fuel
Point(1) = {r_fuel, 0, 0};
Point(2) = {0, r_fuel, 0};
Point(3) = {-r_fuel, 0, 0};
Point(4) = {0, -r_fuel, 0};

// Punti per l'interfaccia fuel-cladding
Point(5) = {r_clad_inner, 0, 0};
Point(6) = {0, r_clad_inner, 0};
Point(7) = {-r_clad_inner, 0, 0};
Point(8) = {0, -r_clad_inner, 0};

// Punti per il cladding esterno
Point(9)  = {r_clad_outer, 0, 0};
Point(10) = {0, r_clad_outer, 0};
Point(11) = {-r_clad_outer, 0, 0};
Point(12) = {0, -r_clad_outer, 0};

// ------------------------
// 3️⃣ DEFINIZIONE DELLE CURVE CHIUSE (CERCHI)
// ------------------------

// Fuel
Circle(1) = {1, 100, 2};  
Circle(2) = {2, 100, 3};
Circle(3) = {3, 100, 4};
Circle(4) = {4, 100, 1};

// Cladding interno
Circle(5) = {5, 100, 6};  
Circle(6) = {6, 100, 7};
Circle(7) = {7, 100, 8};
Circle(8) = {8, 100, 5};

// Cladding esterno
Circle(9)  = {9, 100, 10};  
Circle(10) = {10, 100, 11};
Circle(11) = {11, 100, 12};
Circle(12) = {12, 100, 9};

// Connessioni radiali
Line(13) = {1, 5};
Line(14) = {2, 6};
Line(15) = {3, 7};
Line(16) = {4, 8};

Line(17) = {5, 9};
Line(18) = {6, 10};
Line(19) = {7, 11};
Line(20) = {8, 12};

// ------------------------
// 4️⃣ DEFINIZIONE DELLE SUPERFICI
// ------------------------
Curve Loop(1) = {1, 2, 3, 4}; // Fuel
Curve Loop(2) = {5, 6, 7, 8}; // Cladding interno
Curve Loop(3) = {9, 10, 11, 12}; // Cladding esterno

Curve Loop(4) = {1, 14, -5, -13}; 
Curve Loop(5) = {2, 15, -6, -14}; 
Curve Loop(6) = {3, 16, -7, -15}; 
Curve Loop(7) = {4, 13, -8, -16}; 

Curve Loop(8) = {5, 18, -9, -17}; 
Curve Loop(9) = {6, 19, -10, -18}; 
Curve Loop(10) = {7, 20, -11, -19}; 
Curve Loop(11) = {8, 17, -12, -20}; 

Plane Surface(1) = {1}; // Fuel
Plane Surface(2) = {4}; // Primo settore Cladding
Plane Surface(3) = {5}; 
Plane Surface(4) = {6}; 
Plane Surface(5) = {7}; 

// Cladding esterno
Plane Surface(6) = {8}; 
Plane Surface(7) = {9}; 
Plane Surface(8) = {10}; 
Plane Surface(9) = {11}; 

// ------------------------
// 5️⃣ ESTRUSIONE
// ------------------------
Extrude {0, 0, h} {
    Surface {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Layers {8}; // Numero di layer in altezza
    Recombine;
}

// ------------------------
// 6️⃣ MESH STRUTTURATA
// ------------------------
Transfinite Line {1:20} = 16; 
Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9};  
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9};

// ------------------------
// 7️⃣ DEFINIZIONE FISICA (PER FEM)
// ------------------------
Physical Surface("fuel_bottom") = {1};
Physical Surface("cladding_bottom") = {2, 3, 4, 5, 6, 7, 8, 9};

Physical Surface("fuel_lateral") = {10};
Physical Surface("fuel_top") = {11};
Physical Surface("cladding_inner_lateral") = {12, 13, 14, 15};
Physical Surface("cladding_outer_lateral") = {16, 17, 18, 19};
Physical Surface("cladding_top") = {20};

Physical Volume ("fuel", 1) = {1};
Physical Volume ("cladding", 2) = {2, 3, 4, 5, 6, 7, 8, 9};

// ------------------------
// 8️⃣ GENERAZIONE MESH
// ------------------------
Mesh.ElementOrder = 1;
Mesh 3;

// **Salvataggio**
Save "fuel_cladding_3d_structured.msh";
