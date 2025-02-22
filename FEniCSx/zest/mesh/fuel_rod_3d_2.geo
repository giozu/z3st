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
r_clad = 2.0;      // Raggio esterno Cladding

// Punti per il fuel
Point(1) = {r_fuel, 0, 0};
Point(2) = {0, r_fuel, 0};
Point(3) = {-r_fuel, 0, 0};
Point(4) = {0, -r_fuel, 0};

// Punti per il cladding esterno
Point(5)  = {r_clad, 0, 0};
Point(6)  = {0, r_clad, 0};
Point(7)  = {-r_clad, 0, 0};
Point(8)  = {0, -r_clad, 0};

// ------------------------
// 3️⃣ DEFINIZIONE DELLE CURVE CHIUSE (CERCHI)
// ------------------------

// Fuel
Circle(1) = {1, 100, 2};  
Circle(2) = {2, 100, 3};
Circle(3) = {3, 100, 4};
Circle(4) = {4, 100, 1};

// Cladding esterno
Circle(5)  = {5, 100, 6};  
Circle(6)  = {6, 100, 7};
Circle(7)  = {7, 100, 8};
Circle(8)  = {8, 100, 5};

// Connessioni radiali
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// ------------------------
// 4️⃣ DEFINIZIONE DELLE SUPERFICI
// ------------------------
Curve Loop(1) = {1, 2, 3, 4}; // Fuel
Curve Loop(2) = {5, 6, 7, 8}; // Cladding esterno

Curve Loop(3) = {1, 10, -5, -9}; 
Curve Loop(4) = {2, 11, -6, -10}; 
Curve Loop(5) = {3, 12, -7, -11}; 
Curve Loop(6) = {4, 9, -8, -12}; 

Plane Surface(1) = {1}; // Fuel
Plane Surface(2) = {3}; // Primo settore Cladding
Plane Surface(3) = {4}; 
Plane Surface(4) = {5}; 
Plane Surface(5) = {6}; 

// ------------------------
// 5️⃣ ESTRUSIONE
// ------------------------
Extrude {0, 0, h} {
    Surface {1, 2, 3, 4, 5};
    Layers {8}; // Numero di layer in altezza
    Recombine;
}

// ------------------------
// 6️⃣ MESH STRUTTURATA
// ------------------------
Transfinite Line {1:12} = 16; 
Transfinite Surface {1, 2, 3, 4, 5};  
Recombine Surface {1, 2, 3, 4, 5};

// ------------------------
// 7️⃣ DEFINIZIONE FISICA (PER FEM)
// ------------------------
Physical Surface("fuel_bottom") = {1};
Physical Surface("cladding_bottom") = {2, 3, 4, 5};

Physical Surface("fuel_lateral") = {6};
Physical Surface("fuel_top") = {7};
Physical Surface("cladding_outer_lateral") = {8, 9, 10, 11};
Physical Surface("cladding_top") = {12};

Physical Volume ("fuel", 1) = {1};
Physical Volume ("cladding", 2) = {2, 3, 4, 5};

// ------------------------
// 8️⃣ GENERAZIONE MESH
// ------------------------
Mesh.ElementOrder = 1;
Mesh 3;

// **Salvataggio**
Save "fuel_cladding_3d_structured.msh";
