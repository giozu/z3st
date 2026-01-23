// -----------------------------------------------------------------------------
//  4-Block Structured Mesh for Notched Plate (Z3ST Damage Model)
//  Garantisce 4 angoli per superficie e matching dei nodi
// -----------------------------------------------------------------------------
SetFactory("Built-in");

L = 0.1; W = 0.1; a = 0.005; b = 0.002;
y1 = (W - b)/2; y2 = y1 + b;

// Punti (Griglia 3x4)
Point(1) = {0, 0, 0};   Point(2) = {a, 0, 0};   Point(3) = {L, 0, 0};
Point(4) = {0, y1, 0};  Point(5) = {a, y1, 0};  Point(6) = {L, y1, 0};
Point(7) = {0, y2, 0};  Point(8) = {a, y2, 0};  Point(9) = {L, y2, 0};
Point(10)= {0, W, 0};   Point(11)= {a, W, 0};   Point(12)= {L, W, 0};

// Linee Orizzontali
Line(1) = {1, 2}; Line(2) = {2, 3};
Line(3) = {4, 5}; Line(4) = {5, 6};
Line(5) = {7, 8}; Line(6) = {8, 9};
Line(7) = {10, 11}; Line(8) = {11, 12};

// Linee Verticali
Line(9)  = {1, 4};  Line(10) = {2, 5};  Line(11) = {3, 6};
Line(12) = {5, 8};  Line(13) = {6, 9};
Line(14) = {7, 10}; Line(15) = {8, 11}; Line(16) = {9, 12};

// Superfici (I 4 blocchi attorno all'intaglio vuoto)
Curve Loop(1) = {1, 10, -3, -9};   Plane Surface(1) = {1}; // Basso-SX
Curve Loop(2) = {2, 11, -4, -10};  Plane Surface(2) = {2}; // Basso-DX
Curve Loop(3) = {4, 13, -6, -12};  Plane Surface(3) = {3}; // Centro-DX (fronte intaglio)
Curve Loop(4) = {5, 15, -7, -14};  Plane Surface(4) = {4}; // Alto-SX
Curve Loop(5) = {6, 16, -8, -15};  Plane Surface(5) = {5}; // Alto-DX

// --- Transfinite Setup ---
NxA = 15;  // Lunghezza intaglio
NxB = 40;  // Resto del provino
Ny  = 25;  // Zone esterne
NyN = 10;  // Zona intaglio

// Assegnazione nodi (matching obbligatorio)
Transfinite Curve {1, 3, 7, 5} = NxA;
Transfinite Curve {2, 4, 6, 8} = NxB;
Transfinite Curve {9, 10, 11, 14, 15, 16} = Ny;
Transfinite Curve {12, 13} = NyN;

Transfinite Surface "*";
Recombine Surface "*";

// --- Physical Groups per Z3ST ---
Physical Surface("steel") = {1, 2, 3, 4, 5};
Physical Curve("ymin") = {1, 2};
Physical Curve("ymax") = {7, 8};
Physical Curve("xmax") = {11, 13, 16};
Physical Curve("xmin") = {9, 14};
Physical Curve("notch") = {3, 12, 5}; // Il profilo che "morde" il materiale

Mesh 2;
Save "mesh.msh";