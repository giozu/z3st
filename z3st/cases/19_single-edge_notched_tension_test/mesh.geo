// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a single-edge notched test
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 0.001;
Ly = 0.001;
Dn = 0.0005;

lc_damage = 0.000004;
h_fine = lc_damage / 4.0; //instead of /7.0 -> h/l=0.25, still at the valid limit
h_coarse = Lx / 75;

// Corner points
Point(1) = {0,  0,  0, h_coarse};
Point(2) = {Lx, 0,  0, h_coarse};
Point(3) = {Lx, Ly, 0, h_coarse};
Point(4) = {0,  Ly, 0, h_coarse};
Point(5) = {Dn, Ly/2, 0, h_fine};
Point(6) = {0,  Ly/2, 0, h_fine};

// Nouveau point : bord droit à mi-hauteur, extrémité du trajet de fissure attendu
Point(7) = {Lx - 20e-6, Ly/2, 0, h_fine}; //20 um before the right edge, not over it

// Edges
Line(1) = {1, 2}; // y-min
Line(2) = {2, 3}; // x-max
Line(3) = {3, 4}; // y-max
Line(4) = {4, 6}; // x-min
Line(5) = {6, 1}; // x-min
Line(6) = {6, 5}; // entaille (crack)

// Nouvelle ligne : trajet de fissure attendu, de la pointe (5) au bord droit (7)
// Purement géométrique -- sert de guide au maillage, ne referme pas de contour.
Line(7) = {5, 7};

// Surface
Line Loop(1) = {1, 2, 3, 4, 6, -6, 5};
Plane Surface(1) = {1};

// Incorpore la ligne 7 dans la surface sans la découper en deux surfaces
Line{7} In Surface{1};

// Physical groups
Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3};
Physical Curve("xmin") = {4,5};
Physical Curve("crack") = {6};
Physical Surface("uo2") = {1};

// --- Raffinement gradué le long de l'entaille + du trajet de fissure attendu ---
Field[1] = Distance;
Field[1].CurvesList = {6, 7};
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_fine;
Field[2].SizeMax = h_coarse;
Field[2].DistMin = 0;
Field[2].DistMax = 30e-6;

Background Field = 2;


// Mesh 2;
// Save "mesh.msh";
