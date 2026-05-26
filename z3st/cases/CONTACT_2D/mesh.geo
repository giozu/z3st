
SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
Lx_ext = 0.20;  // Longueur extérieure (m)
Ly_ext = 0.30;  // Largeur extérieure (m)
Lx_int = 0.19;  // Longueur intérieure (m)
Ly_int = 0.29;  // Largeur intérieure (m)

// === Points du rectangle extérieur ===
Point(1) = {0, 0, 0};
Point(2) = {Lx_ext, 0, 0};
Point(3) = {Lx_ext, Ly_ext, 0};
Point(4) = {0, Ly_ext, 0};

// === Points du rectangle intérieur ===
Point(5) = {(Lx_ext - Lx_int)/2, (Ly_ext - Ly_int)/2, 0};
Point(6) = {(Lx_ext + Lx_int)/2, (Ly_ext - Ly_int)/2, 0};
Point(7) = {(Lx_ext + Lx_int)/2, (Ly_ext + Ly_int)/2, 0};
Point(8) = {(Lx_ext - Lx_int)/2, (Ly_ext + Ly_int)/2, 0};

// === Lignes du rectangle extérieur ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// === Lignes du rectangle intérieur ===
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};  // Boucle extérieure
Curve Loop(2) = {5, 6, 7, 8};  // Boucle intérieure

// === Surfaces ===
// Surface entre les deux rectangles (domaine extérieur - trou intérieur)
Plane Surface(1) = {1, 2};

// Surface du rectangle intérieur (à mailler aussi)
Plane Surface(2) = {2};

// === Maillage structuré ===
nx = 81;  // Divisions en X 
ny = 81;  // Divisions en Y

// Transfinite pour les lignes extérieures
Transfinite Line{1, 3} = nx;
Transfinite Line{2, 4} = ny;

// Transfinite pour les lignes intérieures
Transfinite Line{5, 7} = nx - 20;  
Transfinite Line{6, 8} = ny - 40;  

// Transfinite pour les surfaces
// La surface 1 contient un trou (boucle intérieure) et ne peut pas
// être marquée `Transfinite Surface` (elle a 8 "coins").
// On garde le transfinite uniquement pour la surface intérieure.
//Transfinite Surface{1};
Transfinite Surface{2};

// Recombine pour obtenir des quadrilatères
Recombine Surface{1};
Recombine Surface{2};

// === Groupes physiques ===
Physical Surface("steel_o", 10) = {1};  // Zone entre les rectangles
Physical Surface("steel_i", 11) = {2}; // Rectangle intérieur maillé
Physical Curve("inner_radius_o", 1) = {4};
Physical Curve("outer_radius_o", 2) = {2};
Physical Curve("bottom_o", 3) = {1};
Physical Curve("top_o", 4) = {3};
Physical Curve("inner_radius_in", 5) = {8};
Physical Curve("outer_radius_in", 6) = {6};
Physical Curve("bottom_in", 7) = {5};
Physical Curve("top_in", 8) = {7};

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";