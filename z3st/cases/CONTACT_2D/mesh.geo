
SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
Lx_ext = 0.02;  // Longueur extérieure (m)
Ly_ext = 0.30;  // Largeur extérieure (m)
Lx_int = 0.019;  // Longueur intérieure (m)
Ly_int = 0.299;  // Largeur intérieure (m)
Lx_mid = 0.015; // Longueur du rectangle central (m)
Ly_mid = 0.25; // Largeur du rectangle central (m)

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

// === Points du rectangle central ===
Point(9) = {(Lx_ext - Lx_mid)/2, (Ly_ext - Ly_mid)/2, 0};
Point(10) = {(Lx_ext + Lx_mid)/2, (Ly_ext - Ly_mid)/2, 0};
Point(11) = {(Lx_ext + Lx_mid)/2, (Ly_ext + Ly_mid)/2, 0};
Point(12) = {(Lx_ext - Lx_mid)/2, (Ly_ext + Ly_mid)/2, 0};

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

// === Lignes du rectangle central ===
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};  // Boucle extérieure
Curve Loop(2) = {5, 6, 7, 8};  // Boucle intérieure
Curve Loop(3) = {9, 10, 11, 12}; // Boucle du rectangle central

// === Surfaces ===
// Surface entre les deux rectangles (domaine extérieur - trou intérieur)
Plane Surface(1) = {1, 2};


// Surface du rectangle central
Plane Surface(3) = {3};

// === Maillage structuré ===
nx = 162;  // Divisions en X 
ny = 482;  // Divisions en Y

// Transfinite pour les lignes extérieures
Transfinite Line{1, 3} = nx;
Transfinite Line{2, 4} = ny;

// Transfinite pour les lignes intérieures
Transfinite Line{5, 7} = nx ;  
Transfinite Line{6, 8} = ny ;  

// Transfinite pour les lignes du rectangle central
mid_nx = 41;
mid_ny = 61;
Transfinite Line{9, 11} = mid_nx;
Transfinite Line{10, 12} = mid_ny;

// Transfinite pour la surface centrale
Transfinite Surface{3};

// Recombine pour obtenir des quadrilatères sur la surface centrale
Recombine Surface{3};

// === Groupes physiques ===

Physical Surface("steel_o", 13) = {1}; // Rectangle "extérieur" maillé
Physical Surface("steel_mid", 14) = {3}; // Rectangle central maillé
Physical Curve("left_radius_o", 1) = {4};
Physical Curve("right_radius_o", 2) = {2};
Physical Curve("bottom_o", 3) = {1};
Physical Curve("top_o", 4) = {3};
Physical Curve("left_radius_in", 5) = {8};
Physical Curve("right_radius_in", 6) = {6};
Physical Curve("bottom_in", 7) = {5};
Physical Curve("top_in", 8) = {7};
Physical Curve("left_mid", 9) = {12};
Physical Curve("bottom_mid", 10) = {9};
Physical Curve("right_mid", 11) = {10};
Physical Curve("top_mid", 12) = {11};

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";