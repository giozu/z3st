SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
R_o = 0.05;   // Longueur totale du domaine (m)
Ly_rect = 0.1;   // Hauteur du domaine (m)
R_mid = 0.04;  // Largeur de chaque rectangle (m)
R_int = 0.01;  // Décalage du centre des rectangles par rapport à l'axe central (m)
gap = 0.0001;   // Largeur du gap (m) - comme dans before-contact

// === Points des trois rectangles adjacents ===
Point(1) = {R_int, 0, 0};
Point(2) = {R_mid - gap/2, 0, 0};
Point(3) = {R_mid - gap/2, Ly_rect, 0};
Point(4) = {R_int, Ly_rect, 0};
Point(5) = {R_mid + gap/2, 0, 0};
Point(6) = {R_mid + gap/2, Ly_rect, 0};
Point(7) = {R_o, 0, 0};
Point(8) = {R_o, Ly_rect, 0};

// === Lignes du rectangle gauche ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// === Lignes du rectangle du milieu (gap) ===
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

// === Lignes du rectangle droit ===
Line(8) = {5, 7};
Line(9) = {7, 8};
Line(10) = {8, 6};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};        // Rectangle gauche
Curve Loop(2) = {5, 6, 7, -2};       // Rectangle milieu (gap) rempli
Curve Loop(3) = {8, 9, 10, -6};      // Rectangle droit

// === Surfaces ===
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// === Maillage structuré ===
nx = 40;   // divisions en X pour chaque rectangle
ny = 50;  // divisions en Y
nx_gap = 10; // divisions en X pour la zone du gap

Transfinite Line{1, 3, 8, 10} = nx;
Transfinite Line{5, 7} = nx_gap;
Transfinite Line{2, 4, 6, 9} = ny;

Transfinite Surface{1, 2, 3};

// === Groupes physiques ===
Physical Surface("steel_in", 11) = {1};
Physical Surface("gap_region", 12) = {2};
Physical Surface("steel_o", 13) = {3};

Physical Curve("bottom_in", 1) = {1};
Physical Curve("top_in", 2) = {3};
Physical Curve("bottom_o", 6) = {8};
Physical Curve("top_o", 7) = {10};
Physical Curve("bottom_gap", 8) = {5};
Physical Curve("top_gap", 9) = {7};
Physical Curve("inner_radius", 3) = {4};
Physical Curve("outer_radius", 4) = {9};
Physical Curve("interface_in", 5) = {2};
Physical Curve("interface_out", 10) = {6};

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";
