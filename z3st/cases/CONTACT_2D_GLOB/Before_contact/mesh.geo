SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
R_o = 0.05;   // Longueur totale du domaine (m)
Ly_rect = 0.1;   // Hauteur du domaine (m)
R_mid = 0.03;  // Largeur de rectangle du milieu (m)
R_int = 0.01;  // Décalage du centre des rectangles par rapport à l'axe central (m)

// === Points des deux rectangles adjacents ===
Point(1) = {R_int, 0, 0};
Point(2) = {R_mid, 0, 0};
Point(3) = {R_mid, Ly_rect, 0};
Point(4) = {R_int, Ly_rect, 0};
Point(5) = {R_o, 0, 0};
Point(6) = {R_o, Ly_rect, 0};
Point(7) = {R_o-R_int, Ly_rect, 0};
Point(8) = {R_o-R_int, 0, 0};

// === Lignes du rectangle gauche ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// === Lignes du rectangle droit ===
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};


// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};        // Rectangle gauche
Curve Loop(2) = {5, 6, 7, 8};       // Rectangle droit

// === Surfaces ===
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// === Maillage structuré ===
nx = 30;   // divisions en X pour chaque rectangle
ny = 50;  // divisions en Y

Transfinite Line{1, 5} = nx;
Transfinite Line{3, 7} = nx;
Transfinite Line{2, 4, 6, 8} = ny;

Transfinite Surface{1, 2};

// === Groupes physiques ===
Physical Surface("steel_in", 9) = {1};
Physical Surface("steel_o", 10) = {2};
Physical Curve("bottom_in", 1) = {1};
Physical Curve("top_in", 2) = {3};
Physical Curve("left_radius_in", 3) = {4};
Physical Curve("left_radius_o", 4) = {2};
Physical Curve("bottom_o", 5) = {5};
Physical Curve("top_o", 6) = {7};
Physical Curve("right_radius_in", 7) = {8};
Physical Curve("right_radius_o", 8) = {6};


// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";