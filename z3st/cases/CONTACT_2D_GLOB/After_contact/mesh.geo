SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
R_o = 0.005;   // Longueur totale du domaine (m)
Ly_rect = 0.01;   // Hauteur du domaine (m)
R_mid = 0.004;  // Largeur de chaque rectangle (m)
R_int = 0.001;  // Décalage du centre des rectangles par rapport à l'axe central (m)

// === Points des deux rectangles adjacents ===
Point(1) = {R_int, 0, 0};
Point(2) = {R_mid, 0, 0};
Point(3) = {R_mid, Ly_rect, 0};
Point(4) = {R_int, Ly_rect, 0};
Point(5) = {R_o, 0, 0};
Point(6) = {R_o, Ly_rect, 0};

// === Lignes du rectangle gauche ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// === Lignes du rectangle droit ===
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};        // Rectangle gauche
Curve Loop(2) = {5, 6, 7, -2};       // Rectangle droit

// === Surfaces ===
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// === Maillage structuré ===
nx = 40;   // divisions en X pour chaque rectangle
ny = 50;  // divisions en Y

Transfinite Line{1, 5} = nx;
Transfinite Line{3, 7} = nx;
Transfinite Line{2, 4, 6} = ny;

Transfinite Surface{1, 2};

// === Groupes physiques ===
Physical Surface("steel_in", 6) = {1};
Physical Surface("steel_o", 7) = {2};
Physical Curve("bottom", 1) = {1, 5};
Physical Curve("top", 2) = {3, 7};
Physical Curve("inner_radius", 3) = {4};
Physical Curve("outer_radius", 4) = {6};
Physical Curve("interface", 5) = {2};

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";
