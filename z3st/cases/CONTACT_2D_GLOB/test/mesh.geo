SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
R_int = 0.01;   // x-coordinate of the left boundary of steel_in
R_mid = 0.04;   // x-coordinate of the interface (with a small gap)
R_o = 0.05;     // x-coordinate of the right boundary of steel_o
Ly_rect = 0.1;  // Height of the rectangles
gap = 0.001;    // Small gap between the two rectangles to avoid shared nodes

// === Points des deux rectangles ===
// Rectangle steel_in (gauche)
Point(1) = {R_int, 0, 0};
Point(2) = {R_mid - gap/2, 0, 0};  // Décalé pour éviter le contact initial
Point(3) = {R_mid - gap/2, Ly_rect, 0};
Point(4) = {R_int, Ly_rect, 0};

// Rectangle steel_o (droit)
Point(5) = {R_mid + gap/2, 0, 0};  // Décalé pour éviter le contact initial
Point(6) = {R_o, 0, 0};
Point(7) = {R_o, Ly_rect, 0};
Point(8) = {R_mid + gap/2, Ly_rect, 0};

// === Lignes des rectangles ===
// Rectangle steel_in
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Rectangle steel_o
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};  // steel_in
Curve Loop(2) = {5, 6, 7, 8};  // steel_o

// === Surfaces ===
Plane Surface(1) = {1};  // steel_in
Plane Surface(2) = {2};  // steel_o

// === Maillage structuré ===
nx = 40;
ny = 50;

Transfinite Line{1, 3, 5, 7} = nx;
Transfinite Line{2, 4, 6, 8} = ny;

Transfinite Surface{1};
Transfinite Surface{2};

// === Groupes physiques ===
// steel_in
Physical Surface("steel_in", 1) = {1};
Physical Curve("steel_in_bottom", 1) = {1};
Physical Curve("steel_in_top", 2) = {3};
Physical Curve("steel_in_left", 3) = {4};
Physical Curve("steel_in_right", 4) = {2};  // Surface de contact côté steel_in

// steel_o
Physical Surface("steel_o", 2) = {2};
Physical Curve("steel_o_bottom", 5) = {5};
Physical Curve("steel_o_top", 6) = {7};
Physical Curve("steel_o_right", 7) = {6};
Physical Curve("steel_o_left", 8) = {8};   // Surface de contact côté steel_o

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde des maillages séparés ===
// Sauvegarder steel_in
Save "mesh.msh";