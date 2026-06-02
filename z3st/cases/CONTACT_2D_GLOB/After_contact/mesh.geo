SetFactory("OpenCASCADE");

// === Paramètres géométriques ===
Lx_ext = 0.03;  // Longueur extérieure (m)
Ly_ext = 0.50;  // Largeur extérieure (m)
Lx_int = 0.029; // Largeur du rectangle interne (m)
Ly_int = 0.499; // Hauteur du rectangle interne (m)

// === Points du grand rectangle extérieur ===
Point(1) = {0, 0, 0};
Point(2) = {Lx_ext, 0, 0};
Point(3) = {Lx_ext, Ly_ext, 0};
Point(4) = {0, Ly_ext, 0};

// === Points du rectangle interne en contact ===
Point(5) = {(Lx_ext - Lx_int)/2, (Ly_ext - Ly_int)/2, 0};
Point(6) = {(Lx_ext + Lx_int)/2, (Ly_ext - Ly_int)/2, 0};
Point(7) = {(Lx_ext + Lx_int)/2, (Ly_ext + Ly_int)/2, 0};
Point(8) = {(Lx_ext - Lx_int)/2, (Ly_ext + Ly_int)/2, 0};

// === Lignes du grand rectangle extérieur ===
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// === Lignes du rectangle interne ===
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// === Boucles de courbes ===
Curve Loop(1) = {1, 2, 3, 4};  // Surface extérieure
Curve Loop(2) = {5, 6, 7, 8};  // Surface intérieure en contact


Plane Surface(1) = {2};
Plane Surface(2) = {1,2};

// === Maillage structuré ===
nx = 162;  // divisions en X
ny = 4820;  // divisions en Y

Transfinite Line{1, 3} = nx*10;
Transfinite Line{2, 4} = ny;
Transfinite Line{5, 7} = nx;
Transfinite Line{6, 8} = ny;

Transfinite Surface{1, 2};  // surfaces structurées en triangles


// === Groupes physiques ===
Physical Surface("steel_o", 1) = {1};
Physical Surface("steel_mid", 2) = {2};
Physical Curve("left_o", 1) = {4};
Physical Curve("right_o", 2) = {2};
Physical Curve("bottom_o", 3) = {1};
Physical Curve("top_o", 4) = {3};
Physical Curve("left_mid", 5) = {8};
Physical Curve("right_mid", 6) = {6};
Physical Curve("bottom_mid", 7) = {5};
Physical Curve("top_mid", 8) = {7};

// === Génération du maillage ===
Mesh 2;

// === Sauvegarde du maillage ===
Save "mesh.msh";
