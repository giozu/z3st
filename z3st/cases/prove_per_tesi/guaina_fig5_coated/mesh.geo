SetFactory("OpenCASCADE");

Ri    = 4.18e-3;    // Raggio interno Zircaloy (m)
Ro    = 4.75e-3;    // Raggio esterno Zircaloy / interfaccia Cr (m)
Rcr   = 4.765e-3;   // Raggio esterno Cromo (coating = 15 um) (m)
Rplug = 2.50e-3;    // Raggio interno tappi terminali (m)
Lplug = 5.0e-3;     // Altezza tappi terminali (m)
Lz    = 77.0e-3;    // Altezza totale provino (m)

// 1. Corpo Zircaloy (tappo inferiore + parete centrale + tappo superiore)
Rectangle(1) = {Rplug, 0.0,   0.0, Ro - Rplug, Lplug};
Rectangle(2) = {Ri,    Lplug, 0.0, Ro - Ri,    Lz - 2.0 * Lplug};
Rectangle(3) = {Rplug, Lz - Lplug, 0.0, Ro - Rplug, Lplug};
BooleanUnion{ Surface{1}; Delete; }{ Surface{2, 3}; Delete; }

// 2. Liner esterno in Cromo lungo tutta la barretta (compresi i tappi)
Rectangle(4) = {Ro, 0.0, 0.0, Rcr - Ro, Lz};

// Mesh conforme tra Zircaloy e Cromo
BooleanFragments{ Surface{1}; Delete; }{ Surface{4}; Delete; }

Mesh.CharacteristicLengthMin = 1.5e-4;
Mesh.CharacteristicLengthMax = 4.0e-4;

eps = 1e-6;
c_bot[]   = Curve In BoundingBox{Rplug-eps, -eps,       -eps, Rcr+eps, eps,          eps};
c_top[]   = Curve In BoundingBox{Rplug-eps, Lz-eps,     -eps, Rcr+eps, Lz+eps,       eps};
c_outer[] = Curve In BoundingBox{Rcr-eps,   -eps,       -eps, Rcr+eps, Lz+eps,       eps};
c_inner[] = Curve In BoundingBox{Ri-eps,    Lplug-eps,  -eps, Ri+eps,  Lz-Lplug+eps, eps};

Physical Surface("zircaloy", 1) = {1};
Physical Surface("chromium", 2) = {4};

Physical Curve("bottom", 11) = {c_bot[]};
Physical Curve("outer",  12) = {c_outer[]};
Physical Curve("top",    13) = {c_top[]};
Physical Curve("inner",  14) = {c_inner[]};