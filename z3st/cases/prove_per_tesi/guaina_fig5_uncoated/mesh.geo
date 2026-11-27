SetFactory("OpenCASCADE");

Ri    = 4.18e-3;
Ro    = 4.75e-3;
Rplug = 2.50e-3;
Lplug = 5.0e-3;
Lz    = 77.0e-3;

Rectangle(1) = {Rplug, 0.0,        0.0, Ro - Rplug, Lplug};
Rectangle(2) = {Ri,    Lplug,      0.0, Ro - Ri,    Lz - 2.0 * Lplug};
Rectangle(3) = {Rplug, Lz - Lplug, 0.0, Ro - Rplug, Lplug};

BooleanUnion{ Surface{1}; Delete; }{ Surface{2, 3}; Delete; }

Mesh.CharacteristicLengthMin = 2.0e-4;
Mesh.CharacteristicLengthMax = 5.0e-4;

eps = 1e-6;
c_bot[]   = Curve In BoundingBox{Rplug-eps, -eps,      -eps, Ro+eps, eps,         eps};
c_top[]   = Curve In BoundingBox{Rplug-eps, Lz-eps,    -eps, Ro+eps, Lz+eps,      eps};
c_outer[] = Curve In BoundingBox{Ro-eps,    -eps,      -eps, Ro+eps, Lz+eps,      eps};
c_inner[] = Curve In BoundingBox{Ri-eps,    Lplug-eps, -eps, Ri+eps, Lz-Lplug+eps, eps};

Physical Surface("zircaloy", 1) = {1};
Physical Curve("bottom",   11) = {c_bot[]};
Physical Curve("outer",    12) = {c_outer[]};
Physical Curve("top",      13) = {c_top[]};
Physical Curve("inner",    14) = {c_inner[]};
