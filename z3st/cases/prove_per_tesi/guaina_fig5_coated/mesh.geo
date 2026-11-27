SetFactory("OpenCASCADE");

Ri    = 4.18e-3;    // Raggio interno Zircaloy (m)
Ro    = 4.75e-3;    // Raggio esterno Zircaloy / interfaccia Cr (m)
Rcr   = 4.765e-3;   // Raggio esterno Cromo (coating = 15 um) (m)
Rplug = 2.50e-3;    // Raggio interno tappi terminali (m)
Lplug = 5.0e-3;     // Altezza tappi terminali (m)
Lz    = 77.0e-3;    // Altezza totale provino (m)

n_C = 11;     // nodes across the cladding
n_L = 5;      // nodes across the coating
n_z = 2;      // nodes along z (solution is z-independent)

Point(1) = {Ri, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {Rcr, 0, 0};
Point(4) = {Rcr, Lz, 0};
Point(5) = {Ro, Lz, 0};
Point(6) = {Ri, Lz, 0};

Line(1) = {1, 2};  // bottom, cladding
Line(2) = {2, 3};  // bottom, coating
Line(3) = {3, 4};  // outer radius of coating
Line(4) = {4, 5};  // top, coating
Line(5) = {5, 6};  // top, cladding
Line(6) = {6, 1};  // outer radius of cladding
Line(7) = {2, 5};  // inner radius of cladding

Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Transfinite Line {1, 5} = n_C;
Transfinite Line {2, 4} = n_L;
Transfinite Line {3, 6, 7} = n_z;
Transfinite Surface {1, 2};
Recombine Surface {1, 2};

Physical Surface("zircaloy", 10) = {1};
Physical Surface("chromium", 11) = {2};
Physical Curve("inner", 1) = {6};
Physical Curve("outer", 2) = {3};
Physical Curve("bottom", 3) = {1, 2};
Physical Curve("top", 4) = {4, 5};