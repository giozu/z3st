// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a box Lx Ly with a structured mesh
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

Lx = 1.000;
Ly = 0.500;

W_notch = 0.02;
D_notch = 0.10;

X_left  = (Lx/2) - (W_notch/2);
X_right = (Lx/2) + (W_notch/2);
Y_tip   = Ly - D_notch;

lc_fine = 0.00055;
lc_coarse = 0.0075;

Point(1) = {0,  0,  0, lc_coarse};
Point(2) = {Lx, 0,  0, lc_coarse};
Point(3) = {Lx, Ly, 0, lc_coarse};
Point(4) = {X_right, Ly,    0, lc_fine}; 
Point(5) = {Lx/2,    Y_tip, 0, lc_fine};
Point(6) = {X_left,  Ly,    0, lc_fine};

Point(7) = {0,  Ly, 0, lc_coarse};

Line(1) = {1, 2}; // ymin
Line(2) = {2, 3}; // xmax
Line(3) = {3, 4}; // ymax
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7}; // ymax
Line(7) = {7, 1}; // xmin

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

Field[1] = Distance;
Field[1].PointsList = {5};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_fine;
Field[2].LcMax = lc_coarse;
Field[2].DistMin = 0.02;
Field[2].DistMax = 0.2;
Background Field = 2;

Physical Curve("ymin") = {1};
Physical Curve("xmax") = {2};
Physical Curve("ymax") = {3, 6};
Physical Curve("xmin") = {7};
Physical Surface("steel") = {1};

// Mesh 2;
// Save "mesh.msh";
