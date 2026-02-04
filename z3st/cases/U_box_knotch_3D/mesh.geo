// -----------------------------------------------------------------------------
//
//  Structured mesh for box LxWxH with a rectangular notch at x = 0
//  Multi-block 2D (5 rectangles) + extrude to 3D with Nz layers -> hexahedra
//
//  Author: Giovanni Zullo
//
// -----------------------------------------------------------------------------

SetFactory("Built-in");

// ---------------- Parameters ----------------
L  = 0.1;      // length  (x)
W  = 0.1;      // width   (y)
H  = 0.004;    // height  (z)

// Notch (carved out from x in [0,a], y in [y1,y2])
a  = 0.005;     // notch depth along x
b  = 0.002;     // notch width along y

y1 = (W - b)/2;
y2 = y1 + b;

// Number of points along each direction per block (points, not elements)
// Ensure all are >= 2
NxA = 11;       // along x on [0, a]
NxB = 31;       // along x on [a, L]

// Split y into three bands: [0,y1], [y1,y2], [y2,W]
Ny0 = 21;        // along y on [0,  y1]
Ny1 = 9;        // along y on [y1, y2] (the notch band)
Ny2 = 21;        // along y on [y2, W]

Nz  = 11;       // along z on [0,  H]

// Basic checks to avoid invalid settings
If (NxA < 2 || NxB < 2 || Ny0 < 2 || Ny1 < 2 || Ny2 < 2 || Nz < 2)
  Error("All transfinite counts must be >= 2");
EndIf

// ---------------- 2D Multi-block layout ----------------
// Blocks in (x,y):
//  A : [0,a] x [0,y1]
//  B : [0,a] x [y2,W]
//  C0: [a,L] x [0,y1]
//  C1: [a,L] x [y1,y2]
//  C2: [a,L] x [y2,W]

// --- Points ---
// Lower band
Point(1)  = {0,  0,  0};
Point(2)  = {a,  0,  0};
Point(3)  = {L,  0,  0};
Point(4)  = {0,  y1, 0};
Point(5)  = {a,  y1, 0};
Point(6)  = {L,  y1, 0};

// Middle band
Point(7)  = {0,  y2, 0};
Point(8)  = {a,  y2, 0};
Point(9)  = {L,  y2, 0};

// Upper band
Point(10) = {0,  W,  0};
Point(11) = {a,  W,  0};
Point(12) = {L,  W,  0};

// --- Rectangles A, C0 on lower band ---
// A: [0,a] x [0,y1]
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 1};
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(10) = {10}; // sA

// C0: [a,L] x [0,y1]
Line(5) = {2, 3};
Line(6) = {3, 6};
Line(7) = {6, 5};
Line(8) = {5, 2};
Line Loop(11) = {5, 6, 7, 8};
Plane Surface(11) = {11}; // sC0

// --- Middle band rectangle C1: [a,L] x [y1,y2]
Line(9)  = {5, 8};
Line(10) = {8, 9};
Line(11) = {9, 6};
Line(12) = {6, 5};
Line Loop(12) = {9, 10, 11, 12};
Plane Surface(12) = {12}; // sC1

// --- Upper band rectangles (B, C2) ---
// B: [0,a] x [y2,W]
Line(13) = {7, 8};
Line(14) = {8, 11};
Line(15) = {11, 10};
Line(16) = {10, 7};
Line Loop(13) = {13, 14, 15, 16};
Plane Surface(13) = {13}; // sB

// C2: [a,L] x [y2,W]
Line(17) = {8, 9};
Line(18) = {9, 12};
Line(19) = {12, 11};
Line(20) = {11, 8};
Line Loop(14) = {17, 18, 19, 20};
Plane Surface(14) = {14}; // sC2

// ---------------- Transfinite (2D) ----------------
// A
Transfinite Curve {1, 3} = NxA Using Progression 1; // along x
Transfinite Curve {2, 4} = Ny0 Using Progression 1; // along y
Transfinite Surface {10};
// C0
Transfinite Curve {5, 7} = NxB Using Progression 1;
Transfinite Curve {6, 8} = Ny0 Using Progression 1;
Transfinite Surface {11};
// C1
Transfinite Curve {10, 12} = NxB Using Progression 1;
Transfinite Curve {9, 11}  = Ny1 Using Progression 1;
Transfinite Surface {12};
// B
Transfinite Curve {13, 15} = NxA Using Progression 1;
Transfinite Curve {14, 16} = Ny2 Using Progression 1;
Transfinite Surface {13};
// C2
Transfinite Curve {17, 19} = NxB Using Progression 1;
Transfinite Curve {18, 20} = Ny2 Using Progression 1;
Transfinite Surface {14};

// Make quad mesh in 2D
Recombine Surface {10, 11, 12, 13, 14};

// ---------------- Extrude to 3D (Nz layers) ----------------
outA[]  = Extrude {0, 0, H} { Surface{10}; Layers{Nz}; Recombine; };
outC0[] = Extrude {0, 0, H} { Surface{11}; Layers{Nz}; Recombine; };
outC1[] = Extrude {0, 0, H} { Surface{12}; Layers{Nz}; Recombine; };
outB[]  = Extrude {0, 0, H} { Surface{13}; Layers{Nz}; Recombine; };
outC2[] = Extrude {0, 0, H} { Surface{14}; Layers{Nz}; Recombine; };

// ---------------- Physical groups ----------------
// Volumes
Physical Volume("steel") = {outA[1], outC0[1], outC1[1], outB[1], outC2[1]};

// Bottom / Top surfaces
Physical Surface("zmin") = {10, 11, 12, 13, 14}; // Original 2D surfaces are the bottom
Physical Surface("zmax")    = {outA[0], outC0[0], outC1[0], outB[0], outC2[0]};

// x faces
Physical Surface("xmin") = {outA[5], outB[5]};

Physical Surface("xmax") = {outC0[3], outC1[4], outC2[3]};

// y faces
Physical Surface("ymin") = {outA[2], outC0[2]};

// face_yy_1 (y=W plane)
Physical Surface("ymax") = {outB[4], outC2[4]};

// Notch faces (optional groups)
Physical Surface("notch") = {outC1[2]};
Physical Surface("notch_sides_y1") = {outA[4]};
Physical Surface("notch_sides_y2") = {outB[2]};

// Save
Mesh.MshFileVersion = 2.2;
Mesh 3;
Save "mesh.msh";
