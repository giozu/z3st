// -----------------------------------------------------------------------------
// Structured HEX mesh for box LxWxH with a rectangular notch at x = 0
// Multi-block 2D (5 rectangles) + extrude to 3D with Nz layers -> hexahedra
// -----------------------------------------------------------------------------
SetFactory("Built-in");

// ---------------- Parameters ----------------
L  = 0.1;      // length  (x)
W  = 0.1;      // width   (y)
H  = 0.004;    // height  (z)

// Notch (carved out from x in [0,a], y in [y1,y2])
a  = 0.02;     // notch depth along x
b  = 0.02;     // notch width  along y
y1 = (W - b)/2;
y2 = y1 + b;

// ---------------- Structured counts (EDIT ME) ----------------
// Number of points along each direction per block (points, not elements)
// Ensure all are >= 2
NxA = 9;        // along x on [0, a]
NxB = 12;       // along x on [a, L]   (total along x ~ NxA + NxB - 1 when joined)
// Split y into three bands: [0,y1], [y1,y2], [y2,W]
Ny0 = 7;        // along y on [0,  y1]
Ny1 = 8;        // along y on [y1, y2] (the notch band)
Ny2 = 7;        // along y on [y2, W]
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
// Lower band [0, y1]
Point(1)  = {0,  0,  0};
Point(2)  = {a,  0,  0};
Point(3)  = {L,  0,  0};
Point(4)  = {0,  y1, 0};
Point(5)  = {a,  y1, 0};
Point(6)  = {L,  y1, 0};

// Middle band [y1, y2]
Point(7)  = {0,  y2, 0};
Point(8)  = {a,  y2, 0};
Point(9)  = {L,  y2, 0};

// Upper band [y2, W]
Point(10) = {0,  W,  0};
Point(11) = {a,  W,  0};
Point(12) = {L,  W,  0};

// --- Rectangles (A, C0) on lower band ---
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
// Capture extruded entities.
// outX[0] = Top surface
// outX[1] = Volume
// outX[2 to 2+N-1] = Side surfaces (N is number of lines in base surface)
// The order of side surfaces in outX[] corresponds to the order of lines
// in the Line Loop of the base surface.

// Volumes
Physical Volume("fuel") = {outA[1], outC0[1], outC1[1], outB[1], outC2[1]};

// Bottom / Top surfaces
Physical Surface("bottom") = {10, 11, 12, 13, 14}; // Original 2D surfaces are the bottom
Physical Surface("top")    = {outA[0], outC0[0], outC1[0], outB[0], outC2[0]};

// Side surfaces
// A: Line Loop(10) = {1, 2, 3, 4} -> {P1-P2, P2-P5, P5-P4, P4-P1}
// C0: Line Loop(11) = {5, 6, 7, 8} -> {P2-P3, P3-P6, P6-P5, P5-P2}
// C1: Line Loop(12) = {9, 10, 11, 12} -> {P5-P8, P8-P9, P9-P6, P6-P5}
// B: Line Loop(13) = {13, 14, 15, 16} -> {P7-P8, P8-P11, P11-P10, P10-P7}
// C2: Line Loop(14) = {17, 18, 19, 20} -> {P8-P9, P9-P12, P12-P11, P11-P8}

// x faces
// face_xx_0 (x=0 plane)
// This comes from Line 4 (P4-P1) of surface 10 (A) and Line 16 (P10-P7) of surface 13 (B)
Physical Surface("face_xx_0") = {outA[5], outB[5]};

// face_xx_1 (x=L plane)
// This comes from Line 6 (P3-P6) of surface 11 (C0), Line 11 (P9-P6) of surface 12 (C1),
// and Line 18 (P9-P12) of surface 14 (C2)
Physical Surface("face_xx_1") = {outC0[3], outC1[4], outC2[3]};

// y faces
// face_yy_0 (y=0 plane)
// This comes from Line 1 (P1-P2) of surface 10 (A) and Line 5 (P2-P3) of surface 11 (C0)
Physical Surface("face_yy_0") = {outA[2], outC0[2]};

// face_yy_1 (y=W plane)
// This comes from Line 15 (P11-P10) of surface 13 (B) and Line 19 (P12-P11) of surface 14 (C2)
Physical Surface("face_yy_1") = {outB[4], outC2[4]};

// Notch faces (optional groups)
// Back wall at x=a (y in [y1,y2]) - this is the "back" of the notch, facing +x
// This surface is generated from Line 9 (P5-P8) of surface 12 (C1) AND Line 14 (P8-P11) of surface 13 (B)
Physical Surface("notch") = {outC1[2], outB[3]}; // (outC1[2] for L9, outB[3] for L14 - assuming the line order in loop for B is L13,L14,L15,L16)

// Side walls at y=y1 and y=y2 for x in [0,a]
// y=y1 for x in [0,a]:
// This is the interior boundary between A and C0 (Line 3: P5-P4, Line 8: P5-P2, reverse direction)
// More accurately, it's the extruded surfaces from Line 3 of surface 10 and Line 12 of surface 12.
Physical Surface("notch_sides") = {outA[4], outC1[5]}; // outA[4] for line 3 (P5-P4), outC1[5] for line 12 (P6-P5)

// y=y2 for x in [0,a]:
// This is the interior boundary between B and C1 (Line 13: P7-P8, Line 9: P5-P8, reverse direction)
// More accurately, it's the extruded surfaces from Line 13 of surface 13 and Line 9 of surface 12.
Physical Surface("notch_sides") += {outB[2], outC1[2]}; // outB[2] for line 13 (P7-P8), outC1[2] for line 9 (P5-P8) - Note: L9 is also used for the "notch" back wall, so be careful.
                                                      // This one is tricky - Line 9 is shared between B and C1. When extruded it forms two distinct surfaces.
                                                      // The "notch" surface from Line 9 is outC1[2] (part of the x=a plane).
                                                      // The "side wall" at y=y2 is the extruded surface from Line 13 (P7-P8) of surface B.

// Let's refine the "notch_sides" carefully based on the definitions:
// "Side walls at y=y1 and y=y2 for x in [0,a]"
// The definition "x in [0,a]" is critical.
// At y=y1, x in [0,a]: This is the extruded surface from Line 3 (P5-P4) of sA.
// At y=y2, x in [0,a]: This is the extruded surface from Line 13 (P7-P8) of sB.
Physical Surface("notch_sides_y1") = {outA[4]};
Physical Surface("notch_sides_y2") = {outB[2]};


// Save
Mesh.MshFileVersion = 2.2;
Mesh 3;
Save "mesh.msh";