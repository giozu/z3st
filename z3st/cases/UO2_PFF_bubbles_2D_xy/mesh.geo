// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  2D elliptical cavity in a rectangular plate (cross-section)
//  Coherent with a given semi-dihedral angle
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// --- Parameters (overridable via -setnumber) ---
// Work on a quarter of the system
If(!Exists(Rp))
  Rp = 11.2e-6;        // Projected radius of the cavity (m)
EndIf
If(!Exists(theta_deg))
  theta_deg = 48.1;    // Semi-dihedral angle (deg)
EndIf
If(!Exists(Lx))
  Lx = 60e-6;
EndIf
If(!Exists(Ly))
  Ly = 60e-6;
EndIf
If(!Exists(h_cavity))
  h_cavity = 0.5e-6;   // fine mesh size at the cavity boundary
EndIf
If(!Exists(h_plate))
  h_plate = 4.0e-6;    // coarse mesh size
EndIf

// Rectangular plate (Quart de symétrie : coin inférieur gauche en 0,0)
Rectangle(1) = {0, 0, 0, Lx, Ly};

//----------Lenticular shape------------------
xc = 0;  // Centre de la bulle (dans le quart)

Macro MakeLenticularbubble
  p1 = newp; Point(p1) = {xc - Rp, 0, 0};
  p2 = newp; Point(p2) = {xc + Rp, 0, 0};
  p5 = newp; Point(p5) = {xc, yc, 0};
  p6 = newp; Point(p6) = {xc, -yc, 0};

  c1 = newc; Circle(c1) = {p1, p5, p2};
  c2 = newc; Circle(c2) = {p2, p6, p1};

  ll = newll; Curve Loop(ll) = {c1, c2};
  s  = news;  Plane Surface(s) = {ll};

  lens_surface = s;
Return

// theta_deg = 90 is the analytical-verification case: the two arc
// centers (p5, p6) both fall on the origin. Building them as a single
// shared center point keeps the geometry a clean circle instead of two
// coincident but distinct points.
If (theta_deg == 90)
  ay = Rp;
  yc = 0;

  p1 = newp; Point(p1) = {xc - Rp, 0, 0};
  p2 = newp; Point(p2) = {xc + Rp, 0, 0};
  pc = newp; Point(pc) = {xc, 0, 0};

  c1 = newc; Circle(c1) = {p1, pc, p2};
  c2 = newc; Circle(c2) = {p2, pc, p1};

  ll = newll; Curve Loop(ll) = {c1, c2};
  s  = news;  Plane Surface(s) = {ll};

  lens_surface = s;
Else
  ay = Rp * Tan(theta_deg * Pi / 360.0);
  yc = (ay*ay - Rp*Rp) / (2*ay);

  Call MakeLenticularbubble;
EndIf

// Subtract the lenticular bubble from the rectangle
BooleanDifference{ Surface{1}; Delete; }{ Surface{lens_surface}; Delete; }

//---------------------------------------------------------
// --- Groupes Physiques ---
Physical Surface("uo2") = {1};

// Sélection dynamique des frontières par boîte englobante pour éviter les erreurs d'IDs post-opération booléenne
// eps must stay well above OpenCASCADE's built-in confusion tolerance
// (Precision::Confusion() = 1e-7): boolean-cut curves get their bbox
// padded by that amount, so eps <= ~1e-7 makes the BoundingBox tests
// flip unpredictably and silently drop the "cavity" group.
eps = 0.5e-6;
c_ymin() = Curve In BoundingBox{Rp-eps, -eps, -eps, Lx+eps, eps, eps};
c_xmin() = Curve In BoundingBox{-eps, ay-eps, -eps, eps, Ly+eps, eps};
c_xmax() = Curve In BoundingBox{Lx-eps, -eps, -eps, Lx+eps, Ly+eps, eps};
c_ymax() = Curve In BoundingBox{-eps, Ly-eps, -eps, Lx+eps, Ly+eps, eps};
c_cavity() = Curve In BoundingBox{-eps, -eps, -eps, Rp+eps, ay+eps, eps};

Printf("Boundary group ymin: %g curve(s) found", #c_ymin());
Printf("Boundary group xmin: %g curve(s) found", #c_xmin());
Printf("Boundary group xmax: %g curve(s) found", #c_xmax());
Printf("Boundary group ymax: %g curve(s) found", #c_ymax());
Printf("Boundary group cavity: %g curve(s) found", #c_cavity());

Physical Curve("ymin") = {c_ymin()};
Physical Curve("xmin") = {c_xmin()};
Physical Curve("xmax") = {c_xmax()};
Physical Curve("ymax") = {c_ymax()};
Physical Curve("cavity") = {c_cavity()};

// Mesh Refinement
Field[1] = Distance;
Field[1].CurvesList = {c_cavity(), c_ymin()};
Field[1].NumPointsPerCurve = 400;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_cavity;
Field[2].SizeMax = h_plate;
Field[2].DistMin = 0.5e-6;   // Taille minimale jusqu'à 0.5 µm des courbes
Field[2].DistMax = 10.0e-6;  // Transition progressive vers la taille grossière sur 10 µm

Background Field = 2;

// Mesh generation options
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; // Frontal-Delaunay for better quality in 2D
