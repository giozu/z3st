// ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---
//
//  Gmsh GEO for a 3D hexagonal shell (assembly wrapper), hexahedral mesh
//
//  Author: Bianca Funaro
//
//  Parameters can be overridden from the command line, e.g.
//    gmsh -setnumber D 0.25 -setnumber H 1.5 mesh.geo -3
//
// ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---

SetFactory("Built-in");

If (!Exists(D))  D  = 0.200;  EndIf  // outer circumscribed diameter (m)
If (!Exists(H))  H  = 1.000;  EndIf  // height (m)
If (!Exists(t))  t  = 0.0045; EndIf  // wall thickness (m)
If (!Exists(nf)) nf = 20;     EndIf  // elements along each flat
If (!Exists(nt)) nt = 3;      EndIf  // elements through the thickness
If (!Exists(nz)) nz = 40;     EndIf  // elements along the height

Ro = D / 2;                          // outer circumradius
Ri = Ro - 2 * t / Sqrt(3);           // inner circumradius (flats offset by t)

If (!Exists(orientation)) orientation = "flat"; EndIf // hexagone orientation 
// it can be modified with -setstring orientation flat | vertical

// flat:     vertices at 0, 60, ..., 300 deg  -> flats parallel to the x-axis
// vertical: vertices at 90, 150, ..., 30 deg -> flats parallel to the y-axis
If (StrCmp(orientation, "vertical") == 0)
  a0 = Pi / 2;
ElseIf (StrCmp(orientation, "flat") == 0)
  a0 = 0;
Else
  Error("orientation must be 'flat' or 'vertical'");
  Abort;
EndIf

// Vertices at 0, 60, ..., 300 deg: flats parallel to the x-axis
For i In {0:5}
  a = a0+ i * Pi / 3;
  Point(1 + i) = {Ri * Cos(a), Ri * Sin(a), 0, 1.0};  // inner: 1-6
  Point(7 + i) = {Ro * Cos(a), Ro * Sin(a), 0, 1.0};  // outer: 7-12
EndFor

For i In {0:5}
  j = (i + 1) % 6;
  Line(1 + i)  = {1 + i, 1 + j};  // inner flats: 1-6
  Line(7 + i)  = {7 + i, 7 + j};  // outer flats: 7-12
  Line(13 + i) = {1 + i, 7 + i};  // radial: 13-18
EndFor

Transfinite Curve {1:12}  = nf + 1;
Transfinite Curve {13:18} = nt + 1;

// Six trapezoidal sectors, one per flat, nf x nt quadrilaterals each
For i In {0:5}
  j = (i + 1) % 6;
  // inner i->j, radial j, outer j->i, radial i->inner (clockwise from +z)
  Curve Loop(1 + i) = {1 + i, 13 + j, -(7 + i), -(13 + i)};
  Plane Surface(1 + i) = {1 + i};
  Transfinite Surface {1 + i};
  Recombine Surface {1 + i};
EndFor

// Extrude each sector into nz hexahedral layers: out[0] top, out[1] volume,
// out[2..5] lateral faces in curve-loop order (inner, radial, outer, radial)
// Radial faces out[3], out[5] are interior and shared between neighbouring
// sectors (one surface per radial curve): not tagged
top[] = {}; vol[] = {}; inner[] = {}; outer[] = {};
For i In {0:5}
  out[] = Extrude {0, 0, H} { Surface{1 + i}; Layers{nz}; Recombine; };
  top[] += out[0];
  vol[] += out[1];
  inner[] += out[2];
  outer[] += out[4];
EndFor

// Physical groups: tags match geometry.yaml (volumes 1-9, surfaces 10+)
Physical Surface("bottom", 10) = {1:6};
Physical Surface("top", 11)    = top[];
Physical Surface("inner", 12)  = inner[];
Physical Surface("outer", 13)  = outer[];
Physical Volume("steel", 1)    = vol[];