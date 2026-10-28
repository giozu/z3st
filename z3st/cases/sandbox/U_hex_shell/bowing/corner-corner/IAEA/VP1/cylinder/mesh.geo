// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 3D circular tube with the same mid-wall perimeter as the
//  hexagonal shell of ../mesh.geo, hexahedral mesh with the same layout:
//  six 60-degree sectors, nf elements along each sector, nt through the
//  wall, nz layers along the height. Reference for the corner effect in VP1.
//
//  Author: Bianca Funaro
//  Version: 0.4.1 (2026)
//
//  Parameters from geometry.yaml through z3st.utils.geo_args, e.g.
//    gmsh $(python3 -m z3st.utils.geo_args) mesh.geo -3
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

If (!Exists(Rm)) Rm = 0.0716197244; EndIf // mid-wall radius (m): 6 L_mid / (2 pi), L_mid = 0.075 m
If (!Exists(H))  H  = 4.000;  EndIf        // height (m)
If (!Exists(t))  t  = 0.003;  EndIf        // wall thickness (m)
If (!Exists(nf)) nf = 8;      EndIf        // elements along each 60-degree sector
If (!Exists(nt)) nt = 1;      EndIf        // elements through the thickness
If (!Exists(nz)) nz = 400;    EndIf        // elements along the height

Ri = Rm - t / 2;
Ro = Rm + t / 2;

If (!Exists(orientation)) orientation = "flat"; EndIf
// same sector boundaries as the hexagon vertices:
// flat: 0, 60, ..., 300 deg; vertical: 90, 150, ..., 30 deg
If (StrCmp(orientation, "vertical") == 0)
  a0 = Pi / 2;
ElseIf (StrCmp(orientation, "flat") == 0)
  a0 = 0;
Else
  Error("orientation must be 'flat' or 'vertical'");
  Abort;
EndIf

Point(100) = {0, 0, 0, 1.0};  // centre of the arcs
For i In {0:5}
  a = a0 + i * Pi / 3;
  Point(1 + i) = {Ri * Cos(a), Ri * Sin(a), 0, 1.0};  // inner: 1-6
  Point(7 + i) = {Ro * Cos(a), Ro * Sin(a), 0, 1.0};  // outer: 7-12
EndFor

For i In {0:5}
  j = (i + 1) % 6;
  Circle(1 + i)  = {1 + i, 100, 1 + j};  // inner arcs: 1-6
  Circle(7 + i)  = {7 + i, 100, 7 + j};  // outer arcs: 7-12
  Line(13 + i)   = {1 + i, 7 + i};       // radial: 13-18
EndFor

Transfinite Curve {1:12}  = nf + 1;
Transfinite Curve {13:18} = nt + 1;

// Six sectors, nf x nt quadrilaterals each (same loop as the hexagon)
For i In {0:5}
  j = (i + 1) % 6;
  Curve Loop(1 + i) = {1 + i, 13 + j, -(7 + i), -(13 + i)};
  Plane Surface(1 + i) = {1 + i};
  Transfinite Surface {1 + i};
  Recombine Surface {1 + i};
EndFor

// out[0] top, out[1] volume, out[2..5] lateral faces (inner, radial, outer, radial)
top[] = {}; vol[] = {}; inner[] = {}; outer[] = {};
For i In {0:5}
  out[] = Extrude {0, 0, H} { Surface{1 + i}; Layers{nz}; Recombine; };
  top[] += out[0];
  vol[] += out[1];
  inner[] += out[2];
  outer[] += out[4];
EndFor

// Physical groups: tags match geometry.yaml
Physical Surface("bottom", 10) = {1:6};
Physical Surface("top", 11)    = top[];
Physical Surface("inner", 12)  = inner[];
Physical Surface("outer", 13)  = outer[];
Physical Volume("steel", 1)    = vol[];
