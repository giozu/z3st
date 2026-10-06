// ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---
//
//  Gmsh GEO for a 2D hexagonal shell (assembly wrapper), hexahedral mesh
//
//  Author: Bianca Funaro
//
//  Parameters can be overridden from the command line, e.g.
//    gmsh -setnumber D 0.25 -setnumber H 1.5 mesh.geo -3
//
// ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---

SetFactory("Built-in");
// gmsh cannot read yaml, so py helper converts the yaml into command line flags
// .geo uses the if(!Exist(x)) to default
// z3st/utils/geo_args.py: reads the top-level scalars in geometry.yaml and prints them as gmsh flags. It skips name and nested blocks

If (!Exists(D))  D  = 0.200;  EndIf  // outer circumscribed diameter (m)
If (!Exists(t))  t  = 0.0045; EndIf  // wall thickness (m)
If (!Exists(nf)) nf = 20;     EndIf  // elements along each flat
If (!Exists(nt)) nt = 3;      EndIf  // elements through the thickness

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

For i In {0:5}
  a = a0 + i * Pi / 3;
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

// Physical groups: tags match geometry.yaml (volumes 1-3, surfaces 4+)
Physical Curve("inner", 4) = {1:6};
Physical Curve("outer", 5) = {7:12};
Physical Surface("steel", 2)  = {1:6};
