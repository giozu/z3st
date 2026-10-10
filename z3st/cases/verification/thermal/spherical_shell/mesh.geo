// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
// Spherical shell, one octant (x, y, z >= 0), structured hexahedra.
//
// The octant of each sphere is split into three quadrilateral patches that
// meet at (1,1,1)/sqrt(3) ("cube-sphere"), so the shell is three hexahedral
// blocks. Patches are true spherical surfaces, lateral faces lie on the planes
// of great circles, and the radial lines carry a geometric progression so the
// cells are thin next to the inner surface, where the gamma source decays.
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

Ri = DefineNumber[ 2.0, Name "Parameters/Ri" ];
Ro = DefineNumber[ 2.5, Name "Parameters/Ro" ];
Na = DefineNumber[ 8, Name "Parameters/Na" ];     // cells along each patch edge
Nr = DefineNumber[ 20, Name "Parameters/Nr" ];    // cells across the thickness
q  = DefineNumber[ 1.140175, Name "Parameters/q" ];   // radial growth ratio, sqrt(1.3), inner to outer

Point(1) = {0, 0, 0};

// p = 10 (inner) or 20 (outer) + index: X Y Z XY YZ XZ C
For s In {0:1}
  R = (s == 0) ? Ri : Ro;
  b = 10 + 10*s;
  a = R/Sqrt(2); c = R/Sqrt(3);
  Point(b+1) = {R, 0, 0};
  Point(b+2) = {0, R, 0};
  Point(b+3) = {0, 0, R};
  Point(b+4) = {a, a, 0};
  Point(b+5) = {0, a, a};
  Point(b+6) = {a, 0, a};
  Point(b+7) = {c, c, c};
  // arcs, numbered b+1 .. b+9
  Circle(b+1) = {b+1, 1, b+4};   // X  -> XY   (z = 0)
  Circle(b+2) = {b+4, 1, b+2};   // XY -> Y    (z = 0)
  Circle(b+3) = {b+2, 1, b+5};   // Y  -> YZ   (x = 0)
  Circle(b+4) = {b+5, 1, b+3};   // YZ -> Z    (x = 0)
  Circle(b+5) = {b+3, 1, b+6};   // Z  -> XZ   (y = 0)
  Circle(b+6) = {b+6, 1, b+1};   // XZ -> X    (y = 0)
  Circle(b+7) = {b+4, 1, b+7};   // XY -> C
  Circle(b+8) = {b+5, 1, b+7};   // YZ -> C
  Circle(b+9) = {b+6, 1, b+7};   // XZ -> C
  // spherical patches, numbered b+1 .. b+3 (near X, Y, Z)
  Curve Loop(b+1) = {b+1, b+7, -(b+9), b+6};
  Surface(b+1) = {b+1} In Sphere {1};
  Curve Loop(b+2) = {b+3, b+8, -(b+7), b+2};
  Surface(b+2) = {b+2} In Sphere {1};
  Curve Loop(b+3) = {b+5, b+9, -(b+8), b+4};
  Surface(b+3) = {b+3} In Sphere {1};
EndFor

// radial lines inner -> outer, numbered 31 .. 37 (X Y Z XY YZ XZ C)
For i In {1:7}
  Line(30+i) = {10+i, 20+i};
EndFor

// lateral planar faces, one per arc: numbered 41 .. 49
// arc k joins points (pa[k], pb[k]); radial line of point j is 30+j
pa[] = {1, 4, 2, 5, 3, 6, 4, 5, 6};
pb[] = {4, 2, 5, 3, 6, 1, 7, 7, 7};
For k In {1:9}
  Curve Loop(40+k) = {10+k, 30+pb[k-1], -(20+k), -(30+pa[k-1])};
  Plane Surface(40+k) = {40+k};
EndFor

// three blocks: near X (arcs 1, 7, 9, 6), Y (3, 8, 7, 2), Z (5, 9, 8, 4)
Surface Loop(1) = {11, 21, 41, 47, 49, 46};
Volume(1) = {1};
Surface Loop(2) = {12, 22, 43, 48, 47, 42};
Volume(2) = {2};
Surface Loop(3) = {13, 23, 45, 49, 48, 44};
Volume(3) = {3};

// structured meshing
Transfinite Curve {11:19, 21:29} = Na + 1;
Transfinite Curve {31:37} = Nr + 1 Using Progression q;
Transfinite Surface {11:13, 21:23, 41:49};
Recombine Surface {11:13, 21:23, 41:49};
Transfinite Volume {1:3};
Recombine Volume {1:3};

// physical groups (tags must match geometry.yaml)
Physical Surface("outer", 1) = {21, 22, 23};
Physical Surface("inner", 2) = {11, 12, 13};
Physical Volume("mat0", 3)   = {1, 2, 3};
Physical Surface("xmin", 4)  = {43, 44};   // x = 0 plane (arcs 3, 4)
Physical Surface("ymin", 5)  = {45, 46};   // y = 0 plane (arcs 5, 6)
Physical Surface("zmin", 6)  = {41, 42};   // z = 0 plane (arcs 1, 2)
