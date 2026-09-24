// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO: axisymmetric wall cut in two at r = Rm, structured mesh.
//  split = 1 tags the two parts as two materials; split = 0 tags them as
//  one. The geometry and the node set are the same either way, which is
//  what lets non-regression.py compare the two runs node by node.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

DefineConstant[ split = {1, Name "split"} ];

Ri = 2.000;  // inner radius (m)
Rm = 2.050;  // interface radius (m): near the inner face, where the source is high
Ro = 2.400;  // outer radius (m)
Lz = 0.200;  // height (m)

// Coarse on purpose: the doubled interface source scales with the element
// width, so a fine radial mesh would hide the bug the case exists to catch.
n_in  = 6;   // nodes across Ri..Rm
n_out = 36;  // nodes across Rm..Ro
n_z   = 5;   // nodes along z

Point(1) = {Ri, 0, 0};
Point(2) = {Rm, 0, 0};
Point(3) = {Ro, 0, 0};
Point(4) = {Ro, Lz, 0};
Point(5) = {Rm, Lz, 0};
Point(6) = {Ri, Lz, 0};

Line(1) = {1, 2};  // bottom, inner part
Line(2) = {2, 3};  // bottom, outer part
Line(3) = {3, 4};  // outer radius
Line(4) = {4, 5};  // top, outer part
Line(5) = {5, 6};  // top, inner part
Line(6) = {6, 1};  // inner radius
Line(7) = {2, 5};  // interface r = Rm

Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Transfinite Line {1, 5} = n_in;
Transfinite Line {2, 4} = n_out;
Transfinite Line {3, 6, 7} = n_z;
Transfinite Surface {1, 2};

If (split == 1)
  Physical Surface("steel_in", 10) = {1};
  Physical Surface("steel_out", 11) = {2};
Else
  Physical Surface("steel", 10) = {1, 2};
EndIf
Physical Curve("inner_radius", 1) = {6};
Physical Curve("outer_radius", 2) = {3};
Physical Curve("bottom", 3) = {1, 2};
Physical Curve("top", 4) = {4, 5};
