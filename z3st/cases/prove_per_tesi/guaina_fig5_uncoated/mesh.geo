SetFactory("OpenCASCADE");

// Same tube as guaina_fig5_coated, without the Cr layer. End plugs are not
// meshed: the closed end enters as the end-cap load on "top"
// (boundary_conditions.yaml), "bottom" is the mid-plane symmetry.

Ri = 4.18e-3;    // inner radius (m)
Ro = 4.75e-3;    // outer radius (m)
Lz = 77.0e-3;    // height (m); the solution does not depend on z

n_C = 11;        // nodes across the cladding
n_z = 2;         // nodes along z (solution is z-independent)

Point(1) = {Ri, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {Ro, Lz, 0};
Point(4) = {Ri, Lz, 0};

Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // outer radius
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // inner radius

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = n_C;
Transfinite Line {2, 4} = n_z;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface("zircaloy", 1) = {1};
Physical Curve("bottom", 11) = {1};
Physical Curve("outer", 12) = {2};
Physical Curve("top", 13) = {3};
Physical Curve("inner", 14) = {4};
