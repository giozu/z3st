// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D single-edge notched tension (SENT) test.
//
//  Reproduces the Ambati et al. 2015 (Comput Mech 55:383-405) §4.1
//  benchmark geometry: a 1 mm x 1 mm square plate with a horizontal
//  notch running from the left edge to the centre, at y = Ly/2.
//
//  The notch is a zero-width slit: curve 6 is traversed twice (forward
//  and backward) in the curve loop, so the surface mesh covers both
//  sides and the slit interior is part of the boundary. A Dirichlet
//  D = 1 BC is applied on "crack" in boundary_conditions.yaml.
//
//  IMPORTANT: this is a *2D* mesh. Generated with `gmsh -2`, the result
//  is a triangulation in the z=0 plane with no z thickness, as required
//  by `regime: 2d` (plane strain) in input.yaml.
//
//  Mesh refinement: a Distance + Threshold field is anchored on the
//  notch slit (curve 6) and on an auxiliary embedded line that runs
//  from the notch tip (5) to a midpoint on the right edge (7). The
//  embedded line is NOT part of the surface boundary -- it's marked
//  `Line { 7 } In Surface { 1 };` so gmsh only sees it for sizing.
//  This keeps h ~ lc/5 along the expected Mode-I crack path; without
//  it, the AT2 bandwidth would be undersampled past the notch tip and
//  the crack would smear into a triangular diffuse zone rather than
//  forming a sharp line.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

// Geometry
Lx = 1.0e-3;                          // plate width  (m) = 1 mm
Ly = 1.0e-3;                          // plate height (m) = 1 mm
Dn = 0.5e-3;                          // notch length (m) = 0.5 mm (half the width)

// Mesh sizing
lc_damage = 4.0e-6;                   // phase-field regularisation length (m)
h_fine    = lc_damage / 5.0;          // 0.8 um at the notch tip AND along the expected
                                      // Mode-I crack path (h/lc = 0.2 -- finer than the
                                      // Borden/Miehe ideal of lc/4, ensuring AT2 bandwidth
                                      // localisation).
h_coarse  = Lx / 75.0;                // ~13.3 um in the bulk (far from the crack).

// Corner / notch / crack-path points
Point(1) = {0,  0,    0, h_coarse};   // (0, 0)         bottom-left
Point(2) = {Lx, 0,    0, h_coarse};   // (Lx, 0)        bottom-right
Point(3) = {Lx, Ly,   0, h_coarse};   // (Lx, Ly)       top-right
Point(4) = {0,  Ly,   0, h_coarse};   // (0, Ly)        top-left
Point(5) = {Dn, Ly/2, 0, h_fine};     // (Dn, Ly/2)     notch tip (interior)
Point(6) = {0,  Ly/2, 0, h_fine};     // (0, Ly/2)      notch mouth (on left edge)
Point(7) = {Lx, Ly/2, 0, h_fine};     // (Lx, Ly/2)     right-edge midpoint (boundary)

// Boundary edges (counterclockwise)
Line(1) = {1, 2};                     // y-min : bottom edge
Line(2) = {2, 7};                     // x-max : right edge, lower half
Line(3) = {7, 3};                     // x-max : right edge, upper half
Line(4) = {3, 4};                     // y-max : top edge
Line(5) = {4, 6};                     // x-min : upper half of left edge
Line(6) = {6, 1};                     // x-min : lower half of left edge
Line(7) = {6, 5};                     // notch slit (boundary): from mouth (6) to tip (5)

// Auxiliary embedded line (NOT part of the surface boundary): runs along the
// expected Mode-I crack path from the notch tip (5) to the right-edge
// midpoint (7). Anchors the Distance refinement field so h stays at h_fine
// through the propagation corridor.
Line(8) = {5, 7};

// Curve loop: enter the notch at Point 6 going inward (forward 7), reach the
// tip, come back out (backward -7), then continue along the boundary.
// This creates a one-sided slit of zero width.
Curve Loop(1) = {1, 2, 3, 4, 5, 7, -7, 6};
Plane Surface(1) = {1};

// Embed the auxiliary line in the surface so gmsh meshes around it but does
// NOT treat it as a boundary that splits the domain.
Line { 8 } In Surface { 1 };

// Physical groups (tags must match those in geometry.yaml)
Physical Surface("steel", 6) = {1};
Physical Curve("ymin",  1)   = {1};
Physical Curve("xmax",  2)   = {2, 3};      // both halves of the right edge
Physical Curve("ymax",  3)   = {4};
Physical Curve("xmin",  4)   = {5, 6};      // both halves of the left edge
Physical Curve("crack", 5)   = {7};

// Refinement field: keep h = h_fine through a 25-um band around the notch
// slit (curve 7) and along the auxiliary crack-path line (curve 8); ramp up
// to h_coarse beyond.
Field[1] = Distance;
Field[1].CurvesList = {7, 8};

Field[2] = Threshold;
Field[2].InField   = 1;
Field[2].SizeMin   = h_fine;
Field[2].SizeMax   = h_coarse;
Field[2].DistMin   = 0.0;
Field[2].DistMax   = 1.0e-4;          // 100 um (= 25 lc) transition band

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm                  = 6;

// Usage:
//   gmsh -2 mesh.geo -format msh2
