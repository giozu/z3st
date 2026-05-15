// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 2D circular UO2 pellet cross-section (transverse cut),
//  with a 60-degree cold-contact arc on the right side. Modeled as the
//  upper semicircle (y >= 0) with mirror symmetry on the y=0 diameter,
//  matching the McClenny et al. JNM 565 (2022) experiment geometry.
//
//  A radial pre-crack ("crack_seed") is embedded at theta = 15 deg, the
//  middle of the upper-half contact arc. By mirror symmetry this represents
//  two cracks (at +/-15 deg) in the full 60-deg wedge -- consistent with
//  McClenny's "two major (longer) radial cracks" observation. The pre-crack
//  is a zero-width geometric slit (case-19 SENT/SENS style: the slit line
//  is traversed forward and backward in the curve loop, so the surface mesh
//  is on both sides and the slit interior is part of the boundary). A
//  Dirichlet D = 1 BC is applied along the slit in boundary_conditions.yaml,
//  giving the staggered scheme a well-posed problem from step 0.
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

// Geometry
R = 10.0e-3;                          // pellet radius (m) = 10 mm
half_contact_deg = 30.0;              // half of the 60-degree contact arc
half_contact = half_contact_deg * Pi / 180.0;

// Pre-crack: radial slit at theta = 15 deg (middle of upper-half contact arc)
crack_theta_deg  = 15.0;
crack_theta      = crack_theta_deg * Pi / 180.0;
crack_length     = 2.5e-4;            // 250 um = 5 * lc; long enough to be well-resolved,
                                      // short enough to leave room for crack propagation
crack_R_inner    = R - crack_length;  // r-coordinate of the crack tip

// Mesh sizing
lc_outer  = 2.5e-5;                   // 25 um  (lc/2 with phase-field lc = 50 um)
lc_center = 2.0e-4;                   // 200 um (4 lc; coarsens away from rim)
pin_size  = 5.0e-5;                   // 50 um pin segment at (-R, 0) to fix x-translation

// Points
Point(1) = {0, 0, 0, lc_center};                                                    // disc center
Point(2) = {R, 0, 0, lc_outer};                                                      // (+R, 0): start of contact arc (mirror-symmetric edge)
Point(3) = {R*Cos(crack_theta), R*Sin(crack_theta), 0, lc_outer};                    // crack mouth on the rim (theta = 15 deg)
Point(4) = {R*Cos(half_contact), R*Sin(half_contact), 0, lc_outer};                  // end of contact arc on upper half (theta = 30 deg)
Point(5) = {-R, 0, 0, lc_outer};                                                     // (-R, 0): leftmost diameter endpoint
Point(6) = {-R + pin_size, 0, 0, lc_outer};                                          // (-R + delta, 0): inner end of the pin segment
Point(7) = {crack_R_inner*Cos(crack_theta), crack_R_inner*Sin(crack_theta), 0, lc_outer};   // crack tip (interior)

// Boundary curves (counterclockwise around the upper-half disc, with the
// pre-crack slit traversed twice -- forward and backward -- in the curve
// loop so it forms a zero-width interior slit).
Circle(1) = {2, 1, 3};                // 0 -> +15 deg : COLD CONTACT (lower half of the upper-half arc)
Circle(2) = {3, 1, 4};                // +15 -> +30 deg : COLD CONTACT (upper half of the upper-half arc)
Circle(3) = {4, 1, 5};                // +30 -> +180 deg : insulated (zero-flux)
Line(4)   = {5, 6};                   // (-R, 0) -> (-R + pin_size, 0) : pin segment (Clamp_x)
Line(5)   = {6, 1};                   // (-R + pin_size, 0) -> (0, 0) : symmetry diameter, left
Line(6)   = {1, 2};                   // (0, 0) -> (+R, 0) : symmetry diameter, right
Line(7)   = {3, 7};                   // pre-crack slit: from rim mouth (Point 3) to interior tip (Point 7)

// Curve loop: enter the slit at Point 3 going inward (forward 7), reach the
// tip, come back out (backward -7), then continue along the contact arc.
// This creates a one-sided slit of zero width; the surface mesh covers both
// sides and the line "7" is part of the boundary.
Curve Loop(1) = {1, 7, -7, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Physical groups (tags must match those in geometry.yaml)
Physical Surface("uo2", 10) = {1};
Physical Curve("contact_wall", 20)   = {1, 2};     // both halves of the cold contact arc
Physical Curve("insulated_wall", 21) = {3};        // remaining 5/6 of the perimeter
Physical Curve("pin", 30)            = {4};        // 50-um Clamp_x segment
Physical Curve("symmetry", 40)       = {5, 6};     // y=0 mirror symmetry plane (Clamp_y)
Physical Curve("crack_seed", 50)     = {7};        // pre-crack slit; Dirichlet D = 1 applied here

// Refinement field: fine near the contact arc AND near the pre-crack slit,
// coarse toward the centre.
Field[1] = Distance;
Field[1].CurvesList = {1, 2, 7};

Field[2] = Threshold;
Field[2].InField   = 1;
Field[2].SizeMin   = lc_outer;
Field[2].SizeMax   = lc_center;
Field[2].DistMin   = 0.0;
Field[2].DistMax   = 8.0e-3;          // mesh stays fine through the expected crack-growth zone

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm                  = 6;

// Usage:
//   gmsh -2 mesh.geo -format msh2
