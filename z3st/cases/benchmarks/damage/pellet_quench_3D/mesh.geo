// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a full cylinder (3D extruded circle)
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

// Parameters
R = 10.0e-3;             // Radius (m) = 10 mm
H = 10.0e-3;             // Height (m) = 10 mm
lc_outer = 2.5e-5;       // Mesh size at outer edge (m)  = lc/2 (lc_phase_field = 50 um)
lc_center = 2.0e-4;      // Mesh size at center (m)      = 4*lc (was 5e-4 = 10*lc, too coarse for crack path)
n_layers = 20;           // Number of layers in extrusion

// Contact region: from 0° to 60° (1/6 of circumference)
Point(1) = {0, 0, 0, lc_center};                          // Center
Point(2) = {R, 0, 0, lc_outer};                           // 0°
Point(3) = {R*Cos(Pi/3), R*Sin(Pi/3), 0, lc_outer};       // 60°
Point(4) = {0, R, 0, lc_center};                          // 90°
Point(5) = {-R, 0, 0, lc_center};                         // 180°
Point(6) = {0, -R, 0, lc_center};                         // 270°

// Arcs
Circle(1) = {2, 1, 3};    // 0° → 60°  (contact, 1/6)
Circle(2) = {3, 1, 4};    // 60° → 90°
Circle(3) = {4, 1, 5};    // 90° → 180°
Circle(4) = {5, 1, 6};    // 180° → 270°
Circle(5) = {6, 1, 2};    // 270° → 360°

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// Extrude to 3D
out[] = Extrude {0, 0, H} {
    Surface{1};
    Layers{n_layers};
};

// Physical groups
Physical Volume("uo2", 10) = {out[1]};
Physical Surface("contact_wall", 20) = {out[2]};                            // 1/6 quench
Physical Surface("insulated_wall", 21) = {out[3], out[4], out[5], out[6]};  // 5/6 insulated
Physical Surface("bottom", 30) = {1};
Physical Surface("top", 40) = {out[0]};

// Mesh settings
Mesh.Algorithm = 6;

// Refine near the outer boundary
Field[1] = Distance;
Field[1].CurvesList = {1};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_outer;
Field[2].SizeMax = lc_center;
Field[2].DistMin = 0.0;
Field[2].DistMax = 8.0e-3;

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Usage:
//   gmsh -3 mesh.geo -format msh2
