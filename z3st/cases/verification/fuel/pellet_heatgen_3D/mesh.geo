// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a solid 3D pellet
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ro = 0.005;        // Outer radius (m)
Lz = 0.005;        // Modelled height (m)

s  = 0.45 * Ro;          // half-side of the central square (m)
q  = Ro / Sqrt(2.0);     // outer block corners sit at +/-45 deg on the circle

n_t = 18;          // divisions along each square side / arc
n_r = 26;          // radial divisions in the outer blocks
n_z = 12;          // divisions along the height

// --- points ---
Point(0) = {0, 0, 0};       // arc centre
Point(1) = { s,  s, 0};     // inner square corners
Point(2) = {-s,  s, 0};
Point(3) = {-s, -s, 0};
Point(4) = { s, -s, 0};
Point(5) = { q,  q, 0};     // outer (circle) corners, 45/135/225/315 deg
Point(6) = {-q,  q, 0};
Point(7) = {-q, -q, 0};
Point(8) = { q, -q, 0};

// --- lines: inner square ---
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// --- lines: radials (inner corner -> outer corner) ---
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};

// --- arcs: outer circle ---
Circle(9)  = {5, 0, 6};   // top
Circle(10) = {6, 0, 7};   // left
Circle(11) = {7, 0, 8};   // bottom
Circle(12) = {8, 0, 5};   // right

// --- surfaces: central square + 4 outer blocks ---
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 9, -6, -1};      // top block
Plane Surface(2) = {2};
Curve Loop(3) = {6, 10, -7, -2};     // left block
Plane Surface(3) = {3};
Curve Loop(4) = {7, 11, -8, -3};     // bottom block
Plane Surface(4) = {4};
Curve Loop(5) = {8, 12, -5, -4};     // right block
Plane Surface(5) = {5};

// --- structured (hex) base mesh ---
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12} = n_t + 1;
Transfinite Curve {5, 6, 7, 8} = n_r + 1;

Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Surface {2} = {1, 5, 6, 2};
Transfinite Surface {3} = {2, 6, 7, 3};
Transfinite Surface {4} = {3, 7, 8, 4};
Transfinite Surface {5} = {4, 8, 5, 1};
Recombine Surface {1, 2, 3, 4, 5};

// --- 3D extrusion ---
out[] = Extrude {0, 0, Lz} { Surface{1, 2, 3, 4, 5}; Layers{n_z}; Recombine; };

// --- physical groups (tags must match geometry.yaml::labels) ---
Physical Surface("bottom", 3) = {1, 2, 3, 4, 5};
Physical Surface("top", 4)    = {out[0], out[6], out[12], out[18], out[24]};
Physical Surface("outer", 2)  = {out[9], out[15], out[21], out[27]};
Physical Volume("uo2", 10)    = {out[1], out[7], out[13], out[19], out[25]};

// --- mesh generation ---
Mesh.ElementOrder = 1;
Mesh 3;

// Usage:
//   gmsh -3 mesh.geo -format msh2
