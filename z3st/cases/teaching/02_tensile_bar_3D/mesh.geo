// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a round tensile bar, axis along z
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

R  = 0.005;   // bar radius    (m)
Lz = 0.050;   // gauge length  (m)
h  = 0.001;   // element size in the cross-section (m)
nz = 25;      // element layers along the bar

// The cross-section, a disk
Point(1) = {0, 0, 0, h};
Point(2) = {R, 0, 0, h};
Point(3) = {0, R, 0, h};
Point(4) = {-R, 0, 0, h};
Point(5) = {0, -R, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point{1} In Surface{1};

// Extruded, not meshed as a free solid: every facet on the curved surface then lies
// in a vertical plane. A free tetrahedral mesh puts triangles across different heights,
// whose planes tilt, and sigma_zz then pulls on the side of the bar and spoils the
// uniform solution at the 1e-4 level.
out[] = Extrude {0, 0, Lz} { Surface{1}; Layers{nz}; };
// out[0] top, out[1] volume, out[2..5] the faces swept by the four arcs

// Labels match geometry.yaml.
Physical Surface("bottom", 1) = {1};
Physical Surface("top", 2)    = {out[0]};
Physical Surface("outer", 3)  = {out[2], out[3], out[4], out[5]};
Physical Volume("steel", 4)   = {out[1]};
