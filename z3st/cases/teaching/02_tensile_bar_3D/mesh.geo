// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a quarter of a round tensile bar, axis along z
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

R  = 0.005;   // bar radius    (m)
Lz = 0.050;   // gauge length  (m)
h  = 0.001;   // element size in the cross-section (m)
nz = 25;      // element layers along the bar

// Quarter of the cross-section, in the plane z = 0
Point(1) = {0, 0, 0, h};
Point(2) = {R, 0, 0, h};
Point(3) = {0, R, 0, h};
Line(1)   = {1, 2};
Circle(2) = {2, 1, 3};
Line(3)   = {3, 1};
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

// Extruded, not meshed as a solid: every facet on the curved surface then lies
// in a vertical plane. A free tetrahedral mesh of the cylinder puts triangles
// across different heights, whose planes tilt, and sigma_zz then pulls on the
// side of the bar and spoils the uniform solution at the 1e-4 level.
out[] = Extrude {0, 0, Lz} { Surface{1}; Layers{nz}; };
// out[0] top, out[1] volume, out[2..4] the faces swept by curves 1, 2, 3

// Labels match geometry.yaml.
Physical Surface("symx", 1)   = {out[4]};   // swept by curve 3, on x = 0
Physical Surface("symy", 2)   = {out[2]};   // swept by curve 1, on y = 0
Physical Surface("bottom", 3) = {1};
Physical Surface("top", 4)    = {out[0]};
Physical Surface("outer", 5)  = {out[3]};
Physical Volume("steel", 6)   = {out[1]};
