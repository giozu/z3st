// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a whole round tensile bar, axis along z
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

R  = 0.005;   // bar radius    (m)
Lz = 0.050;   // gauge length  (m)
h  = 0.001;   // element size in the cross-section (m)
nz = 25;      // element layers along the bar

// The cross-section as four quarter disks, so that the planes x = 0 and y = 0 exist
// as interior surfaces of the mesh: the bar is held on them (see boundary_conditions.yaml).
Point(1) = {0, 0, 0, h};
Point(2) = {R, 0, 0, h};
Point(3) = {0, R, 0, h};
Point(4) = {-R, 0, 0, h};
Point(5) = {0, -R, 0, h};
Line(1) = {1, 2}; Line(2) = {1, 3}; Line(3) = {1, 4}; Line(4) = {1, 5};
Circle(5) = {2, 1, 3}; Circle(6) = {3, 1, 4}; Circle(7) = {4, 1, 5}; Circle(8) = {5, 1, 2};
Curve Loop(1) = {1, 5, -2}; Plane Surface(1) = {1};
Curve Loop(2) = {2, 6, -3}; Plane Surface(2) = {2};
Curve Loop(3) = {3, 7, -4}; Plane Surface(3) = {3};
Curve Loop(4) = {4, 8, -1}; Plane Surface(4) = {4};

// Extruded, not meshed as a free solid: every facet on the curved surface then lies
// in a vertical plane. A free tetrahedral mesh puts triangles across different heights,
// whose planes tilt, and sigma_zz then pulls on the side of the bar and spoils the
// uniform solution at the 1e-4 level.
Extrude {0, 0, Lz} { Surface{1, 2, 3, 4}; Layers{nz}; }
Coherence;

e = 1e-6;
symx[]   = Surface In BoundingBox{-e, -R-e, -e, e, R+e, Lz+e};
symy[]   = Surface In BoundingBox{-R-e, -e, -e, R+e, e, Lz+e};
bottom[] = Surface In BoundingBox{-R-e, -R-e, -e, R+e, R+e, e};
top[]    = Surface In BoundingBox{-R-e, -R-e, Lz-e, R+e, R+e, Lz+e};
all[]    = Surface{:};
outer[]  = all[];
outer[] -= {symx[], symy[], bottom[], top[]};

// Labels match geometry.yaml.
Physical Surface("symx", 1)   = {symx[]};
Physical Surface("symy", 2)   = {symy[]};
Physical Surface("bottom", 3) = {bottom[]};
Physical Surface("top", 4)    = {top[]};
Physical Surface("outer", 5)  = {outer[]};
Physical Volume("steel", 6)   = Volume{:};
