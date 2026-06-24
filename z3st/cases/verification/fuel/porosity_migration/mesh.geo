// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a 22.5 degree sector (2D Cartesian)
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

Ro = 0.002675; // Pellet radius (2.675 mm)
angle = 22.5 * Pi / 180.0; // 22.5 degrees

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {Ro, 0.0, 0.0};
Point(3) = {Ro * Cos(angle), Ro * Sin(angle), 0.0};

Line(1) = {1, 2};             // Bottom boundary (theta = 0)
Circle(2) = {2, 1, 3};        // Outer boundary (r = Ro)
Line(3) = {3, 1};             // Top boundary (theta = 22.5 deg)

Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

// Structured mesh: 100 radial elements, 6 angular elements
Transfinite Line {1, 3} = 101;
Transfinite Line {2} = 7;
Transfinite Surface {1};

Physical Surface("fuel", 10) = {1};
Physical Curve("bottom", 3) = {1};
Physical Curve("outer", 2) = {2};
Physical Curve("top", 4) = {3};
