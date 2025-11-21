// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for a spherical shell, with specified radii
//
//  Author: Giovanni Zullo
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("OpenCASCADE");

// parameters
Ri = DefineNumber[ 2.0, Name "Parameters/Ri" ];
Ro = DefineNumber[ 2.5, Name "Parameters/Ro" ];
ndens = DefineNumber[ 0.08, Name "Parameters/ndens" ];

// define inner and outer spheres (volumes)
Sphere(1) = {0, 0, 0, Ro};
Sphere(2) = {0, 0, 0, Ri};

// boolean difference: shell = outer - inner
v() = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// get surfaces of the resulting shell
s[] = Boundary{ Volume{v()}; };

// you now have two spherical surfaces in s[]
// assign physical groups
Physical Surface("outer") = {s[0]};   // outer sphere
Physical Surface("inner") = {s[1]};   // inner sphere
Physical Volume("mat0")   = {v()};

// mesh settings
Mesh.CharacteristicLengthMin = ndens;
Mesh.CharacteristicLengthMax = ndens;
Mesh 3;
