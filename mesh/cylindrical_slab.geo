//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {0, 10, 0, 1.0};
//+
Point(4) = {4, 0., 0, 1.0};
//+
Point(5) = {0., 4, 0, 1.0};
//+
Point(6) = {4, 4, 0, 1.0};
//+
Line(1) = {5, 6};
//+
Line(2) = {6, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 5};
//+
Line(5) = {4, 2};
//+
Line(6) = {5, 3};
//+
Point(7) = {7.0710678118654755, 7.0710678118654755, 0, 1.0};
//+
Line(7) = {6, 7};
//+
Circle(8) = {2, 1, 7};
//+
Circle(9) = {7, 1, 3};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
Transfinite Line {3} = 10 Using Progression 1.;
Transfinite Line {4} = 10 Using Progression 1.;
Transfinite Line {1} = 10;
Transfinite Line {2} = 10;
Transfinite Surface {1} = {1, 4, 6, 5};
Recombine Surface {1};

//+
Curve Loop(2) = {7, 9, -6, 1};
//+
Plane Surface(2) = {2};
Transfinite Line {7} = 10 Using Progression 1.;
Transfinite Line {9} = 10 Using Progression 1.;
Transfinite Line {6} = 10;
Transfinite Line {1} = 10;
Transfinite Surface {2} = {5, 6, 7, 3};
Recombine Surface {2};
//+
Curve Loop(3) = {5, 8, -7, 2};
//+
Plane Surface(3) = {3};
Transfinite Line {5} = 10 Using Progression 1.;
Transfinite Line {7} = 10 Using Progression 1.;
Transfinite Line {2} = 10;
Transfinite Line {8} = 10;
//Transfinite Surface {3} = {3,6,4,2};

Transfinite Surface {3} = {4, 2, 7, 6};
Recombine Surface {3};

//+
//Physical Surface(3) = {3};
//+
//+
Transfinite Surface {2} = {5, 6, 7, 3};
//+
//Physical Surface(4) = {1};
//+
//Physical Surface(5) = {3};
//+
//Physical Surface(6) = {2};
//+
Extrude {0, 0, 5} {
  Surface{1}; Layers {10}; Recombine;
}
Extrude {0, 0, 5} {
  Surface{2}; Layers {10}; Recombine;
}
Extrude {0, 0, 5} {
  Surface{3}; Layers {10}; Recombine;
}
//+
//+
Physical Volume(76) = {1, 2, 3};
//+
//Transfinite Volume(76);
//Recombine Volume(76);

Mesh 3; // 3 = volumes are meshed

Save "cylindrical_slab.msh";
