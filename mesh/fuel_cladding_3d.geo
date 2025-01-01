SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1.8, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1.5, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 2, 0, 2*Pi};
//+
Curve Loop(1) = {3};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {2};
//+
Plane Surface(2) = {3};
//+
Physical Surface("cladding_bottom", 2) = {1};
//+
Physical Surface("fuel_bottom", 1) = {2};
//+
Extrude {0, 0, 4} {
  Curve{3}; Curve{1}; Curve{2}; Surface{1}; Surface{2}; 
}
//+
Physical Surface("cladding_outer", 16) = {3};
//+
Physical Surface("cladding_inner", 17) = {4};
//+
Physical Surface("fuel_outer", 18) = {5};
//+
Physical Surface("cladding_top", 19) = {8};
//+
Physical Surface("fuel_top", 20) = {10};
//+
Recombine Surface {5, 10, 2};
//+
Recombine Surface {3, 4};
//+
Transfinite Surface {2};
//+
Transfinite Surface {2};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {1};
//+
Recombine Surface {2, 1, 3, 4, 5, 8, 10};
//+
Recombine Surface {3, 4, 1, 8};
//+
Recombine Surface {2, 5, 10};
//+
Recombine Surface {2};
//+
Recombine Surface {2, 10, 5};
//+
Recombine Surface {8, 4, 3, 1};
//+
Recombine Surface {10, 5, 2};
