// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
//
//  Gmsh GEO for the 1D beam of IAEA VP1 (beam_fx.py)
//
//  gmsh $(python -m z3st.utils.geo_args) mesh.geo -1 -o mesh.msh
//
// --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

SetFactory("Built-in");

// H and nz from geometry.yaml (geo_args); the stations block is nested, so
// its values are repeated here and checked against geometry.yaml by beam_fx.py
If (!Exists(H))  H = 4.000; EndIf    // height (m)
If (!Exists(nz)) nz = 80;   EndIf    // elements along the height, h = H / nz

z_cb  = 1.5;   // core_bottom (m)
z_ct  = 2.5;   // core_top (m)
z_lrp = 3.0;   // LRP (m)

// Stations: the beam axis is x
Point(1) = {0.0,   0, 0, 1.0};
Point(2) = {z_cb,  0, 0, 1.0};
Point(3) = {z_ct,  0, 0, 1.0};
Point(4) = {z_lrp, 0, 0, 1.0};
Point(5) = {H,     0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

// Structured mesh, uniform h = H / nz when the stations are multiples of h
Transfinite Curve {1} = Round(z_cb / H * nz) + 1;
Transfinite Curve {2} = Round((z_ct - z_cb) / H * nz) + 1;
Transfinite Curve {3} = Round((z_lrp - z_ct) / H * nz) + 1;
Transfinite Curve {4} = Round((H - z_lrp) / H * nz) + 1;

// Physical groups: the stations are the regions of boundary_conditions.yaml
Physical Point("bottom")      = {1};
Physical Point("core_bottom") = {2};
Physical Point("core_top")    = {3};
Physical Point("LRP")         = {4};
Physical Point("top")         = {5};
Physical Curve("steel")       = {1, 2, 3, 4};
