# UO2 pellet thermal-shock fracture - phase-field simulation
## z3st application case

### Reference
McClenny, Butt, Abdoelatef et al.,
"Experimentally validated multiphysics modeling of fracture induced by
thermal shocks in sintered UO2 pellets",
Journal of Nuclear Materials 565 (2022) 153719.

### Problem description
A sintered UO2 fuel pellet (R = 10 mm, H = 10 mm) is preheated to 750 deg C in a
molten salt bath and then quenched against a cold contact surface at -10 deg C.
Only 1/6 of the outer circumference (a 60 deg wedge) is in contact with the
cold bath; the remaining 5/6 is insulated, mirroring the experimental setup
of McClenny et al. (Fig. 4 of the reference). The rapid surface cooling drives
a steep radial temperature gradient: the cooled periphery contracts while the
hot interior resists, generating a tensile hoop stress at the periphery and
compressive stress at the centre. When the tensile stress exceeds the local
UO2 fracture strength, radial cracks nucleate at the cold surface and
propagate inwards.

This is the same physical mechanism that drives the first-power-ramp pellet
cracking pattern observed in light-water-reactor fuel.

### Simulation setup
| Parameter                     | Value                  | Source                                      |
|-------------------------------|------------------------|---------------------------------------------|
| Pellet radius R               | 10 mm                  | McClenny et al. (Section 6)                 |
| Axial height H                | 10 mm                  | McClenny et al. (Section 6)                 |
| Contact area                  | 1/6 of circumference   | McClenny et al. Fig. 4                      |
| Young's modulus E             | 358 GPa                | McClenny et al. Table 3                     |
| Poisson's ratio nu            | 0.23                   | McClenny et al. Table 3                     |
| Thermal conductivity k        | 5 W/(m K)              | McClenny et al. Table 3                     |
| Specific heat C               | 280 J/(kg K)           | McClenny et al. Table 3 (typical UO2 value) |
| Density rho                   | 10 970 kg/m^3          | McClenny et al. Table 3 (~95% TD)           |
| Thermal expansion alpha       | 1e-5 K^-1              | McClenny et al. Table 3                     |
| Energy release rate Gc        | 80 000 J/m^2           | McClenny et al. Table 3 (80 MPa mm)         |
| Phase-field length scale lc   | 50 um                  | This work (coarsened from McClenny's 1 um)  |

The phase-field length scale lc has been coarsened from McClenny's lc = 1 um
to 50 um in this case so that the simulation remains tractable in 3D on
commodity hardware. The mesh refinement near the contact wall (lc_outer = 25 um)
satisfies the standard h <= lc/2 phase-field resolution criterion. The crack
pattern recovered is therefore qualitatively faithful to the McClenny reference
(radial cracks initiating at the cold contact, propagating inwards) but with
diffuse cracks ~100 um wide rather than the ~2 um of the reference simulation.

### Boundary conditions
- Initial: T = 1023.15 K (750 deg C) uniform
- Thermal: Dirichlet T = 263.15 K (-10 deg C) on the 60 deg contact_wall;
  zero-flux (natural) on the 300 deg insulated_wall
- Mechanical: Dirichlet u = 0 on the bottom face (region 30); traction-free
  elsewhere
- Phase-field: zero-flux (natural) on the entire boundary

### Expected results
1. Thermal: a steep radial temperature gradient develops, with the cooled
   periphery dropping to T_quench within milliseconds of the quench while the
   centre remains near T_initial; the field equilibrates over ~100 s on the
   thermal-diffusion timescale tau = R^2 rho cp / k.
2. Stress: tensile hoop stress at the cold periphery (peaking at ~150-200 MPa
   in the early transient), compressive at the centre.
3. Fracture: radial cracks nucleate at the cold contact and propagate inwards.
   The number and depth of cracks depend on Gc and lc; with the present
   parameters, two to four major radial cracks are expected, consistent with
   the experimental observations of McClenny et al. (Fig. 8).

### Verification steps
1. Thermal only: compare T(r, t) on a wedge through the contact region against
   the analytical solution for a cylinder with a sudden surface-temperature
   change (Bessel-function series). The agreement is qualitative because the
   analytical solution assumes uniform Dirichlet on the entire outer surface,
   while this case applies it on only 1/6 of the circumference.
2. Thermo-elastic (no fracture): verify that the early-time hoop-stress
   distribution along the contact wedge matches the analytical Kingery
   solution for a thermally-stressed disk.
3. Fully coupled: compare the recovered crack pattern with Fig. 8 of
   McClenny et al.

### Running the case
```bash
cd 14_full_cylinder_cracking/
gmsh -3 mesh.geo -format msh2
python3 -m z3st > log.z3st
python3 non-regression.py
```

### Variants
- Full circumference quench: change `contact_wall` to cover the entire outer
  surface and remove the `insulated_wall`. The corresponding axisymmetric
  variant lives in `14_full_cylinder_cracking_2D_rz/`.
- Reactor power ramp: replace the thermal-shock BC with a volumetric source
  q''' = LHR/(pi R^2) and a coolant Robin BC on the outer surface. This gives
  the classic in-reactor pellet cracking pattern.
- Gc sensitivity: sweep Gc in the range 1.5-150 MPa mm to reproduce the
  parametric study in Table 4 / Fig. A.14 of the reference.
