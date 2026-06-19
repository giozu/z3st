# UO2 pellet thermal-shock fracture - phase-field simulation (3D)
## z3st application case (McClenny reproducer)

### Role

This is the **3D** reproduction of the McClenny et al.
(JNM 565, 2022) UO2 pellet thermal-shock experiment. It is the only
case-14 variant that captures all dimensions of the experiment:

- A 60-deg azimuthal contact wedge on the lateral surface (1/6 of the
  total pellet surface area, matching McClenny Fig. 4).
- Zero-flux thermal BCs on the top and bottom faces, consistent with
  McClenny's alumina-spacer design which deliberately suppressed
  axial cooling (their Section 3) so that heat transfer was confined
  to the radial direction.
- Full 3D mechanical equilibrium (no plane-strain or plane-stress
  approximation).

The 2D variants serve more limited roles:

- `benchmarks/damage/pellet_quench_2D_xy/` is the 2D plane-strain analogue
  (transverse cross-section with the 60-deg arc). It captures the
  Fig. 8 crack pattern at much lower compute cost. Justified by the
  alumina-spacer's suppression of axial gradients.
- An axisymmetric (r-z) variant is deliberately absent: axisymmetric
  mode cannot represent the 60-deg wedge, so it is unsuitable for the
  cracking benchmark (a thermal-only verification variant existed
  until 2026-06-11).

### Reference
McClenny, Butt, Abdoelatef et al.,
"Experimentally validated multiphysics modeling of fracture induced by
thermal shocks in sintered UO2 pellets",
Journal of Nuclear Materials 565 (2022) 153719.

The phase-field formulation used here follows Ambati, Gerasimov & De Lorenzis,
"A review on phase-field models of brittle fracture and a new fast hybrid
formulation", Computational Mechanics 55 (2015) 383-405. The case is
positioned as a benchmark of the Ambati hybrid model against the Miehe
anisotropic + viscous Allen-Cahn formulation used in McClenny et al. -
see "Phase-field formulation" below.

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

### Phase-field formulation
McClenny et al. solve the Miehe *anisotropic* model with a rate-dependent
Allen-Cahn evolution (their Eq. 10):
- stress: sigma = (1-c)^2 dPsi+/de + dPsi-/de  (mechanical equation is non-linear),
- evolution: eta dc/dt = 2(1-c) Psi+ - (Gc/l)(c - l^2 Lap c), with viscosity
  eta = 1e-8 s/mm.

Z3ST solves the Ambati *hybrid* formulation (Ambati et al. Eq. 27):
- stress: sigma = (1-D)^2 dPsi0/de  (full energy degraded isotropically;
  mechanical equation stays *linear*),
- damage: -l^2 Lap D + D = (2 l / Gc)(1-D) H+, with H+ the spectral-split
  positive elastic-energy history (Miehe split + irreversibility),
- hybrid constraint: in cells where Psi- > Psi+, the driving force H is set
  to zero (z3st softens Ambati's hard "D := 0" projection to "H = 0", which
  is equivalent under monotonic loading like a thermal shock and more
  physical under unloading).

The implementation lives in `z3st/models/damage_model.py` and
`z3st/core/solver.py::_damage_step`. The single-edge-notched benchmark
(`benchmarks/notched_plate_2D`) verifies the formulation against the
Miehe-Ambati canonical results.

The methodological point of running this case is therefore not just
reproducing McClenny's UO2 thermal-shock crack pattern, but doing so with
a *different and cheaper* phase-field formulation: the linear mechanical
block of the hybrid model has roughly an order-of-magnitude lower
per-iteration cost than the anisotropic + viscous formulation (Ambati
Section 3.1), and Z3ST's elliptic AT2 minimisation removes the
viscosity-tuned crack kinetics entirely.

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
   The McClenny reference (Fig. 8 and Fig. 9) reports **two long radial
   cracks** plus a fan of shorter surface cracks distributed along the 60 deg
   contact wedge. The two long cracks do not reach the pellet centre: they
   are arrested by the central compression zone. In Z3ST this arrest
   emerges from two compounding mechanisms: the Miehe spectral split sends
   compressive eigenvalue contributions into Psi-, and the hybrid
   constraint then sets the driving force H = 0 in cells where Psi- > Psi+.
   Z3ST is expected to reproduce the same topology with diffuser crack
   bands (~100 um vs ~2 um), as a consequence of the lc coarsening
   discussed above.

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
cd benchmarks/damage/pellet_quench_3D/
gmsh -3 mesh.geo -format msh2
python3 -m z3st > log.z3st
python3 non-regression.py
```

### Variants
- Full circumference quench: change `contact_wall` to cover the entire outer
  surface and remove the `insulated_wall`.
- Reactor power ramp: replace the thermal-shock BC with a volumetric source
  q''' = LHR/(pi R^2) and a coolant Robin BC on the outer surface. This gives
  the classic in-reactor pellet cracking pattern.
- Gc sensitivity: sweep Gc in the range 1.5-150 MPa mm to reproduce the
  parametric study in Table 4 / Fig. A.14 of the reference.
