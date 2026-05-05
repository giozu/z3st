# UO2 Pellet Thermal Shock Fracture — Phase-Field Simulation
## Z3ST Application Case

### Reference
McClenny, Butt, Abdoelatef et al.,
*"Experimentally validated multiphysics modeling of fracture induced by
thermal shocks in sintered UO2 pellets"*,
JNM **565** (2022) 153719.

### Problem Description
A sintered UO2 fuel pellet (R = 5 mm) is heated to 750 °C in a molten salt
bath, then quenched by contact with a cold bath at −10 °C. The rapid surface
cooling creates a steep radial temperature gradient: the outer surface contracts
while the hot interior resists, generating **tensile hoop stress at the periphery**
and compressive stress at the center. When the tensile stress exceeds the UO2
fracture strength (~150 MPa), radial cracks nucleate at the outer surface and
propagate inward.

This is the same physical mechanism that causes fuel pellet cracking during
the first power ramp in a power reactor — making it a directly relevant
engineering demonstration.

### Simulation Setup
| Parameter | Value | Source |
|-----------|-------|--------|
| Pellet radius | 5 mm | McClenny et al. |
| Young's modulus | 358 GPa | Govers et al. (2007) |
| Poisson's ratio | 0.23 | Govers et al. (2007) |
| Thermal conductivity | ~5 W/(m·K) | Badry et al. (2019) |
| Specific heat | ~280 J/(kg·K) | Kavazauri et al. (2016) |
| Density | 10,970 kg/m³ | ~95% TD |
| CTE | 1×10⁻⁵ K⁻¹ | McClenny et al. |
| Energy release rate Gc | 80 MPa·mm | McClenny et al. (fitted) |
| Phase-field length scale ℓ | 0.1 mm | This work |

### Boundary Conditions
- **Initial**: T = 750 °C uniform
- **Thermal**: Dirichlet T = −10 °C on outer circumference (full quench)
- **Mechanical**: Traction-free (natural BC)
- **Phase-field**: Zero flux (natural BC)

### Expected Results
1. **Thermal**: Radial temperature gradient develops from outside in, reaching
   steady state in ~300 s.
2. **Stress**: Tensile hoop stress at periphery (~150–200 MPa), compressive at center.
3. **Fracture**: Radial cracks nucleate at outer surface within ~0.01 s of quench,
   propagate inward. Number of cracks depends on Gc and ℓ.

### Verification Steps
1. **Thermal only**: Compare T(r,t) against analytical solution for a cylinder
   with sudden surface temperature change (Bessel function series).
2. **Thermo-elastic** (no fracture): Verify hoop stress matches the analytical
   Kingery solution for a thermally stressed disk.
3. **Full coupled**: Compare crack patterns against Fig. 8 of McClenny et al.

### Running the Case
```bash
cd z3st_uo2_fracture/
gmsh -2 mesh.geo -format msh2
python3 -m z3st > log.z3st
```

### Variants
- **Partial contact**: Apply cold BC on 1/6 of outer circumference only
  (matches the experiment copper tube setup). Compare to Fig. 8.
- **Reactor power ramp**: Replace thermal shock BC with volumetric heating
  q''' = LHR/(π·R²) and coolant BC at outer surface. This gives the
  classic in-reactor pellet cracking pattern.
- **Gc sensitivity**: Sweep Gc = 10, 30, 50, 80, 100, 150 MPa·mm to reproduce
  the parametric study in Fig. A.14 of the paper.
