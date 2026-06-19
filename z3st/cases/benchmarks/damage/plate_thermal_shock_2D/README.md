# Plate under thermal shock — pure crack nucleation

## z3st damage benchmark case (Kamagate et al. 2025 replica)

Reproduces the quenching-plate benchmark of Kamagate, Cheng, Abdelmoula,
Danho & Kondo (2025), *An incremental variational method to the coupling
between gradient damage, thermoelasticity and heat conduction*,
C. R. Mecanique **353**:1063-1084 (Fig. 1-3).

### Scope

There is **no pre-crack and no damage seed** anywhere (see
`boundary_conditions.yaml`: there is no `damage:` block). Cracks must
**nucleate** from the quenched edges. This is the
gradient-damage instability, not seeded propagation: the AT1 model has an
elastic threshold `w1 = 3*Gc/(8*lc)` (strength
`sigma_c = sqrt(3*E*Gc/(8*lc)) ~ 243 MPa`), and where the transient
thermoelastic tension at the cooled surface exceeds it, damage localises
into bands of width ~`lc`. A continuous damage front along the edge is
itself unstable to a periodic perturbation, so it breaks into the
discrete array of edge cracks.

### Set-up

| Item | Value |
|---|---|
| Plate | `L = 25 mm` x `H = 9.8 mm` |
| Quenched edges | bottom, top, left -> Dirichlet `T_B = 300 K` |
| Right edge | `u_x = 0` symmetry plane, adiabatic |
| Pin | 50-um `Clamp_y` segment (removes rigid y-translation) |
| Initial / stress-free T | `T_0 = T_ref = 550 K` (so `dT = 250 K`) |
| Model | AT1, `lc = 0.092 mm`, Amor split, hybrid constraint |
| Material | `plate_ceramic.yaml` (paper Table 1: E=340 GPa, nu=0.22, Gc=42.47 J/m2, alpha=8e-6) |

### Run

```
./Allrun
```

then open `output/fields.xdmf` in ParaView and colour by `damage`. You
should see short, roughly parallel cracks appear at the bottom and top
edges and grow inward, with smaller cracks interleaved between the
longer ones (Kamagate Fig. 2d).

### Knobs

- **Shock amplitude:** in `plate_ceramic.yaml` set `T_initial = T_ref = 880.0`
  for the `dT = 580 K` case — expect more and deeper cracks (paper Fig. 3).
- **Mesh:** `lc_fine` in `mesh.geo` is `lc/3` for a first look; drop it to
  ~`2e-5` (`lc/4.5`) once nucleation is confirmed, for crack counts that
  no longer move with refinement.

### Caveats

- The **number** of cracks is sensitive to mesh / `lc` / heterogeneity —
  a known feature of this benchmark. Match the pattern and the
  `dT`-trend (more cracks at higher `dT`), not an exact count.
- This reproduces the cracking **physics/pattern** via the staggered
  thermo -> mech -> damage loop. It does **not** reproduce the paper's
  kinetic-entropy incremental-variational *formulation*; for this
  benchmark (constant k,c, moderate coupling, damage heat neglected) the
  temperature field is effectively one-way and the pattern is governed by
  the thermoelastic stress + AT1, both of which z3st has.
