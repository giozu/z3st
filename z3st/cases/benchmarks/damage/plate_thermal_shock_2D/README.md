# Plate under thermal shock: crack nucleation

## z3st damage benchmark case, after Kamagate et al. (2025)

Configured after the quenching-plate benchmark of Kamagate, Cheng, Abdelmoula,
Danho and Kondo (2025), *An incremental variational method to the coupling
between gradient damage, thermoelasticity and heat conduction*,
C. R. Mecanique 353:1063-1084 (Fig. 1-3).

Excluded from the local suite (`suite_exclude.txt`, about 2 h). Run it by hand.

### Scope

There is no pre-crack and no damage seed (`boundary_conditions.yaml` has no
`damage:` block). Cracks nucleate from the quenched edges. The AT1 model has an
elastic threshold `w1 = 3*Gc/(8*lc)` (strength
`sigma_c = sqrt(3*E*Gc/(8*lc)) ~ 243 MPa`). Where the transient
thermoelastic tension at the cooled surface exceeds it, damage localises
into bands of width about `lc`. A continuous damage front along the edge is
unstable to a periodic perturbation and breaks into an array of edge cracks.

### Set-up

| Item | Value |
|---|---|
| Plate | `L = 25 mm` x `H = 9.8 mm` |
| Quenched edges | bottom, top, left -> Dirichlet `T_B = 300 K` |
| Right edge | `u_x = 0` symmetry plane, adiabatic |
| Pin | 50-um `Clamp_y` segment (removes rigid y-translation) |
| Initial / stress-free T | `T_0 = T_ref = 550 K` (so `dT = 250 K`) |
| Model | AT1, `lc = 0.092 mm`, `split: amor`, `hybrid_constraint: true` |
| Material | `plate_ceramic.yaml` (paper Table 1: E=340 GPa, nu=0.22, Gc=42.47 J/m2, alpha=8e-6) |
| Time | 0 to 5 ms, 50 steps of 100 us |
| Output | XDMF (`output/fields.xdmf`) |

### Run

```
./Allrun
```

Open `output/fields.xdmf` in ParaView and colour by `Damage`. Kamagate Fig. 2d
shows short, roughly parallel cracks at the bottom and top edges, growing
inward, with shorter cracks between the longer ones.

`Allrun` ends with `plot_damage.py`, which writes `output/damage_field.png` from
the XDMF file. `./Allrun_mpi [NP]` runs the solver on `NP` MPI ranks (default 4),
then `plot_damage.py`.

### Parameters to vary

- Shock amplitude: in `plate_ceramic.yaml` set `T_initial = T_ref = 880.0`
  for the `dT = 580 K` case (paper Fig. 3: more and deeper cracks).
- Mesh: `lc_fine` in `mesh.geo` is `3e-5` (about `lc/3`). The `mesh.geo` comment
  gives `2e-5` (`lc/4.5`) for converged crack counts.

### Comparison with Kamagate et al.

Counted by hand on the damage field of the stored run at `t = 5 ms` (no script
computes these yet): 13 cracks on each quenched edge, spacing about 1.9 mm,
penetrating 0.2 to 1.6 mm (mean 0.87 mm), alternating deep and shallow. An
automatic count in `non-regression.py` is planned.
Differences with the reference:

- Time. The gold state is at `t = 5 ms`. Kamagate Fig. 2 and 3 are at
  `t = 10 us`. The first time step is 100 us, so there is no snapshot at the
  reference time.
- Crack count. Fig. 2d has roughly 20 to 25 cracks per edge, about twice as
  many as here. The `lc/4.5` mesh has not been run.
- The `dT = 580 K` case has not been run.

Gold values at the last step: `D_max_final = 1.0`, `D_mean_final = 0.0338`,
`E_el_final = 3.549 J`, `E_frac_final = 1.910 J`, `T_mean_final = 403.9 K`.

### Limits

- The number of cracks depends on mesh, `lc` and heterogeneity. The checks to make
  are the pattern and the `dT` trend (more cracks at higher `dT`), not an exact count.
- The coupling is the staggered thermal -> mechanical -> damage loop. The
  kinetic-entropy incremental variational formulation of the paper is not
  implemented. In this benchmark (constant k and c, damage heat neglected) the
  temperature field is one-way coupled.
