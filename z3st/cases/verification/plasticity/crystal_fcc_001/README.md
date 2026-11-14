# verification/plasticity/crystal_fcc_001: multi-slip crystal plasticity, FCC crystal along [001]

Uniaxial tension of a FCC single crystal along a cube axis, with the library
law `materials/crystal_plasticity.py` on the 12 {111}<110> systems. The stress-strain curve and the
cumulative slip are compared with the same scheme integrated at a material
point in numpy (`uniaxial_reference`).

In the local suite (carries a gold). Not in `cases_ci.txt`: it takes about 80 s.

## Model

```
tau_a       = sigma_tr : P_a,   sigma_tr = C : (eps - ep_n)
dgamma_a/dt = gamma0 (tau_a/g0) |tau_a/g0|^(n-1)
ep          = ep_n + dt sum_a dgamma_a/dt P_a
sigma       = C : (eps - ep)
p           = p_n + dt sum_a |dgamma_a/dt|
```

Small strain, isotropic elasticity, same slip resistance on every system, no
hardening. The slip rates are evaluated on the trial stress, so the update is
explicit in the plastic strain and first-order in time. The Newton tangent
comes from UFL differentiation.

Card (`fcc_crystal.yaml`): `crystal: fcc`, `euler_angles: [0, 0, 0]` (Bunge, degrees),
E = 200 GPa, nu = 0.3, g0 = 200 MPa, gamma0 = 0.001 1/s, n = 5.

## Setup

- Unit cube, 8 hexahedra. Symmetry planes on `xmin`, `ymin`, `zmin`; `xmax` and
  `ymax` free; `zmax` pulled to 1 % strain in 2 s (41 steps, dt = 0.05 s).
- With the crystal axes on the sample axes the response is homogeneous and the
  stress is uniaxial along z.
- Active systems along [001]: 8 of 12, Schmid factor 1/sqrt(6) = 0.408.

## Checks (tolerance 1e-5)

| Metric | Reference | Result |
|---|---|---|
| `sigma_zz_curve` | material-point curve, maximum deviation over the 41 steps | 6.3e-16 |
| `cumulative_slip_final` | material-point cumulative slip | 3.7e-16 |
| `other_stress_components` | zero | 6.8e-16 |
| `saturation_gap` | tracked, not checked | 10.8 % |

`saturation_gap` is the distance of the final stress (475.8 MPa) from the
exact saturation stress of the rate law,

```
sigma_sat = g0 (strain_rate / (gamma0 sum_a |m_a|^(n+1)))^(1/n) = 533.5 MPa
```

It is the time-discretisation error of the scheme, about E times the strain
increment per step (50 MPa here), and it halves with the time step (see the
self-check, `python materials/crystal_plasticity.py`).

![Stress-strain curve](output/stress_strain_curve.png)

Verified against the blessed gold with Z3ST 0.4.1 and FEniCSx 0.11.0.
