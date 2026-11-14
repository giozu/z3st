# verification/plasticity/crystal_fcc_001_hardening: crystal plasticity with slip-system hardening

The problem of `../crystal_fcc_001_implicit` (FCC single crystal, 12 {111}<110>
systems, uniaxial tension along [001], implicit law through an external
operator) with a slip resistance that evolves on every system.

In the local suite (carries a gold). Not in `cases_ci.txt`: it takes about
140 s and needs `dolfinx-external-operator`.

## Model

```
dgamma_a/dt = gamma0 (tau_a/g_a) |tau_a/g_a|^(n-1)
dg_a/dt     = sum_b q_ab h0 (1 - g_b/g_sat)^a_hard |dgamma_b/dt|,   g_a(0) = g0
q_ab        = 1 (a = b), q_latent (a != b)
```

Backward Euler on the plastic strain and on the 12 slip resistances together:
21 unknowns per quadrature point, solved by Newton in numpy; the consistent
tangent follows from the Jacobian of the local problem. The slip resistances
are kept per material in `model.slip_state_n` and are rolled back with the rest of
the state when a time step is rejected.

Card (`fcc_crystal.yaml`): `hardening: saturation`, g0 = 200 MPa, h0 = 5 GPa, g_sat = 300 MPa,
a_hard = 2.25, q_latent = 1.4; the rest as in `../crystal_fcc_001`.

## Checks (tolerance 1e-5)

The reference is the material-point integration of the same scheme with
`scipy.optimize.root` on strain, plastic strain and slip resistances
(`uniaxial_reference(..., implicit=True)`), which does not use the local
Newton solve of the finite-element law.

| Metric | Reference | Result |
|---|---|---|
| `sigma_zz_curve` | material-point curve, maximum deviation over the 41 steps | 2.5e-12 |
| `cumulative_slip_final` | material-point cumulative slip | 1.2e-12 |
| `other_stress_components` | zero | 5.1e-16 |
| `saturation_gap` | tracked, not checked | 30.5 % |

The final stress is 555.9 MPa, against 533.5 MPa without hardening. The
saturation stress at g_sat is 800.2 MPa; with a_hard = 2.25 the slip
resistance approaches g_sat slowly and 1 % strain is far from it, so
`saturation_gap` is only tracked here.

![Stress-strain curve](output/stress_strain_curve.png)

Verified against the blessed gold with Z3ST 0.4.1, FEniCSx 0.11.0 and
dolfinx-external-operator 0.10.1.
