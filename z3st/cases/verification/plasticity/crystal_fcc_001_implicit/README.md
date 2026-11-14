# verification/plasticity/crystal_fcc_001_implicit: implicit crystal plasticity through an external operator

The same problem as `../crystal_fcc_001` (FCC single crystal, 12 {111}<110>
systems, uniaxial tension along [001]) with the backward-Euler time
integration of `materials/crystal_plasticity.py`. The stress is an external
operator of the strain (`dolfinx-external-operator`): at every quadrature point
a Newton solve in numpy finds the end-of-step plastic strain, and the weak form
receives the stress and its consistent tangent.

In the local suite (carries a gold). Not in `cases_ci.txt`: it takes about 70 s
and needs `dolfinx-external-operator`.

## Model

```
tau_a       = sigma : P_a,   sigma = C : (eps - ep)          (end-of-step stress)
dgamma_a/dt = gamma0 (tau_a/g0) |tau_a/g0|^(n-1)
ep          = ep_n + dt sum_a dgamma_a/dt P_a                (local Newton, 9 unknowns)
d sigma / d eps = lambda I(x)I + 2 mu J^-1,   J = I + 2 mu dt sum_a (d rate_a / d tau_a) P_a (x) P_a
```

Selected on the card by one more key next to `stress_function`:

```yaml
stress_function: "materials.crystal_plasticity.crystal_stress"            # output stress
stress_operator: "materials.crystal_plasticity.crystal_stress_operator"   # stress in the weak form
```

## Setup

As `../crystal_fcc_001`: unit cube, 8 hexahedra, symmetry planes, 1 % strain in
2 s in 41 steps (dt = 0.05 s).

## Checks (tolerance 1e-5)

The reference is the material-point integration of the same backward-Euler
scheme, solved with `scipy.optimize.root` on strain and plastic strain together
(`uniaxial_reference(..., implicit=True)`); it does not share the local Newton
solve or the tangent with the finite-element law.

| Metric | Reference | Result |
|---|---|---|
| `sigma_zz_curve` | material-point curve, maximum deviation over the 41 steps | 2.3e-12 |
| `cumulative_slip_final` | material-point cumulative slip | 1.1e-12 |
| `other_stress_components` | zero | 4.3e-16 |
| `saturation_gap` | tracked, not checked | 2.1e-6 |

The final stress is 533.453 MPa against an exact saturation stress of
533.454 MPa. The explicit scheme of `../crystal_fcc_001` stops at 475.8 MPa
with the same time step.

![Stress-strain curve](output/stress_strain_curve.png)

Verified against the blessed gold with Z3ST 0.4.1, FEniCSx 0.11.0 and
dolfinx-external-operator 0.10.1.
