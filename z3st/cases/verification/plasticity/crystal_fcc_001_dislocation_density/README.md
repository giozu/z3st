# verification/plasticity/crystal_fcc_001_dislocation_density: crystal plasticity with dislocation-density hardening

The problem of `../crystal_fcc_001_implicit` (FCC single crystal, 12 {111}<110>
systems, uniaxial tension along [001], implicit law through an external
operator) with a dislocation density on every slip system.

In the local suite (carries a gold). Not in `cases_ci.txt`: it takes about
90 s and needs `dolfinx-external-operator`.

## Model

```
dgamma_a/dt = gamma0 (tau_a/g_a) |tau_a/g_a|^(n-1)
F_a         = sum_b A_ab rho_b,      A_ab = 1 (a = b), a_latent (a != b)
g_a         = tau0 + alpha_taylor G burgers sqrt(F_a)        (Taylor)
drho_a/dt   = (k1 sqrt(F_a) - k2 rho_a) |dgamma_a/dt|        (Kocks-Mecking)
```

Backward Euler on the plastic strain and on the 12 densities together (21
unknowns per quadrature point, Newton in numpy, consistent tangent from the
Jacobian of the local problem). When Newton from the last converged state
fails, because the trial stress is far above the slip resistance, the solution
for half the time step is used as the guess (continuation in dt).

Card (`fcc_crystal.yaml`): tau0 = 20 MPa, alpha_taylor = 0.3,
burgers = 0.25 nm, rho0 = 1e11 1/m2 per system, k1 = 1.3e8 1/m, k2 = 10,
a_latent = 1; E = 200 GPa, nu = 0.3, gamma0 = 0.001 1/s, n = 5. The values are
chosen for the verification, not fitted to a material.

## Checks (tolerance 1e-5)

The reference is the material-point integration of the same scheme with
`scipy.optimize.root` on strain, plastic strain and densities
(`uniaxial_reference(..., implicit=True)`); it writes the hardening law again
and does not use the local solve of the finite-element law.

| Metric | Reference | Result |
|---|---|---|
| `sigma_zz_curve` | material-point curve, maximum deviation over the 41 steps | 5.7e-13 |
| `cumulative_slip_final` | material-point cumulative slip | 4.1e-13 |
| `other_stress_components` | zero | 5.5e-15 |
| `saturation_gap` | tracked, not checked | 94.4 % |

The initial slip resistance is 26.3 MPa. The final stress is 93.2 MPa; with
the resistance frozen at its initial value the stress would saturate at
70.2 MPa. The steady state of the densities (k1 sqrt(F) = k2 rho, stress
1653.8 MPa) needs a strain far beyond 1 % and beyond small strain, so
`saturation_gap` is only tracked; the self-check of the module
(`python materials/crystal_plasticity.py`) verifies that steady state on the
material-point integration.

![Stress-strain curve](output/stress_strain_curve.png)

Verified against the blessed gold with Z3ST 0.4.1, FEniCSx 0.11.0 and
dolfinx-external-operator 0.10.1.
