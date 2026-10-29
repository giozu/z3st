# Thermoelastic cylinder: deviatoric/hydrostatic split of von Mises and Tresca

## Overview

Open-ended thick cylinder under internal pressure and a radial temperature
gradient, axisymmetric. The total stress is checked against the analytic
solution, then split into its hydrostatic and deviatoric parts, and von Mises
and Tresca are computed on each part.

In the local suite (carries a gold).

- Lamé (open end) superposed on the free-end thermal stresses
- von Mises and Tresca on the total, deviatoric and hydrostatic stress
- Tresca from the principal stresses (`Tresca` field of the writer)

## Physical model

### Split of the stress tensor

$$
\boldsymbol\sigma = \mathbf{s} + p\mathbf{I},
\qquad
p = \frac{\mathrm{tr}\ \boldsymbol\sigma}{3},
\qquad
\mathrm{tr}\ \mathbf{s} = 0
$$

The hydrostatic part $p\mathbf{I}$ changes the volume, the deviator
$\mathbf{s}$ changes the shape.

### Equivalent stresses

$$
\sigma_{\mathrm{vM}}(\boldsymbol\sigma) = \sqrt{\frac{3}{2}\ \mathbf{s}:\mathbf{s}}
$$

$$
\sigma_{\mathrm{T}}(\boldsymbol\sigma) = \max_{i,j} |\sigma_i - \sigma_j| = \sigma_1 - \sigma_3,
\qquad \sigma_1 \ge \sigma_2 \ge \sigma_3
$$

with $\sigma_1, \sigma_2, \sigma_3$ the principal stresses, i.e. the
eigenvalues of $\boldsymbol\sigma$. Both criteria depend only on the
differences of the principal stresses, so neither sees the hydrostatic part.
Adding $p\mathbf{I}$ shifts every principal stress by $p$ and leaves the
differences unchanged:

$$
\sigma_{\mathrm{vM}}(\mathbf{s}) = \sigma_{\mathrm{vM}}(\boldsymbol\sigma),
\qquad
\sigma_{\mathrm{T}}(\mathbf{s}) = \sigma_{\mathrm{T}}(\boldsymbol\sigma),
\qquad
\sigma_{\mathrm{vM}}(p\mathbf{I}) = \sigma_{\mathrm{T}}(p\mathbf{I}) = 0
$$

The two criteria are bounded by each other:

$$
\frac{\sqrt{3}}{2}\ \sigma_{\mathrm{T}} \le \sigma_{\mathrm{vM}} \le \sigma_{\mathrm{T}}
$$

with $\sigma_{\mathrm{vM}} = \sigma_{\mathrm{T}}$ when two principal stresses
coincide.

### Analytic solution

Temperature, imposed on both faces:

$$
T(r) = T_i + (T_o - T_i)\ \frac{\ln(r/R_i)}{\ln(R_o/R_i)}
$$

Lamé, open end (no axial load), with $P_o = 0$:

$$
\sigma_{rr} = A - \frac{B}{r^2},
\qquad
\sigma_{\theta\theta} = A + \frac{B}{r^2},
\qquad
\sigma_{zz} = 0,
\qquad
A = \frac{P_i R_i^2}{R_o^2 - R_i^2},
\qquad
B = \frac{P_i R_i^2 R_o^2}{R_o^2 - R_i^2}
$$

Thermal stresses, free end (zero axial force), any $T(r)$:

$$
\sigma_{rr} = \frac{\alpha E}{1-\nu}\ \frac{1}{r^2}
\left[ \frac{r^2 - R_i^2}{R_o^2 - R_i^2} \int_{R_i}^{R_o} T r\ dr - \int_{R_i}^{r} T r\ dr \right]
$$

$$
\sigma_{\theta\theta} = \frac{\alpha E}{1-\nu}
\left[ \frac{1 + R_i^2/r^2}{R_o^2 - R_i^2} \int_{R_i}^{R_o} T r\ dr + \frac{1}{r^2} \int_{R_i}^{r} T r\ dr - T(r) \right]
$$

$$
\sigma_{zz} = \frac{\alpha E}{1-\nu}
\left[ \frac{2}{R_o^2 - R_i^2} \int_{R_i}^{R_o} T r\ dr - T(r) \right]
$$

The reference is the sum of the two.

## Geometry and loading

| Parameter | Value | Description |
|-----------|-------|-------------|
| $R_i$, $R_o$ | 20 mm, 30 mm | inner and outer radius |
| $L_z$ | 500 mm | height |
| $P_i$ | 50 MPa | internal pressure |
| $T_i$, $T_o$ | 500 K, 400 K | inner and outer temperature |
| steel | E = 200 GPa, ν = 0.3, α = 1e-5 1/K | `materials/steel.yaml` |

- Axisymmetric, structured quadrilateral mesh, 40 × 80 cells
- Bottom: axial displacement blocked, radial free (`Clamp_y`)
- Top: free, so the cylinder is open-ended

## Checks

1. Total stress on the mid-plane against Lamé + thermal (relative L2 error):
   this is the only check that tests the FE solution.
2. The identities above, in every cell, to round-off ($10^{-12}$). They are
   algebraic, so they hold whatever the FE does: they test the tensor layout
   and the writer's `VonMises` and `Tresca` fields, not the physics.
3. The bounds between von Mises and Tresca, in every cell.

## Results

| Check | Value |
|-------|-------|
| $\sigma_{rr}$, $\sigma_{\theta\theta}$, $\sigma_{zz}$ vs analytic | 9.0e-4, 8.4e-5, 1.1e-4 |
| $\sigma_{\mathrm{vM}}$, $\sigma_{\mathrm{T}}$ vs analytic | 1.6e-4, 1.7e-4 |
| $\max \lvert\sigma_{\mathrm{vM}}(\mathbf{s}) - \sigma_{\mathrm{vM}}(\boldsymbol\sigma)\rvert$ / max | 1.5e-16 |
| $\max \lvert\sigma_{\mathrm{T}}(\mathbf{s}) - \sigma_{\mathrm{T}}(\boldsymbol\sigma)\rvert$ / max | 3.1e-16 |
| $\sigma_{\mathrm{vM}}(p\mathbf{I})$, $\sigma_{\mathrm{T}}(p\mathbf{I})$ | 0 |

What the plot `output/invariants.png` shows:

- The hydrostatic stress is large, from -80 MPa to +108 MPa through the wall,
  yet removing it changes neither criterion.
- Heating from the inside unloads the inner wall and loads the outer one: the
  thermal hoop stress is compressive inside and tensile outside, so the
  maximum moves from the inner face (pressure alone) to the outer face.
- Tresca has a kink at $r \approx 24$ mm, where $\sigma_{zz}$ crosses
  $\sigma_{rr}$ and the smallest principal stress changes. There two principal
  stresses coincide and $\sigma_{\mathrm{vM}} = \sigma_{\mathrm{T}}$.

Near the free top face the analytic solution no longer holds (end effect,
Saint-Venant): $\sigma_{rz}$ reaches 16 MPa and the
largest equivalent stresses of the domain, 195 MPa von Mises and 208 MPa
Tresca, sit there. With shear present, Tresca from the diagonal components
alone would be wrong by up to 12 MPa: that is why the writer takes the
eigenvalues.

## Run

```bash
./Allrun
```

See also `thermoelastic_cylinder_superposition_2D`, which splits the same
problem into its thermal and mechanical parts.
