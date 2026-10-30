# Thermoelastic cylinder: superposition of thermal and mechanical stresses

## Overview

The cylinder of `thermoelastic_cylinder_invariants_2D` (internal pressure plus
a radial temperature gradient) solved as a whole and as its thermal and
mechanical subproblems, on the same mesh. The stress tensors superpose; von
Mises and Tresca do not. The case then repeats the exercise with a hydrostatic
mechanical load, the one load for which the equivalent stresses do add up
whatever the thermal field.

In the local suite (carries a gold).

- Five runs on one mesh, one subfolder each
- Superposition of the tensors, to solver precision
- Triangle inequality for von Mises and Tresca, and when it becomes an equality

## Physical model

### Superposition of the tensors

With linear elasticity and temperature-independent properties, the problem is
linear in its loads, so the stress of the whole is the sum of the stresses of
the parts:

$$
\boldsymbol\sigma_{\mathrm{tot}} = \boldsymbol\sigma_{\mathrm{th}} + \boldsymbol\sigma_{\mathrm{mech}}
$$

$\boldsymbol\sigma_{\mathrm{th}}$ needs a solve of its own: the thermal stress depends on how the structure restrains the expansion, i.e. on compatibility. It is not the eigenstress $-(3\lambda + 2G) \alpha (T - T_{\mathrm{ref}}) \mathbf{I}$, which is hydrostatic and is the stress of a body that cannot expand at all.

### The equivalent stresses do not superpose

Write $f$ for either $\sigma_{\mathrm{vM}}$ or $\sigma_{\mathrm{T}}$. Both are seminorms of the stress tensor: $f \ge 0$, $f(k\boldsymbol\sigma) = |k| f(\boldsymbol\sigma)$, and

$$
f(\boldsymbol\sigma_{\mathrm{th}} + \boldsymbol\sigma_{\mathrm{mech}})
\le f(\boldsymbol\sigma_{\mathrm{th}}) + f(\boldsymbol\sigma_{\mathrm{mech}})
$$

from which also

$$
f(\boldsymbol\sigma_{\mathrm{tot}}) \ge
\left| f(\boldsymbol\sigma_{\mathrm{th}}) - f(\boldsymbol\sigma_{\mathrm{mech}}) \right|
$$

So $f(\boldsymbol\sigma_{\mathrm{tot}})$ lies in a band, and adding the equivalent stresses of the parts gives its upper edge: an estimate that is never unsafe, and can be far too high. The tensors are added first, the equivalent stress is taken after. ASME III NB-3200 does the same: the stress components of each category are summed ($P_m + P_b + Q$) and the stress intensity is computed on the sum.

### When the inequality becomes an equality

With a hydrostatic load, $\boldsymbol\sigma_{\mathrm{mech}} = -p\mathbf{I}$, one has $f(\boldsymbol\sigma_{\mathrm{mech}}) = 0$ and $f(\boldsymbol\sigma_{\mathrm{th}} - p\mathbf{I}) = f(\boldsymbol\sigma_{\mathrm{th}})$, so

$$
f(\boldsymbol\sigma_{\mathrm{tot}}) = f(\boldsymbol\sigma_{\mathrm{th}}) + f(\boldsymbol\sigma_{\mathrm{mech}})
$$

for any thermal field. It is the only load with this property for every $\boldsymbol\sigma_{\mathrm{th}}$: taking $\boldsymbol\sigma_{\mathrm{th}} = -\boldsymbol\sigma_{\mathrm{mech}}$ gives $0 = 2 f(\boldsymbol\sigma_{\mathrm{mech}})$, so the deviator of $\boldsymbol\sigma_{\mathrm{mech}}$ vanishes.

It is not the only case of equality, though:

- von Mises is, up to a factor, the norm of the deviator, so equality holds whenever the two deviators are proportional, $\mathbf s_{\mathrm{mech}} = k \mathbf s_{\mathrm{th}}$ with $k \ge 0$; the hydrostatic load is the limit $k = 0$.
- Tresca is more permissive: equality holds wherever the two states share their principal directions and the ordering $\sigma_1 \ge \sigma_2 \ge \sigma_3$, because then $\sigma_1 - \sigma_3$ is linear.

## Subproblems

All five share `mesh.geo` and `geometry.yaml` of the case root.

| Subfolder | Mechanical load | Temperature |
|-----------|-----------------|-------------|
| `total/` | $P_i$ = 50 MPa on the inner face | 500/400 K |
| `thermal/` | none | 500/400 K |
| `mechanical/` | $P_i$ = 50 MPa on the inner face | $T_{\mathrm{ref}}$ = 300 K |
| `total_hydrostatic/` | $p$ = 50 MPa on the inner, outer and top faces | 500/400 K |
| `mechanical_hydrostatic/` | $p$ = 50 MPa on the inner, outer and top faces | $T_{\mathrm{ref}}$ = 300 K |

At $T = T_{\mathrm{ref}}$ there is no thermal strain. The pressure on every face, with the bottom free to slide radially, gives $\boldsymbol\sigma = -p\mathbf{I}$ everywhere. The analytic references are those of `thermoelastic_cylinder_invariants_2D`.

## Checks

1. Each subproblem on the mid-plane against its analytic solution.
2. Superposition of the tensors in every cell, for both mechanical loads.
3. The triangle inequality in every cell, for both criteria and both loads.
4. The deviatoric/hydrostatic identities on each subproblem.
5. The equality under the hydrostatic load; the fraction of cells where
   Tresca adds up exactly under the internal pressure is tracked.

## Results

| Check | Value |
|-------|-------|
| subproblems vs analytic (relative L2) | 6e-6 to 2.9e-3 |
| $\max \lvert\boldsymbol\sigma_{\mathrm{tot}} - \boldsymbol\sigma_{\mathrm{th}} - \boldsymbol\sigma_{\mathrm{mech}}\rvert$ / max | 6.7e-13 (pressure), 5.8e-13 (hydrostatic) |
| $\max \lvert\boldsymbol\sigma_{\mathrm{mech}} + p\mathbf{I}\rvert / p$, hydrostatic load | 9.3e-9 |
| internal pressure: $f(\mathrm{th}) + f(\mathrm{mech}) - f(\mathrm{tot})$, up to | 116% (vM), 123% (Tresca) of max $f(\mathrm{tot})$ |
| hydrostatic load: same quantity | 6.5e-12 (vM), 4.3e-12 (Tresca) |
| internal pressure: cells where Tresca adds up exactly | 44.5% |

The plot `output/superposition.png` has von Mises on the left, Tresca on the
right; the internal pressure on the top row, the hydrostatic load on the
bottom row. Markers are FE, lines analytic. Black is the total, orange the
thermal part, blue the mechanical part; the dashed and dotted green lines are
the upper and lower edges of the band.

- On the top row, at the inner face, thermal and mechanical von Mises are both about 158 MPa, their sum 316 MPa, the total only 120 MPa: the thermal hoop stress is compressive (-157 MPa, the hot inside is restrained), the pressure hoop stress tensile (+130 MPa), and they cancel in the tensor.
- On the top row, Tresca of the total coincides with the sum for $r > 24.5$ mm. Both states have $\sigma_1 = \sigma_{\theta\theta}$ and $\sigma_3 = \sigma_{rr}$, so Tresca is $\sigma_{\theta\theta} - \sigma_{rr}$ in both and adds linearly: on the outer face 122 + 80 = 202 MPa. Inside, the thermal state has $\sigma_1 = \sigma_{rr}$ and the mechanical one $\sigma_1 = \sigma_{\theta\theta}$: the ordering differs and the sum overestimates, 337 MPa against 130 MPa.
- On the bottom row the mechanical part is zero for both criteria, and total, thermal and both edges of the band collapse onto one curve.

## Run

```bash
./Allrun
```

`Allrun` meshes once, runs z3st in each subfolder, then runs
`non-regression.py` on the five results.
