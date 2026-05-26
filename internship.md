# Internship

## Overview

| Field                       | Value                                                                                            |
|-----------------------------|--------------------------------------------------------------------------------------------------|
| **Intern**                  | Baptiste Touron                                                                                  |
| **Supervisor**              | Giovanni Zullo                                                                                   |
| **Tool**                    | Z3ST (https://github.com/giozu/z3st)                                                             |
| **Status**                  | Planned                                                                                          |
| **Period**                  | 2026-05-26 → 2026-09-11 (~16 weeks)                                                              |
| **Project title**           | Phase-field fracture modelling and analysis of fission gas bubbles in UO₂                        |
| **Fork of the repository**  | `feature/PFF` off `develop`                                                      |

## Abstract

The internship focuses on phase-field modelling of fracture in UO₂ nuclear
fuel driven by fission gas bubble populations. Using Z3ST (FEniCSx-based,
AT1/AT2 already available), the student will apply the existing phase-field
formulation to representative pellet geometries containing explicit pore /
bubble distributions, to investigate crack initiation and propagation from
bubble-induced stress concentrations under thermo-mechanical loading.

## Plan

The work consists of a review and calibration phase on existing Z3ST
phase-field cases, followed by a new case applying the framework to UO₂
with explicit porosity, plus a supporting material card update.

### 1. Review

- Axisymmetric UO₂ pellet under thermal shock (extends the existing
`14_full_cylinder_cracking_2D_xy` case)

$$
\sigma_c \;=\; \sqrt{\dfrac{27\,E\,G_c}{256\,\ell}} \;\;\text{(AT2)}, \qquad
\sigma_c \;=\; \sqrt{\dfrac{3\,E\,G_c}{8\,\ell}} \;\;\text{(AT1)}
$$

- Single-edge notched tension/shear benchmarks (cases
  `19_single-edge_notched_*`).

Acceptance criteria:

- damage activates at physically meaningful loadings on the thermal-shock
  case;
- $(G_c, \ell)$ pair yields $\sigma_c$ in the
  $[100, 200]\,\text{MPa}$ range;
- both AT1 and AT2 paths reproduced and benchmarked.

### 2. Bubble-driven fracture case — `z3st/cases/UO2_PFF_bubbles_2D_rz/`

Axisymmetric UO₂ pellet geometry containing $N_b$ explicit circular holes
representing porosity / fission gas bubbles, meshed directly in Gmsh.
Bubble loading prescribed either as internal gas pressure (Neumann boundary
condition on the hole surfaces) or as a far-field thermo-mechanical load
(thermal gradient and/or PCI-induced hoop loading, interfacing with the
contact infrastructure developed in Romain's `feature/viscoplasticity`
branch, when merged).

The student investigates:

- crack initiation site as a function of bubble spacing, size distribution,
  and applied loading;
- transition from isolated bubble cracking to bubble-coalescence-driven
  fracture;
- comparison with refs. [1] and [3].

Optional stretch goal — stress-based phase-field variant (ref. [2]) added
to `damage_model.py` alongside the existing AT1/AT2 routes, controlled by
`damage.type: stress_based` in `input.yaml`. This would address the
AT2-only sensitivity to $\sigma_c$ via $(G_c, \ell)$ that motivated case 1.

Acceptance criteria:

- crack initiation pattern reproducible across mesh refinements;
- qualitative agreement with bubble-driven fracture morphologies in refs.
  [1, 3];
- (stretch) stress-based variant validated on the same SENT benchmark as
  AT2 in case 1.

### 3. Material cards — update

- `z3st/materials/uo2.yaml`: review and recalibrate phase-field entries
  (`Gc`, `sigma_c`, `lc`) based on the calibration outcome of the Review.

## Key references

1. Three-dimensional phase-field modeling of porosity-dependent intergranular fracture in UO₂ — https://www.sciencedirect.com/science/article/pii/S0927025619305683
2. Damage mechanics challenge: Predictions from an adaptive finite element implementation of the stress-based phase-field fracture model — https://www.sciencedirect.com/science/article/pii/S0013794424004156
3. Aagesen et al., *Phase-field simulations of fission gas bubbles in high burnup UO₂ during steady-state and LOCA transient conditions* — https://www.sciencedirect.com/science/article/pii/S0022311522002094

## Constraints

- Backward-compatibility with existing Z3ST cases.