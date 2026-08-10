# creep_shrink_fit_2D — contact-pressure relaxation by Norton creep

2D axisymmetric (r,z) UO₂ pellet inside a Zircaloy-4 cladding. The pellet
expands thermally, closes a 10 μm cold radial gap, and loads the cladding
through the penalty contact model. The cladding then creeps under sustained
load and the contact pressure **relaxes** in time.

Reference for the closed form: Esposito, Bruno, Bertocco, *Int. J. Pressure
Vessels and Piping* **185** (2020) 104126 — eq. (21), shrink-fit assembly
relaxing by Norton creep:

    Pk(t) = Δu_el(t) · f,   Δu_el(t) = Δ · (1 + φ(1-n)·t / Δ^(1-n))^(1/(1-n))
    φ = -(A/b) · f^n · (k1 + k2)

`plots.py::plot_contact_pressure_evolution` evaluates this law and overlays it
on the simulated `output/history.csv`. Δ is the interference of the assembled
joint at temperature, held fixed across the transient — the law relaxes it
itself — and `f` is corrected for the penalty spring in series with the joint,
`1/(1/f + 1/K_PEN)`. Both matter: see "Result" and "Convergence" below.

Author of the case and of the analysis scripts: Romain Turgis (ENSTA Paris),
internship 2026-05-25 → 2026-07-31, supervised by G. Zullo.

## Parameter handling

`case_params.py` is the single source of truth for the post-processing: it
reads `geometry.yaml`, `input.yaml`, `clad.yaml` and `fuel.yaml` and exposes
the resolved constants to `plots.py` and `non-regression.py`. No
geometric or material value is restated as a literal in a script — an
analytical overlay computed with parameters that differ from the ones the
solver ran is not a verification.

`case_params.check_consistency()`, called from `diagnostics.py` before step 0,
aborts the run on the traps that invalidate
the comparison silently: a cold gap from `geometry.yaml` whose magnitude
disagrees with `models.contact.initial_gap`, an `initial_gap` that is a
clearance rather than an interference, and irradiation creep left switched on.
`plots.py` calls it on every run and warns before plotting.

The gap check compares magnitudes, not signed values. A conforming mesh has to
keep the pellet and the cladding disjoint, so `geometry.yaml` can only carry a
positive clearance, while the interference is expressed by a negative
`initial_gap` that the contact model uses in place of the mesh-derived value.
The two therefore differ in sign by construction, and only their magnitudes are
required to agree.

## Result: the elastic point verifies, the relaxation trend does not yet

From a clean serial `./Allrun` on the committed configuration (non-fissile
pellet, cracking off, traction-free clad OD). Only the t = 0 row is asserted by
`non-regression.py`; the rest of the table is transcribed by hand and is not
re-checked by any run.

| t [days] | Z3ST [MPa] | eq. (21) [MPa] | ratio |
|---|---|---|---|
| 0 | 24.70 | 25.47 | 0.970 |
| 600 | 7.26 | 8.54 | 0.850 |
| 1240 | 5.14 | 6.12 | 0.840 |
| 2500 | 3.65 | 4.37 | 0.834 |

The elastic starting point agrees to **3.0%** against a closed form that shares
nothing with the solver output — see "Elastic point at t = 0".

The relaxation trend does not. The eq. (21) column above holds Δ at the
assembled interference, 14.73 μm, which is what the law assumes.
`plots.py::plot_contact_pressure_evolution` instead passes `Δ = -gap_um` **per
step** (`plots.py:432`) into a law that relaxes Δ itself: the interference
handed to eq. (21) has already relaxed in the simulation, so the relaxation is
applied twice.

That changes the reading. The Tresca argument below predicts **1.540**: the
paper's formulation overestimates the Pk decrease, so the closed form should
sit *below* Z3ST and the ratio should exceed 1. At 0.83 Z3ST relaxes ~17 %
*faster* than the closed form (0.834 at 2500 days), which is the opposite
direction, and the
criterion difference cannot account for it.

Checked two independent ways, because `K_PEN = 3.5e12` sits within 3 % of the
joint stiffness `f = 3.42e12` and the penalty contact therefore carries half
the compliance of the problem:

- soft penalty, closed form corrected for the two springs in series
  (`dPk/dt = -f_ser·(A/b)·(k2-k1)·Pk^n`, `f_ser = 1.73e12`): ratio **0.855** at
  2500 days;
- `K_PEN` raised to 3.5e14, penalty compliance negligible (`f_ser/f = 0.990`,
  0.13 μm of penetration against 13.75 μm across the joint), closed form
  applied directly: ratio **0.875** at 2500 days. That run converged in 150 s,
  so the soft penalty is not a robustness requirement.

Both routes agree to a couple of percent, so the ratio is not an artefact of
how the penalty spring is handled, and the convergence study below rules out
the discretisation. Open, and the reason this case does not assert the trend.

### Convergence

The ratio is not a discretisation artefact. Three independent refinements, each
against the committed configuration:

| refinement | p(2500 d) [MPa] | ratio |
|---|---|---|
| committed (`n_r1/n_r2 = 25/7`, P1, `n_steps 56/80`) | 3.648 | 0.834 |
| mesh `n_r1/n_r2 = 49/21` | 3.731* | 0.853* |
| `mechanical.order: 2` (P2 displacement) | 3.730* | 0.853* |
| `n_steps 112/160` | 3.633 | 0.831 |

\* measured at `n_steps 14/20`, against a 3.738 MPa baseline at the same
stepping: mesh and element order move the pressure by 0.2 %.

Space is converged: doubling the radial divisions and raising the displacement
field to P2 land on the same 0.853 from two different directions.

Time is the only knob that moves anything, and it moves *away* from 1.540. The
original `n_steps 14/20` carried an 8.4 % error at 200 days, where the pressure
halves in two hundred days, against a 112/160 reference; the case now runs
56/80, whose residual is 1.3 % at 200 days and 0.4 % at 2500. That refinement
alone took the ratio from 0.855 to 0.834.

`Mesh.ElementOrder = 2` is not usable here: the output writer requires the
output functions to match the mesh degree (`RuntimeError: Degree of output
Function must be same as mesh degree`), which is the limitation noted at
`utils/writer.py:358`. `mechanical.order` (`core/finite_element_setup.py:25`)
raises the displacement field instead and needs no quadratic mesh.

The Tresca argument itself, for the record. Z3ST's creep is J2 (von Mises);
eq. (11) of the paper adopts Tresca for the hub, and the paper's own
Discussion says the formulation is "marginally conservative (overestimate the
Pk decrease) compared to fem results ... due to the adoption of Tresca's
equivalent criterion ... which slightly overestimates the equivalent von Mises
stress adopted by the FEM solver". Tresca exceeds von Mises by 2/√3 on the
equivalent stress here, and the creep rate goes as σ_eq^n, so the predicted
bias is (2/√3)^n = **1.540** at n = 3. It remains the expected magnitude and
sign; nothing in the current numbers reproduces it.

Also open, and a candidate: with the stiff penalty the elastic point falls to
**-5.7 %** (47.00 MPa against 49.86), against -3.0 % with the soft one, and
0.85 μm of the 14.73 μm interference is unaccounted for by the joint-plus-
penalty split. An elastic deficit that grows as the contact stiffens points at
the contact model rather than at the creep law.

### What it took to get here

Three configuration errors had to be removed, each producing a different and
individually plausible wrong answer:

| Configuration | Behaviour | Cause |
|---|---|---|
| fissile pellet, clamped clad OD | P_c **rises** to 160 MPa | burnup swelling grows Δ; eq. (21) assumes Δ fixed |
| non-fissile, clamped clad OD | P_c **freezes** at 12.9 MPa | u_r = 0 at r = c lets creep only redistribute stress until the deviatoric part vanishes, then it stops |
| non-fissile, free clad OD, `initial_gap: +10 µm` | **no contact at all** | the clamp was what closed the gap; nothing else does |
| non-fissile, free clad OD, `initial_gap: -10 µm` | relaxes correctly | a negative gap is a true interference, contact active from t = 0 |

The last row is the paper's shrink fit, and is what the case now carries.

## Configuration constraints

The case carries `output/non-regression_gold.json` and is a member of the local
suite. Five settings are required for the comparison with eq. (21) to mean
anything; `case_params.check_consistency()`, called from `diagnostics.py` before step 0,
aborts the run on the first three.

1. **Irradiation creep off.** `creep_irr_B` and `fast_flux` are commented out in
   `clad.yaml`. Esposito eq. (21) is a pure Norton derivation with no in-pile
   term.

2. **Every constant read from the YAML**, through `case_params.py`. An
   analytical overlay computed with parameters that differ from the ones the
   solver ran is not a verification, and the two drift apart silently.

3. **`creep_Q: 0` — the case is isothermal.** Esposito eq. (9) is
   `eps_eq_dot = A*sigma_eq^n` with `A` a constant `[h^-1]` (their
   nomenclature): there is no Arrhenius factor in the derivation, and
   temperature enters only through the interference via `alpha*dT`.

   The trap: `creep_Q = 0` makes `A(T) = A0` identically, so **`A0` must be the
   Arrhenius factor already folded in**. The Zircaloy card value
   `2.82e-24 Pa^-n s^-1` is the prefactor that belongs with `creep_Q = 1.2e5`;
   used with `creep_Q = 0` it runs the cladding ~1e10 times too fast, and
   nothing in the output says so. `clad.yaml` therefore carries the re-based
   `A0 = 2.82e-24 * exp(-1.2e5/(R*600 K)) = 1.0083e-34 Pa^-n s^-1`, the same
   `A(T)` that `verification/fuel/creep` reports at 600 K.

4. **`fissile: false`, and cracking off.** Eq. (21) is integrated under a
   *fixed* interference assembled once and left to relax. A fissile pellet grows
   the interference through the swelling eigenstrain and the contact pressure
   rises instead of relaxing (see the table above). Zero burnup switches off both
   the swelling and the densification terms, so the interference changes only by
   creep; the eigenstrain entry stays in the card, inert. Cracking would degrade
   the pellet to `E_iso/E = 0.112` while the analytical overlay uses the card's
   nominal `E = 200 GPa`, and it fires on the nominal LHR even with
   `fissile: false`, i.e. on power that is never deposited. Esposito's shaft is
   an elastic solid.

5. **Clad outer surface traction-free.** `outer_2` carries no restraint in
   `boundary_conditions.yaml`, which is Esposito's hub boundary condition
   `sigma_r(c) = 0`. Under `u_r = 0` at `r = c` creep can only redistribute
   stress until the deviatoric part vanishes and then stops, freezing the
   contact pressure at a plateau instead of decaying as a power law.

## Open point

`A0` is re-based at a round 600 K stand-in for the clad mean temperature. A hand
estimate at LHR = 20 kW/m (`h_conv = 3.5e4`, `T_ext = 580 K`) puts the clad
between ~597 K outer and ~628 K inner, so ~610 K mean. `A` enters phi linearly,
so a 10 K error is a few percent on the relaxation rate: re-derive `A0` from the
computed clad mean temperature and re-bless the gold.

### Elastic point at t = 0

`contact_pressure_elastic_MPa` is asserted against eq. (4) with the joint and
the penalty contact treated as springs in series, within a 5 % band:

    p = Delta / (1/f + 1/K_PEN)     25.47 MPa predicted, 24.70 MPa simulated, +3.0 %

Nothing on the right-hand side is read back from the solver: `f` comes from the
geometry and the moduli, `K_PEN` from `input.yaml`, `Delta` from the material
cards and the run temperature. The measured split at t = 0 is 7.23 um across
the joint and 7.06 um of penalty penetration, which sum to `Delta`.

`Delta` is the interference of the assembled joint **at temperature**, not the
cold 10 um of `models.contact.initial_gap`. At 580 K the pellet outgrows the
cladding by 4.73 um (`r*alpha*(T - T_ref)` on each body, exact here because the
field is uniform), so the assembled interference is 14.73 um. Using the cold
value instead predicts 17.29 MPa and reads as a 43 % error — a variant of the
same trap as the table above.

The case tolerance is 5e-2, set by this
metric; the burnup assertions are unaffected, they compare against exactly zero.

The relaxation *trend* is likewise gold-only. Gating on the full Pk(t) would mean
asserting the (2/√3)^n Tresca bias itself, which is a claim about the two
equivalent-stress criteria rather than about this solver.

## Run

    conda activate z3st
    ./Allrun      # gmsh → z3st → non-regression → plots
    ./Allclean
