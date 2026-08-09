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
on the simulated `output/history.csv`.

Author of the case and of the analysis scripts: Romain Turgis (ENSTA Paris),
internship 2026-05-25 → 2026-07-31, supervised by G. Zullo.

## Provenance

Adapted from `regression/pwr_rod_2D` into a different physical problem — a
zero-to-10 um gap, a traction-free clad outer surface, isothermal, over a
2500-day transient — and kept as a case of its own so that `pwr_rod_2D` stays a
PWR rod with its own gold.

## Parameter handling

`case_params.py` is the single source of truth for the post-processing: it
reads `geometry.yaml`, `input.yaml`, `clad.yaml` and `fuel.yaml` and exposes
the resolved constants to `plots.py` and `non-regression.py`. No
geometric or material value is restated as a literal in a script — an
analytical overlay computed with parameters that differ from the ones the
solver ran is not a verification.

`case_params.check_consistency()` fails loudly on the traps that invalidate
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

## Result: the case reproduces Esposito to within the paper's own Tresca bias

From a clean serial `./Allrun` on the committed configuration (non-fissile
pellet, cracking off, traction-free clad OD). Only the t = 0 row is asserted by
`non-regression.py`; the rest of the table is transcribed by hand and is not
re-checked by any run.

| t [days] | Z3ST [MPa] | eq. (21) [MPa] | ratio |
|---|---|---|---|
| 0 | 24.70 | 24.13 | 1.024 |
| 500 | 8.31 | 5.33 | 1.560 |
| 1240 | 5.33 | 3.40 | 1.569 |
| 2500 | 3.74 | 2.39 | 1.565 |

The elastic starting point agrees to **2.4%**, which verifies the elastic
factor `f`, the geometry and the material data. Both curves then relax as
power laws and the ratio settles at ~1.57, flat from day 500 onward.

That residual is **expected and quantified**. Z3ST's creep is J2 (von Mises);
eq. (11) of the paper adopts Tresca for the hub, and the paper's own
Discussion says the formulation is "marginally conservative (overestimate the
Pk decrease) compared to fem results ... due to the adoption of Tresca's
equivalent criterion ... which slightly overestimates the equivalent von Mises
stress adopted by the FEM solver". Tresca exceeds von Mises by 2/√3 on the
equivalent stress here, and the creep rate goes as σ_eq^n, so the predicted
bias is (2/√3)^n = **1.540** at n = 3 against **1.565** observed at 2500 days —
agreement within 2% on the bias itself.

So the disagreement is not a defect: it is a known criterion difference,
reproduced at the right magnitude. A gold blessed on this case should gate on
that, not on a naive overlay.

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
anything; `case_params.check_consistency()` fails loudly on the first three.

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

`non-regression.py` gates on the **elastic starting point**,
`contact_pressure_elastic_MPa`: at t = 0 no creep has accumulated, so
Δu_el = Δ and eq. (4) gives Pk = f·Δ. That point is purely elastic Lamé, so it
carries none of the Tresca/von Mises bias, and it is asserted against a 5 %
band (measured deviation +2.4 %, set by the finite penalty stiffness).

The relaxation *trend* is still gold-only. Gating on the full Pk(t) would mean
asserting the (2/√3)^n Tresca bias itself, which is a claim about the two
equivalent-stress criteria rather than about this solver.

## Run

    conda activate z3st
    ./Allrun      # gmsh → z3st → non-regression → plots
    ./Allclean
