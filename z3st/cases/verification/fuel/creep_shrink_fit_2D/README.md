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

`pressure_creep_analysis.py` evaluates this law and overlays it on the
simulated `output/history.csv`.

Author of the case and of the analysis scripts: Romain Turgis (ENSTA Paris),
internship 2026-05-25 → 2026-07-31, supervised by G. Zullo.

## Provenance

The case was developed on the branch `merging/contact` **on top of**
`regression/pwr_rod_2D`, which it progressively mutated into a different
physical problem (zero-to-10 μm gap, radially clamped clad outer surface,
isothermal variants, 2500-day transient). It has been re-planted here as a
case of its own so that `regression/pwr_rod_2D` stays a PWR rod and keeps its
blessed gold. Snapshot taken from commit `15b4b4e` ("free boundary
conditions", 2026-07-16), the most advanced committed state.

## Parameter handling

`case_params.py` is the single source of truth for the post-processing: it
reads `geometry.yaml`, `input.yaml`, `clad.yaml` and `fuel.yaml` and exposes
the resolved constants to `plots.py` and `pressure_creep_analysis.py`. No
geometric or material value is restated as a literal in a script — an
analytical overlay computed with parameters that differ from the ones the
solver ran is not a verification.

`case_params.check_consistency()` fails loudly on the two traps that
invalidate the comparison silently: a cold gap from `geometry.yaml` that
disagrees with `models.contact.initial_gap`, and irradiation creep left
switched on. `pressure_creep_analysis.py` calls it on every run and warns
before plotting.

## Result (2026-07-28): the case reproduces Esposito to within the paper's own Tresca bias

| t [days] | Z3ST [MPa] | eq. (21) [MPa] | ratio |
|---|---|---|---|
| 0 | 24.70 | 24.13 | 1.024 |
| 500 | 7.30 | 5.14 | 1.421 |
| 1240 | 5.20 | 3.33 | 1.559 |
| 2500 | 3.88 | 2.45 | 1.586 |

The elastic starting point agrees to **2.4%**, which verifies the elastic
factor `f`, the geometry and the material data. Both curves then relax as
power laws and the ratio settles at ~1.58.

That residual is **expected and quantified**. Z3ST's creep is J2 (von Mises);
eq. (11) of the paper adopts Tresca for the hub, and the paper's own
Discussion says the formulation is "marginally conservative (overestimate the
Pk decrease) compared to fem results ... due to the adoption of Tresca's
equivalent criterion ... which slightly overestimates the equivalent von Mises
stress adopted by the FEM solver". Tresca exceeds von Mises by 2/√3 on the
equivalent stress here, and the creep rate goes as σ_eq^n, so the predicted
bias is (2/√3)^n = **1.540** at n = 3 against **1.586** observed — agreement
within 3% on the bias itself.

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

## Remaining open points

There is no `non-regression_gold.json` in this directory on purpose.
**The figure Romain reported on 2026-07-20 (excellent agreement with the
analytical law out to 2500 days) is not reproducible from any commit** — it
was produced with hand-edited values that were never committed. Two of the
four causes have been fixed here; two remain open.

1. ~~**Irradiation creep left active.**~~ Fixed: `creep_irr_B` and
   `fast_flux` are commented out in `clad.yaml`, and `check_consistency()`
   now refuses to let them come back unnoticed. Esposito's eq. (21) is a pure
   Norton derivation with no in-pile term.

2. ~~**Constants hardcoded in the analysis scripts.**~~ Fixed: they are read
   from the YAML via `case_params.py`. For the record, the drift that was
   there — `R_CLAD_I` 0.00453 vs 0.00451, `K_PEN` 5.0e13 vs 3.5e12,
   `T_TOTAL` 1800 d vs 2500 d, `NB_STEPS` 45/50 vs 34 — and the reported
   figure is titled `k_pen=3,47e12`, matching none of them.

3. ~~**`creep_Q`.**~~ Resolved: the case is **isothermal**, `creep_Q: 0`.
   Esposito's eq. (9) is `ε̇_eq = A·σ_eq^n` with `A` a constant `[h⁻¹]` (see
   the paper's nomenclature) — there is no Arrhenius factor anywhere in the
   derivation, and temperature enters only through the interference Δ via
   `α·ΔT`. This is the path the `models/creep_model.py` fix unblocks.

   The trap this opens: `creep_Q = 0` makes `A(T) = A0` identically, so `A0`
   must be the Arrhenius factor *already folded in*. The Zircaloy card value
   `2.82e-24 Pa⁻ⁿ s⁻¹` is the prefactor that belongs with `creep_Q = 1.2e5`;
   used with `creep_Q = 0` it runs the cladding ~10¹⁰ times too fast, and
   nothing in the output says so. `clad.yaml` therefore carries the re-based
   value `A0 = 2.82e-24 · exp(-1.2e5/(R·600 K)) = 1.0083e-34 Pa⁻ⁿ s⁻¹`, which
   is the same `A(T)` that `verification/fuel/creep` reports at 600 K.
   `case_params.check_consistency()` now flags an un-folded prefactor.

   Still to confirm from a run: 600 K is a round stand-in for the clad mean
   temperature. A hand estimate at LHR = 20 kW/m (`h_conv = 3.5e4`,
   `T_ext = 580 K`) puts the clad between ~597 K outer and ~628 K inner, so
   ~610 K mean. Re-derive `A0` from the computed clad mean temperature before
   blessing a gold; `A` enters φ linearly, so a 10 K error is a few percent
   on the relaxation rate.

4. **OPEN — the interference is not fixed: burnup swelling drives it up.**
   This is the big one, confirmed by running the case (2026-07-28). The
   simulated contact pressure does **not** relax onto Esposito's curve — it
   dips to ~9 MPa around day 300 and then climbs monotonically to **160 MPa
   at 2500 days**, while eq. (21) decays to ~3 MPa. Total divergence.

   The cause is not the creep model. `fuel.yaml` has `fissile: true` and
   `eigenstrain: materials.fuel_swelling.solid_gas_densification` with
   `swelling_rate: 7.0e-4` per MWd/kgU, and the run reaches **70.8 MWd/kgU**
   at 20 kW/m over 2500 days. The pellet keeps growing, so the interference
   keeps growing, so `P_k` rises. Esposito eq. (21) is integrated under a
   *fixed* interference Δ assembled once and left to relax — the two loadings
   are not the same problem. This is exactly the mismatch flagged in the
   supervisor's 2026-07-10 e-mail ("the interference is not fixed and your
   P_k rises"), and the committed snapshot never resolved it.

   To actually reproduce the paper the burnup-driven interference growth has
   to be removed — `fissile: false`, or dropping the swelling eigenstrain, so
   that Δ is set once and only creep changes it. That is a physics decision
   about what the case is for, so it is left open rather than patched in.

5. **OPEN — modelling choice to document.** `boundary_conditions.yaml`
   applies `Clamp_x` on `outer_2`, i.e. the cladding outer surface is
   radially restrained rather than loaded by coolant pressure. That is what
   makes the configuration behave like Esposito's rigid hub; it is a
   deliberate idealisation, not a PWR boundary condition, and the report
   should say so.

Once 3 is settled and the case runs clean, `non-regression.py` should gate on
the analytical Pk(t) (a real verification), not only on a gold snapshot.

## What was deliberately NOT carried over

- The `core/solver.py` edit from `15b4b4e`. It re-computes the `_mech_cache`
  rebuild condition ~20 lines above the existing block in
  `_mechanical_step`, which then overwrites it — so the added
  `_force_creep_rebuild` term never reaches the assembly decision and the
  change is a no-op. It also moves the `nonlinear_constitutive` assignment
  into an `else` branch, leaving it conditionally unbound. Dropping it does
  not change any result.
- The `__main__.py` edit casting `time` / `lhr` to float arrays. Harmless but
  unnecessary; both are only forwarded to `generate_power_history`.
- The `coaxial_contact` rework from `15b4b4e` (mesh, BCs, geometry,
  `non-regression.py`, **and its gold**). That is a separate piece of work on
  an existing verification case and needs its own review — it must not ride
  in on this branch.

## Run

    conda activate z3st
    ./Allrun      # gmsh → z3st → non-regression → plots → pressure_creep_analysis
    ./Allclean
