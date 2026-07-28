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

## Status: NOT yet verified — do not bless a gold from this as-is

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

3. **OPEN — `creep_Q`.** This snapshot runs with `creep_Q: 1.2e5` (Arrhenius
   active, non-isothermal). Esposito's derivation is isothermal, which is why
   the `creep_Q: 0` path was needed — see the companion fix in
   `models/creep_model.py`. Decide explicitly which of the two the
   verification claims, and state it. `pressure_creep_analysis.py` evaluates
   φ with `exp(-Q/RT)` at a clad temperature read from the solution, so the
   two choices are not equivalent.

4. **OPEN — modelling choice to document.** `boundary_conditions.yaml`
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
