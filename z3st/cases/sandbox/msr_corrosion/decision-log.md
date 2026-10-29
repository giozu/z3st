# z3st MSR corrosion: decision log

Chronological record of the choices made on branch `feature/msr-corrosion`, with the reasoning behind each, so the process can be followed and cited in the papers. Newest entries at the bottom. Related files in this folder: `bhave2023_reference.md` (benchmark data), `design-reduced-model.md` (phase 1 design).

---

## 2026-10-29 · Findings from Bhave et al. 2023 (J. Nucl. Mater. 574, 154147)

- **Temperature settled.** All Bhave simulations and the Zheng 2015 experiments are at 973 K (700 °C), the same temperature as the Sankar & Singh tests.
- **In 1D their model reduces to ordinary diffusion.** Their sensitivity study (Fig. 3) shows the Cr diffusivity in the alloy dominates by more than 20x; the electric field, exchange current and interface energy hardly matter. The full phase-field result matches a simple diffusion solution with a fixed surface Cr content to within 5 % (Sec. 4, Fig. 4).
- **The salt redox state enters in one place.** It shifts the free energy of Cr2+ in the salt (Eqs. 30-31), which sets the Cr content at the surface. This is the link the impurity work builds on: impurity level -> redox potential -> surface Cr -> depletion depth.
- **Check of the reduction against their data.** Using their Hart effective diffusivity at 973 K in the analytic diffusion solution (Eq. 50), the 1000 h mass loss is 0.60 mg/cm2 against 0.42 measured for Ni-5Cr, and 4.5 against 7.1 for Ni-20Cr. This is the same over/under pattern they report for their own 1D runs, so the parameter set reads correctly.

| Alloy | D_eff (um2/s) | Eq. 50, 50 / 500 / 1000 h (mg/cm2) | Exp., Table 3 | 2 sqrt(D t) at 1000 h |
| --- | --- | --- | --- | --- |
| Ni-5Cr | 4.40e-5 | 0.13 / 0.42 / 0.60 | 0.09 / 0.27 / 0.42 | 25 um |
| Ni-20Cr | 1.53e-4 | 1.01 / 3.18 / 4.50 | 0.67 / 3.06 / 7.12 | 47 um |

## 2026-10-29 · D1: reduced model first, full phase field later

**Options considered.**
- (a) Rebuild Bhave's full model in z3st: chemical potentials, electric potential, order parameters. Needed later for grain-resolved runs, voids and a moving interface; this is what makes a methods paper.
- (b) Start with the reduced model: Cr diffusion in the alloy plus a redox-dependent surface condition. Enough for the 1D benchmark and much faster to results.

**Decision.** (b) first, then (a) built on top of it. The reduced model then serves as the 1D check of the full phase-field model.

**Rationale.** In 1D the two are equivalent within 5 % (Bhave Sec. 4), so (b) reaches the benchmark with a fraction of the effort, and it keeps a verified reference for (a).

**Correction (same day).** The recommendation originally said (b) would also be enough for the Sankar & Singh impurity series. The finding below shows this is unlikely; (b) remains the right first step, but Paper 2 needs more (see "Open for Paper 2").

## 2026-10-29 · Phase 1 design (`design-reduced-model.md`)

- New `CorrosionModel` in `z3st/models/corrosion_model.py`, switched on with `models: corrosion: true`, wired like `PorosityMigrationModel`.
- Unknown: Cr atomic fraction in the alloy only; the salt is not meshed. The salt-facing boundary carries the surface condition.
- Diffusivity: Arrhenius bulk + grain boundary, combined with the Hart rule (Bhave Eq. 46), evaluated at the local temperature (hook for thermal coupling).
- Surface condition: `fixed` (surface Cr given) or `redox` (computed from the fluoride potential E_F2/F- and the salt Cr2+/Ni2+ fractions). `redox` is the entry point for the impurity work.
- Outputs: Cr field, mass loss (mg/cm2), depletion depth (threshold fraction of the initial Cr, default 0.843 as in Bhave Fig. 2), surface Cr flux.
- Acceptance: V1, analytic erfc solution within 1 %; V2, Bhave 1D Ni-5Cr and Ni-20Cr within 2 % of the analytic mass loss, with the comparison to their experiments reported.

**Calls made in the design (to be checked by the user).**
1. The benchmark uses the paper's surface value, 1.8e-3. The `redox` value computed from their Table 1 and Appendix A data is 1.86e-4 (mass fraction 1.65e-4), about 11x lower. It is reported alongside and listed as a known item to chase (start from the H_xs/S_xs units and the kB T term in Eqs. B.8-B.9). Effect on the 1000 h mass loss: about 3 %, because the initial Cr content is much larger than the surface value.
2. Eq. 31 as printed omits E0_F2/F- and, taken literally, gives a fully Cr-saturated surface (x_surf = 1). The Eq. 30 form (shift E_F2/F- - E0_F2/F-) is used for both Cr and Ni.
3. Fe and multi-element depletion are left out of phase 1 (needed later for SS316H).
4. Design documents live in this case folder, not in `docs/`. Nothing is committed.

## 2026-10-29 · Open for Paper 2: the reduced model alone is unlikely to reproduce the Sankar trend

- Under Bhave's conditions the surface Cr is close to zero, so corrosion is limited by diffusion in the alloy. Changing the redox potential then barely changes the depletion depth.
- Sankar's 200 um attack in 100 h implies an effective diffusivity of order 0.03 um2/s, about 100x the Hart value for Ni-Cr at the same temperature. 316H (Fe-based) diffuses faster, but probably not by that much.
- So the 10 -> 200 um trend with the added impurity most likely needs an extra mechanism, to be chosen after re-reading Sankar & Singh 2022/2023:
  - limited oxidant supply in the salt (closed capsule, finite impurity inventory);
  - reaction-limited kinetics at the surface (flux condition instead of a fixed surface value);
  - intergranular penetration (grain-resolved, ties in with the full phase-field model).
- Also to check: whether Sankar's "attack depth" is the depth of intergranular attack rather than an averaged depletion depth. If so, it is not the same observable as the 1D depletion depth.
- Phase 1 is not affected.

## Earlier decisions (2026-10-27 to 2026-10-29, from the conversation)

- **Scope of the branch.** Study how salt impurities change corrosion of Cr-bearing structural alloys and carry it through to the component in z3st. Salt and alloy stay as inputs (not chosen yet).
- **Novelty.** One chain in one tool: salt chemistry and impurity level -> corrosion front (depletion, voids, oxide) -> local property change -> stress and damage in the wall. Existing phase-field work stops at the microstructure.
- **Scale.** Both: continuum wall first (effective diffusivity, full chain end to end), grain-resolved second (checks against the phase-field papers and calibrates the effective numbers for the continuum model).
- **Output.** Papers (P1 method + verification, P2 impurities, P3 mechanics coupling) forming a thesis chapter. No experimental data of our own; validation on published and partner data.
- **AIMD.** Only if the literature leaves key parameters unknown; decide once the parameter list is fixed.
- **Layout.** Physics in `z3st/models/corrosion_model.py`; cases in `cases/sandbox/msr_corrosion/` during development, moved to `cases/verification/` or `cases/benchmarks/` with non-regression once they pass. PDFs in `0.biblio/` stay untracked (`*.pdf` is gitignored).

## 2026-10-29 · Design approved; implementation plan written

- `design-reduced-model.md` approved by the user as written.
- Implementation plan (kept locally), 6 tasks: (1) pure functions, (2) model wiring + V1 erfc case, (3) diagnostics and `output/corrosion.tsv`, (4) V2 Bhave 1D benchmark, (5) temperature from the thermal model, (6) no-regression check and records.
- Choices made while planning:
  - Unknown is the Cr atomic fraction; Bhave's c_surf = 1.8e-3 is a mass fraction, so V2 uses its atomic equivalent x_surf = 2.0315e-3.
  - Two references for mass loss: the exact integral of the erfc profile (the model must match within 1 %) and Bhave Eq. 50 (design acceptance, 2 %). They differ by 0.7 % (Ni-5Cr) and 1.2 % (Ni-20Cr) at 1000 h because Eq. 50 treats mass fraction as linear in atomic fraction.
  - Checks are plain assert scripts run in the `z3st` conda env (pytest is not installed).
  - V2 uses one sub-folder per alloy because z3st always reads `./input.yaml`.
  - Axisymmetric is rejected for corrosion in phase 1.

## 2026-10-29 · Phase 1 implemented and verified

Implementation: `z3st/models/corrosion_model.py` (new), hooks in `core/config.py`, `core/finite_element_setup.py`, `core/solver.py`, `core/spine.py`, `utils/writer.py`, `__main__.py`. Switched on with `models: corrosion: true`; output field `CrFraction`, history `output/corrosion.tsv` (t, mass loss in kg/m2 and mg/cm2, depletion depth, surface flux).

**Verification.**
- `check_corrosion_functions.py`: Hart D, mass/atomic fraction, redox x_surf, analytic references. Pass.
- V1 (`V1_erfc`), erfc solution, D = 1e-16 m2/s, 1000 h: profile 0.51 %, mass loss 0.46 %, depth 0.97 %, surface flux 1.4 %. Errors scale as first order in dt (60 steps gave 1.67 %, 200 give 0.51 %).
- V1 with the thermal model on (`V1_erfc_thermal`), uniform 973 K: matches V1 to 2e-9. With the surface end at 1000 K the depth grows from 37.6 to 50.8 um.
- 2D variant, 1 vs 2 MPI ranks: identical mass loss, depth and flux.
- V2 (`V2_bhave1D`), Bhave 1D benchmark:

| Alloy | t (h) | z3st | Bhave Eq. 50 | Experiment (mg/cm2) |
| --- | --- | --- | --- | --- |
| Ni-5Cr | 50 / 500 / 1000 | 0.133 / 0.421 / 0.595 | 0.133 / 0.421 / 0.596 | 0.09 / 0.27 / 0.42 |
| Ni-20Cr | 50 / 500 / 1000 | 1.014 / 3.205 / 4.533 | 1.007 / 3.184 / 4.503 | 0.67 / 3.06 / 7.12 |

  Within 0.46 % of the exact solution and 0.66 % of Eq. 50. Against the experiments, the same over/under pattern as Bhave's own 1D Hart runs. Ni-5Cr depth at 500 h: 17.4 um (Bhave Fig. 2: 18 um). Plot: `V2_bhave1D/mass_loss.png`.
- Existing CI case list: 23/23 pass, gold regressions included.

**Choices made during implementation** (full list in the execution ledger):
- Gas constant R = 8.314462618 J/mol/K; the planning values used 8.314 (0.13 % in D).
- Every material in a corrosion run needs a `corrosion` card. Thermal-coupled runs also need E, nu, alpha on each material (existing z3st requirement of `get_results`).
- The D unit-slip warning is evaluated at `corrosion.T`, or at the first solved temperature when the thermal model supplies T.
- A salt-facing region with no boundary facets stops the run with an error.

**Finding for phase 2 (corrects an earlier note).** In `redox_surface_fraction`, applying the same E_F shift to Cr and Ni with the salt Ni2+ fraction held fixed makes E_F cancel when the surface is Ni-rich (Bhave's conditions, E_F above about 3 V): x_surf then depends only on x_Cr_salt / x_Ni_salt. At lower E_F it does depend on E_F (x_Cr = 9.7e-6 at 2.6 V against 1.86e-4 at 3.36 V). For the redox potential to act as the impurity lever in phase 2, the salt Ni2+ (or the Ni chemical potential) has to follow E_F rather than stay at Bhave's fixed value.

**Known limits left as they are:** depth assumes one planar salt face along `depth_axis`; molar masses are fixed to Ni-Cr; the diffusion step uses the thermal field while the redox x_surf uses `corrosion.T`; no case exercises adaptive time-step rollback (the code path was reviewed, not run).

## Pending

- [x] User approval of `design-reduced-model.md`, then the step-by-step implementation plan.
- [x] User review of the implementation plan, then implementation.
- [ ] Move V1 and V2 to `cases/verification/` and `cases/benchmarks/` with a gold file (until then the suite does not protect the corrosion model).
- [ ] Resolve the 1.8e-3 vs 1.86e-4 surface value discrepancy.
- [ ] Re-read Sankar & Singh 2022/2023 for the attack-depth definition and capsule/impurity inventory, before planning phase 2.
