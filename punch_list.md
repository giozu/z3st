# z3st punch list

Audit performed 2026-05-18 against develop @ dbcee1d, generated from four parallel reviewers (code, cases, README, paper). Re-audited 2026-06-10 against develop @ 131bc73 (four parallel reviewers: core, models, utils/cases, docs-vs-code) — verdicts folded into the items below. Item IDs are stable — reference them in commits / PRs.

## Start here

All P0 items resolved (CODE-P0-5 fixed 2026-06-10). A full four-agent re-audit was run on 2026-06-10 against develop @ 131bc73; verdicts and new items are folded in below. The structural items now worth attention: **[CODE-P1-14]** (orphaned `heat_flux` diagnostic — keep/rewire/delete decision), **[CASES-FOLLOWUP-2]** (the `plane_stress` regime has **zero** exercising cases since `U_thick_cylindrical_shell_plane_stress` was deleted in f1bb70b), and **[CASES-P1-7]** (`U_pwr_rod_2D` not wired into any regression suite).

## Resolved

- **2026-06-10 — [CODE-P0-5](#code-p0-5) `plane_stress` missing from solver regime lists.** Added `"plane_stress"` to the traction-normal branch (solver.py:343) and the body-force branch (solver.py:366), matching the existing handling in `mechanical_model.py`. NOTE: no case currently exercises this regime (see [CASES-FOLLOWUP-2](#cases-followup-2)), so the fix is read-verified only.
- **2026-06-10 — regime validation added (new, no prior ID).** `config.py` now validates `regime ∈ {2d, 3d, axisymmetric, plane_stress}` after lowercasing and raises a clear `ValueError` on anything else (previously a typo silently fell through to undefined downstream branches).
- **2026-06-10 — [CODE-P2-1](#) debug scaffolding removed.** spine.py (function-space introspection prints) and finite_element_setup.py (mixed-element scaffolding) cleaned. gap_model.py alt-distance comment also removed.
- **2026-06-10 — [CODE-P2-2](#) Italian error message translated** (mechanical_model.py: "Formato displacement non valido" → "Invalid displacement format").
- **2026-06-10 — [CODE-P2-3](#) side-effect list-comprehension** at solver.py (bcs_m population) converted to a for-loop. (Only one instance remained; the second cited occurrence no longer exists.)
- **2026-06-10 — [CODE-P2-4](#) material dump gated behind `--debug`** (spine.py::load_materials; checks `"--debug" in sys.argv`, same convention as `__main__.py`).
- **2026-06-10 — [CODE-P2-5](#) stag_tol defaults unified** at 1e-4: `solve_staggered` signature defaults now mirror `Spine.solve` (the only caller, which always passes explicitly), with a sync comment.
- **2026-06-10 — export_vtu doc stragglers removed** (follow-up to CODE-P2-7): README.md directory tree + post-processing table, docs/source/getting_started.rst (stale `python3 -m z3st.utils.export_vtu` command), docs/source/index.rst module list. api.rst gained `z3st.utils.writer`, `contact_model`, `plasticity_model` automodule entries. README case tree now shows `U_pwr_rod_2D` and `V_*`.
- **2026-06-10 — contact_conductance silent degradation now warns** (gap_model.py): when contact coupling is enabled but fewer than two materials carry a *numeric* `k` (e.g. symbolic-k cards), the Ross-Stoute term silently returned 0.0; it now prints a `[WARNING]` before doing so. Full symbolic-k support in the harmonic mean is deferred (relates to [CODE-P2-6](#)).
- **2026-06-10 — traction 'raw' error message enriched** (solver.py): now reports the offending type/value and the expected format.
- **2026-06-10 — [INFRA-P0-1] gold-regression verdict was silently dropped — fixed + golds re-blessed.** Discovered while smoke-testing: `utils_verification.py::regression_check`'s PASS/FAIL was printed to stdout but written nowhere — `non-regression.sh` reports only the analytic `"summary"` JSON key (written *before* the gold check), and CI (`non-regression_github.sh`) gates only on Allrun's exit code while the per-case `non-regression.py` always exits 0. **Net effect: neither a gold regression nor even an analytic-tolerance FAIL could ever fail the local summary or CI.** Three patches: (i) `regression_check` now persists `"regression": "PASS"/"FAIL"` into `non-regression.json`; (ii) `non-regression.sh` prints a "Gold regression" line per case (and the no-JSON fallback now says "NOT VALIDATED" — addresses the misleading-summary half of [CASES-P0-3](#cases-p0-3)); (iii) `non-regression_github.sh` fails a case (and CI) when `summary` or `regression` is FAIL, with the reason in the summary line. Verified end-to-end on `1_thin_slab_2D` including the negative path (injected 5% gold drift → `regression: FAIL` → shell parse → CI gate). **Consequence of the blind spot:** 8 cases had silently drifted past the 1e-3 gold tolerance (worst: 90% rel diff on `3_thin_slab_adiabatic_2D` L2_error_T). Audit of all drifted metrics vs the *analytic* truth showed the drift is almost entirely accuracy **improvements** (case 8 errors ↓ ~10×, case 17 ↓ ~4 orders — consistent with the CODE-P1-11 quadrature bump and the solver/BC refactors since blessing); only 3 metrics negligibly worse (e.g. 15_spherical σ_tt 0.04463 → 0.04467). Golds re-blessed for: `1_thin_slab_2D`, `2_thin_cylindrical_shell_2D`, `3_thick_slab_adiabatic_3D`, `3_thin_slab_adiabatic_2D`, `8_thick_cylindrical_shell_plane_strain_2D`, `12_cylindrical_shell_thermal_gradient_2D`, `15_spherical_pressurized_cavity`, `17_stress_strain_curve_stress`. Full 28-case local suite green (analytic) on 2026-06-10 with all of today's code changes in.

- **2026-05-19 — [CODE-P0-3](#code-p0-3) Missing MPI allreduce + rank-0 guard on diagnostics.** Three patches: `damage_model.py::compute_energy_balance` allreduces `E_el`/`E_frac`; `thermal_model.py::heat_flux` allreduces all five scalars; `__main__.py:208-218` guards the `energies.txt` write with `if comm.rank == 0`. Verification owed under `mpirun -n 2`.
- **2026-05-19 — [CODE-P0-2](#code-p0-2) `plane_stress` damage degradation + tensor-shape mismatch.** Single patch in `mechanical_model.py:537-558`: drop 3×3 padding and early `return` from the plane_stress branch in `sigma_mech`. σ is now 2×2 (matching ε(u), ε(v), and σ_th in this regime) and falls through to the damage block. Verification owed via `cases/U_thick_cylindrical_shell_plane_stress/` (no-damage smoke test) and a damage-enabled plane-stress case to be added.
- **2026-05-19 — [CODE-P0-1](#code-p0-1) `gap_model.py` AttributeError.** `self.locateFacetDofs(...)` → `self.mgr.locate_facets_dofs(...)` at gap_model.py:59-60. `Gas` gap path is now invokable; exercised by `cases/U_coaxial_contact_2D/`.
- **2026-05-19 — [README-P0-1](#readme-p0-1) wrong cases path in Quick installation.** Two `cd ~/z3st/cases/00_example/` lines updated to `cd ~/z3st/z3st/cases/00_example/`.
- **2026-05-18 — [PAPER-P0-1](#paper-p0-1) AT1 `Gc`/`ψc` factor-of-4 error.** Option A taken: code is correct; paper, `CONTEXT.md` §9c, and `uo2.yaml:18-19` updated to `Gc ≈ 372 J/m²`, `ψ_c ≈ 1.4 MJ/m³`. Calibration paragraph reframed: σ_c=1 GPa is the AT1 sweet spot under the rigid `Gc = (8/3) ℓc σ_c²/E` identity; McClenny's macro-tuned `Gc ≈ 80 kJ/m²` is acknowledged but is not directly reachable in strict AT1. Root cause of the error: `σ_c = 500 MPa` was inadvertently used in the arithmetic (the 2 GPa values were always correct because they were derived independently). No re-run needed.

---

## P0 — load-bearing bugs and errors

### Code

<a id="code-p0-1"></a>
- [x] **CODE-P0-1** — *Resolved 2026-05-19 (verification deferred).* `gap_model.py:59-60` updated to `self.mgr.locate_facets_dofs(...)`. Fix is read against the API; **end-to-end Gas-gap verification is still owed** — `cases/U_coaxial_contact_2D/` (the only Gas-gap case) needs three edits before it actually exercises this code path: (i) add `gap_conductance: { type: Gas, value: <prefactor> }` under `models:` in `input.yaml`; (ii) introduce an actual gap in geometry (`inner_radius_2: 0.056` in `geometry.yaml` to leave 1 mm) and mirror in `mesh.geo`; (iii) re-mesh and run. Tracked separately as **[CASES-FOLLOWUP-1](#cases-followup-1)** below. Also: `set_gap_temperature` lines 71-72 still `break/break` after the first Robin BC of the first material — single-gap-pair-only by design; worth revisiting if multi-pair gap support is ever needed.

<a id="code-p0-2"></a>
- [x] **CODE-P0-2** — *Resolved 2026-05-19 (verification deferred).* `mechanical_model.py:537-558` updated: dropped the 3×3 padding and the early `return` in the `plane_stress` branch, so σ is now built 2×2 (matching `epsilon(u)`, `epsilon(v)`, and `sigma_th`, which were already 2×2 in plane_stress) and falls through to the damage `g(D)` block. Both sub-bugs fixed in one patch (Bug A: damage degradation silently skipped; Bug B: 3×3 vs 2×2 shape mismatch with `sigma_th` / `epsilon(v)` that would fail UFL form compilation). End-to-end verification owed by running `cases/U_thick_cylindrical_shell_plane_stress/` (no-damage path) and a damage-enabled plane-stress case (none exists yet — would need to be added).

<a id="code-p0-3"></a>
- [x] **CODE-P0-3** — *Resolved 2026-05-19 (verification deferred).* Three patches: (i) `damage_model.py::compute_energy_balance` now `allreduce`s `E_el` and `E_frac` over `self.mesh.comm` before returning (with `from mpi4py import MPI` added); (ii) `thermal_model.py::heat_flux` now `allreduce`s `q_integral`, `qx/qy/qz_integral`, and `volume` before averaging; (iii) `__main__.py:208-218` guards the `energies.txt` write with `if problem.mesh.comm.rank == 0:`. The `print` statements remain unguarded — they will duplicate under MPI with N ranks (cosmetic) but never corrupt. Verification owed: run any damage case under `mpirun -n 2` and confirm `energies.txt` is clean and `E_el` matches the serial value.

<a id="code-p0-4"></a>
- [x] **CODE-P0-4** — *Resolved 2026-05-19.* Fixed incidentally while refactoring the output section: `__main__.py:149` (was :157) now reads `from z3st.core.mesh.plotter import MeshPlotter`. `--mesh_plot` no longer ImportError's.

<a id="code-p0-5"></a>
- [x] **CODE-P0-5** — *Resolved 2026-06-10.* `"plane_stress"` added to both regime branches in `solver.py::_mechanical_step` (traction-normal at :343, body-force at :366), consistent with `mechanical_model.py`. Read-verified only — no case exercises the regime (see [CASES-FOLLOWUP-2](#cases-followup-2)).

### Cases

<a id="cases-p0-1"></a>
- [ ] **CASES-P0-1** — `z3st/cases/U_cluster_dynamics_test/` and `z3st/cases/U_quarter_block/` are missing `non-regression.py`. They silently "PASS" under `cases/non-regression.sh` via the file-absence fallback.

<a id="cases-p0-2"></a>
- [ ] **CASES-P0-2** — `z3st/cases/II_attenuation_map/` has neither `Allrun`/`Allclean` nor `non-regression.py` (uses custom `run_map.py`, `plot_map.py`). It is listed in CONTEXT.md §6.1 as a regular case but cannot be driven by `non-regression.sh`. Either add wrappers or document it as a non-standard utility.

<a id="cases-p0-3"></a>
- [ ] **CASES-P0-3** — `z3st/cases/I_mesh_sensitivity_2D/` is a custom mesh-sweep driver by design, but `non-regression_summary.txt` records it as "OK" through the same no-JSON fallback. Make this explicit so the summary is not misleading.

### README

<a id="readme-p0-1"></a>
- [x] **README-P0-1** — *Resolved 2026-05-19.* Two `cd ~/z3st/cases/00_example/` lines (README.md:51 and :60) updated to `cd ~/z3st/z3st/cases/00_example/`. The relative path `../../utils/plot_convergence.py` already works from there.

### Paper

<a id="paper-p0-1"></a>
- [x] **PAPER-P0-1** — *Resolved 2026-05-18.* AT1 `Gc`/`ψc` were off by 4× (root cause: σ_c=500 MPa used in arithmetic instead of 1 GPa). Code is correct: `core/spine.py:147-159` uses the canonical AT1 identity `Gc = (8/3) ℓc σ_c²/E`. With (σ_c=1 GPa, ℓ_c=50 μm, E=358 GPa) the correct values are `Gc ≈ 372 J/m²` and `ψ_c = σ_c²/(2E) ≈ 1.4 MJ/m³`. Updated in three locations: `main.tex:317`, `CONTEXT.md:650`, `uo2.yaml:18-19`. Calibration paragraph reframed (McClenny's macro-tuned `Gc ≈ 80 kJ/m²` is acknowledged as unreachable in strict AT1).

<a id="paper-p0-2"></a>
- [ ] **PAPER-P0-2** — `main.tex:181` says case 17 plate is "1.0 × 0.5 m"; `cases/17_stress_strain_curve_displacement/geometry.yaml` is `Lx=1.0, Ly=0.001` (1 m × 1 mm).

<a id="paper-p0-3"></a>
- [ ] **PAPER-P0-3** — `main.tex:193` says case 20 has hardening `H = 20 GPa`; `materials/steel.yaml` has `hardening_modulus: 10.0e+9`.

<a id="paper-p0-4"></a>
- [ ] **PAPER-P0-4** *(refined 2026-05-19)* — `main.tex:202` describes case 9 as "2D generalised plane strain ... out-of-plane strain treated as a single global degree of freedom". The case actually uses `regime: axisymmetric` with the GPS condition `ε_zz = const` **enforced via a prescribed axial-end displacement** (`Clamp_y` at `top` with value `-1.359155e-6 m = ε_zz·Lz`), not via a global Lagrange-multiplier DoF. The *physics* is a valid GPS configuration; only the implementation description is misleading. **Fix**: tighten the wording to "the out-of-plane strain is held constant across the section through a prescribed axial-end displacement". No section rename or new feature needed. Numerical side resolved separately — see [CASES-FOLLOWUP-3](#cases-followup-3).

<a id="paper-p0-5"></a>
- [ ] **PAPER-P0-5** — `main.tex:321` and `:366` quote `E_frac` ramps "0.56 J → 3.74 J"; actual `cases/14_full_cylinder_cracking_2D_xy/energies.txt` step 0 is **0.316 J**, step 1 is **0.77 J**, end is 3.737 J. `CONTEXT.md:684` (§9c) carries the same stale 0.56 J. Either re-run case 14 and update both, or fix the paper to "0.32 → 3.74 J".

---

## P1 — correctness / clarity

### Code

<a id="code-p1-1"></a>
- [x] **CODE-P1-1** — *Resolved 2026-05-19.* `DamageModel.__init__` (`damage_model.py`) now validates `damage.type ∈ {"AT1", "AT2"}` up-front, raising a clear `ValueError` on any other value (lowercase typo, hyphenated form, etc.). Also reordered: the missing-block check now runs *before* the `setdefault("linear_solver", ...)` so an entirely absent damage block still surfaces clearly rather than being silently filled in.

<a id="code-p1-2"></a>
- [x] **CODE-P1-2** — *Resolved 2026-05-19.* `spine.py:160` now uses `isinstance(Gc, (int, float, np.floating, np.integer))` so `Gc` provided as a numpy scalar (or Python `int`) is recognised and the σc-from-Gc conversion runs as intended.

<a id="code-p1-3"></a>
- [x] **CODE-P1-3** — *Resolved 2026-05-19.* `solver.py:559` diagnostic print now formats `Gc` / `sigma_c` conditionally: scalars use `{:.2e}` as before; non-scalars (UFL expressions, when the material's `Gc` comes from a Python callable such as `materials/oxide.py::Gc`) display `<ClassName>` instead. No more `TypeError` on AT1 runs with symbolic `Gc`.

<a id="code-p1-4"></a>
- [x] **CODE-P1-4** — *Resolved 2026-05-19 (false alarm, comment added).* Re-derivation shows the current ordering `α·D_new + (1-α)·D_old → max(·, D_old) → clip` is algebraically identical to the audit's proposed "relax-the-increment" formulation: `α·D_new + (1-α)·D_old = D_old + α·(D_new - D_old)`, and `max(D_old + α·X, D_old) = D_old + α·max(X, 0)` for α ≥ 0, so the composite equals `D_old + α·ΔD⁺` with `ΔD⁺ = max(0, D_new - D_old)`. The `max` projection isn't "clobbering" relaxation; it's correctly enforcing irreversibility when the raw `D_new` from the linear solve happens to dip below `D_old` (which the AT1 elliptic problem doesn't prevent by construction). Added an inline comment in `solver.py:582-591` explaining the equivalence so a future reader doesn't reorder it wrong.

<a id="code-p1-5"></a>
- [x] **CODE-P1-5** — *Resolved 2026-05-19.* `__main__.py` now sets `problem.n_steps = len(times)` immediately after `generate_power_history` and before `set_boundary_conditions`, so BC setters validate against the actual time-loop length, not the raw YAML `n_steps`.

<a id="code-p1-6"></a>
- [x] **CODE-P1-6** — *Resolved 2026-05-19.* `mechanical_model.py` Neumann branch now validates `len(traction_list) == self.n_steps` and exits with a clear `[ERROR]` message on mismatch — mirrors the existing Dirichlet check at L100.

<a id="code-p1-7"></a>
- [x] **CODE-P1-7** — *Resolved 2026-05-19 (documented, not extended).* `thermal_model.py` Dirichlet branch now detects list inputs and emits a clear `[ERROR]` with a pointer to the mechanical-block pattern, instead of silently fixing T to `temperature[0]` or crashing inside dolfinx. Extending thermal Dirichlet to step-dependent lists is deferred — would require wiring per-step BC updates into `_thermal_step` and is non-trivial. Logged as [CODE-FOLLOWUP-1](#code-followup-1) below.

<a id="code-p1-8"></a>
- [x] **CODE-P1-8** — *Resolved 2026-05-19.* Moved the `lame → plasticity` promotion (for materials with `yield_strength` when global plasticity is on) from inside `MechanicalModel.sigma_mech` to `Spine.load_materials`. The material dict is now deterministic at load time; `sigma_mech` just reads `material["constitutive_mode"]` instead of side-effect-mutating it on the first call.

<a id="code-p1-9"></a>
- [x] **CODE-P1-9** — *Resolved 2026-05-19.* `set_power` now guards against the singular case: `gamma_heating > 0` with `geometry_type ∈ {cyl, cylinder, sphere}` and `inner_radius == 0.0` raises a clear `ValueError` instead of silently producing zero (cylindrical: `K_0(0) = +∞` in the denominator) or NaN (spherical: `inner_radius/r` with `r → 0`) heating.

<a id="code-p1-10"></a>
- [x] **CODE-P1-10** — *Resolved 2026-05-19.* `set_power` now `+=` into `q_third` for both the fissile branch and the gamma-heating branch (previously both used `=`, so a material with both flags lost the fissile contribution). Physically correct: the two heating sources superpose on a material that has both. No behavioural change for cases that set only one of the two (the `q_third` array is zero-initialised, so `+= v` == `= v` on first write).

<a id="code-p1-11"></a>
- [x] **CODE-P1-11** — *Resolved 2026-05-19.* `q_degree = mech_degree + 1` → `q_degree = 2 * mech_degree + 1` in `finite_element_setup.py`. Aligns with the standard rule for full integration of the J2 return-mapping form (whose integrand contains nonlinear functions of `u` through `n = s_trial / |s_trial|_eq`). For `mech_degree = 1` (the default), `q_degree` goes from 2 → 3. **Owed verification:** rerun `cases/20_plasticity_2D/` and `demo_CP_single_grain/`. If their non-regression numbers shift (very plausible — small precision delta in the plastic-strain integration), re-bless their gold JSONs. Tracked in [CASES-FOLLOWUP-4](#cases-followup-4) below.

<a id="code-p1-12"></a>
- [x] **CODE-P1-12** — *Resolved 2026-05-19.* Added explicit `else: raise ValueError(...)` after the AT1 branch in `crack_driving_force` (damage_model.py) **and** in the parallel `gamma_density` function (same shape-of-bug). A typo in `damage.type` now fails with a clear message instead of returning `None` and crashing later inside `fem.Expression`.

<a id="code-p1-13"></a>
- [x] **CODE-P1-13** — *Resolved 2026-05-19.* Each of `mech_cfg`, `th_cfg`, `dmg_cfg` now `.setdefault("linear_solver", "iterative_hypre")` right after `self.input_file.get(...)`, so the solver no longer `KeyError`s when the user omits the key. Choice matches the README example. Per-case overrides still respected.

### Cases

<a id="cases-p1-1"></a>
- [x] **CASES-P1-1** — *Resolved 2026-05-19.* CONTEXT.md per-case structure now reads `mesh.msh    generated mesh (NOT tracked in git per .gitignore — Allrun regenerates it with gmsh on each run)`.

<a id="cases-p1-2"></a>
- [x] **CASES-P1-2** — *Resolved 2026-05-19.* CONTEXT.md §6.1 case-14 line now lists all three damage-active variants (`14_full_cylinder_cracking` 3D, `14_full_cylinder_cracking_2D_xy` plane-strain workhorse, `14_full_cylinder_thermal_2D_rz` thermal verification) with a back-pointer to §9c for the three-variant rationale.

<a id="cases-p1-3"></a>
- [x] **CASES-P1-3** — *Resolved 2026-05-19.* Removed the stale `# "15_box_elliptical_cavity_2D"` line from `non-regression.sh` (the name doesn't exist on disk) and replaced with commented entries for the actual gold-less cases (`15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`, `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch`, `18_box_knotch_2D`) with a pointer to [CASES-FOLLOWUP-5](#cases-followup-5). The summary text file is regenerated on every `./non-regression.sh` run, so it'll self-heal next pass.

<a id="cases-p1-4"></a>
- [x] **CASES-P1-4** — *Resolved 2026-05-19.* `non-regression_github.py` now resolves the shell driver via `Path(__file__).resolve().parent / "non-regression_github.sh"` (the script sits in the same directory as this Python wrapper). Replaces the stale `parent.parent / "tests" / "non-regression_github.sh"` path.

<a id="cases-p1-5"></a>
- [x] **CASES-P1-5** — *Resolved 2026-05-19.* Added a comment block at the head of `non-regression_github.sh` documenting the inclusion policy: tight 9-case subset for fast CI turnaround, one case per orthogonal physics path, blessed-gold requirement; explicit exclusion list (active case-14 cracking variants, case-19 SENT/SENS) with the reason (too long for CI, golds being blessed) and the local-suite escape hatch (`cases/non-regression.sh`).

<a id="cases-p1-6"></a>
- [x] **CASES-P1-6** — *Triaged 2026-05-19 → tracked as [CASES-FOLLOWUP-5](#cases-followup-5).* Inventory is correct: the three case-14 cracking variants (WIP, expected), plus 7 cases (`15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`, `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch`, `U_box_knotch_3D`, `U_coaxial_contact_2D`, `U_slab_contact`) that have a `non-regression.py` but no gold. Resolution requires running each and copying `output/non-regression.json` → `output/non-regression_gold.json` — substantial compute and per-case sanity-check work, deferred to the follow-up.

### README

<a id="readme-p1-1"></a>
- [x] **README-P1-1** — *Resolved 2026-05-19.* Key Features section rewritten to cover the full capability set: 4 kinematic regimes, 5 constitutive routes (lame/voigt/hyperelastic/plasticity/custom), AT1/AT2 + Miehe/Amor splits + Ambati-hybrid constraint, gap-conductance modes, full BC inventory, Python material modules, PETSc backends, the new unified VTU/XDMF writer, full material card list.

<a id="readme-p1-2"></a>
- [x] **README-P1-2** — *Resolved 2026-05-19.* Example YAML now uses canonical lowercase `regime: 2d` (with the full enumeration as an inline comment), `true`/`false` instead of `True`/`False`. Also added the missing `damage:` and `gap_conductance:` blocks as commented placeholders, plus an `output: format:` block. Per-case `input.yaml` files on disk remain mixed-case in places — a future sweep could lowercase them all, but it's lower-priority cleanup.

<a id="readme-p1-3"></a>
- [x] **README-P1-3** — *Resolved 2026-05-19.* Directory tree rewritten from scratch to match `CONTEXT.md` §2.1: correct nesting (everything under `z3st/`), `core/mesh/` subpackage shown, `models/plasticity_model.py` listed, full material card inventory, complete `utils/` listing including the new `writer.py`, `cases/` shows the right `1_thin_slab_2D/` (was `1_thin_thermal_slab/`) plus case-14 and case-19 highlights with one-line annotations. With a pointer at the bottom: "the full case catalogue and per-module details are in CONTEXT.md".

<a id="readme-p1-4"></a>
- [x] **README-P1-4** — *Resolved 2026-05-19.* Author name typo fixed: `Scrogggs` → `Scroggs` in the basix2022a BibTeX entry. The basix2022b entry already had the correct spelling.

<a id="readme-p1-5"></a>
- [x] **README-P1-5** — *Resolved 2026-05-19.* `CITATION.cff` rewritten to the cff-version 1.2.0 schema: `title`, `authors` (with `affiliation: Politecnico di Milano`), `version`, `date-released`, `license`, `repository-code`, `url`, `doi`, `identifiers`, and a `keywords` block. Matches the README BibTeX `Z3ST2025` entry.

### Paper

<a id="paper-p1-1"></a>
- [x] **PAPER-P1-1** — *Resolved 2026-05-19.* Added `\label{sec:plasticity}` at the J2 plasticity subsubsection (main.tex:115). The cross-reference in the plasticity verification paragraph (main.tex:194) is now `Section~\ref{sec:plasticity}` instead of the wrong `Section~\ref{sec:phase-field}`.

<a id="paper-p1-2"></a>
- [x] **PAPER-P1-2** — *Resolved 2026-05-19.* Case-13 description updated to say "a temperature-dependent thermal conductivity $k(T)$ supplied through the Python callable extensibility hook (`materials.ceramic.k`; see Section~2.2)", reflecting what the material card actually does. The misleading scalar `k = 50 W/(m K)` is gone.

<a id="paper-p1-3"></a>
- [x] **PAPER-P1-3** *(optional polish — closed as no-action 2026-05-19)*: PAPER-P0-1 is now resolved (the arithmetic + reframed calibration is in place at main.tex:317). Re-reading the case-14 narrative with fresh eyes, the AT1-sharp-threshold mechanism is **already explicit** in two places: (a) the Results paragraph at main.tex:321 ("$\sigma_{\theta\theta}$ recovers above the AT1 threshold a few elements further along the arc and nucleates the next crack"), and (b) the Discussion at main.tex:366 ("we attribute the difference to the AT1 sharp elastic threshold, which prevents the damage variable from accumulating spuriously in cells whose driving energy is below $\psi_c = 3 G_c / (16 \ell_c)$ ..."). No further tightening warranted.

---

## P1 additions — 2026-06-10 re-audit

<a id="code-p1-14"></a>
- [x] **CODE-P1-14** — *Resolved 2026-06-10 (decision: fix + re-wire, per Giovanni).* `thermal_model.py::heat_flux` rewritten: (i) flux components now follow the mesh's geometric dimension (with `r`/`z` labels in axisymmetric) instead of unconditional x/y/z indexing; (ii) accepts symbolic-k cards — when `material["k"]` is the UFL expression spine resolved at init, it is used directly instead of being wrapped in `fem.Constant` (this also closes [CODE-P2-6](#)). Re-wired in `__main__.py` after each step's `get_results()`, gated behind `--debug` so default logs are unchanged. Smoke-tested on `1_thin_slab_2D --debug`: |q| = q_x = 4810 W/m² (matches the imposed flux), q_y = q_z = 0. The CODE-P0-3 MPI allreduce is preserved.

<a id="cases-p1-7"></a>
- [x] **CASES-P1-7** — *Resolved 2026-06-10.* `U_pwr_rod_2D` now has a `non-regression.py` reading end-state PCMI scalars from `output/history.csv` (format-independent — the run writes XDMF only): final mean/peak burnup, gap, contact pressure, final/peak T_max. One analytic anchor: the nodal-mean burnup matches the closed form `Σ lhr_k·Δt_k/(area·ρ·HM·8.64e10)` to 7e-8 (validates the state-bus arithmetic end-to-end); the PCMI scalars are gold-protected (`rel_error = 0` so the analytic gate ignores them). Allrun extended to chain `non-regression.py` + `plots.py`; wired into `non-regression.sh`. **BC fix (same day, Giovanni's request):** 2 MPa He fill-gas pressure added on all gap/plenum-facing surfaces (`lateral_1`, `top_1`, `inner_2`) — previously the clad saw 15.5 MPa outside and vacuum inside. Effect, in the physically expected direction: PCMI onset 28.8 → 30.4 MWd/kgU (later — gas opposes closure), final contact pressure 40.3 → 38.8 MPa, final gap −0.81 → −0.78 µm; burnup unchanged (closed-form anchor still 7e-8). Gold blessed from the fill-gas run. Known residual simplifications, deliberate: Dirichlet coolant (no film drop — [CODE-FEATURE-3](#code-feature-3)), no end-cap axial load / GPS top (needs `dolfinx_mpc`), gas + contact tractions coexist after closure (2 MPa ≪ 39 MPa, stands in for interface-roughness gas). Also added the four `V_*` verification cases to the local suite (they were CI-only — local is now a superset of CI) and dropped the stale `15_box_elliptical_cavity_2D` comment. Still CI-excluded (run is ~3 min — Giovanni's call whether to add it to the tight CI subset).

<a id="code-p2-9"></a>
- [x] **CODE-P2-9** — *Resolved 2026-06-10.* `core/mesh/reader.py` XDMF fallback: the bare `except:` around the "Grid"-name `read_meshtags` retry (which would even swallow `KeyboardInterrupt`) narrowed to `except RuntimeError`, matching the outer handler.</a>

---

## P2 — cleanup / style

### Code

- [x] **CODE-P2-1** — *Resolved 2026-06-10.* Commented scaffolding removed in spine.py, finite_element_setup.py, gap_model.py.
- [x] **CODE-P2-2** — *Resolved 2026-06-10.* Italian error message translated (mechanical_model.py).
- [x] **CODE-P2-3** — *Resolved 2026-06-10.* The one surviving side-effect list-comprehension (solver.py, bcs_m) converted to a for-loop.
- [x] **CODE-P2-4** — *Resolved 2026-06-10.* Material dump gated behind `--debug` (spine.py::load_materials).
- [x] **CODE-P2-5** — *Resolved 2026-06-10.* `solve_staggered` signature defaults aligned to `Spine.solve`'s 1e-4 (the defaults were dead code — spine always passes explicitly — but represented intent mismatch).
- [x] **CODE-P2-6** — *Resolved 2026-06-10 together with [CODE-P1-14](#code-p1-14).* `heat_flux` now uses the UFL expression directly when `material["k"]` is symbolic (spine resolves `k(T)` at field init) and only wraps genuine scalars in `fem.Constant`.
- [x] **CODE-P2-7** — *Resolved 2026-06-10. Dead-code removal (verified zero callers).* Deleted `z3st/utils/export_vtu.py` (290 LoC, fully superseded by `writer.py`; its only reference was a — itself broken — lazy-facade entry in `__init__.py`, now removed) and four unused functions in `utils_extract_vtu.py` (`extract_VonMises`, `_detect_VonMises`, `extract_spherical_stresses`, `save_csv_principal_stress`). Full suite re-verified green. **Deliberately KEPT (false positives):** `solid_swelling`/`rim_peaking` (dispatched via YAML strings), intra-file helpers (`_detect_*`, `_data_coords`, `_hydrostatic`, `_radius`), and standalone entry-point scripts (`plot_convergence.py`, `z-gui.py`).
- [ ] **CODE-P2-8** — *Further dead-function audit. Verification status updated by the 2026-06-10 re-audit:*
  - (i) `core/diagnostic.py` matrix-debug helpers (`analyze_matrix`, `assemble_and_analyze_matrix`, `check_fixed_dofs`, `check_mechanical_constraints`, `debug_dirichlet_mechanical`, `dense_condition`, `estimate_condition_number_sparse`, `is_symmetric`) — **CONFIRMED zero callers** across z3st/, cases/, docs/, notebooks (grep over .py + .ipynb). Safe to delete; Giovanni may want to keep some as solver-debugging aids — his call.
  - (ii) `utils/output.py` — **dead confirmed**: `export_xdmf`, `plot_residuals`, `plot_scalar`, `plot_radial_displacement`, `thermal_shield_slab_analytical_*`. **Live, keep**: `lame_solutions` (notebook), `plot3d_*`/`plot1d_*` (notebook + case diagnostics).
  - (iii) `utils/utils_plot.py` orphans (`plot_field_along_r_xy`, `plot_field_along_x`, `plot_sigma_cyl`, `plot_sigma_principal`, `rescale_axis`) — still need a full notebook audit.
  - (iv) `utils/utils_load.py` (`clear_log_file`, `get_ksp_reason_name`, `read_power_history`) — **no grep match anywhere**; likely dead, final notebook check before removal.
  - Also: `thermal_model.py::heat_flux` joined the orphan set — tracked separately as [CODE-P1-14](#code-p1-14) because of its latent bugs. The thin `extract_*` wrappers in `utils_extract_vtu.py` are used but overlap with `extract_field` — could be collapsed in the same pass.

### Cases

- [ ] **CASES-P2-1** — Allrun convention drift in `U_cluster_dynamics_test/` (uses `log.txt` instead of `log_z3st.md`, `gmsh -1`, no nr step) and `U_quarter_block/` (truncated pipeline).
- [ ] **CASES-P2-2** — `Allclean` is only called from `non-regression.sh`, never from inside individual Allrun scripts. Confirm design intent.
- [ ] **CASES-P2-3** — *Case-naming taxonomy cleanup (deferred — deliberate one-time refactor, not urgent).* The descriptive part of every case name is good (geometry + thickness + physics + dimension, e.g. `thick_cylindrical_shell_adiabatic_2D`) and should be preserved. The problem is the *prefix/role* axis: there are ~7 conventions for "what kind of case is this" — bare number (39 cases), `V_` (4, verification, added 2026-06-09), `U_` (5, demos/utility — ambiguous letter), `I_` (1, study), `II_` (1, study with two I's), `demo_` (1), and a `teaching/` folder. Key finding that complicates the obvious "numbered = benchmark, V_ = verification" split: **~34 of the 39 numbered cases already verify against a closed-form solution** (Mariotte, Lamé, slab/shell analytics, GPS, plasticity), so *verification is the norm, not a distinct category* — only a handful (`14_full_cylinder_cracking_2D_xy`, `15_single_elliptical_cavity_2D`, `18_box_crack_2D`, `18_box_knotch_2D`) are genuinely regression-only (no analytical truth). This means the new `V_` prefix overlaps in kind with the numbered cases rather than naming a new category. **Open design decision (Giovanni's call) before any rename:** either (a) treat analytical verification as the unmarked default — numbered progression *is* the verification suite, fold the 4 `V_` cases into the numbering, and only prefix the exceptions (demo / study / no-closed-form); or (b) split by *domain* — numbered = classical thermo-mechanical verifications, `V_` = the fuel-performance verifications (burnup/swelling/PCMI), the nuclear layer (maps onto the z3st-for-nuclear framing). Either way, pick a small fixed vocabulary (e.g. V/B/D/S) and document it in a new `cases/README.md`. **Objective nits to fix regardless of the taxonomy choice:** `knotch` → `notch` typo (`17_stress_strain_curve_knotch`, `18_box_knotch_2D`); hyphen → underscore in `19_single-edge_notched_*`; zero-pad the leading numbers (`1_` → `01_`) so they sort correctly (`1_,11_,12_,…,2_` currently sorts wrong). Any rename must sweep all references (`non-regression.sh`, `non-regression_github.sh`, `docs/source/physics_models.rst`, `main.tex`, `CONTEXT.md`) — the `U_→V_` rename on 2026-06-09 left dangling pointers in the CI script + docstrings + rst that had to be chased down.

### README

- [ ] **README-P2-1** — No public Sphinx URL listed. Either link `https://giozu.github.io/z3st/` (the URL mentioned in `main.tex:387`) or drop the "fully documented API" claim.
- [ ] **README-P2-2** — Cited BibTeX Zenodo DOIs are unverified manually.

### Paper

- [ ] **PAPER-P2-1** — 10 declared labels are never `\ref`d anywhere: `fig:case17, fig:case11, fig:case12, fig:case13, fig:case16, fig:case15_single, fig:case15_two, fig:plasticity_curve, fig:case9, fig:caseI, sec:cluster`. Either inline-reference each figure in its paragraph or drop the labels.

---

## Verified clean (no findings)

These were audited and found in order — recording so we don't re-audit:

- BC step-list lengths match `input.yaml::n_steps` exactly in `15_two_elliptical_cavities_2D` (10/10), `17_stress_strain_curve_double_crack` (5/5), `17_stress_strain_curve_knotch` (6/6), `17_stress_strain_curve_stress` (5/5), `20_plasticity_2D` (21/21), `demo_CP_single_grain` (41/41). The case-19 SENT failure mode (`701` BCs vs `n_steps: 141`) is not replicated elsewhere.
- All `Allrun` / `Allclean` files have the executable bit set.
- All `input.yaml::materials:` relative paths resolve to existing YAML cards.
- All 16 figures referenced in `main.tex` exist in `figures/`.
- All citation keys cited in `main.tex` are present in `references.bib`.
- Case-14 paper claims verified ✓ (other than [PAPER-P0-1](#paper-p0-1) and [PAPER-P0-5](#paper-p0-5)): `σc = 1 GPa`, `ℓc = 50 μm`, `h = 12.5 μm`, `T₀ = 1023.15 K`, `T_cold = 263.15 K`, `R = 10 mm`, 60° contact wedge, 100 steps × 1 ms, `E = 358 GPa`, `ν = 0.23`, `ρ = 10 970 kg/m³`, `Cp = 280 J/(kg K)`, `k = 5 W/(m K)`, `α = 1e-5 K⁻¹`, "~5 discrete radial cracks".

---

## Batch H — Literature-driven enhancements

Items spawned by the 2026-05-19 literature pass on `~/Bibliography - PFF/`
(Vicentini 2024 star-convex, Kumar 2020 / Kamarei 2024 strength critique,
Fajardo Lacave 2026 multi-cohesive). These are forward-looking improvements,
not bug fixes — orthogonal to the P0/P1/P2 tracks.

<a id="lit-1"></a>
- [x] **LIT-1 (code)** — *Resolved 2026-05-19; γ⋆ scope corrected 2026-05-19.* Added `psi_star_convex(u, material, T)` to `damage_model.py` implementing Vicentini et al. 2024 Eqs. (38)-(39): `ψ⁺ = G|dev ε|² + (λ/2)[⟨tr ε⟩₊² − γ⋆⟨tr ε⟩₋²]`, `ψ⁻ = (1+γ⋆)(λ/2)⟨tr ε⟩₋²`. Reads `γ⋆` from `self.dmg_cfg["gamma_star"]` (i.e. the `damage:` block of input.yaml, alongside `lc` and `hybrid_constraint`); defaults to 0 → reduces to Amor exactly, including the λ-not-κ convention. **γ⋆ is intentionally a model parameter, not a material property** — same scope as `lc`. Extended `psi_split` to honour an explicit `damage.split: amor | miehe | star_convex` config key; defaults preserve the historical AT1→Amor / AT2→Miehe pairing when the key is absent, so all existing cases are unaffected. CONTEXT.md §4.3 and README example YAML updated to document the new option. **Sanity check 2026-05-19**: case-14 2D-xy with `damage.split: star_convex` + (implicit) `gamma_star = 0` matches the Amor baseline — confirmed by the user. **Owed verification**: γ⋆ sweep over `{0.5, 1.0, 5.0}` on the same case — tracked as [CASES-FOLLOWUP-6](#cases-followup-6).

<a id="lit-2"></a>
- [ ] **LIT-2 (paper)** — Acknowledge the **Kumar 2020 / Kamarei 2024 nucleation critique** in `main.tex` §4 (Discussion). Per Kumar, Bourdin, Francfort, Lopez-Pamies (2020), *"Revisiting nucleation in the phase-field approach to brittle fracture"*, JMPS 142:104027, classical variational phase-field — including AT1 and AT2 — *cannot* predict crack nucleation independently of elasticity and toughness; it lacks the material strength as an independent ingredient. Kamarei, Dolbow, Lopez-Pamies (2024), J. Appl. Mech. 92:014502, makes this concrete in the first octant: AT1's predicted equi-biaxial strength `s_bs = √(3 Gc E / (16(1-ν) ℓ))` and hydrostatic strength `s_hs = √(Gc E / (8(1-2ν) ℓ))` depend unphysically on Poisson's ratio (`s_hs → ∞` as `ν → 1/2`). Add a single paragraph that (i) acknowledges this critique, (ii) frames the present case-14 contribution as a *performance demonstration* of the AT1 + Ambati-hybrid framework against McClenny's published crack topology, **not** as a contribution to the nucleation theory itself, (iii) flags Kumar's strength-based amendment as a future direction. Doesn't weaken the paper; strengthens its intellectual honesty.

<a id="lit-3"></a>
- [ ] **LIT-3 (paper)** — Expand §1 (Introduction) lit survey with **recent UO2 PFF work** that postdates McClenny: Gencturk et al. (2025), *"Thermo-Mechanical Phase-Field Modeling of Fracture in High-Burnup UO2 Fuels Under Transient Conditions"*, Materials 18(5):1162 (UO2 fracture as a *stochastic phase transition*); Xiong et al. (2024), *"Three-dimensional fracture of UO2 ceramic pellets by phase field modeling"*, Sci. China Phys. Mech. Astron. (cohesive-PFF for 3D UO2). One or two sentences situating z3st's case-14 alongside these contemporaries.

<a id="lit-4"></a>
- [ ] **LIT-4 (paper)** — Add **Fajardo Lacave, Vicentini, Welschinger, De Lorenzis (2026)**, *"A variational phase-field model for anisotropic fracture accounting for multiple cohesive lengths"*, JMPS 212:106585 — the newest paper in the bibliography — to the roadmap section of `main.tex`. Single bullet: future direction for grain-boundary-resolved polycrystalline UO2 with directional toughness. Z3st's current implementation does not handle anisotropic strength; this paper's multi-cohesive framework (one damage variable + directional cohesive lengths) would be the reference for any future extension.

<a id="lit-5"></a>
- [x] **LIT-5 (refs)** — *Resolved 2026-05-19 (content) — one small `main.tex` wiring step deferred.* User added `~/research-manuscripts/z3st_paper/z3st - PFF.bib` containing all the literature-driven references needed for LIT-1 through LIT-4: `vicentini_energy_2024`, `kumar_revisiting_2020`, `kamarei_nucleation_2024`, `xiong_three-dimensional_2024`, `gencturk_thermo-mechanical_2025`, `fajardo_lacave_variational_2026`, plus the foundational entries (Ambati 2015, Miehe 2015, Gerasimov 2019, McClenny 2022, MOOSE, BISON, JAX-FEM, OpenFOAM, Aagesen ×2, Simon 2024, Williamson ×2, Hales) matching keys already cited in `main.tex`. **Owed**: tiny `main.tex:27` update — `\addbibresource{../references.bib}` currently does not include the new file. Either (a) add `\addbibresource{"z3st - PFF.bib"}` (filename has a space; quotes for safety) or (b) rename to `z3st-PFF.bib` and add without quotes, or (c) merge the new entries into `../references.bib` and keep the existing single `\addbibresource`. Note: the new entries use **snake_case** keys (e.g. `vicentini_energy_2024`) while many existing citations in `main.tex` use **camelCase** (e.g. `ambatiReviewPhasefieldModels2015`) — when writing LIT-2/3/4 paragraphs, cite the new entries using their exact snake_case keys as declared in the new bib.

<a id="lit-6"></a>
- [ ] **LIT-6 (refs expansion)** — `z3st - PFF.bib` currently has 19 entries and covers the LIT-1/2/3/4 needs. Several PDFs in `~/Bibliography - PFF/` and a few foundational / very-recent papers aren't yet represented and should be added so future citations have a key to reach for:
  - **Foundational fracture mechanics**: Francfort & Marigo 1998 (variational fracture), Bourdin, Francfort & Marigo 2000/2008 (numerical regularization), Griffith 1921, Pham, Amor, Marigo, Maurini 2011 (gradient damage models), Sicsic & Marigo 2013 (sharp brittle limit), Borden et al. 2012 (mesh / regularization convergence in PFF).
  - **Original tension/compression splits** (z3st implements both; citing the originals would strengthen `damage_model.py`'s docstrings + the paper's §3.4.3): Amor, Marigo & Maurini 2009 (the Amor split paper — PDF already in `~/Bibliography - PFF/Amor2009_VariationalBrittleFracture.pdf`), Miehe, Welschinger & Hofacker 2010 (the original Miehe spectral split — PDF in `Miehe2010_CrackPropagation.pdf`).
  - **Gerasimov / De Lorenzis prior work** (PDFs already in the folder, not yet in `z3st - PFF.bib`): Gerasimov & De Lorenzis 2016 (monolithic Newton-Krylov), Gerasimov 2022 (second-order PFF). The 2019 penalization paper is already in the bib.
  - **Ductile / coupled extensions** (PDFs in the folder, candidates for the §4 roadmap or §1 lit-survey): Ambati, Gerasimov & De Lorenzis 2015 — the ductile PFF paper (`Ambati2015_PFFductile.pdf`), distinct from the 2015 review already cited; Sur 2025 — ductile + hydrogen (`Sur2025_DuctileFractureHyd.pdf`); Jiang 2020 — porosity-dependent intergranular fracture in UO2 (`Jiang2020_PFF_UO2.pdf`).
  - **Strength-based PFF amendment** (companion of LIT-2 and the latest in the Lopez-Pamies critique line): Lopez-Pamies, Dolbow, Francfort & Larsen 2025 — *"Classical variational phase-field models cannot predict fracture nucleation"*, CMAME 2025; Kumar et al. — *"A comprehensive macroscopic phase-field theory for nucleation and propagation of fracture"* (the strength-driving-force proposal). Neither PDF is in `~/Bibliography - PFF/` yet.
  - **ML-side** (the folder has one): Kiyani 2025 — DeepONet for crack nucleation (`Kiyani2025_PredictingCrackNucleation_DON.pdf`); add if the §4 roadmap mentions ML-aided PFF.
  - **Solver-side recent**: arXiv 2511.23064 (Nov 2025) — *"Iterative convergence in phase-field brittle fracture computations: exact line search is all you need"* — directly relevant to z3st's staggered loop.

  Acceptance criterion: each cited paper in `main.tex` has a matching key in the bib (no `??` in the rendered PDF). Use the same snake_case convention as the new entries when adding (e.g. `amor_regularized_2009`, `miehe_phase_2010`, `bourdin_numerical_2000`, `lopez-pamies_classical_2025`, `kumar_comprehensive_2025`).

## Batch I — Paper expansion

The current `main.tex` (~390 lines) describes most physics well but is light on
figures, tables, and concrete numbers, and several recent cases / refactors are
not represented. Items below scope the missing content. None of these are
blockers for the existing narrative; they are additions that strengthen the
paper's empirical and pedagogical content.

<a id="paper-exp-1"></a>
- [ ] **PAPER-EXP-1 — workflow schematic.** Add a figure to §2 (Methods) sketching the Spine driver: input YAMLs → Config / FE setup / mesh → per-physics mixins (Thermal / Mechanical / Damage / Plasticity / Gap / Cluster) → staggered iteration → OutputWriter. *User to draw; figure goes in `figures/workflow_schematic.png`*. One paragraph of caption + cross-reference from the existing §2.1 finite-element description.

<a id="paper-exp-2"></a>
- [ ] **PAPER-EXP-2 — mesh figures for representative cases.** Add a 2×2 or 3×2 panel showing the actual meshes that drive the verification: case 9 GPS axisymmetric, case 14 2D-xy (showing the contact-arc refinement + pre-crack seed), case 16 multi-material coaxial 3D, case 19 SENT/SENS, case 20 plasticity, `demo_CP_single_grain`. Captions should give node/element counts (the paper currently describes geometries only in text). Goes into `figures/case_meshes.png` (or multiple files).

<a id="paper-exp-3"></a>
- [ ] **PAPER-EXP-3 — crystal-plasticity demo content.** `demo_CP_single_grain` is currently a one-sentence mention at the end of §3.1.3 ("A demonstration of the `plasticity.mode: custom` hook..."). Promote it to a proper subsubsection with: (i) the single-grain geometry and orientation, (ii) the Schmid-factor verification (the case's `non-regression.py` already computes this — μ ≈ 0.408 for (111)[01̄1] slip under z-tension), (iii) the resolved-shear-stress vs accumulated-slip plot, (iv) a note that this exercises the `material.stress_function = "pkg.mod.func"` extensibility hook end-to-end.

<a id="paper-exp-4"></a>
- [ ] **PAPER-EXP-4 — table of cases.** Add a tabular catalogue of all suite cases (one row each): case name, regime, active physics, geometry, mesh size, wall-time on a reference machine, status. Mirrors the CONTEXT.md §6.1 catalogue but in a paper-suitable form. Could go in §3 (Results) as an opening overview before the per-case subsections, or as an appendix table.

<a id="paper-exp-5"></a>
- [ ] **PAPER-EXP-5 — performance / scaling table.** Add a small table reporting wall-time, peak memory, number of staggered iterations, and PETSc backend for ~6 representative cases (small linear, large 3D, transient damage, plasticity, cluster). Lets a reader gauge what the code actually costs. Goes in §4 (Discussion) or as an appendix.

<a id="paper-exp-6"></a>
- [ ] **PAPER-EXP-6 — OutputWriter section.** §2 (Methods) currently doesn't describe the I/O infrastructure. Add a short subsubsection on the unified `OutputWriter` (pre-allocated FE spaces, pre-compiled UFL Expressions, same field set in VTU and XDMF, the cumulative-plastic-strain projection fallback for `custom` plasticity). Half a page; situates the engineering choice that makes long damage runs tractable.

<a id="paper-exp-7"></a>
- [ ] **PAPER-EXP-7 — experimental-vs-numerical comparison for case 14.** Currently main.tex:321 *cites* McClenny's Fig. 8 ("two major radial cracks plus a fan of shorter surface cracks") but doesn't reproduce it alongside the z3st result. Add a side-by-side panel: z3st damage field at t = 0.1 s (already in `figures/case14_damage_field.png`) next to McClenny's experimental optical micrograph (their Fig. 8 top) or numerical crack field (their Fig. A.13/A.14). Permissions check needed before reproducing the McClenny figure; alternative is a re-traced schematic.

<a id="paper-exp-8"></a>
- [ ] **PAPER-EXP-8 — energy-balance and convergence plots.** Add two diagnostic plots for case 14: (i) `E_el(t)` and `E_frac(t)` from `energies.txt` over the full 100-step window (the current `case14_thermal_shock_results.png` shows T(r) profiles, not energies); (ii) the staggered convergence history (residuals(u), residuals(D)) read from `log_z3st.md` via `utils/plot_convergence.py`. These directly visualise the claims made in lines 321 and 366.

<a id="paper-exp-9"></a>
- [ ] **PAPER-EXP-9 — concrete numbers throughout.** The paper currently uses qualitative language in places where the implementation supplies concrete numbers. Examples: §2.5.1 says "adaptive relaxation" with no growth/shrink factors stated (defaults are 1.1 / 0.9 / 0.95 max / 0.1 min); §2.5.2 doesn't quote default rtol / stag_tol; §3.4.3 case-14 paragraph could state mesh size (currently 17 965 nodes per CONTEXT.md §9c), wall-time (311 s per CONTEXT.md §9c). A surgical sweep that promotes these from CONTEXT.md into the paper would help reproducibility.

## Batch J — June 2026 session (docs rendering + per-material γ-heating)

Work done in the 2026-06-05 session: a live-docs math-rendering fix, a docs
accuracy + figures overhaul, and a new per-material γ-attenuation reference
radius. Recorded here for traceability.

<a id="docs-1"></a>
- [x] **DOCS-1** — *Resolved 2026-06-05.* Math on the live Sphinx site (`giozu.github.io/z3st`) rendered as **raw LaTeX**. Root cause: `sphinx.ext.imgmath` was listed *after* `sphinx.ext.mathjax` in `docs/source/conf.py` and became the active math renderer; imgmath rasterises equations via `latex` + `dvipng`, which are absent in the Pages CI, so every formula fell back to source (and produced the `dvipng cannot be run` warning). Fix: removed `imgmath` from `extensions` and pinned `html_math_renderer = "mathjax"` (client-side, no toolchain). Verified locally — MathJax script injected, equations wrapped in `\(...\)`, dvipng warning gone. Deployed via PR #26. *Note:* the docs deploy is `main`-only (`static.yml`); `develop` work must reach `main` to publish.

<a id="docs-2"></a>
- [x] **DOCS-2** — *Resolved 2026-06-05.* `physics_models.rst` rewritten for accuracy against `main.tex`: hyperelasticity corrected from "(planned)" to **implemented**; a **crystal-plasticity** section added; the consistent **g(D) degradation of the thermal stress** documented (with the `1/K` argument); gap-conductance (gas correlation) and cluster-dynamics (advection–diffusion / DG1 / SIPG) equations aligned to the implementation. `examples.rst` gained the **SEN-shear (Miehe 2010)** and **UO₂ thermal-shock (McClenny 2022)** cases with result figures, and all inline math converted to proper `:math:` directives (also fixed a malformed `\dot{\gamma}` accent). New figures tracked under `docs/source/images/` via a `.gitignore` negation; `conf.py` excludes stray `images/**/*.md` from the build. Relates to [README-P2-1](#readme-p2-1) (public Sphinx URL now demonstrably useful).

<a id="code-feature-3"></a>
- [ ] **CODE-FEATURE-3** — *Added 2026-06-10 (Giovanni's request).* **Coolant heat-transfer module.** Today the coolant is a Dirichlet temperature on the clad outer surface (`U_pwr_rod_2D`: 580 K), i.e. an infinite heat-transfer coefficient — no convective film drop (~15–25 K at 25 kW/m with PWR-typical h ≈ 30–50 kW/m²K), so clad/gap/fuel temperatures are all biased low. Scope, in two stages:
  - **Stage 1 — convective film:** a `models/coolant_model.py` (or an extension of the thermal Robin BC) that computes `h_conv` from a forced-convection correlation (Dittus–Boelter / Gnielinski; Todreas & Kazimi Ch. 10) given coolant properties (a `h2o.yaml`-style card: ρ, cp, μ, k, plus mass flux and hydraulic diameter from the case YAML), and applies `Robin(h_conv, T_bulk)` on the labelled clad-outer surface. The existing `Robin: h_conv + T_ext` machinery is the natural plug-in point — stage 1 is essentially "compute h instead of hard-coding it".
  - **Stage 2 — axial coolant channel:** a 1D axial energy balance `ṁ·cp·dT_bulk/dz = q'(z)·per` coupled to the clad-outer flux, giving `T_bulk(z)` along the rod (the classic 1.5D thermo-hydraulic boundary; enables multi-slice / full-rod cases). Two-phase / subcooled-boiling regimes out of scope for now.
  - **Interim:** `U_pwr_rod_2D` keeps the Dirichlet coolant BC until stage 1 lands; the fill-gas pressure fix (2026-06-10, see [CASES-P1-7](#cases-p1-7)) is independent and already applied. Fits the z3st-for-nuclear track — this is the missing boundary block between Z3ST and a FRAPCON-class rod problem.

<a id="code-feature-2"></a>
- [ ] **CODE-FEATURE-2** — *Added 2026-06-05 (staged, not yet committed).* Per-material γ-attenuation reference radius. `spine.py::set_power` (the `f(x)` heating closure) now reads `gamma_inner_radius` from each material card (`mat.get("gamma_inner_radius", self.inner_radius)`) and normalises the cylindrical `K₀(μr)/K₀(μ·R_ref)` and spherical `(R_ref/r)·exp(−μ(r−R_ref))` profiles at that per-material surface instead of the single geometry `inner_radius`. Backward-compatible (defaults to `inner_radius` → existing cases unchanged). Motivation: layered γ-heating where an inboard layer (e.g. a thermal shield ahead of the pressure vessel) must be normalised at its own inner surface. Documented in `CONTEXT.md` §2.2/§5. **Owed:** no material card or case sets `gamma_inner_radius` yet — wire it into a layered shield+vessel attenuation case (likely a variant of `II_attenuation_map`) and add a non-regression check. Commit with a `feat(thermal):` message.

## Follow-ups (work spawned by other fixes)

<a id="code-feature-1"></a>
- [x] **CODE-FEATURE-1** — *Resolved 2026-05-20; relax_* propagation fixed 2026-05-20.* Implemented hot-reload of allow-listed `input.yaml` parameters: at the start of each time step, `__main__.py::_reload_hot_params` re-reads `input.yaml`, compares against the in-memory config dicts, and propagates any allowlisted changes in-place (shared by reference with `problem.dmg_cfg` / `mech_cfg` / `th_cfg` so the solver sees the new value on the next step). User can edit `input.yaml` mid-run to retune tolerances, iteration limits, relaxation factors, hybrid constraint, or γ⋆; changes apply at the next step boundary (latency ≤ one step). Robust to mid-edit reads (transient `yaml.YAMLError` / `FileNotFoundError` → silent skip). Documented in `CONTEXT.md` §2.3. Mesh / geometry / regime / model toggles / `damage.type`/`lc`/`split` / time history / `n_steps` are *not* hot-reloadable and edits to them are silently ignored mid-run. **Follow-up patch (same day):** `solver_settings.relax_*` keys were initially broken — they're cached as `self.relax_*` plain attributes on the Spine at `Solver.__init__` (not re-read from the dict per step), so the dict-side update didn't reach the solver. Fixed by also `setattr(problem, key, caster(new_val))` for solver_settings keys, with a `_SOLVER_SETTINGS_CASTS` map matching the type cast in `Solver.__init__`. The fix takes effect on the next `python -m z3st` launch; running processes that imported the old `__main__.py` are unaffected and need to be restarted to pick it up.

<a id="cases-followup-6"></a>
- [ ] **CASES-FOLLOWUP-6** — Verify the star-convex split ([LIT-1](#lit-1)) on `cases/14_full_cylinder_cracking_2D_xy/`. (i) ✅ Baseline (Amor) preserved as `output_amor_baseline/`. (ii) ✅ Sanity: `damage.split: star_convex` with `damage.gamma_star = 0` (or absent) matches the Amor baseline — confirmed 2026-05-19. (iii) ⏳ Sweep `damage.gamma_star ∈ {0.5, 1.0, 5.0}` (set in the `damage:` block of `input.yaml`, *not* in `uo2.yaml` — γ⋆ is a model parameter, like `lc`). For each run, record damage field, crack count, E_el / E_frac trajectories. Side-by-side comparison vs. Amor + hybrid is a paper-worthy data point (feeds the LIT-2 narrative).

<a id="cases-followup-5"></a>
- [ ] **CASES-FOLLOWUP-5** *(list corrected 2026-06-10)* — Bless `output/non-regression_gold.json` for the cases currently without one: `15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`, `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch` (these 4 have a `non-regression.py` but no gold). `U_box_knotch_3D` and `U_slab_contact` were **deleted from the repo in f1bb70b** — dropped from this list. `U_coaxial_contact_2D` has **no `non-regression.py` at all** (not just no gold) — needs the script first; folds into [CASES-FOLLOWUP-1](#cases-followup-1). For each: `./Allclean && ./Allrun`, sanity-check the resulting `output/non-regression.json`, then `cp output/non-regression.json output/non-regression_gold.json`. After blessing, uncomment the matching entry in `cases/non-regression.sh` so the case rejoins the local suite. Excluded from this follow-up: the three `14_full_cylinder_cracking*` variants — those need the calibration argument settled first (see `[CONTEXT.md §9c]` and `[PAPER-P0-1]`).

<a id="cases-followup-4"></a>
- [ ] **CASES-FOLLOWUP-4** — Re-run `cases/20_plasticity_2D/` and `cases/demo_CP_single_grain/` after the **[CODE-P1-11](#code-p1-11)** quadrature-degree bump (`mech_degree+1` → `2·mech_degree+1`). Plastic-strain integration may shift in the 4th–6th significant figure. If the non-regression script reports FAIL on either case, inspect — the new run is the more accurate one, so re-bless the gold (`cp output/non-regression.json output/non-regression_gold.json`) once you've sanity-checked the numbers. Standard cases (mechanical/thermal without plasticity) are unaffected (`q_degree` is unused outside the plasticity quadrature path).

<a id="code-followup-1"></a>
- [ ] **CODE-FOLLOWUP-1** — Extend thermal Dirichlet BCs to step-dependent temperature lists, mirroring the `raw` mechanism in `mechanical_model.py`. Requires (i) storing the list on `self.dirichlet_thermal[label]` analogously to mechanical, (ii) updating the `Constant` value per step inside `_thermal_step` (similar to how `_mechanical_step` currently calls `traction_const.value = ...` per step), (iii) validating list length against `self.n_steps`. Use case: time-varying temperature ramps on a boundary (e.g. annealing schedules).

<a id="utils-followup-2"></a>
- [x] **UTILS-FOLLOWUP-2** — *Resolved 2026-05-19.* Unified output backend: new `z3st/utils/writer.py::OutputWriter`. Both backends (VTU, XDMF) share the same field set (T, u, strain, per-material stress + von Mises + hydrostatic + strain-energy density, per-material heat flux, damage + crack-driving-force, cluster density, cumulative plastic strain). FE function spaces, Function targets, and `dolfinx.fem.Expression` objects are all pre-allocated / pre-compiled in `__init__`; per-step `write()` is just `interpolate` + I/O — no UFL JIT compilation in the time loop. `__main__.py` output section shrank from ~70 LoC to ~10 LoC of orchestration. `export_vtu.py` left intact for backward compat (no longer imported from `__main__.py`); user may deprecate in a follow-up. Per-format filename defaults now correct (was the `"fields.xdmf"` dead default for VTU mode). Patched 2026-05-19 to add a graceful interpolation-vs-projection fallback (`_make_interp_or_proj`): Expression construction is tried first; on `ValueError: Mismatch of tabulation points and element points.` (quadrature-sourced UFL, e.g. `plasticity.mode: custom` in `demo_CP_single_grain`), the writer builds a cached `LinearProblem` for L2 projection with the right quadrature metadata and dispatches per-field at write time. Confirmed working on `demo_CP_single_grain` (the only case that exercises the fallback path); standard cases stay on the fast interpolation path.

<a id="utils-followup-3"></a>
- [ ] **UTILS-FOLLOWUP-3** — Revisit VTKHDF as a third writer backend once dolfinx exposes a class-based `VTKHDFFile` with per-field naming. The dolfinx 0.10 `dolfinx.io.vtkhdf` submodule has a functional API (`write_mesh / write_point_data / write_cell_data`) but the write_* wrappers take no `name` parameter (the underlying nanobind binding `write_vtkhdf_data`'s second `str` argument is the data-location selector "PointData"/"CellData", not a field name). z3st's multi-field output (T, u, strain, per-material stress, von Mises, hydrostatic, heat flux, damage, cluster) is therefore not currently expressible in a single VTKHDF file without a forest of one-field-per-file workarounds. A prototype branch was added then reverted on 2026-05-19; see `writer.py` module docstring.

<a id="utils-followup-1"></a>
- [x] **UTILS-FOLLOWUP-1** — *Resolved 2026-05-19.* Convergence plot was blank because the markdown stdout filter in `__main__.py` rewrites `[STEP NN/MM]` → `## Step NN/MM:` and `--- Staggering iteration N/M ---` → `#### Iteration N/M` whenever stdout is redirected (non-TTY), but `plot_convergence.py:35,42` regex hunted for the pre-filter markers. Patched `plot_convergence.py` to normalize either form back to the raw form at parse-time, and changed default filename `log.z3st` → `log_z3st.md`. Companion change: standardized log-file naming across the suite. Allrun now does `gmsh mesh.geo -<dim> > log_mesh.md` and `python -m z3st > log_z3st.md` in every case (43 cases edited by an audit agent; case 19 + the three case-14 cracking variants were already correct or skipped; `II_attenuation_map` and `I_mesh_sensitivity_2D` have no per-case Allrun/z3st invocation). Allclean uniformly sweeps `rm -f log_*.md`. README.md and CONTEXT.md updated to match.

<a id="cases-followup-1"></a>
- [ ] **CASES-FOLLOWUP-1** — End-to-end verification of the **[CODE-P0-1](#code-p0-1)** fix in the Gas-gap path. Touch `cases/U_coaxial_contact_2D/`: add `gap_conductance: { type: Gas, value: <prefactor> }` under `models:` in `input.yaml`; introduce a real gap (`inner_radius_2: 0.056` in `geometry.yaml` and the mirroring change in `mesh.geo`); re-mesh and run. Also write a `non-regression_gold.json` once a good run exists (folds into [CASES-P1-6](#cases-p1-6)).

<a id="cases-followup-2"></a>
- [ ] **CASES-FOLLOWUP-2** *(rescoped 2026-06-10; decision recorded 2026-06-10)* — The `plane_stress` regime has **zero exercising cases**: `U_thick_cylindrical_shell_plane_stress` was deleted in f1bb70b, so the verification owed for [CODE-P0-2](#code-p0-2) (and the 2026-06-10 [CODE-P0-5](#code-p0-5) fix) has no vehicle. **Giovanni's decision (2026-06-10): keep the regime — do NOT deprecate** ("maybe one day we will do something on that"). The regime stays fully plumbed (config validation, `epsilon`, `sigma_mech` λ_ps, BC setters, solver branches); remaining work, unscheduled: recreate a plane-stress Lamé pressure-vessel case (analytical reference in `utils_verification.py`) and optionally a damage-enabled variant for the `g(D)` path.

<a id="cases-followup-3"></a>
- [x] **CASES-FOLLOWUP-3** — *Resolved 2026-05-19.* Case 9 reverted to small-strain Lamé to match the analytical GPS reference; added `mechanical.linear_solver: direct_mumps` (after the framework KeyError'd on the missing key — tracked separately as [CODE-P1-13](#code-p1-13)). All six metrics now PASS vs analytical (σ_zz error 0.315% ↓ from 1.92%, well inside 0.5% tolerance). Gold re-blessed.

## Cross-references

- **[CODE-P1-5](#code-p1-5)** ↔ case-19 SENT failure on 2026-05-18 (same shape-of-bug: list length not validated up front).
- **[PAPER-P0-1](#paper-p0-1)** ↔ `CONTEXT.md:650` ↔ `core/spine.py:147-149` (AT1 `Gc` identity must be the single source of truth).
- **[PAPER-P0-5](#paper-p0-5)** ↔ `CONTEXT.md:684` (case-14 energy numbers).
- **[CODE-P0-2](#code-p0-2)** ↔ **[CODE-P0-5](#code-p0-5)** (both `plane_stress` regime gaps).
- **[CASES-P0-1](#cases-p0-1)**, **[CASES-P0-2](#cases-p0-2)**, **[CASES-P0-3](#cases-p0-3)** all show the same file-absence-fallback design issue in `non-regression.sh` — fix the script's fallback to be loud, then the three cases become individually addressable.
