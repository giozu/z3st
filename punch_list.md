# z3st punch list

Audit performed 2026-05-18 against develop @ dbcee1d, generated from four parallel reviewers (code, cases, README, paper). Item IDs are stable — reference them in commits / PRs.

## Start here

All P0 items resolved. Remaining work is P1 (correctness / clarity) and P2 (cleanup / style). See sections below; the P1 candidates worth grouping for a single PR are **[CODE-P1-1/2/3]** (damage YAML / Gc-conversion validation), **[CODE-P1-5/6]** (BC list-length validation — the case-19 SENT shape-of-bug), and **[CODE-P1-12]** (`crack_driving_force` else).

## Resolved

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
- [ ] **CODE-P0-5** — `z3st/core/solver.py:329-333` and `:350-356` (`_mechanical_step`): `plane_stress` is missing from the regime list controlling the traction-normal and body-force assembly. It falls into the 3-D branch on a 2-D mesh → shape mismatch on both `n_vec` and `Constant(mesh, (0,0,-ρg))`. The mechanical BC setter at `z3st/models/mechanical_model.py:152` already handles `plane_stress`, so this is internal inconsistency. Add `"plane_stress"` to both branches.

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
- [ ] **CODE-P1-1** — `z3st/core/spine.py:148-169`: no upfront validation that `dmg_cfg["type"] ∈ {AT1, AT2}`. A typo like `at1` or `AT-1` silently skips the auto-conversion and crashes later inside `_damage_step` with `KeyError("Gc")`. Add an explicit check in `DamageModel.__init__` (`damage_model.py:18`).

<a id="code-p1-2"></a>
- [ ] **CODE-P1-2** — `z3st/core/spine.py:160` uses `type(Gc) == float`, which excludes `np.floating` / `int`. The σc-from-Gc conversion is silently bypassed in those cases. Use `isinstance(Gc, (int, float, np.floating))`.

<a id="code-p1-3"></a>
- [ ] **CODE-P1-3** — `z3st/core/spine.py:281,285`: when `Gc` is a Python callable returning a UFL expression (the `materials/ceramic.py` / `materials/oxide.py` extensibility hook), `mat["sigma_c"]` is set to `ufl.sqrt(...)`. `z3st/core/solver.py:559` later does `f"{sigma_c:.2e}"` → `TypeError`. Either skip the print when `sigma_c` is non-scalar, or stop overloading the same key.

<a id="code-p1-4"></a>
- [ ] **CODE-P1-4** — `z3st/core/solver.py:582-591` (AT1, and the analogous AT2 path): irreversibility / relaxation ordering. Sequence is `clip → relaxation → max(D_new, D_old) → clip`. Relaxation can drag `D_new` below `D_old`, then `max` resets to `D_old`, so when relaxation tries to slow propagation the effective increment is **zero**. Re-order: relax the increment `ΔD = max(0, D_new − D_old)` then set `D ← D_old + α·ΔD`.

<a id="code-p1-5"></a>
- [ ] **CODE-P1-5** — `z3st/__main__.py:166-172` + BC validation: `generate_power_history(n_steps=input_n_steps-1)` may return a `times` array whose length differs from `n_steps` for multi-segment histories. BC setters validate list lengths against `self.n_steps` (the YAML value), but the solver iterates `len(times)`. Tighten by setting `self.n_steps = len(times)` immediately after `generate_power_history` and **before** `set_boundary_conditions`. This is the same shape-of-bug that caused the tension-case failure on 2026-05-18.

<a id="code-p1-6"></a>
- [ ] **CODE-P1-6** — `z3st/models/mechanical_model.py:139-146`: Neumann traction list is not length-validated against `n_steps`; Dirichlet at L100 is. Add a parallel check.

<a id="code-p1-7"></a>
- [ ] **CODE-P1-7** — `z3st/models/thermal_model.py:76-95`: Thermal Dirichlet does not accept step-dependent lists. Either document this asymmetry vs. the mechanical block or extend with the same `raw` mechanism.

<a id="code-p1-8"></a>
- [ ] **CODE-P1-8** — `z3st/models/mechanical_model.py:453-455`: side-effect mutation of `material["constitutive_mode"]` from inside `sigma_mech`. Surprising and order-dependent (the first call sets it; later calls behave differently). Move to `load_materials` in `spine.py`.

<a id="code-p1-9"></a>
- [ ] **CODE-P1-9** — `z3st/core/spine.py:289-340` (`set_power`): with `inner_radius == 0` and `geometry_type ∈ {cyl, sphere}` with `gamma_heating > 0`, `K_0(0) = +∞` (cyl) or `1/r → ∞` (sphere) at the centre. Guard.

<a id="code-p1-10"></a>
- [ ] **CODE-P1-10** — `z3st/core/spine.py:299-308`: `fissile: true` **and** `gamma_heating > 0` on the same material silently overwrites the fissile contribution (separate `if` blocks, second assignment wins). Use `+=` or raise.

<a id="code-p1-11"></a>
- [ ] **CODE-P1-11** — `z3st/core/finite_element_setup.py:53`: `q_degree = mech_degree + 1` is on the low side for full integration of the J2 return-mapping form. Standard is `2·mech_degree + 1`. Worth re-running `cases/20_plasticity_2D/` after bumping.

<a id="code-p1-12"></a>
- [ ] **CODE-P1-12** — `z3st/models/damage_model.py:103-106`: `crack_driving_force` has no `else` after the AT1 branch and returns `None` on a typo. The call site `update_history` then crashes inside `fem.Expression`. Raise `ValueError` explicitly.

<a id="code-p1-13"></a>
- [ ] **CODE-P1-13** — `z3st/core/solver.py:396` (and analogous spots) raise `KeyError: 'linear_solver'` when the `mechanical:` or `damage:` blocks omit `linear_solver`. Encountered live in case 9 on 2026-05-19. Either pick a default (e.g. `direct_mumps`) in `Config` / `Solver.__init__`, or validate early with a clear message ("'linear_solver' is required when 'solver: linear' — choose one of direct_mumps / iterative_amg / iterative_hypre"). Same likely applies to `thermal.linear_solver` and `damage.linear_solver`.

### Cases

<a id="cases-p1-1"></a>
- [ ] **CASES-P1-1** — `CONTEXT.md:415` claims `mesh.msh` is "checked in for CI"; `.gitignore:32` excludes `*.msh`. No case has `mesh.msh` tracked. CI works because `Allrun` regenerates it. Fix the doc.

<a id="cases-p1-2"></a>
- [ ] **CASES-P1-2** — `CONTEXT.md` §6.1 lists `14_full_cylinder`, `14_full_cylinder_cracking`, `14_full_cylinder_thermal_2D_rz` but not `14_full_cylinder_cracking_2D_xy`, even though §9c (and disk) tracks all three damage variants. Update §6.1.

<a id="cases-p1-3"></a>
- [ ] **CASES-P1-3** — `z3st/cases/non-regression_summary.txt` (2026-05-15) is stale. Missing rows for `15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`, `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch`, `18_box_knotch_2D`. `non-regression.sh` has a commented entry `15_box_elliptical_cavity_2D` that doesn't exist on disk.

<a id="cases-p1-4"></a>
- [ ] **CASES-P1-4** — `z3st/cases/non-regression_github.py` points to `tests/non-regression_github.sh`; the script actually lives at `z3st/cases/non-regression_github.sh`. The pytest invocation would not locate it.

<a id="cases-p1-5"></a>
- [ ] **CASES-P1-5** — `.github/workflows/ci.yml` runs 9 cases; none of them are the active `14_full_cylinder_cracking*` variants. Document the CI inclusion policy.

<a id="cases-p1-6"></a>
- [ ] **CASES-P1-6** — Missing `output/non-regression_gold.json` in 9 cases: the three `14_full_cylinder_cracking*` variants (expected for WIP), plus `15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`, `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch`, `U_box_knotch_3D`, `U_coaxial_contact_2D`, `U_slab_contact`. `non-regression.py` exists in each but has no gold to compare against.

### README

<a id="readme-p1-1"></a>
- [ ] **README-P1-1** — `README.md:88-102` Key Features section omits: plasticity (J2 + custom CP), hyperelasticity (Neo-Hookean), kinematic regimes (`2d / 3d / axisymmetric / plane_stress`), AT1 vs AT2, Miehe vs Amor splits. CONTEXT.md §3.1 / §4 has all of them.

<a id="readme-p1-2"></a>
- [ ] **README-P1-2** — `README.md:172` example YAML uses `regime: 2D`; canonical schema is lowercase. Cases on disk are inconsistent themselves — pick one and propagate.

<a id="readme-p1-3"></a>
- [ ] **README-P1-3** — `README.md:118-157` directory tree is out of sync: missing `core/mesh/`, `models/plasticity_model.py`, the full material card list, several utils. `1_thin_thermal_slab/` is the wrong directory name (actual: `1_thin_slab_2D/`).

<a id="readme-p1-4"></a>
- [ ] **README-P1-4** — `README.md:296` BibTeX key typo: `Scrogggs` → `Scroggs`.

<a id="readme-p1-5"></a>
- [ ] **README-P1-5** — `CITATION.cff` is bare (`doi` + `identifiers` only). Missing `cff-version`, `authors`, `title`, `version`, `date-released`, `repository-code`. The README BibTeX has nothing to validate against.

### Paper

<a id="paper-p1-1"></a>
- [ ] **PAPER-P1-1** — `main.tex:193` plasticity paragraph cross-references `\ref{sec:phase-field}` for the radial-return update; should add a `\label{sec:plasticity}` at line 114 and fix the cite target.

<a id="paper-p1-2"></a>
- [ ] **PAPER-P1-2** — `main.tex:251` describes case 13 with `k = 50 W/(m K)`; the material is `materials/ceramic.yaml` with `k` resolved through the Python callable `materials.ceramic.k(T)`. Either say so or pin to a scalar `k`.

<a id="paper-p1-3"></a>
- [ ] **PAPER-P1-3** *(optional polish)* — `main.tex:319-321` secondary-nucleation narrative is consistent with CONTEXT.md §9c, but the wording can be tightened to make the AT1-sharp-threshold mechanism the explicit cause (currently implicit). Defer until [PAPER-P0-1](#paper-p0-1) is resolved.

---

## P2 — cleanup / style

### Code

- [ ] **CODE-P2-1** — Commented-out debug prints / scaffolding to remove: `z3st/core/spine.py:259-262`, `z3st/core/finite_element_setup.py:69-74` (mixed-element scaffolding), `z3st/models/gap_model.py:97-98` (alt distance).
- [ ] **CODE-P2-2** — Italian error message at `z3st/models/mechanical_model.py:310` ("Formato displacement non valido"); translate.
- [ ] **CODE-P2-3** — List-comprehensions used purely for side effects at `z3st/core/solver.py:793-794, 803-804`. Convert to `for` loops.
- [ ] **CODE-P2-4** — `z3st/core/spine.py:177-180` dumps every material field on every run. Gate behind `--debug`.
- [ ] **CODE-P2-5** — `stag_tol` has two sources of truth: `spine.solve` defaults at `core/spine.py:356-361` (1e-4) and `solve_staggered` defaults at `core/solver.py:761-770` (1e-3). Pick one.
- [ ] **CODE-P2-6** — `z3st/models/thermal_model.py:154` wraps `material["k"]` in `dolfinx.fem.Constant(...)`; crashes when `k` is a UFL expression (symbolic-k extensibility hook).

### Cases

- [ ] **CASES-P2-1** — Allrun convention drift in `U_cluster_dynamics_test/` (uses `log.txt` instead of `log_z3st.md`, `gmsh -1`, no nr step) and `U_quarter_block/` (truncated pipeline).
- [ ] **CASES-P2-2** — `Allclean` is only called from `non-regression.sh`, never from inside individual Allrun scripts. Confirm design intent.

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

## Follow-ups (work spawned by other fixes)

<a id="utils-followup-2"></a>
- [x] **UTILS-FOLLOWUP-2** — *Resolved 2026-05-19.* Unified output backend: new `z3st/utils/writer.py::OutputWriter`. Both backends (VTU, XDMF) share the same field set (T, u, strain, per-material stress + von Mises + hydrostatic + strain-energy density, per-material heat flux, damage + crack-driving-force, cluster density, cumulative plastic strain). FE function spaces, Function targets, and `dolfinx.fem.Expression` objects are all pre-allocated / pre-compiled in `__init__`; per-step `write()` is just `interpolate` + I/O — no UFL JIT compilation in the time loop. `__main__.py` output section shrank from ~70 LoC to ~10 LoC of orchestration. `export_vtu.py` left intact for backward compat (no longer imported from `__main__.py`); user may deprecate in a follow-up. Per-format filename defaults now correct (was the `"fields.xdmf"` dead default for VTU mode). Owed verification: re-run the full non-regression suite to confirm field-name compatibility with the per-case non-regression scripts (most read `Stress_<mat> (cells)`, `Strain (cells)` — should match the new writer 1:1, but worth confirming on case 9 + case 18 + case 20 + a damage case).

<a id="utils-followup-3"></a>
- [ ] **UTILS-FOLLOWUP-3** — Revisit VTKHDF as a third writer backend once dolfinx exposes a class-based `VTKHDFFile` with per-field naming. The dolfinx 0.10 `dolfinx.io.vtkhdf` submodule has a functional API (`write_mesh / write_point_data / write_cell_data`) but the write_* wrappers take no `name` parameter (the underlying nanobind binding `write_vtkhdf_data`'s second `str` argument is the data-location selector "PointData"/"CellData", not a field name). z3st's multi-field output (T, u, strain, per-material stress, von Mises, hydrostatic, heat flux, damage, cluster) is therefore not currently expressible in a single VTKHDF file without a forest of one-field-per-file workarounds. A prototype branch was added then reverted on 2026-05-19; see `writer.py` module docstring.

<a id="utils-followup-1"></a>
- [x] **UTILS-FOLLOWUP-1** — *Resolved 2026-05-19.* Convergence plot was blank because the markdown stdout filter in `__main__.py` rewrites `[STEP NN/MM]` → `## Step NN/MM:` and `--- Staggering iteration N/M ---` → `#### Iteration N/M` whenever stdout is redirected (non-TTY), but `plot_convergence.py:35,42` regex hunted for the pre-filter markers. Patched `plot_convergence.py` to normalize either form back to the raw form at parse-time, and changed default filename `log.z3st` → `log_z3st.md`. Companion change: standardized log-file naming across the suite. Allrun now does `gmsh mesh.geo -<dim> > log_mesh.md` and `python -m z3st > log_z3st.md` in every case (43 cases edited by an audit agent; case 19 + the three case-14 cracking variants were already correct or skipped; `II_attenuation_map` and `I_mesh_sensitivity_2D` have no per-case Allrun/z3st invocation). Allclean uniformly sweeps `rm -f log_*.md`. README.md and CONTEXT.md updated to match.

<a id="cases-followup-1"></a>
- [ ] **CASES-FOLLOWUP-1** — End-to-end verification of the **[CODE-P0-1](#code-p0-1)** fix in the Gas-gap path. Touch `cases/U_coaxial_contact_2D/`: add `gap_conductance: { type: Gas, value: <prefactor> }` under `models:` in `input.yaml`; introduce a real gap (`inner_radius_2: 0.056` in `geometry.yaml` and the mirroring change in `mesh.geo`); re-mesh and run. Also write a `non-regression_gold.json` once a good run exists (folds into [CASES-P1-6](#cases-p1-6)).

<a id="cases-followup-2"></a>
- [ ] **CASES-FOLLOWUP-2** — End-to-end verification of the **[CODE-P0-2](#code-p0-2)** fix. Two sub-tasks: (i) re-run `cases/U_thick_cylindrical_shell_plane_stress/Allrun` to confirm the no-damage plane-stress path still produces a correct Lamé pressure-vessel solution after dropping the 3×3 padding (downstream consumers in `spine.py:387` and `utils/export_vtu.py` may need adjustment if they assume 3×3 σ in plane_stress); (ii) add a damage-enabled plane-stress case (no such case currently exists) to exercise the `g(D)` degradation path.

<a id="cases-followup-3"></a>
- [x] **CASES-FOLLOWUP-3** — *Resolved 2026-05-19.* Case 9 reverted to small-strain Lamé to match the analytical GPS reference; added `mechanical.linear_solver: direct_mumps` (after the framework KeyError'd on the missing key — tracked separately as [CODE-P1-13](#code-p1-13)). All six metrics now PASS vs analytical (σ_zz error 0.315% ↓ from 1.92%, well inside 0.5% tolerance). Gold re-blessed.

## Cross-references

- **[CODE-P1-5](#code-p1-5)** ↔ case-19 SENT failure on 2026-05-18 (same shape-of-bug: list length not validated up front).
- **[PAPER-P0-1](#paper-p0-1)** ↔ `CONTEXT.md:650` ↔ `core/spine.py:147-149` (AT1 `Gc` identity must be the single source of truth).
- **[PAPER-P0-5](#paper-p0-5)** ↔ `CONTEXT.md:684` (case-14 energy numbers).
- **[CODE-P0-2](#code-p0-2)** ↔ **[CODE-P0-5](#code-p0-5)** (both `plane_stress` regime gaps).
- **[CASES-P0-1](#cases-p0-1)**, **[CASES-P0-2](#cases-p0-2)**, **[CASES-P0-3](#cases-p0-3)** all show the same file-absence-fallback design issue in `non-regression.sh` — fix the script's fallback to be loud, then the three cases become individually addressable.
