# Phase-field calibration envelope for the UO₂ cavity / bubble cases

These cases do **not** target a single fitted `(G_c, ℓ)` pair. Too many effects are
still unmodelled — no burnup dependence, no grain-boundary vs transgranular
distinction, no temperature dependence of the fracture properties, no porosity
history — so a single "final" number would be false precision. Instead we fix a
**physically defensible interval on σ_c**, let every case derive its own `G_c`
from its own `ℓ`, and require results to be reported with the interval.

## The one knob: σ_c

`σ_c`, `G_c` and `ℓ` are not independent. Given a formulation and a mesh, fixing
two fixes the third:

$$
\sigma_c^{\mathrm{AT1}} = \sqrt{\frac{3E G_c}{8\ell}}, \qquad
\sigma_c^{\mathrm{AT2}} = \sqrt{\frac{27 E G_c}{256\ell}}
$$

`spine.py::load_materials` (see `z3st/core/spine.py:208`) implements this: when a
material card declares **`sigma_c`** and the case declares **`lc`**, `G_c` is
derived from the case's own formulation (AT1 or AT2) and **overrides** any `Gc`
written on the card. So `σ_c` is the only quantity we calibrate; `G_c` follows.

**Adopted envelope:**

| Quantity | Value | Basis |
|---|---|---|
| `σ_c` | **150 MPa nominal, [100, 200] MPa envelope** | UO₂ tensile strength; the acceptance criterion in the internship plan |
| `E` | 358–385 GPa | 358 = z3st `uo2.yaml`; 385 = Jiang 2020 |
| `ℓ` | per case, subject to `ℓ ≥ 4h` | mesh resolution floor |
| `G_c` | **derived, never set by hand** | AT1/AT2 identity above |

Report every headline number (rupture stress, critical cracking pressure) as the
value at σ_c = 150 MPa **plus the spread over [100, 200] MPa**. A result that
only holds at one σ_c is not yet a result.

## Coherence of the three existing sources

The internship update deck and the material cards looked contradictory. They are
not — two of the three agree exactly, and the third was never a σ_c calibration.

| Source | Model | E | G_c | ℓ | ⇒ implied σ_c |
|---|---|---|---|---|---|
| Deck, case 14 (thermal shock) | AT1 | 358 GPa | 8.38 J/m² | 50 µm | **150.0 MPa** ✅ |
| Deck, case 19 (SENT) | AT1 | 358 GPa | 0.67 J/m² | 4 µm | **150.0 MPa** ✅ |
| `uo2_jiang.yaml`, as written | AT2 | 385 GPa | 2.0 J/m² | 0.5 µm | 403 MPa ⚠️ |
| `uo2_jiang.yaml`, as written | AT1 | 385 GPa | 2.0 J/m² | 2 µm | 380 MPa ⚠️ |

The deck's two calibrations are the *same* calibration at two mesh resolutions —
σ_c = 150 MPa exactly, in both. There is no inconsistency there.

`Gc: 2.0 J/m²` in `uo2_jiang.yaml` is Jiang's **grain-boundary fracture energy**,
a material datum — not a phase-field regularisation choice. Combined with the ℓ
our meshes can afford, it implies a σ_c of 380–400 MPa, i.e. 2.5× UO₂'s real
tensile strength: those cases would crack far too late. Reconciling it requires
an ℓ we cannot mesh:

| To reach σ_c with G_c = 2.0 J/m² | AT2 | AT1 |
|---|---|---|
| σ_c = 100 MPa | ℓ = 8.12 µm | ℓ = 28.9 µm |
| σ_c = 150 MPa | ℓ = 3.61 µm | ℓ = 12.8 µm |
| σ_c = 200 MPa | ℓ = 2.03 µm | ℓ = 7.2 µm |

An ℓ of 3.6–12.8 µm against a cavity semi-axis of 11.2 µm would smear the crack
over the entire defect — the geometry would stop being resolved. This is the same
trap as the case-14 / McClenny `G_c ≈ 80 kJ/m²` story: **a literature G_c is not
transferable to a phase-field model at a different ℓ.**

Resolution: `uo2_jiang.yaml` now declares `sigma_c`, keeps Jiang's `Gc: 2.0` as a
documented comment, and lets each case derive its own G_c.

### Two studies, two cards (decided 2026-08-03)

There are two calibrations in circulation both informally called "case 14". They are
not in conflict -- they belong to different studies and now live on different cards:

| Card | sigma_c | Used by | Study |
|---|---|---|---|
| `uo2_jiang.yaml` | **150 MPa** (envelope [100, 200]) | the five cavity / bubble cases | isothermal micromechanics, after Jiang 2020 |
| `uo2.yaml` (shared) | 1.0 GPa | `pellet_quench_2D_xy`, `pellet_quench_3D` | thermal-shock pellet cracking, after McClenny 2022 |

**150 MPa is canonical for the cavity/bubble study** -- it is UO2's physical tensile
strength. The 1 GPa on the shared card is a calibration choice specific to the
thermal-shock geometry (it sits below the cold-rim driving energy and above the bulk)
and must not be carried over. Neither number is "the" UO2 value; each belongs to its
own case family.

## Formulation: why AT2 here and AT1 there (decided 2026-08-03)

The five cases deliberately run **two** formulations. This is not drift:

| Case | type | split | why |
|---|---|---|---|
| `elliptical_cavity_tension_2D` | AT2 | miehe | reproduces Jiang 2020; AT2 has no elastic threshold, so damage onset is smooth and the full stress-strain curve is resolved |
| `elliptical_cavity_pressurized_2D` | AT2 | miehe | same geometry as above, load case swapped -- must share its formulation to be comparable |
| `two_elliptical_cavities_2D` | AT2 | miehe | cavity-interaction study in the same family |
| `bubble_fracture_2D` | AT1 | amor | AT1's sharp elastic threshold is the right model for asking *does this crack at all* below a bound |
| `spherical_void_tension_3D` | AT1 | amor | 3D RVE; AT1 is better behaved in bulk/compression |

`damage_model.py:328` defaults AT2 to Miehe and AT1 to Amor. Two cases previously
relied on that default; the split is now written out explicitly in every case so the
formulation is readable from `input.yaml` alone.

Note this differs from the AT1 + `star_convex` path used by `pellet_quench_2D_xy`,
which is a thermal-shock case with different requirements. Cross-comparisons between
the cavity cases and `pellet_quench` are therefore *not* like-for-like.

### One `lc` for the 2D family

All four 2D cases run **`lc = 0.5 um`**, so each derives the same `G_c` from `sigma_c`
and cavity-shape / load-case comparisons are like-for-like. `bubble_fracture_2D` was
moved from 2.0 to 0.5 um for this, which required refining `h_cavity` 0.5 -> 0.125 um
to keep `lc >= 4h`.

`spherical_void_tension_3D` keeps **`lc = 1.0 um`** as a documented exception: matching
0.5 um would need `h_sph <= 0.125 um`, roughly 8x the elements in 3D. Its `G_c` differs
from the 2D family accordingly, so its absolute strength is not directly comparable
with them -- only its trend is.

### Temperature

`uo2_jiang.yaml` sets `T_ref = T_initial = 293.15 K`. The thermal eigenstrain is then
zero **by construction**, which makes these cases structurally immune to the
`T_ref != T_initial` defect Baptiste diagnosed in May 2026, and matches Jiang's
isothermal setup. The shared `uo2.yaml` deliberately keeps `T_ref = 298.15` with
`T_initial = 1023.15` because that mismatch *is* the quench eigenstrain that drives
`pellet_quench` -- the same card cannot serve both intents, which is why these cases
have their own.

## Per-case derived G_c across the envelope

Each case's own `(model, ℓ)` with E = 385 GPa. The four 2D cases share ℓ = 0.5 µm, so
the three AT2 cases share a G_c exactly; `bubble_fracture_2D` differs only because AT1
and AT2 imply different G_c for the same σ_c.

| Case | Model | ℓ | G_c @ 100 MPa | G_c @ 150 MPa | G_c @ 200 MPa |
|---|---|---|---|---|---|
| `elliptical_cavity_tension_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `elliptical_cavity_pressurized_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `two_elliptical_cavities_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `bubble_fracture_2D` | AT1 | 0.5 µm | 0.035 | **0.078** | 0.139 J/m² |
| `spherical_void_tension_3D` | AT1 | 1.0 µm | 0.069 | **0.156** | 0.277 J/m² |

## Mesh floor — all cases now comply

`ℓ ≥ 4h` at the cavity boundary, so the damage band is resolved:

| Case | h_cavity | 4h | ℓ | |
|---|---|---|---|---|
| `elliptical_cavity_tension_2D` | 0.050 µm | 0.20 µm | 0.5 µm | OK |
| `elliptical_cavity_pressurized_2D` | 0.050 µm | 0.20 µm | 0.5 µm | OK |
| `two_elliptical_cavities_2D` | 0.125 µm | 0.50 µm | 0.5 µm | OK |
| `bubble_fracture_2D` | 0.125 µm | 0.50 µm | 0.5 µm | OK |
| `spherical_void_tension_3D` | 0.250 µm | 1.00 µm | 1.0 µm | OK |

Three of these sit exactly at `ℓ = 4h`, which satisfies the rule with no margin. A
convergence check at `ℓ = 6h` or `8h` is worth running before blessing anything whose
value moved — see the note on `two_elliptical_cavities_2D` below.

## Measured results (run 2026-08-03, sigma_c = 150 MPa)

| Case | Metric | Value | Baptiste's earlier figure |
|---|---|---|---|
| `elliptical_cavity_pressurized_2D` | critical cracking pressure | **179.06 MPa** (step 294, max D = 0.586) | 177 MPa — **1.2% apart** |
| `spherical_void_tension_3D` | peak macroscopic σ_zz | **116.6 MPa** at E_zz = 3.40e-4, softening → 0.0005 | (case had no metric) |
| `elliptical_cavity_tension_2D` | remote stress at initiation | **59.03 MPa** (peak σ_yy 82.66, informational) | 188 MPa — **reconciled**: local tip stress, see below |
| `two_elliptical_cavities_2D` | peak macroscopic σ_yy | **53.65 MPa**, max D = 1.0 | **116.52 MPa** |
| `bubble_fracture_2D` | max damage attained | **1.0000** — cracks at 77.4 MPa (D ≥ 0.5) | — |

### The mesh floor was inflating the two-cavity strength by 2.2x

`two_elliptical_cavities_2D` dropped from 116.52 to 53.65 MPa. Two things changed at
once — the mesh (`h_cavity` 0.15 → 0.125 µm, fixing the `ℓ ≥ 4h` violation) and the
material repoint (shared `uo2.yaml` → `uo2_jiang.yaml`, E 358 → 385 GPa with `G_c`
re-derived) — so the drop was initially ambiguous.

`elliptical_cavity_pressurized_2D` separates them. It took **only** the material
repoint; its mesh was already inside the floor (`h_cavity` = 0.05 µm → 4h = 0.2 µm vs
ℓ = 0.5 µm) and was not touched. It reproduces Baptiste's 177 MPa to within 1.2%. Both
cases run AT2 at the same `lc = 0.5 µm`, so the comparison is like-for-like.

**Conclusion: the material repoint is neutral; the entire two-cavity change is the mesh
floor.** Under-resolving the damage band makes the crack artificially hard to form, so
the specimen carries more load before failing. The original 116.52 MPa was a mesh
artifact, not a material property — which is precisely why `ℓ ≥ 4h` is a blessing
precondition and not a style preference.

Sanity of the remaining numbers at σ_c = 150 MPa: a macroscopic peak of 82.66 MPa for a
single cavity implies a stress-concentration factor of ~1.8, and the pressurised number
is a *pressure*, not a stress, so it is not required to track σ_c at all (see
`elliptical_cavity_pressurized_2D/compare_tip_stress.py`).

### ℓ is not just a numerical knob — it changed the bubble result

`bubble_fracture_2D` at ℓ = 2.0 µm (h = 0.5 µm) showed **`max(Damage)` ≡ 0** over the
whole 80 MPa ramp. At ℓ = 0.5 µm (h = 0.125 µm) the same case, same σ_c, same load,
**cracks completely**:

| ℓ | h_cavity | first damage | D ≥ 0.5 | max D |
|---|---|---|---|---|
| 2.0 µm | 0.5 µm | never | never | 0.0000 |
| 0.5 µm | 0.125 µm | 75.0 MPa | **77.4 MPa** | 1.0000 at 78.4 MPa |

σ_c is 150 MPa in both, so the *initiation stress* is unchanged by construction. What
changed is the resolution of the tip stress field: ℓ = 2 µm is comparable to the cavity
tip radius (ρ = ay²/ax ≈ 2.2 µm), so it averaged the concentration away; ℓ = 0.5 µm
resolves it. **For a stress-concentration problem ℓ acts as the averaging length of the
Point/Line Method — it sets the answer, not just the crack width.** This is the same
effect Baptiste identified for the tip-stress probe in `compare_tip_stress.py`.

Consequence: the earlier finding "fission-gas pressure alone does not crack the matrix
in the realistic range" was an ℓ artefact and does **not** hold. The bubble cracks at
77.4 MPa — between high-burnup (~50 MPa) and LOCA (~100 MPa) pressures.

Recorded as-is for Baptiste to take forward. Note the crossing sits only 3% below the
80 MPa ramp ceiling, so confirming it with a wider ramp should come before any blessing.

### The 188 MPa is reconciled — it is a local tip stress

`compare_tip_stress.py`, run 2026-08-03 against both cases' fresh output, probing each
case at *its own* max(Damage) location at standoffs in multiples of ℓ:

| | tension | pressure |
|---|---|---|
| critical step | 32 | 294 |
| applied load at initiation | 9.6 nm on ymax | 179.06 MPa cavity pressure |
| initiation point | (11.200, 0.000) µm, **0.0°** | (6.608, 3.471) µm, **27.7°** |
| remote σ_yy | 59.03 MPa | n/a (pressure-driven) |
| local σ_yy @ 0.5 ℓ | **189.99 MPa** | 92.57 MPa |
| local σ_yy @ 1.0 ℓ | 173.31 MPa | 94.15 MPa |
| raw nodal value | 302.69 MPa (singular) | 93.87 MPa |

**Baptiste's 188 MPa is the local σ_yy at roughly a half-ℓ standoff from the tip
(189.99 MPa), not a remote stress.** With the remote initiation stress at 59.03 MPa that
is a stress-concentration factor of 3.22. The two figures describe the same event
measured in different places; neither is wrong, and both should be quoted with their
location.

Two of his claims are confirmed exactly:

- tension initiates **at the major-axis tip** (0.0°), as Inglis predicts;
- internal pressure initiates **~28° off the major-axis tip** (27.7° measured) — pure
  pressurisation of an elongated cavity shifts the tensile concentration away from the
  tip. This is a real and non-obvious result.

His corner-singularity caution also holds: the raw nodal value (302.69 MPa) is 60% above
the half-ℓ value over a distance of 250 nm, which is why the standoff (Point/Line) method
is the right probe and a single nodal reading is not.

Open question for Baptiste: the two load cases initiate at very different local σ_yy
(190 vs 93 MPa), so "the same local stress triggers fracture regardless of load path"
does **not** hold for σ_yy. Von Mises at 1 ℓ is closer (99.7 vs 110.7 MPa, 11% apart),
but neither is the AT2/Miehe driving force, which keys off positive principal strains.
The right comparator is ψ⁺ itself.

## Case status

| Case | Verdict mechanism | Note |
|---|---|---|
| `elliptical_cavity_tension_2D` | `pass_fail_check` present, reference = 0.0 | reference is a placeholder, not a gold |
| `elliptical_cavity_pressurized_2D` | critical cracking pressure at max(D) ≥ 0.5 | |
| `two_elliptical_cavities_2D` | `pass_fail_check` present | blocked on the mesh floor |
| `bubble_fracture_2D` | `pass_fail_check` present | |
| `spherical_void_tension_3D` | **none — emits no `non-regression.json`** | see below |

`spherical_void_tension_3D` was renamed from `spherical_cavity_pressurized`: per
Jiang's methodology ("the pressure inside bubbles is not considered") the internal
Neumann pressure was removed, and `make_bcs.py` now applies only symmetry clamps
plus a `Dirichlet_z` displacement ramp on `zmax`. It is a remote-tension RVE with
an unpressurised void — the old name described a load case that no longer exists.

Its post-processing was rewritten at the same time to extract the macroscopic
Σ_zz–E_zz curve, replacing the original analytic σ_rr/σ_tt/σ_vm check. The old
`non-regression_gold.json` (σ_rr/σ_tt/σ_vm) was orphaned by that rewrite and has
been deleted; the case's `non-regression.py` never calls `pass_fail_check`, so the
suite reports `MISSING non-regression.json` and the case cannot pass. **A metric
must be chosen before it can be blessed.**

## Parametric sweeps

`study_theta.py`, `study_jiang.py` and `study_pressure_sweep.py` run each sample in
an isolated copy of the case under `.sweep/` — see `z3st/utils/case_sweep.py`. They
never write the tracked `mesh.geo` / `geometry.yaml` / `input.yaml` /
`boundary_conditions.yaml`. Keep it that way: the in-place variant silently
corrupted `bubble_fracture_2D`'s angle parameterisation once already, by
overwriting the line that *computed* `ay` from `theta_deg` with a literal.

## Decisions taken 2026-08-03 (supervisor)

Recorded so they are not re-litigated. Each is implemented in the cases.

1. **Formulation** — keep AT2/Miehe for the three Jiang-reproduction cavity cases and
   AT1/Amor for the bubble and 3D cases; document the reason rather than unify. Splits
   are now explicit in every `input.yaml`.
2. **ℓ** — one value (0.5 µm) for the whole 2D family so `G_c` is shared and comparisons
   are like-for-like; 3D keeps 1.0 µm as a meshing-cost exception.
3. **Calibration** — σ_c = 150 MPa is canonical for this study, on its own card
   (`uo2_jiang.yaml`). `pellet_quench`'s 1 GPa stays on the shared card. Two studies,
   two cards.
4. **`elliptical_cavity_tension_2D` headline** — the guarded metric is the **remote
   reaction stress at crack initiation** (first step with max(D) ≥ 0.5), which pairs
   with the pressurised sibling's critical cracking pressure. Peak macroscopic σ_yy is
   reported for information only.
5. **`bubble_fracture_2D` scope** — keep the 80 MPa ceiling, which spans the physical
   range (BOL a few MPa, high burnup ~50, LOCA ~100). The metric is **max damage
   attained**, and the case's finding is a *bound*: fission-gas pressure alone does not
   crack the matrix in the realistic range, so coalescence or a far-field load is
   required. Zero damage is the expected result, not a failure.
6. **Temperature** — isothermal at 293.15 K (`T_ref = T_initial`), eigenstrain-free by
   construction, matching Jiang.
7. **Naming** — Baptiste's case renamed `elliptical_cavity_tension_2D` to state its load
   case and stop colliding with develop's elastic `regression/elliptical_cavity_2D`.

## Still open

- **Convergence margin.** Three cases sit exactly at `ℓ = 4h`. Given how far
  `two_elliptical_cavities_2D` moved when its floor was fixed (116.52 → 53.65 MPa), a
  check at `ℓ = 6h`/`8h` should precede blessing that case.
- **Blessing.** No golds written. `elliptical_cavity_pressurized_2D` (179.06 MPa,
  reproducing Baptiste's 177 MPa to 1.2%) is the strongest candidate.

## Rules

1. Never write `Gc` as an active key on a card these cases use — set `sigma_c`.
2. Never quote a rupture stress or cracking pressure without the σ_c it came from.
3. Changing `lc` silently changes `G_c`. Re-run, do not re-fit.
4. Any case violating `ℓ ≥ 4h` is not eligible for a blessed gold.
