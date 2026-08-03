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

## Per-case derived G_c across the envelope

Each case's own `(model, ℓ)` with E = 385 GPa:

| Case | Model | ℓ | G_c @ 100 MPa | G_c @ 150 MPa | G_c @ 200 MPa |
|---|---|---|---|---|---|
| `elliptical_cavity_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `elliptical_cavity_pressurized_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `two_elliptical_cavities_2D` | AT2 | 0.5 µm | 0.123 | **0.277** | 0.493 J/m² |
| `bubble_fracture_2D` | AT1 | 2.0 µm | 0.139 | **0.312** | 0.554 J/m² |
| `spherical_void_tension_3D` | AT1 | 0.5 µm | 0.035 | **0.078** | 0.139 J/m² |

## Mesh floor — two cases currently violate it

`ℓ ≥ 4h` at the cavity boundary, so the damage band is resolved:

| Case | h_cavity | 4h | ℓ | |
|---|---|---|---|---|
| `elliptical_cavity_2D` / `_pressurized` | 0.05 µm | 0.20 µm | 0.5 µm | OK |
| `two_elliptical_cavities_2D` | 0.15 µm | 0.60 µm | 0.5 µm | **VIOLATED** |
| `bubble_fracture_2D` | 0.5 µm | 2.0 µm | 2.0 µm | OK (at the limit) |
| `spherical_void_tension_3D` | 0.25 µm | 1.0 µm | 0.5 µm | **VIOLATED** |

`two_elliptical_cavities_2D` needs either `h_cavity ≤ 0.125 µm` or `lc ≥ 0.6 µm`
before its result (peak σ_yy = 116.5 MPa at ε = 7.17e-4) can be blessed.
`spherical_void_tension_3D` needs either `h_sph ≤ 0.125 µm` or `lc ≥ 1.0 µm`.

## Measured results (run 2026-08-03, sigma_c = 150 MPa)

| Case | Metric | Value | Baptiste's earlier figure |
|---|---|---|---|
| `elliptical_cavity_pressurized_2D` | critical cracking pressure | **179.06 MPa** (step 294, max D = 0.586) | 177 MPa — **1.2% apart** |
| `spherical_void_tension_3D` | peak macroscopic σ_zz | **116.6 MPa** at E_zz = 3.40e-4, softening → 0.0005 | (case had no metric) |
| `elliptical_cavity_2D` | peak macroscopic σ_yy | **82.66 MPa** | 188 MPa — *different quantity* |
| `two_elliptical_cavities_2D` | peak macroscopic σ_yy | **53.65 MPa**, max D = 1.0 | **116.52 MPa** |
| `bubble_fracture_2D` | critical cracking pressure | **none** — D ≡ 0 up to 80 MPa | — |

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

## Case status

| Case | Verdict mechanism | Note |
|---|---|---|
| `elliptical_cavity_2D` | `pass_fail_check` present, reference = 0.0 | reference is a placeholder, not a gold |
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

## Open — next steps for Baptiste

Two decisions are deliberately left open. Both are physics calls, not maintenance.

### 1. `bubble_fracture_2D` — what should the ramp answer?

As configured the case ramps internal cavity pressure to **80 MPa over 401 steps and
`max(Damage)` stays at exactly 0.0000** — correct AT1 behaviour below the elastic
threshold, not a numerical failure. The comparable `elliptical_cavity_pressurized_2D`
needs ~180–240 MPa to crack.

80 MPa is not an arbitrary ceiling: high-burnup bubble pressure is ~50 MPa and LOCA
transients reach ~100 MPa. So the case currently states *fission-gas pressure alone,
in the realistic range, does not crack the matrix under this calibration* — plausibly
the interesting result, implying coalescence or a superimposed far-field load is
needed to initiate fracture.

- **Option A — extend the ramp to ~250 MPa.** Yields a critical cracking pressure
  comparable with the elliptical case, but at pressures well above anything physical.
  Treat it as a material property, not an operating condition.
- **Option B — keep the 80 MPa ceiling** and record the negative result deliberately.
  Then the metric should become *max damage attained*, not *critical pressure*, and
  the case documents a bound rather than a crossing.

Doing both is reasonable: keep 80 MPa as the physical statement, add the extended ramp
as a parametric study via `study_pressure_sweep.py`.

Until one is chosen the case cannot be blessed — `non-regression.py` exits non-zero
rather than reporting a false PASS.

### 2. Which number is `elliptical_cavity_2D`'s headline?

(The mesh-vs-material confound noted here earlier is resolved — see
"Measured results" above. What remains is purely which quantity the case guards.)

The run gives **peak macroscopic σ_yy = 82.66 MPa**. Baptiste's notes quote **188 MPa**,
but that is a different quantity — the *remote reaction stress at crack initiation*
from `compare_tip_stress.py`, not the peak over the ramp. With σ_c = 150 MPa a
macroscopic peak of ~83 MPa implies a stress concentration factor of ~1.8 at the
cavity, which is sane. Confirm which one the case is meant to guard before blessing.

## Rules

1. Never write `Gc` as an active key on a card these cases use — set `sigma_c`.
2. Never quote a rupture stress or cracking pressure without the σ_c it came from.
3. Changing `lc` silently changes `G_c`. Re-run, do not re-fit.
4. Any case violating `ℓ ≥ 4h` is not eligible for a blessed gold.
