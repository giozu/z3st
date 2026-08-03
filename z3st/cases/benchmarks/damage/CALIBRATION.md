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
| `spherical_cavity_pressurized` | AT1 | 0.5 µm | 0.035 | **0.078** | 0.139 J/m² |

## Mesh floor — one case currently violates it

`ℓ ≥ 4h` at the cavity boundary, so the damage band is resolved:

| Case | h_cavity | 4h | ℓ | |
|---|---|---|---|---|
| `elliptical_cavity_2D` / `_pressurized` | 0.05 µm | 0.20 µm | 0.5 µm | OK |
| `two_elliptical_cavities_2D` | 0.15 µm | 0.60 µm | 0.5 µm | **VIOLATED** |
| `bubble_fracture_2D` | 0.5 µm | 2.0 µm | 2.0 µm | OK (at the limit) |

`two_elliptical_cavities_2D` needs either `h_cavity ≤ 0.125 µm` or `lc ≥ 0.6 µm`
before its result (peak σ_yy = 116.5 MPa at ε = 7.17e-4) can be blessed.

## Consistency check against measured results

At σ_c = 150 MPa the cases produce: 116.5 MPa peak macroscopic stress
(two cavities, remote tension), 188 MPa remote reaction stress and 177 MPa
internal cavity pressure at initiation (single lenticular cavity). All sit inside
or adjacent to the [100, 200] MPa envelope, which is the expected behaviour —
macroscopic rupture stress is reduced below σ_c by the stress concentration, and
the pressurised number is a *pressure*, not a stress, so it is not required to
match σ_c at all (see `elliptical_cavity_pressurized_2D/compare_tip_stress.py`).

## Rules

1. Never write `Gc` as an active key on a card these cases use — set `sigma_c`.
2. Never quote a rupture stress or cracking pressure without the σ_c it came from.
3. Changing `lc` silently changes `G_c`. Re-run, do not re-fit.
4. Any case violating `ℓ ≥ 4h` is not eligible for a blessed gold.
