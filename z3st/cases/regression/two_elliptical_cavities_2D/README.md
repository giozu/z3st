# two_elliptical_cavities_2D

Micromechanical model of **intergranular fracture driven by fission-gas bubbles
on a grain boundary** in oxide fuel. Two lenticular gas bubbles (elliptical
cavities) sit on a grain boundary (GB) at `y = 0`; internal bubble pressure both
concentrates stress at the bubble tips and weakens the load-bearing GB ligament.
An AT2 phase-field damage model, with a GB-weakened toughness `Gc(y)`, lets a
crack nucleate at the tips and link the bubbles along the GB.

The scientific target is the **GB-failure response vs the bubble fractional
coverage `Fc` and the GB toughness `Gc`** — relevant to fuel fragmentation and
fission-gas release.

Units: micron / second / kg / micronewton / MPa / picojoule (a consistent set:
`µN = kg·µm/s²`, `MPa = µN/µm²`, `pJ = µN·µm`; and `pJ/µm² = J/m²` numerically).

## Files

- `parametric_study.py` — **the main deliverable.** Two analytical `p_crit(Fc,
  Gc)` references (a fracture-mechanics SIF estimate, recommended, plus the
  superseded strength estimate, for contrast) and an FEM sweep that regenerates
  the mesh at each `Fc`, ramps the bubble pressure, detects GB percolation, and
  saves a per-`Fc` crack figure. Edit the **CONFIGURATION block** at the top
  (`CASES`, `GC_VALUES`, `CRACK_FACES`, `SIGMA_H`) to set the cases, then run all.
    - `python3 parametric_study.py`            → analytical reference map
    - `python3 parametric_study.py --fem`      → sweep the configured `CASES`
    - `python3 parametric_study.py --fem 0.3 0.5` → override the `Fc` list
- `non-regression.py` — per-run diagnostic: extracts the fields and writes
  `fields_overview.png`, `stress_profile_tip.png`, `gc_profile_check.png`.
- `mesh.geo` — parametrised by `Fc_target` (`gmsh -setnumber Fc_target 0.3 …`),
  refined to `h_cavity = lc/2` so the AT2 damage band can localise.
- `input.yaml` / `boundary_conditions.yaml` — the **committed config**: the
  equilibrium-pressure baseline (bubble pressure ramped to 15 MPa).

## Output state (IMPORTANT — read before trusting the plots)

Two distinct families of figures live in `output/` (all PNGs are gitignored):

- **Baseline diagnostics** — `fields_overview.png`, `stress_profile_tip.png`,
  `gc_profile_check.png`. These are the *committed 15 MPa no-fracture* state,
  written by `non-regression.py`. `./Allrun` (and the suite) regenerate them;
  `Allclean` removes them. Expect `d_max ≈ 0` here — the baseline does not crack.
- **Sweep deliverables** (FRACTURED, ~GPa loads) from `parametric_study.py --fem`,
  preserved by `Allclean`:
    - `pcrit_vs_Fc_sweep.png` — FEM `p_crit` points overlaid on the SIF map
    - `sweep_Fc0.20/0.40/0.60.png` — per-coverage Damage + `σ_yy` (crack
      morphology: isolated tips → bridging → fully linked)
    - `pcrit_vs_Fc_reference.png` — analytical map only (SIF vs strength)

## Analytical reference — strength vs fracture mechanics (SIF)

The original reference `p_crit = σ_c·(1−Fc)/K_t` (AT2 peak stress × tip
concentration / ligament) is **superseded but kept for contrast**. It mixes a
*pointwise* `K_t` with an `lc`-scale `σ_c` and predicts tip *nucleation*, not GB
*percolation*, so it sits **~17× below the FEM** and is systematically
unconservative. The recommended reference is now a fracture-mechanics SIF
(Chakraborty, Tonks & Pastore 2014 — see `paper/`):

      p_crit = K_Ic / (F(Fc)·√(π·R)) + σ_h,   K_Ic = √(E·Gc/(1−ν²)),   R = ax

with `F(Fc)` their non-dimensional Mode-I SIF (`F_sif`, Eqs. 5/9) and `σ_h` a
compressive hydrostatic restraint (Eqs. 6/8 — the external-restraint lever). For
this case: `K_Ic = 0.138 MPa·µm½`, giving `Fc 0.2/0.4/0.6 → 414/365/308 MPa` —
the **right trend and magnitude** (hundreds of MPa), closing the gap to the FEM
from ~17× to **~4–6×**.

The residual ~4–6× is **expected and physical**: Chakraborty assumes a *sharp*
pre-crack (LEFM), but the phase-field bubble tip is a *smooth* ellipse
(`ρ = ay²/ax ≈ 22 nm`) with finite `lc = 4 nm`. The material characteristic
length `ℓ_ch = K_Ic²/σ_c² ≈ 42 nm ≈ ρ`, so the case sits in the
**strength↔toughness transition** — neither pure-strength (old formula) nor
pure-LEFM (Chakraborty) is exact, and the FEM (which resolves a smooth-tip,
finite-`lc`, *percolating* crack) is rightly the highest of the three:
`strength < SIF-LEFM < phase-field-percolation`.

## OPEN DECISION (deferred) — the (Gc, lc) regime

With the current `Gc_gb = 0.1 pJ/µm²` and `lc = 4 nm`, `σ_c ≈ 670 MPa` and the
smooth-tip/finite-`lc` transition above push GB failure to **~1–2.4 GPa** bubble
pressure (validated: `Fc 0.2/0.4/0.6 → p_crit 2368/2053/1105 MPa`).

Equilibrium pressures for these ~0.1 µm bubbles are only **tens of MPa**
(`p_eq ≈ 2γ/r`), so the model fractures only under strong **transient
over-pressurisation** — which the literature supports as the actual trigger
(see references below), so the GPa-scale `p_crit` is qualitatively right, not a
sign the regime is wrong. To model GB cracking at realistic *equilibrium*
pressures instead, move `(Gc, lc)` to a softer regime so `σ_c` is reachable
there; the alternative is to keep the regime and model the over-pressure /
restraint (`σ_h`) state explicitly.

## Non-regression metrics (in the suite, with a gold)

The `non-regression.py` metrics describe the *committed equilibrium baseline*
(15 MPa, no fracture) and are now canonical:

- `max_damage` → reference **0** (the baseline must not crack; percolation is
  0.9). Measured `d_max ≈ 3.3e-5`.
- `max_stress_yy` → the bubble-tip concentration `K_t/(1−Fc)·p_applied`, with
  `p_applied` **read from the ramp** (15 MPa), not hard-coded. The isolated `K_t`
  underpredicts ~2× (two bubbles share load through the ligament); the
  ligament-corrected estimate is accurate to ~20% (`σ_yy ≈ 99` vs `82` MPa).
  Analytic tolerance is loose (0.25) on purpose — the precise drift guard is the
  **gold regression check** (rtol 1e-3).

The case **is in the local regression suite** (`Allrun` + a blessed
`output/non-regression_gold.json`; ~120 s, so kept out of the CI subset). Re-bless
with `cp output/non-regression.json output/non-regression_gold.json` after a
sanity-checked baseline run. `Allclean` preserves the gold and the sweep figures.

## Already fixed this round

- Non-regression metrics reconciled to the no-fracture baseline (`max_damage`
  ref 0, `max_stress_yy` ref = ligament-corrected tip concentration with
  `p_applied` read from the ramp) + a gold + `regression_check`; case joined the
  local suite. `Allclean` no longer wipes `output/` wholesale (gold preserved).
- Analytical reference upgraded: fracture-mechanics SIF (`p_crit_sif`,
  Chakraborty 2014) added alongside the superseded strength estimate; see the
  "Analytical reference" section.
- `parametric_study.py` made config-driven (edit the `CASES` block, run all).
- `Gc` reconciled: `materials/oxide.py` is the single source (`Gc(mesh)` UFL +
  `Gc_numpy(y)`); the script imports it (a stray `*1e-6` in a local copy is gone).
- Output field name corrected (`Stress_solid → Stress`).
- Mesh refined (`h_cavity = lc/2`) and parametrised by `Fc`.
- Cleaner diagnostic plots (field overview, GB stress profile, Gc, fracture map).

## Background literature (GB-bubble overpressurisation → fragmentation / FGR)

PDFs are in `paper/` (named `AuthorYear_Keyword.pdf`). Context for the open
`(Gc, lc)` decision: the literature does **not** expect GB
cracking at equilibrium pressure — it expects it under *transient
over-pressurisation*, modulated by external mechanical restraint. So the GPa-scale
`p_crit` this case finds is qualitatively consistent with over-pressure being the
trigger, not a sign the regime is wrong.

Closest analytical/FEM precedents to swap in for the current first-approximation
`σ_c·(1−Fc)/K_t` reference (which is ~17× below this case's FEM `p_crit`):

- Chakraborty et al. (2014), *J. Nucl. Mater.* — Mode-I non-dimensional **SIF for
  lenticular GB bubbles** under bubble pressure + hydrostatic stress (exact
  geometry of this case; the natural analytical reference).
  https://doi.org/10.1016/j.jnucmat.2014.04.023
- Cappellari et al. (2025), *J. Nucl. Mater.* — PoliMi/**SCIANTIX** physics-based
  GB FGR model applying fracture mechanics to bubble-overpressure micro-cracking;
  ABAQUS FE stress intensification vs bubble density/shape/size (direct sibling /
  validation target).
  https://doi.org/10.1016/j.jnucmat.2025.156116

Mechanism of over-pressurisation (why bubbles sit above `2γ/r`):

- Cooper et al. (2024), *J. Nucl. Mater.* — irradiation-produced U interstitials
  over-pressurise bubbles; high-pressure-at-low-T → low-pressure-at-high-T
  transition; over-pressure builds during steady state.
  https://doi.org/10.1016/j.jnucmat.2024.155452
- Gruber (1982), *J. Nucl. Mater.* — cellular diffusional growth of over-pressured
  intergranular bubbles; swelling threshold on rapid heating, quench on cooling.
  https://doi.org/10.1016/0022-3115(82)90150-7
- Aagesen et al. (2021), *J. Nucl. Mater.* — phase-field HBS bubbles stay above
  equilibrium while growing; bubble-pressure response to a LOCA vs size and
  external restraint.
  https://doi.org/10.1016/j.jnucmat.2021.153267

Coverage / percolation and the restraint coupling:

- Aagesen et al. (2019), *Comput. Mater. Sci.* — GB **fractional coverage** and
  triple-junction saturation vs percolation; high semi-dihedral angle promotes it
  (relevant to the `Fc` sweep and the 50° dihedral).
  https://doi.org/10.1016/j.commatsci.2019.01.019
- Transient-tested BWR fuel SEM observations — under high compressive restraint
  bubbles ripen with little interlinkage; intergranular cracks open at *end of
  transient* when restraint is removed (cracking depends on restraint, not just
  internal pressure).
