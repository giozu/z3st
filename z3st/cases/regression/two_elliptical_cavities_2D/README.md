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

- `parametric_study.py` — **the main deliverable.** Analytical `p_crit(Fc, Gc)`
  reference (ligament + elliptical-tip concentration) and an FEM sweep that
  regenerates the mesh at each `Fc`, ramps the bubble pressure, detects GB
  percolation, and saves a per-`Fc` crack figure.
    - `python3 parametric_study.py`            → analytical reference map
    - `python3 parametric_study.py --fem 0.2 0.4 0.6` → sweep + per-case figures
- `non-regression.py` — per-run diagnostic: extracts the fields and writes
  `fields_overview.png`, `stress_profile_tip.png`, `gc_profile_check.png`.
- `mesh.geo` — parametrised by `Fc_target` (`gmsh -setnumber Fc_target 0.3 …`),
  refined to `h_cavity = lc/2` so the AT2 damage band can localise.
- `input.yaml` / `boundary_conditions.yaml` — the **committed config**: the
  equilibrium-pressure baseline (bubble pressure ramped to 15 MPa).

## Output state (IMPORTANT — read before trusting the plots)

The figures in `output/` are an **illustrative FRACTURED example at `Fc = 0.4`**
produced by a *high-pressure* `parametric_study` run (ramp ~3 GPa), **not** the
committed `input.yaml` (15 MPa). They are kept on purpose to show the crack; a
plain `./Allrun` at the committed 15 MPa load will overwrite them with the
no-fracture equilibrium state. The persistent sweep deliverables are:

- `pcrit_vs_Fc_sweep.png` — `p_crit` vs `Fc` (FEM points on the analytical map)
- `sweep_Fc0.20.png` / `sweep_Fc0.40.png` / `sweep_Fc0.60.png` — per-coverage
  Damage + `σ_yy` (crack morphology: isolated tips → bridging → fully linked)

## OPEN DECISION (deferred) — the (Gc, lc) regime

With the current `Gc_gb = 0.1 pJ/µm²` and `lc = 4 nm`, the AT2 peak stress is
`σ_c ≈ 670 MPa`, and the *sharp bubble-tip* nucleation pushes the effective GB
resistance higher still — so the GB only fractures at **~1–2.4 GPa** bubble
pressure (validated: `Fc 0.2/0.4/0.6 → p_crit 2368/2053/1105 MPa`).

Equilibrium pressures for these ~0.1 µm bubbles are only **tens of MPa**
(`p_eq ≈ 2γ/r`), so the model fractures only under strong **transient
over-pressurisation**. To model GB cracking at realistic *equilibrium*
pressures, move `(Gc, lc)` to a softer regime so `σ_c` is reachable there.

Consequently the `non-regression.py` pass/fail metrics are **not yet
canonical**: `max_damage → 1` is unreachable at 15 MPa, and `max_stress_yy`
compares against a hard-coded `p_applied = 1.0 MPa` that does not match the
applied ramp. Both should be reconciled once the `(Gc, lc)` regime is chosen.
This case is therefore **not in the regression suite** (no gold).

## Already fixed this round

- `Gc` reconciled: `materials/oxide.py` is the single source (`Gc(mesh)` UFL +
  `Gc_numpy(y)`); the script imports it (a stray `*1e-6` in a local copy is gone).
- Output field name corrected (`Stress_solid → Stress`).
- Mesh refined (`h_cavity = lc/2`) and parametrised by `Fc`.
- Cleaner diagnostic plots (field overview, GB stress profile, Gc, fracture map).
