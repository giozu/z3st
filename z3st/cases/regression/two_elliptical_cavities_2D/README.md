# two_elliptical_cavities_2D

Micromechanical model of intergranular fracture driven by fission-gas bubbles
on a grain boundary in oxide fuel. Two lenticular gas bubbles (elliptical
cavities) sit on a grain boundary (GB) at `y = 0`. The internal bubble pressure
concentrates stress at the bubble tips and loads the GB ligament.
An AT2 phase-field damage model, with a GB-weakened toughness `Gc(y)`, lets a
crack nucleate at the tips and link the bubbles along the GB.

The target is the GB failure pressure as a function of the bubble fractional
coverage `Fc` and the GB toughness `Gc` (fuel fragmentation, fission-gas release).

Units: micron / second / kg / micronewton / MPa / picojoule (a consistent set:
`µN = kg·µm/s²`, `MPa = µN/µm²`, `pJ = µN·µm`; and `pJ/µm² = J/m²` numerically).

## Files

- `parametric_study.py`: the sweep driver. Two analytical `p_crit(Fc, Gc)`
  references (a fracture-mechanics SIF estimate and a strength estimate) and an
  FEM sweep that regenerates the mesh at each `Fc`, ramps the bubble pressure,
  detects GB percolation, and saves a per-`Fc` crack figure. The CONFIGURATION
  block at the top (`CASES`, `GC_VALUES`, `CRACK_FACES`, `SIGMA_H`) sets the cases.
    - `python3 parametric_study.py`            → analytical reference map
    - `python3 parametric_study.py --fem`      → sweep the configured `CASES`
    - `python3 parametric_study.py --fem 0.3 0.5` → override the `Fc` list
- `non-regression.py`: per-run diagnostic. Writes `fields_overview.png`,
  `stress_profile_tip.png`, `gc_profile_check.png`.
- `mesh.geo`: parametrised by `Fc_target` (default 0.40, `gmsh -setnumber Fc_target 0.3 …`),
  refined to `h_cavity = lc/2 = 0.002 µm`.
- `input.yaml` / `boundary_conditions.yaml`: the equilibrium-pressure baseline
  (bubble pressure ramped to 15 MPa in 10 steps). AT2, all three solves `iterative_hypre`.

## Figures in `output/`

- Baseline diagnostics: `fields_overview.png`, `stress_profile_tip.png`,
  `gc_profile_check.png`. Written by `non-regression.py` for the 15 MPa
  no-fracture state. `d_max ≈ 0` here.
- Sweep figures (fractured, GPa-scale loads), written by `parametric_study.py`:
    - `pcrit_vs_Fc_sweep.png`: FEM `p_crit` points overlaid on the SIF map
    - `sweep_Fc0.20.png`, `sweep_Fc0.40.png`, `sweep_Fc0.60.png`: per-coverage damage and `σ_yy`
    - `pcrit_vs_Fc_reference.png`: analytical map only (SIF and strength)

`Allclean` removes every `output/*.png`, including the sweep figures, and keeps
`non-regression_gold.json`.

## Analytical references: strength and fracture mechanics (SIF)

Strength estimate: `p_crit = σ_c·(1−Fc)/K_t` (AT2 peak stress × tip
concentration / ligament). It combines a pointwise `K_t` with an `lc`-scale `σ_c`
and predicts tip nucleation, not GB percolation. It lies about 17× below the FEM.

Fracture-mechanics SIF estimate (Chakraborty, Tonks and Pastore 2014):

      p_crit = K_Ic / (F(Fc)·√(π·R)) + σ_h,   K_Ic = √(E·Gc/(1−ν²)),   R = ax

with `F(Fc)` their non-dimensional Mode-I SIF (`F_sif`, Eqs. 5/9) and `σ_h` a
compressive hydrostatic restraint (Eqs. 6/8). For this case
`K_Ic = 0.138 MPa·µm½`, giving `Fc 0.2/0.4/0.6 → 414/365/308 MPa`, about 4 to 6×
below the FEM.

Chakraborty assume a sharp pre-crack (LEFM). The phase-field bubble tip is a
smooth ellipse (`ρ = ay²/ax ≈ 22 nm`) with finite `lc = 4 nm`. The material
characteristic length `ℓ_ch = K_Ic²/σ_c² ≈ 42 nm ≈ ρ`, so the case lies in the
strength-toughness transition. The ordering of the three estimates is
`strength < SIF-LEFM < phase-field percolation`.

## The (Gc, lc) regime

`Gc` has a single source, `z3st/materials/oxide.py` (`Gc(mesh)` as UFL,
`Gc_numpy(y)` for post-processing). The case scripts import it.

With `Gc_gb = 0.1 pJ/µm²` and `lc = 4 nm`, `σ_c ≈ 670 MPa`. The FEM sweep gives GB
failure at `Fc 0.2/0.4/0.6 → p_crit 2368/2053/1105 MPa`.

Equilibrium pressures for these bubbles (about 0.1 µm) are tens of MPa
(`p_eq ≈ 2γ/r`), so the model fractures only under transient
over-pressurisation, the trigger reported in the literature below.
Two options are open:
- move `(Gc, lc)` to a softer regime, with `σ_c` reachable at equilibrium pressure;
- keep the regime and model the over-pressure and restraint (`σ_h`) explicitly.

## Non-regression metrics

`non-regression.py` checks the equilibrium baseline (15 MPa, no fracture):

- `max_damage`: reference 0 (percolation threshold 0.9). Gold `d_max ≈ 3.3e-5`.
- `max_stress_yy`: reference `K_t/(1−Fc)·p_applied`, with `p_applied` read from
  the ramp (15 MPa). The isolated `K_t` underpredicts by about 2× (the two bubbles
  share load through the ligament). The ligament-corrected estimate gives 82 MPa
  against 99 MPa from the FEM (about 20 %).
  Analytic tolerance 0.25. The gold regression check uses rtol 1e-3.

The case is in the local suite (carries a gold) and in `cases_ci.txt` (about 80 s).
To re-bless after a checked baseline run:
`cp output/non-regression.json output/non-regression_gold.json`.

## Background literature (GB-bubble overpressurisation → fragmentation / FGR)

Context for the `(Gc, lc)` regime: the literature reports GB cracking under
transient over-pressurisation, modulated by external mechanical restraint, and
not at equilibrium pressure.

Analytical and FEM precedents:

- Chakraborty et al. (2014), *J. Nucl. Mater.*: Mode-I non-dimensional SIF for
  lenticular GB bubbles under bubble pressure and hydrostatic stress (the
  geometry of this case; source of the SIF reference above).
  https://doi.org/10.1016/j.jnucmat.2014.04.023
- Cappellari et al. (2025), *J. Nucl. Mater.*: PoliMi/SCIANTIX GB FGR model
  applying fracture mechanics to bubble-overpressure micro-cracking, with ABAQUS FE
  stress intensification against bubble density, shape and size.
  https://doi.org/10.1016/j.jnucmat.2025.156116

Mechanism of over-pressurisation (bubble pressure above `2γ/r`):

- Cooper et al. (2024), *J. Nucl. Mater.*: irradiation-produced U interstitials
  over-pressurise bubbles; high-pressure-at-low-T → low-pressure-at-high-T
  transition; over-pressure builds during steady state.
  https://doi.org/10.1016/j.jnucmat.2024.155452
- Gruber (1982), *J. Nucl. Mater.*: cellular diffusional growth of over-pressured
  intergranular bubbles; swelling threshold on rapid heating, quench on cooling.
  https://doi.org/10.1016/0022-3115(82)90150-7
- Aagesen et al. (2021), *J. Nucl. Mater.*: phase-field HBS bubbles stay above
  equilibrium while growing; bubble-pressure response to a LOCA vs size and
  external restraint.
  https://doi.org/10.1016/j.jnucmat.2021.153267

Coverage / percolation and the restraint coupling:

- Aagesen et al. (2019), *Comput. Mater. Sci.*: GB fractional coverage and
  triple-junction saturation vs percolation; high semi-dihedral angle promotes it
  (relevant to the `Fc` sweep and the 50° dihedral).
  https://doi.org/10.1016/j.commatsci.2019.01.019
