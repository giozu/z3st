# creep_shrink_fit_2D: contact-pressure relaxation by Norton creep

2D axisymmetric (r,z) UO₂ pellet inside a Zircaloy-4 cladding. The pellet
expands thermally, closes a 10 μm cold radial gap, and loads the cladding
through the penalty contact model. The cladding then creeps under sustained
load and the contact pressure relaxes in time.

Reference for the closed form: Esposito, Bruno, Bertocco, *Int. J. Pressure
Vessels and Piping* 185 (2020) 104126, eq. (21), shrink-fit assembly
relaxing by Norton creep:

<!-- [TBC] plots.py evaluates phi with (k2 - k1); the formula below has (k1 + k2). Equal here because k1 = 0 (solid shaft); check eq. (22) of the paper -->
    Pk(t) = Δu_el(t) · f,   Δu_el(t) = Δ · (1 + φ(1-n)·t / Δ^(1-n))^(1/(1-n))
    φ = -(A/b) · f^n · (k1 + k2)

`plots.py::plot_contact_pressure_evolution` evaluates this law and overlays it
on the simulated `output/history.csv`. Δ is the interference of the assembled
joint at temperature, held fixed across the transient (the law relaxes it
itself), and `f` is corrected for the penalty spring in series with the joint,
`1/(1/f + 1/K_PEN)`. The relaxation itself is verified against
`reference_1d.py`, an independent radial J2 solution: see "Result" below.

Author of the case and of the analysis scripts: Romain Turgis (ENSTA Paris),
internship 2026-05-25 → 2026-07-31, supervised by G. Zullo.

## Parameter handling

`case_params.py` is the single source of truth for the post-processing: it
reads `geometry.yaml`, `input.yaml`, `clad.yaml` and `fuel.yaml` and exposes
the resolved constants to `plots.py` and `non-regression.py`. No
geometric or material value is restated as a literal in a script.

`case_params.check_consistency()`, called from `diagnostics.py` before step 0,
aborts the run on: a cold gap from `geometry.yaml` whose magnitude
disagrees with `models.contact.initial_gap`, an `initial_gap` that is a
clearance rather than an interference, and irradiation creep left switched on.
`plots.py` calls it on every run and warns before plotting.

The gap check compares magnitudes, not signed values. A conforming mesh has to
keep the pellet and the cladding disjoint, so `geometry.yaml` can only carry a
positive clearance, while the interference is expressed by a negative
`initial_gap` that the contact model uses in place of the mesh-derived value.
The two therefore differ in sign by construction, and only their magnitudes are
required to agree.

## Result

Serial `./Allrun` on the committed configuration (non-fissile pellet, cracking
off, traction-free clad OD). `non-regression.py` asserts the four rows against
the `reference_1d` column (gold values).

| t [days] | Z3ST [MPa] | reference_1d [MPa] | deviation | eq. (21) [MPa] |
|---|---|---|---|---|
| 0 | 24.703 | 24.699 | +0.02 % | 24.699 |
| 600 | 7.256 | 7.120 | +1.90 % | 6.90 |
| 1240 | 5.140 | 5.063 | +1.51 % | 4.90 |
| 2500 | 3.648 | 3.605 | +1.19 % | 3.49 |

The elastic point is the Lame interference pressure of the joint with the
penalty spring in series. The relaxation is checked against
`reference_1d.py`, an independent radial solution of the same joint: elastic
solid shaft, J2 hub in generalised plane strain, penalty spring in series,
explicit creep update converged in time. It shares no code with Z3ST and no
approximation with eq. (21). Doubling `n_steps` takes the 2500-day value from
3.648 to 3.633 MPa (reference 3.605), so the deviation is time-step error.

Esposito eq. (21) lies about 3 % below reference_1d (4-5 % below Z3ST). Its
authors report this direction: their hub uses the Tresca criterion, which
overestimates the von Mises equivalent stress, so the closed form overestimates
the decrease of Pk. It is plotted by `plots.py` and not asserted.

### Notes on the closed form

- Eq. (4) of the paper as printed has the signs of nu swapped between shaft and
  hub (+nu1, -nu2; the Lame interference fit has -nu1, +nu2). With equal
  materials, as in the paper's own validation, the two terms cancel. With UO2
  and Zircaloy the printed signs give f 3 % higher (25.47 MPa elastic point
  instead of 24.70 MPa). `case_params.elastic_factor` uses the Lame signs.
- `k2` in `plots.py` is eq. (20) as published (plane-stress hub, Tresca). The
  von Mises plane-strain shaft form, eq. (16) with c in place of a, is 0.61 times
  eq. (20) for this clad.

### Convergence

Three refinements, each against the committed configuration. Mesh and element
order move the pressure by 0.2 %, the time step by 0.4 %.

| refinement | p(2500 d) [MPa] |
|---|---|
| committed (`n_r1/n_r2 = 25/7`, P1, `n_steps 56/80`) | 3.648 |
| mesh `n_r1/n_r2 = 49/21` | 3.731* |
| `mechanical.order: 2` (P2 displacement) | 3.730* |
| `n_steps 112/160` | 3.633 |

\* measured at `n_steps 14/20`, against a 3.738 MPa baseline at the same
stepping: mesh and element order move the pressure by 0.2 %.

`Mesh.ElementOrder = 2` fails in the output writer (`RuntimeError: Degree of output
Function must be same as mesh degree`, `utils/writer.py`). `mechanical.order`
(`core/finite_element_setup.py`) raises the displacement degree on the linear mesh.

### Configurations that do not reproduce eq. (21)

| Configuration | Behaviour | Cause |
|---|---|---|
| fissile pellet, clamped clad OD | P_c rises to 160 MPa | burnup swelling grows Δ; eq. (21) assumes Δ fixed |
| non-fissile, clamped clad OD | P_c freezes at 12.9 MPa | with u_r = 0 at r = c creep redistributes stress until the deviatoric part vanishes, then stops |
| non-fissile, free clad OD, `initial_gap: +10 µm` | no contact | without the clamp nothing closes the gap |
| non-fissile, free clad OD, `initial_gap: -10 µm` | relaxes | a negative gap is an interference, contact active from t = 0 |

The last row is the paper's shrink fit and is the committed configuration.

## Configuration constraints

The case carries `output/non-regression_gold.json` and is a member of the local
suite. Five settings are required for the comparison with eq. (21).
`case_params.check_consistency()`, called from `diagnostics.py` before step 0,
aborts the run on the first three.

1. Irradiation creep off. `creep_irr_B` and `fast_flux` are commented out in
   `clad.yaml`. Esposito eq. (21) is a pure Norton derivation with no in-pile
   term.

2. Every constant read from the YAML, through `case_params.py`.

3. `creep_Q: 0`: the case is isothermal. Esposito eq. (9) is
   `eps_eq_dot = A*sigma_eq^n` with `A` a constant `[h^-1]` (their
   nomenclature): there is no Arrhenius factor in the derivation, and
   temperature enters only through the interference via `alpha*dT`.

   `creep_Q = 0` makes `A(T) = A0` identically, so `A0` must have the
   Arrhenius factor folded in. The Zircaloy card value
   `2.82e-24 Pa^-n s^-1` is the prefactor that belongs with `creep_Q = 1.2e5`;
   used with `creep_Q = 0` it runs the cladding about 1e10 times too fast.
   `clad.yaml` carries the re-based
   `A0 = 2.82e-24 * exp(-1.2e5/(R*600 K)) = 1.0083e-34 Pa^-n s^-1`, the same
   `A(T)` that `verification/fuel/creep` reports at 600 K.

4. `fissile: false`, and cracking off. Eq. (21) is integrated under a
   fixed interference assembled once and left to relax. A fissile pellet grows
   the interference through the swelling eigenstrain and the contact pressure
   rises instead of relaxing (see the table above). Zero burnup switches off both
   the swelling and the densification terms, so the interference changes only by
   creep; the eigenstrain entry stays in the card, inert. Cracking would degrade
   the pellet to `E_iso/E = 0.112` while the analytical overlay uses the card's
   nominal `E = 200 GPa`, and it fires on the nominal LHR even with
   `fissile: false`, i.e. on power that is never deposited. Esposito's shaft is
   an elastic solid.

5. Clad outer surface traction-free. `outer_2` carries no restraint in
   `boundary_conditions.yaml`, which is Esposito's hub boundary condition
   `sigma_r(c) = 0`. Under `u_r = 0` at `r = c` creep can only redistribute
   stress until the deviatoric part vanishes and then stops, freezing the
   contact pressure at a plateau instead of decaying as a power law.

## Open point

`A0` is re-based at a round 600 K, while the case runs isothermal at 580 K.
Both references use the same `A0`. The value sets how fast this joint relaxes.

### Elastic point at t = 0

`contact_pressure_elastic_MPa` is asserted against the Lame interference
pressure with the joint and the penalty contact treated as springs in series:

    p = Delta / (1/f + 1/K_PEN)     24.6986 MPa predicted, 24.7028 MPa simulated, +0.02 %

Nothing on the right-hand side is read back from the solver: `f` comes from the
geometry and the moduli, `K_PEN` from `input.yaml`, `Delta` from the material
cards and the run temperature.

`Delta` is the interference of the assembled joint at temperature, not the
cold 10 um of `models.contact.initial_gap`. At 580 K the pellet outgrows the
cladding by 4.73 um (`r*alpha*(T - T_ref)` on each body, exact here because the
field is uniform), so the assembled interference is 14.73 um. Using the cold
value instead predicts 17.29 MPa (43 % below).

The case tolerance is 3e-2, set by the time-step error of the relaxation at
600 days (1.9 %); the burnup assertions compare against exactly zero.

## Run

    conda activate z3st
    ./Allrun      # gmsh → z3st → non-regression → plots
    ./Allclean
