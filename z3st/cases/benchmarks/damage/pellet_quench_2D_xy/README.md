# UO2 pellet thermal-shock fracture, 2D transverse cross-section
## z3st test case, configured after McClenny et al. (2022)

### Role

This case is a test case of the phase-field module on a thermal
transient. Its configuration follows the UO2 pellet thermal-shock study of
McClenny et al. (JNM 565, 2022): a transverse cross-section under plane
strain, a circular pellet with a 60-deg cold-contact arc on its perimeter,
producing discrete radial cracks. It is neither a verification nor a
validation, and it does not reproduce their experiment (see Section 3.2 of
the SoftwareX paper): formulation, discretisation, regularisation length and
initial temperature differ from theirs.

It is in the local suite (carries a gold).

An axisymmetric (r-z) model can only produce an annular damage band. This
(x, y) cross-section can break azimuthal symmetry and produce the localised
tensile hoop stress that drives radial cracks.

### Plane-strain idealisation

In the McClenny experiment:
- The cold bath was in contact with 1/6 of the lateral surface
  of the pellet (Fig. 4), a 60-deg azimuthal contact wedge.
- An alumina spacer at the capsule bottom removed axial thermal
  contact (Section 3), so heat transfer was radial.

With no axial gradient, every transverse cross-section sees the same loading,
which is the plane-strain assumption. McClenny use the same 2D representation
for their parametric study (Fig. 8 top row, Fig. A.13, Fig. A.14).

### Geometry: half-disc with mirror symmetry

The contact wedge is centred on the +x axis and spans -30 deg to +30 deg
(60 deg in total). The y = 0 plane is a mirror symmetry plane of loading and
geometry. The mesh covers the upper half of the disc (y >= 0).

`Clamp_y` on the diameter is the symmetry condition (it removes rigid
y-translation and rotation around z). A 50-um pin segment near (-R, 0)
carries a `Clamp_x` to remove rigid x-translation.

### Geometry, material, BCs

| Parameter | Value |
|---|---|
| Pellet radius `R` | 10 mm |
| Domain | upper half-disc (y >= 0); contact arc 0 to +30 deg in the upper half (= 60 deg full) |
| Material | UO2 from `z3st/materials/uo2.yaml`: E = 205 GPa, nu = 0.32, alpha = 1e-5 /K, sigma_c = 1 GPa |
| `T_initial` | 1023.15 K (= 750 deg C), uniform |
| `T_quench` | 263.15 K (= -10 deg C), Dirichlet on the 60-deg contact arc |
| Rest of perimeter | natural Neumann (zero heat flux) |
| Symmetry plane y = 0 | `Clamp_y` (mirror BC) |
| Pin segment at (-R, 0) | `Clamp_x` (50 um, removes rigid x-translation) |
| Damage | AT1, `split: star_convex` with `gamma_star: 0.0` (the Amor split), `lc = 50 um`, `hybrid_constraint: true` |
| Time window | 0.0001 to 0.1 s, n_steps = 100. The 99 intervals are uniform at dt = 1.01 ms, the first step is 0.1 ms |

`sigma_c = 1 GPa` (from `materials/uo2.yaml`) lies above the
plane-strain bulk threshold (about 3.0 MJ/m^3 from the blocked
z thermal expansion, which `damage_model.py` suppresses) and below
the peak surface tensile strain energy at the rim, so damage initiates
only along the cold contact arc. `spine.py` derives `Gc` from it:
Gc = (8/3) lc sigma_c^2 / E = 650 J/m^2 at E = 205 GPa, as printed
in the run log.

<!-- [TBC] "about 3.0 MJ/m^3" and "below the peak surface tensile strain energy" are not checked by any script in this case -->

### Phase-field formulation

Ambati hybrid formulation (Comput. Mech. 55 (2015) 383-405, Eq. 27):
- linear isotropic stress degradation `sigma = (1-D)^2 * dPsi0/de`;
- AT1 surface energy;
- star-convex split with `gamma_star = 0`, identical to the Amor
  (volumetric/deviatoric) split;
- hybrid constraint, `H = 0` in compression-dominated cells.

The damage driving force is evaluated on the elastic strain
`eps_el = eps(u) - alpha (T - T_ref) I`. In 2D plane strain the z-component
of the eigenstrain is suppressed (see the `damage_model.py::_thermal_eigenstrain`
docstring).

McClenny use the Miehe anisotropic formulation with viscous Allen-Cahn
evolution. Agreement with their crack pattern is therefore not a
code-to-code verification.

### Expected results

From McClenny:
- Discrete radial cracks from the cold contact arc within the
  -30 deg to +30 deg wedge.
- "two major (longer) radial cracks" plus a fan of shorter surface
  cracks (p. 7).
- Cracks "immediately appear on the pellet outer surface at about
  10^-2 s right after the instantaneous drop in temperature" (p. 8).

Stored run (`energies.txt`, 100 steps):
- `E_el` is largest at step 0 (364.48 J) and decreases at every step, to
  355.39 J at step 99 (t = 0.1 s).
- `E_frac` is already 0.31 J at step 0 (t = 1e-4 s), 1.28 J at step 1,
  2.52 J at step 10, 4.22 J at step 50 and 5.22 J at step 99.

The regularisation length is 50 um against 1 um in McClenny, so crack bands
are about 50 times wider.

### Running

```bash
cd z3st/cases/benchmarks/damage/pellet_quench_2D_xy
./Allrun
```

`python3 mpi_scaling.py` reruns the case at 1, 2 and 4 MPI ranks and writes
`output/mpi/non-regression_np<N>.json` for each.

### Outputs (in `output/`, written by `non-regression.py`)

- `damage_field.png`: 2D colour map of `D` on the half-disc,
  with the cold contact arc highlighted (compare qualitatively with
  McClenny Fig. 8, top right).
- `temperature_field.png`: `T` at the final time.
- `stress_vm_field.png`: von Mises stress at the final time
  (clipped to the 99th percentile).
- `stress_hoop_field.png`: hoop stress sigma_theta_theta, symmetric
  red/blue colour map (red = tensile).
- `damage_angular.png`: D_max(theta) along the outer ring. Each peak is one radial crack.
- `energy_balance.png`: E_el(t), E_frac(t).
- `non-regression.json`: tracked values `T_final_mean_in_range`, `D_max_final`,
  `crack_count_above_0p5`.

### References

- McClenny et al., JNM 565 (2022) 153719.
- Ambati, Gerasimov, De Lorenzis, Comput. Mech. 55 (2015) 383-405.
