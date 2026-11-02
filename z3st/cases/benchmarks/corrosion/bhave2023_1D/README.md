# Bhave et al. 2023 1D benchmark (corrosion)

C. V. Bhave, G. Zheng, K. Sridharan, D. Schwen, M. R. Tonks, J. Nucl. Mater. 574 (2023) 154147. Experiments: G. Zheng, PhD thesis (2015), via Bhave Table 3. Full parameter sheet: `cases/sandbox/msr_corrosion/bhave2023_reference.md`.

## Set-up

- Ni-5Cr (`ni5cr/`) and Ni-20Cr (`ni20cr/`) in FLiBe, 973 K, 150 um alloy, 1D, salt-facing surface at `xmax`.
- Cr diffusivity: Hart effective value (Bhave Eqs. 44-46) with the constants of the authors' public MOOSE inputs (ref. [64], `simulations/1d_validation/Ni*Cr_hart.i`): Q/R = 34311.274 K (bulk) and 23576.457 K (GB), grain size 16.1122 um (Ni-5Cr) and 4.64 um (Ni-20Cr). D_eff = 4.531e-17 (Ni-5Cr) and 1.568e-16 m2/s (Ni-20Cr), the values of Bhave Table 2. Eq. 45 as printed rounds Q_gb to 196.2 kJ/mol, which gives a D about 2 % lower.
- Initial Cr: 0.0561 and 0.2201 atomic fraction. Density 8700 and 8570 kg/m3.
- Surface: `fixed`, x_surf = 1.7726e-3 atomic. This is the metal-surface fraction in equilibrium with the salt chemical potentials of the authors' inputs (w_Ni = -0.45065 eV, w_Cr = -0.9699 eV, i.e. x_Cr2+ = 2.5e-5 and x_Ni2+ = 3.7e-11 in the salt), quoted in the paper as c_surf = 1.8e-3.
- 202 backward-Euler steps, log-spaced from 1 s to 1000 h, including 50, 500 and 1000 h exactly. Mesh: 600 elements graded to 0.05 um at the surface.

Run with `./Allrun` in the `z3st` env; `./Allclean` removes the outputs and keeps `output/non-regression_gold.json`. The case is in the non-regression suite.

## Results (mass loss, mg/cm2)

| Alloy | t (h) | z3st | Exact erfc integral | Bhave Eq. 50 | Experiment | z3st / exp. |
| --- | --- | --- | --- | --- | --- | --- |
| Ni-5Cr | 50 | 0.136 | 0.136 | 0.135 | 0.09 | 1.51 |
| Ni-5Cr | 500 | 0.429 | 0.431 | 0.428 | 0.27 | 1.59 |
| Ni-5Cr | 1000 | 0.606 | 0.609 | 0.605 | 0.42 | 1.44 |
| Ni-20Cr | 50 | 1.026 | 1.030 | 1.018 | 0.67 | 1.53 |
| Ni-20Cr | 500 | 3.243 | 3.258 | 3.220 | 3.06 | 1.06 |
| Ni-20Cr | 1000 | 4.587 | 4.607 | 4.554 | 7.12 | 0.64 |

- Against the exact solution: max error 0.46 % (backward-Euler lag), the pass criterion (1 %). Against Eq. 50 as published (mass fractions, c_surf = 1.8e-3), computed with the same D: within 0.7 %, printed only.
- Against the experiments: Ni-5Cr overpredicted; Ni-20Cr over at 50 h, close at 500 h, under at 1000 h. This is the pattern Bhave report for their own 1D Hart runs (Sec. 5.3: the Hart approximation suits the Ni-20Cr microstructure better than the Ni-5Cr one).
- Ni-5Cr depletion depth at 500 h (threshold 0.843 x0): 17.7 um; Bhave Fig. 2: 18 um. The authors' scripts do not state a depth definition, so this comparison is indicative.

## Surface value

The `redox` condition with the same salt state (E_F = 3.3607 V, x_Cr2+ = 2.5e-5, x_Ni2+ = 3.696e-11) gives x_surf = 1.7711e-3, within 0.1 % of the fixed value. Both use the SGTE Gibbs energy of BCC Cr from the authors' inputs.

`output/mass_loss.png`: z3st curves, Eq. 50 (dashed) and experiments (points).
