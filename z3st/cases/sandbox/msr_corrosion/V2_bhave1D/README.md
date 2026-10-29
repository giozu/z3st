# V2: Bhave et al. 2023 1D benchmark

C. V. Bhave, G. Zheng, K. Sridharan, D. Schwen, M. R. Tonks, J. Nucl. Mater. 574 (2023) 154147. Experiments: G. Zheng, PhD thesis (2015), via Bhave Table 3. Full parameter sheet: `../bhave2023_reference.md`.

## Set-up

- Ni-5Cr (`ni5cr/`) and Ni-20Cr (`ni20cr/`) in FLiBe, 973 K, 150 um alloy, 1D, salt-facing surface at `xmax`.
- Cr diffusivity: Hart effective value from Bhave Eqs. 44-46, D0 converted from um2/s to m2/s. D_eff = 4.407e-17 (Ni-5Cr) and 1.535e-16 m2/s (Ni-20Cr).
- Initial Cr: 0.0561 and 0.2201 atomic fraction. Density 8700 and 8570 kg/m3.
- Surface: `fixed`, x_surf = 2.0315e-3, the atomic equivalent of Bhave's c_surf = 1.8e-3 mass fraction (Sec. 4).
- 202 backward-Euler steps, log-spaced from 1 s to 1000 h, including 50, 500 and 1000 h exactly. Mesh: 600 elements graded to 0.05 um at the surface.

Run with `./Allrun` in the `z3st` env; `./Allclean` removes the outputs.

## Results (mass loss, mg/cm2)

| Alloy | t (h) | z3st | Exact erfc integral | Bhave Eq. 50 | Experiment | z3st / exp. |
| --- | --- | --- | --- | --- | --- | --- |
| Ni-5Cr | 50 | 0.133 | 0.134 | 0.133 | 0.09 | 1.48 |
| Ni-5Cr | 500 | 0.421 | 0.423 | 0.421 | 0.27 | 1.56 |
| Ni-5Cr | 1000 | 0.595 | 0.598 | 0.596 | 0.42 | 1.42 |
| Ni-20Cr | 50 | 1.014 | 1.018 | 1.007 | 0.67 | 1.51 |
| Ni-20Cr | 500 | 3.205 | 3.220 | 3.184 | 3.06 | 1.05 |
| Ni-20Cr | 1000 | 4.533 | 4.553 | 4.503 | 7.12 | 0.64 |

- Against the exact solution: max error 0.46 % (backward-Euler lag). Against Eq. 50: max 0.66 %. Both pass (1 % and 2 %).
- Against the experiments: Ni-5Cr overpredicted; Ni-20Cr over at 50 h, close at 500 h, under at 1000 h. This is the pattern Bhave report for their own 1D Hart runs (Sec. 5.3: the Hart approximation suits the Ni-20Cr microstructure better than the Ni-5Cr one).
- Ni-5Cr depletion depth at 500 h (threshold 0.843 x0): 17.4 um; Bhave Fig. 2: 18 um, with D = 4.53e-5 um2/s (their Table 2) against 4.41e-5 here.

## Surface value

The `redox` condition computed from Bhave's Table 1 and Appendix A data gives x_surf = 1.858e-4, about 11x below the atomic equivalent of their quoted 1.8e-3. The case uses the quoted value; the mismatch changes the 1000 h mass loss by about 3 % and is listed as open in `../decision-log.md`.

`mass_loss.png`: z3st curves, Eq. 50 (dashed) and experiments (points).
