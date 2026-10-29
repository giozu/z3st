# Bhave et al. 2023: reference sheet for the z3st P1 benchmark

C. V. Bhave, G. Zheng, K. Sridharan, D. Schwen, M. R. Tonks, "An electrochemical mesoscale tool for modeling the corrosion of structural alloys by molten salt", J. Nucl. Mater. 574 (2023) 154147. doi:10.1016/j.jnucmat.2022.154147. PDF: `0.biblio/1-s2.0-S0022311522006262-main.pdf`. Input files: GitHub repository, ref. [64] of the paper.

## System

- Ni-5Cr and Ni-20Cr in FLiBe, **T = 973 K (700 °C)** for all simulations and experiments (Table 1).
- Experiments: G. Zheng PhD thesis (2015), ref. [34]. Mass loss at 50, 500, 1000 h; EDS Cr line scans at 1000 h.
- Driving impurity: H+ from residual HF. The salt redox state enters only through the fluoride potential E_F2/F- (3.36 V vs 2.871 V reference). No impurity is tracked explicitly.

## Model (grand-potential electrochemical phase field, MOOSE)

Unknowns: chemical potentials mu_Ni, mu_Cr; electric potential phi; order parameters eta_s (salt) and eta_mk (grains).

| Eq. | Content |
| --- | --- |
| 7 | d mu_i/dt = (1/chi_i) [div(M_i grad mu_i) + div(z_i F M_i grad phi) - sum_j (d rho_i/d eta_j)(d eta_j/dt)] |
| 8, 10 | Allen-Cahn for eta: d eta/dt = -L d Omega/d eta |
| 12, 13 | h_m = 0.5 - 0.5 tanh(2 pi (eta_s - 0.5)/delta_int), h_s = 1 - h_m |
| 14-18 | Interface: m_S = 6 sigma_S/l, kappa_S = 3 sigma_S l/4, kappa_gb = 3 sigma_gb l/4, gamma = 1.5 |
| 26 | Conserved current condition for phi |
| 28 | Metal: ideal solution in Ni, Cr, vacancy; c_v = 1 - c_Ni - c_Cr |
| 29 | Salt: dilute Ni2+, Cr2+ |
| 30, 31 | G0_Cr2+ = G_bcc_Cr + nF (E0_Cr/Cr2+ + (E_F2/F- - E0_F2/F-)): **this is where the salt redox state, and so the impurities, enter** |
| 41 | L = i0 delta_int / (pi n F kB T l Va^2) |
| 46 | Hart effective diffusivity: D_eff = D_b (1 - d_gb/G) + D_gb d_gb/G (1D only) |
| 47, 48 | 2D/3D: GB switching function h_gb (b = 0.8, delta_gb = 0.01); D = h_m (D_b + h_gb (D_gb - D_b) d_gb/l) + h_s D_s |
| 50 | Analytic check: M_Cr = 2 (c0 - c_surf) rho sqrt(D_Cr t/pi), c_surf = 1.8e-3 (mass fraction) |
| B.1-B.10 | Closed-form rho(mu) and omega for metal and salt phases |

## Parameters

| Quantity | Value | Source in paper |
| --- | --- | --- |
| T | 973 K | Table 1 |
| d_gb | 5e-4 um | Table 1 |
| Va (Ni-5Cr / Ni-20Cr) | 1.113e-11 / 1.109e-11 um^3/atom | Table 1 |
| G0_v | 1.56 - 3.3 kB T eV/atom | Table 1 |
| E0 Ni/Ni2+, Cr/Cr2+ | 0.473 V, -0.39 V | Table 1 |
| E0_F2/F-, E_F2/F- | 2.871 V, 3.36 V | Table 1 |
| i0 | 5.7 A/m^2 | Table 1 |
| sigma_S, sigma_gb | 1.43e7, 6.8e6 eV/um^2 | Table 1 |
| S_e(973 K) | 0.524 S/um | Table 1 |
| D_b Ni | 1.9e8 exp(-2.80e5/RT) um^2/s | Eq. 42 |
| D_gb Ni | 7.0e6 exp(-1.15e5/RT) um^2/s | Eq. 43 |
| D_b Cr | 4.28e8 exp(-2.852e5/RT) um^2/s | Eq. 44 |
| D_gb Cr | 4.85e10 exp(-1.962e5/RT) um^2/s | Eq. 45 |
| D Ni2+, Cr2+ in FLiBe | 500 um^2/s | Sec. 2.3.3 |
| D H+ | 1000 um^2/s (3000 in literature) | Sec. 2.3.3 |
| G_bcc_Cr | Eq. 32 (eV/atom) | Eq. 32 |
| G_fcc_Ni, G_fcc_Cr, L_NiCr | Eqs. A.2, A.3, A.5 (J/mol) | App. A |
| H_xs_Cr, S_xs_Cr | -1.5644e4 J/mol, -1.5601 J/mol/K [TBC unit] | App. A |
| Initial c_Cr (Ni-5Cr / Ni-20Cr) | 0.0561 / 0.2201 at. fr. | Sec. 2.5 |
| Salt c_Cr2+, c_Ni2+ | 2.5e-6, 3.69e-11 | Sec. 2.5 |
| 1D numerics | l = 2 um, 640 linear elements, BDF2, SuperLU, dt0 = 1e-3 s, dt_max = 1e4 s | Sec. 2.4 |
| 1D domain | 150 um metal + 10 um salt | Sec. 2.5 |

## Boundary conditions (1D)

- phi = 0 on the metal end. mu_Ni, mu_Cr fixed at their initial salt values on the salt end (infinite salt). Zero flux elsewhere.

## Validation data (Table 3)

| Alloy | Density (g/cm^3) | Grain size (um) | Mass loss 50 h | 500 h | 1000 h (mg/cm^2) |
| --- | --- | --- | --- | --- | --- |
| Ni-5Cr | 8.70 | 16.22 | 0.09 | 0.27 | 0.42 |
| Ni-20Cr | 8.57 | 4.64 | 0.67 | 3.06 | 7.12 |

Figures with digitisable curves: Fig. 2 (Ni-5Cr 1D profile, 1 and 500 h; depth 18 um at 500 h), Fig. 4 (1D vs Eq. 50), Fig. 5e and 6h (EDS line scans vs 1D/2D/3D).

## Findings that shape the z3st design

1. **The 1D problem is diffusion-limited.** Sensitivity (Fig. 3): D_Cr in the alloy dominates by more than 20x. i0, salt diffusivities, sigma and the electric field have negligible effect. The electromigration term is at least 4 orders of magnitude below the chemical term after a few seconds (Sec. 3).
2. **1D phase field = Fick with a surface Dirichlet value**, to within 5 % (Sec. 4, Fig. 4). The redox state of the salt enters only through c_surf.
3. Our check of Eq. 50 with Hart D_eff at 973 K (D_b = 2.09e-7, D_gb = 1.42 um^2/s):

| Alloy | D_eff (um^2/s) | Eq. 50, 50 / 500 / 1000 h | Exp. | 2 sqrt(D t) at 1000 h |
| --- | --- | --- | --- | --- |
| Ni-5Cr | 4.40e-5 | 0.13 / 0.42 / 0.60 | 0.09 / 0.27 / 0.42 | 25 um |
| Ni-20Cr | 1.53e-4 | 1.01 / 3.18 / 4.50 | 0.67 / 3.06 / 7.12 | 47 um |

   Same pattern as the paper's 1D Hart results: Ni-5Cr overpredicted; Ni-20Cr over at 50 and 500 h, under at 1000 h.
4. The 2D runs underpredict; only 3D matches Ni-5Cr. Mass loss scales with the inverse projected grain size at the surface (Fig. 7).
5. Out of scope in the paper: explicit impurities, Cr3+, oxides, voids, salt infiltration, mechanics.
