# z3st corrosion module, phase 1: reduced Cr-depletion model (design)

2026-10-29 · branch `feature/msr-corrosion` · status: draft for approval

## Intent

First building block of the z3st corrosion module (Paper 1). Selective Cr loss from an alloy into a molten salt is modelled as Cr diffusion in the alloy with a surface condition set by the salt. Bhave et al. 2023 show that their full electrochemical phase-field model reduces to this in 1D (within 5 %, their Sec. 4), so the reduced model reproduces their 1D benchmark at a fraction of the cost. It is also the 1D check for the full phase-field model that comes later.

Why z3st: the Cr field lives on the same mesh as the thermal and mechanical fields, so later phases can use a real temperature field and feed depletion into properties and damage. Phase 1 only has to make that possible, not do it.

Reference data: `bhave2023_reference.md` (same folder).

## Goals

1. A `corrosion` model in z3st that solves transient Cr diffusion in the alloy, in 1D, 2D and 3D, with the regimes z3st already supports.
2. Diffusivity from the material card: Arrhenius bulk and grain-boundary terms, combined with the Hart approximation (Bhave Eq. 46), evaluated at the local temperature.
3. A salt/alloy surface condition in two forms:
   - `fixed`: Cr surface fraction given directly;
   - `redox`: Cr surface fraction computed from equilibrium with the salt (Bhave Eqs. 28-31, B.1-B.4, B.8-B.9) from T, the fluoride potential E_F2/F-, and the Cr2+ and Ni2+ fractions in the salt. This is the entry point for the impurity work in phase 2.
4. Outputs per step: Cr field; Cr mass loss per unit surface area (mg/cm2); depletion depth, defined as the depth where Cr falls below a given fraction of its initial value (input, default 0.843 as in Bhave Fig. 2); Cr flux through the surface.
5. Two verification cases with non-regression scripts.

## Non-goals (this phase)

- Electric potential, order parameters, moving metal/salt interface, voids, salt infiltration.
- Explicit impurity species and the impurity -> redox-potential mapping (phase 2).
- Oxide layer (phase 3).
- Coupling Cr depletion to thermal or mechanical properties, or to damage (later phase). The design must not block it.
- Grain-resolved diffusion (comes with the full phase-field model).
- Fe and multi-component depletion (needed for SS316H later; the model holds one depleting species now).

## Anti-goals

- A benchmark that matches only because c_surf or D was tuned per case. Inputs come from the paper.
- Any change in the results of existing cases. The model is off unless `models: corrosion: true`.

## Physics

Unknown: Cr atomic fraction x in the alloy domain only. The salt is not meshed: Bhave fix the salt chemical potentials at the far boundary ("infinite salt") and the salt diffusivities do not affect the result (their Fig. 3).

    dx/dt = div( D(T) grad x )                       in the alloy
    x = x_surf                                        on the salt-facing boundary
    zero flux                                         elsewhere
    x(t=0) = x0

    D(T) = D_b(T) (1 - d_gb/G) + D_gb(T) d_gb/G       (Hart, Eq. 46)
    D_k(T) = D0_k exp(-Q_k / (R T)),  k = b, gb

`redox` surface value (equilibrium of Cr and Ni between salt and metal):

    G0_Cr2+ = G_bcc_Cr(T) + 2 (E0_Cr/Cr2+ + E_F2/F- - E0_F2/F-)          [eV]
    G0_Ni2+ = G_fcc_Ni(T) + 2 (E0_Ni/Ni2+ + E_F2/F- - E0_F2/F-)
    mu_i    = G0_i2+ + kB T + kB T ln(x_i,salt)                          (B.8, B.9)
    A_i     = exp( (mu_i - (G_i/Ni - G0_v)) / (kB T) )                   (B.3, B.4)
    x_surf  = A_Cr / (1 + A_Ni + A_Cr)                                    (B.2)

The redox potential shift is applied to both Cr and Ni, as Eq. 30 states it; Eq. 31 as printed omits E0_F2/F-, and using E_F2/F- alone gives x_surf = 1, so the shifted form is taken.

Mass loss per unit area: M = rho_alloy * (w0 - mean w over the depth) * depth, i.e. the integral of the drop in Cr mass fraction over the alloy volume divided by the salt-facing area.

## Known discrepancy to resolve (not blocking)

With the Table 1 and Appendix A values, the `redox` formula above gives x_surf = 1.86e-4 (mass fraction 1.65e-4). Bhave quote c_surf = 1.8e-3 (Sec. 4). The difference changes the 1000 h mass loss by about 3 % because x0 >> x_surf. The benchmark case uses `fixed` with 1.8e-3 (paper value) and reports the `redox` value alongside; the cause is chased later, starting from the H_xs/S_xs units and the kB T term in B.8/B.9.

## How it fits z3st (decided)

- New file `z3st/models/corrosion_model.py`, a `CorrosionModel` mixin on `Spine`, initialised when `models: corrosion: true`, following `PorosityMigrationModel` and `ClusterDynamicsModel`.
- Own scalar function space for x; CG1 Lagrange.
- Staggered with the thermal model within each time step, like porosity. With thermal off, T is a constant from the `corrosion` block.
- Inputs:
  - alloy diffusion data (D0, Q for bulk and GB, d_gb, grain size, density, molar masses, initial Cr fraction) in the alloy material YAML;
  - salt and surface data (type `fixed` | `redox`, x_surf or E_F2/F-, E0 values, salt Cr2+/Ni2+ fractions, depletion threshold) in a `corrosion:` block of `input.yaml`;
  - the salt-facing boundary named in `boundary_conditions.yaml` like other BCs.
- Units: SI inside z3st (m, s, K); paper values in um2/s are converted in the case files, not in the model.
- Outputs through the existing writer; scalar histories (mass loss, depth, flux) to the case output like other cases.

## Verification cases (acceptance)

Both under `z3st/cases/sandbox/msr_corrosion/` while in development; moved to `cases/verification/` (V1) and `cases/benchmarks/` (V2) once they pass.

**V1, analytic.** 1D, constant D, `fixed` x_surf, domain long enough to be semi-infinite. Compare with the erfc profile and with M(t) = 2 (w0 - w_s) rho sqrt(D t / pi). Pass: profile and mass-loss error below 1 %, and first-order (backward Euler) or second-order convergence shown in dt and h.

**V2, Bhave 1D.** 973 K, 150 um alloy, Ni-5Cr and Ni-20Cr, Hart D_eff from Eqs. 44-46, x0 = 0.0561 / 0.2201, `fixed` x_surf = 1.8e-3. Compare mass loss at 50, 500, 1000 h with:
- the analytic values (Ni-5Cr 0.133 / 0.421 / 0.596; Ni-20Cr 1.007 / 3.184 / 4.503 mg/cm2): pass below 2 %;
- the experiments in Bhave Table 3, reported (not a pass criterion), expecting the paper's pattern: Ni-5Cr overpredicted, Ni-20Cr over at 50/500 h and under at 1000 h;
- the Ni-5Cr depletion depth at 500 h (18 um in Bhave Fig. 2, with their D): reported.
Also report x_surf from `redox` for both alloys.

## Left to the builder

- Time integrator (backward Euler acceptable; BDF2 if cheap), time-step schedule (geometric growth suits the sqrt(t) start), mesh grading toward the surface.
- Exact names of input keys, as long as they follow z3st conventions.
- Whether the Hart combination is computed in the model or in a material function.

## Constraints

- No git commits on this branch (user rule).
- Code and comments follow the z3st style of the surrounding models.
- Existing non-regression suite unaffected.
