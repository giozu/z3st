# FE101 phase 1: native FIMA → swelling verification

**Verification case only; thermal power is the Z3ST baseline, not OpenMC operational power.**

Coordinate transform: `r_FE_m=r_OpenMC_cm/100`, `z_FE_m=(z_OpenMC_cm−10.20)/100`. The new isolated mesh is 0.3556 m long; the existing 0.356 m regression is unchanged. No stretching, clamping or extrapolation. History: all 21 times through 3500 h, 2100 source records; noncoincident times use linear interpolation of cumulative fissions.

Imported native FIMA and internally accumulated baseline BU are independent. Explicit `fima_source: native_openmc` prevents a silent legacy fallback. The law remains `epsilon_sw=FIMA*I`, `DeltaV/V=3*FIMA`. No q''' import, creep or hydrogen redistribution.

A/B interface tests: PASS (seven tests). Volume relative error 2.442e-15; initial atom error 2.665e-15; final volume-weighted mean FIMA error 2.220e-16; worst fission-sum relative error 2.665e-15.

Native radial/axial trends retained. Peak change 0.0000000%; maximum re-binned radial-profile difference 0.6844338%; maximum individual re-binned domain difference 0.7943175%. Integrals are conserved but mixed radial cells average bin jumps.

C uniform: PASS. Analytic free-swelling displacement error 4.978e-10 m; difference from the existing constant volumetric-swelling channel 1.494e-10 m; tolerance 1e-9 m.

D zero: PASS. Temperature difference from no-swelling baseline 1.603e-11 K; displacement difference 9.216e-19 m. Internal BU remains nonzero, proving imported zero does not fall back to BU-driven swelling.

A–D were PASS before the native solve. Native and uniform-equal-mean solves both completed all 21 times. The uniform control has the same mean at **every saved time**, not merely at the final time.

| Final metric | Native FIMA | Uniform same-mean FIMA |
|---|---:|---:|
| Fuel Tmax [K] | 397.633085051 | 397.633085051 |
| Fuel surface Tmin [K] | 360.611780725 | 360.611780725 |
| Fuel surface Tmean [K] | 360.611780725 | 360.611780725 |
| Fuel surface Tmax [K] | 360.611780725 | 360.611780725 |
| Local gap minimum [um] | 118.815596932 | 118.996070604 |
| Local gap mean [um] | 119.449407185 | 119.449407159 |
| Local gap maximum [um] | 120.385038922 | 119.613310195 |
| Model contact pressure [Pa] | 0 | 0 |
| Fuel maximum radial displacement [um] | 14.1558374457 | 13.9791560294 |
| Clad inner radial displacement min [um] | 2.96989674712 | 2.969896747 |
| Clad inner radial displacement mean [um] | 2.97143437141 | 2.97143437128 |
| Clad inner radial displacement max [um] | 2.97522663374 | 2.97522663361 |
| Imported FIMA volume mean | 0.000153060852823 | 0.000153060852823 |
| Maximum isotropic swelling eigenstrain | 0.000223997107275 | 0.000153060852823 |
| Internal baseline BU [MWd/kgU] | 2.52104000729 | 2.52104000729 |
| Fuel cell von Mises max [MPa] | 13.5661082629 | 16.4395457944 |
| Clad cell von Mises max [MPa] | 0.26927836309 | 0.269278362491 |

Native maximum FIMA: 0.00022399710727541, cell bounds [ri,ro,zlo,zhi] in local metres: [0.01716375, 0.01791, 0.1896533333333332, 0.2015066666666665]; corresponding native source domain r10/z6 (27.98–31.536 cm absolute z).

Stress component extrema are DG0 cell-centre samples of total elastic stress (thermal and swelling eigenstress included), not continuum extrema or nodal axis values.

| Region/component | Native min / max [MPa] | Uniform min / max [MPa] |
|---|---:|---:|
| fuel rr | -9.27894842 / 0.52039616 | -10.2797215 / -0.0187747664 |
| fuel hoop | -9.27894842 / 13.9822775 | -10.2797215 / 16.1910243 |
| fuel zz | -15.1439772 / 13.1557218 | -17.3013678 / 15.9824024 |
| fuel rz | -2.79112833 / 0.312941197 | -3.21942716 / 0.0626807391 |
| clad rr | -0.00360634416 / -0.000566882847 | -0.0036063453 / -0.000566883074 |
| clad hoop | -0.271672753 / 0.266182739 | -0.271672755 / 0.266182738 |
| clad zz | -0.269118491 / 0.262187532 | -0.269118493 / 0.262187531 |
| clad rz | -0.00447013059 / 0.000348445334 | -0.00447013058 / 0.000348445336 |

Contact: native=False, uniform=False; local geometric closure native=False, uniform=False.

Resources before each successful solve (OMP=2, MPI=1, BLAS/NumExpr=1; UCX_TLS=self):

| Solve | RAM available [GiB] | Swap used [GiB] | Disk available [GiB] | CPUs | Runtime [s] |
|---|---:|---:|---:|---:|---:|
| uniform_analytic | 8.2908 | 0.6947 | 923.1067 | 8 | 1.567 |
| uniform_constant_reference | 8.2798 | 0.6947 | 923.1058 | 8 | 2.346 |
| baseline_no_swelling | 8.2902 | 0.6946 | 923.1050 | 8 | 4.917 |
| zero_import | 8.2406 | 0.6946 | 923.0994 | 8 | 4.061 |
| uniform_mean | 8.2398 | 0.6946 | 923.0932 | 8 | 16.433 |
| native | 8.3267 | 0.6945 | 923.0866 | 8 | 18.193 |

Protected manifest: PASS, 538 unchanged files, including all neutron outputs and both original regression trees. Four FE operator files verified identical to the earlier audit hashes.

FE binding tests: PASS. Imported FIMA is independent of BU and q; diagonal swelling equals FIMA, off-diagonal swelling is zero, clad FIMA is zero, adaptive rollback restores coefficient/time state, out-of-range requests fail without state changes, and a missing native field raises instead of falling back to BU. Legacy BU mode remains available.

Approximations/limits:

- DG0 cylindrical overlap averaging on axis-aligned first-order quads; no general unstructured remapper.
- Linear interpolation of cumulative fissions between saved depletion times, no extrapolation.
- Fixed baseline thermal LHR 8780 W/m, insulated ends, outer clad 300 K; no OpenMC power transfer.
- Gap/contact remain mean-interface models; local geometric gaps are diagnostic only.
- U-ZrH swelling correlation unchanged and not independently calibrated by this verification.
- Point stress exports at 31 r=0 axis nodes have u_r/r singularities and are zeroed by the existing writer; reported stress extrema use finite DG0 cell samples.

First uniform solve converged but diagnostics used cells instead of cells0 in FEniCSx 0.11; fixed only diagnostics, retained failed attempt, reran with preflight. The initial postprocessing invocation preceded script creation; resumed postprocessing from saved outputs without rerunning solves.

**Phase 1 verified within the documented data-transfer, baseline thermal and averaged-gap assumptions.** Native heterogeneity produces the recorded mechanical differences from the equal-mean control; these are deterministic FE responses to a Monte Carlo input field, not a new statistical-significance certification. Phase 2 power coupling is not implemented or started.

Spatial maps: phase1_FIMA_swelling_maps.png/pdf. Full cell maps and 21-time fields: runs/native/output/fields_*.vtu. Full diagnostics and comparison: phase1_results.json.

No OpenMC/depletion execution, neutron/regression edits, commit or push.
