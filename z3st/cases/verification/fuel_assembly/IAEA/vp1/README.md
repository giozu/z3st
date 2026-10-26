# VP1: free thermal bow of a single hexagonal duct

IAEA Verification Problem 1 (IWGFR/75, 1990), as stated in ANL/NSE-21/30,
§2.5: one fuel-assembly duct, clamped at the nozzle, free at the top, with a
temperature difference applied **corner to corner** over the top half. The
duct bows away from the hot corner; with nothing to push against, the bow is
stress-free and has a closed-form beam solution. The tip bow is compared with
the IAEA analytical value and with MOOSE.

```bash
conda activate z3st
./Allrun                    # IAEA VP1 (this folder): mesh, solve, plots_vp1.py
./Allrun flats              # same gradient across the flats, free bow of VP3B duct 1
./Allrun flats_1515Ti       # 50 °C across the flats, 15-15Ti card
python view_vp1.py          # interactive pyvista viewer (--case, see --help)
python animate_vp1.py       # gradient applied and removed, output/vp1_bow.gif
python show_mesh.py         # mesh, faces coloured by tag (flats/mesh.msh for a variant)
```

| Folder | Gradient | α (1/K) | T_ref | Tip bow, z3st | Reference |
|---|---|---|---|---|---|
| `.` (VP1) | corner to corner, 50 K / 150 mm | 1.8·10⁻⁵ | 400 °C | −12.13 mm | IAEA −12.25, MOOSE −12.21 mm |
| [`flats/`](flats/README.md) | across the flats, same 333 K/m | 1.8·10⁻⁵ | 400 °C | −12.13 mm | beam −12.25 mm; upper bound for VP3B |
| [`flats_1515Ti/`](flats_1515Ti/README.md) | across the flats, 50 K / 126.9 mm | 1.2·10⁻⁵ | 300 K | −9.56 mm | beam −9.65 mm |

## 1. Data

### Taken from the literature

| Quantity | Value | Source |
|---|---|---|
| Duct corner-to-corner (outer) / wall / length | 150 / 3 / 4000 mm | ANL/NSE-22/38, Table 2.1; ANL/NSE-21/30, §2.5 |
| Core bottom / top | z = 1500 / 2500 mm | ANL/NSE-22/38, Table 2.1 |
| Support | clamped at the bottom, free top | ANL/NSE-21/30, §2.5 |
| T below the core | 400 °C on all corners | ANL/NSE-21/30, §2.5 |
| T rise | linear over the core, constant above | ANL/NSE-21/30, Fig. 2.20 |
| T at the top of the core | corner 4: 550 °C, corners 3, 5: 537.5 °C, corners 2, 6: 512.5 °C, corner 1: 500 °C | ANL/NSE-21/30, §2.5, Fig. 2.21 |
| Tip bow, IAEA analytical | 0.4823 in = **12.25 mm** | ANL/NSE-21/30, §2.5, Table 2.2 |
| Tip bow, MOOSE | 0.4807 in = **12.21 mm** | ANL/NSE-21/30, §2.5, Table 2.2 |

The same values are in `reference/anl_nse-21-30_vp1.csv`.

![IAEA VP1 axial thermal gradient](figures/anl_fig2-2_vp1_thermal_gradient.png)

*ANL/NSE-22/38, Fig. 2-2: VP1 thermal gradient (duct 1 of VP3A, corner to corner).*

### Temperature field

The corner temperatures are linear along the corner 1 → corner 4 axis:
corners 2, 3, 5, 6 sit at ±D/4 along it, and 525 ± 25 · (37.5 / 75) =
537.5 / 512.5 °C. With y along that axis,

```
T(y, z) = 400 + r(z) · (125 + 25 · y / 0.075)     [°C]
r(z)    = clip((z − 1.5) / 1.0, 0, 1)             axial ramp over the core
```

0.075 m is the outer corner radius, so dT/dy = 50 K / 0.150 m = 333.3 K/m.

### Material and beam solution

IWGFR/75 Part I, Appendix 2, Example 1 (p. 108) gives α = 18·10⁻⁶ /°C,
E = 1.7·10⁵ MPa, ν = 0.3. The free bow depends only on α, through the
curvature κ = α · dT/dy. For the clamped cantilever with this ramp,

```
u_tip = ∫₀ᴸ κ(s) (L − s) ds = κ₀ [ ∫₀¹ r (2.5 − r) dr + 1.5² / 2 ] = κ₀ · 2.042 m²
```

with κ₀ = α · 333.3 K/m = 6.00·10⁻³ 1/m: u_tip = 12.25 mm, the spec's
"elementary integration of the curvature" (core top 1 mm, LRP 3.25 mm,
URP 12.25 mm).

The spec puts 150 mm between the **mid-wall** corners, at 20 °C; this case,
as ANL, uses 150 mm **outer** corner to corner. For T linear in y the
curvature, and so the free bow, does not depend on the section size; for the
non-linear profile of VP2 it does (+2 %), and [`../vp2`](../vp2/README.md)
follows the spec.

| Property | Value | Effect on the free bow |
|---|---|---|
| α | 1.8·10⁻⁵ 1/K (IWGFR/75) | proportional |
| T_ref | 673.15 K = 400 °C | none; the unheated duct is strain-free, as ANL takes dimensions at 400 °C |
| E, ν | 170 GPa, 0.3 | none (no restraint) |
| k | 20 W/m·K | none (T imposed on both faces) |

Card: `steel_iaea.yaml`.

## 2. Implementation in z3st

1. **Geometry** (`geometry.yaml`). Dimensions above, `orientation: vertical`:
   `mesh.geo` puts the hexagon vertices at 90°, 150°, …, so corners 4 and 1
   lie on +y and −y and the bow is in the y–z plane.
2. **Mesh.** `Allrun` turns `geometry.yaml` into gmsh flags
   (`python3 -m z3st.utils.geo_args geometry.yaml`) for `mesh.geo`, which
   builds the hexagonal shell as six transfinite sectors extruded along z:
   8-node hexahedra, 10 per flat × 2 through the wall × 160 along z (25 mm),
   19 200 cells, 28 980 nodes. The surface tags (`bottom`, `top`, `inner`,
   `outer`) match `labels` in `geometry.yaml`, which is how z3st finds the
   boundary regions.
3. **Temperature** (`boundary_conditions.yaml`). The field above, in K, as a
   string expression in x, y, z on the `inner` and `outer` faces. z3st
   interpolates a string Dirichlet temperature onto the mesh
   (`models/thermal_model.py`); it is fixed in time. Both faces carry the same
   expression, so T is linear in y through the whole section.
4. **Clamp.** `bottom` (z = 0): u = 0. All other surfaces are traction-free.
5. **Material** (`steel_iaea.yaml`), see §1.
6. **Solve** (`input.yaml`). Steady state (`time: [0, 1]`, one step), thermal
   and mechanical models on, linear solvers with MUMPS. Output
   `output/fields_vp1.{xdmf,h5}`: T, u (nodal); stress, strain, heat flux
   (cell).
7. **Post-processing** (`plots_vp1.py`). Section-mean u_y along z against the
   beam formula and the two published tip values:
   `output/vp1_profiles.png`, `output/vp1_3d.png`.

## 3. Results

### Deflection

| Location | z (m) | z3st section mean | z3st corner mean | Beam | Diff. |
|---|---|---|---|---|---|
| core top | 2.5 | −0.979 mm | −0.979 mm | −1.000 mm | −2.1 % |
| ACLP | 3.0 | −3.211 mm | −3.211 mm | −3.250 mm | −1.2 % |
| TLP (top) | 4.0 | **−12.134 mm** | −12.136 mm | **−12.250 mm** | **−0.9 %** |

| Tip bow, z = 4 m | u_y (mm) | z3st vs. |
|---|---|---|
| IAEA analytical | −12.25 | −0.9 % |
| MOOSE | −12.21 | −0.6 % |

*Corner mean* = average of the 6 outer corners, the centreline estimate used
by ANL for shell models. MOOSE is 0.34 % below the analytical value, z3st
0.9 %. The remaining ~1 % is attributed to the thin-walled cross-section,
which the plane-section beam formula does not represent (not verified by a
mesh study). Only the tip value is published, not a profile along z.

The 12.25 mm is the beam solution for the IAEA α; the comparison with MOOSE
is a code-to-code check on the same problem. IWGFR/75 (Data Table DT1BOW,
p. 51) gives the 11 participants' tip bows, 12.015–12.333 mm, average
12.24 mm: z3st, −12.13 mm, is inside the range.

The section mean of u_x is zero (symmetric about x = 0); locally u_x reaches
0.16 mm from the thermal expansion of the section. The top rises by
u_z = 4.50 mm (mean 525 °C above the core against T_ref = 400 °C).

![profiles](output/vp1_profiles.png)
![3d](output/vp1_3d.png)

### Stress

A free thermal bow with T linear in y produces no bending stress: the median
von Mises stress is 0.18 MPa, the maximum 2.6 MPa at the core bottom
(z = 1.51 m), where the axial gradient of T has a kink. With T_ref = 400 °C
the clamp carries no stress.

### Animation (`animate_vp1.py`)

The problem is linear in T, so a fraction s of the gradient gives
T = 400 + s · (T − 400) and u = s · u (u = 0 at uniform 400 °C = T_ref).
The script ramps s 0 → 1 → 0, draws the bow magnified and writes
`output/vp1_bow.gif`. The bow returns to zero because the model is elastic,
with no creep or swelling.

## 4. Literature

### Available in `sandbox/U_hex_bowing/0.biblio`

**IWGFR/75 Part I** (`22012674.pdf`): *Verification and Validation of LMFBR
Static Core Mechanics Codes, Part I*, IAEA, Vienna, 1990. Example 1 = VP1:
specification and material in Appendix 2 (p. 107–108), temperatures in
Fig. 1 (p. 12), participants' results in Data Table DT1BOW (p. 51).

**ANL/NSE-21/30** (`169472.pdf`): N. Wozniak, E. Shemon, J. Grudzinski,
B. Spencer, *Assessment of MOOSE-Based Tools for Calculating Radial Core
Expansion*, Argonne, 2021, doi:10.2172/1808314. §2.5 (Example 5, report
pp. 18–20) states VP1 with the corner temperatures; Table 2.2 (p. 24) gives
the IAEA analytical and MOOSE tip bow.

**ANL/NSE-22/38** (`ANL.pdf`): *Continued Verification of MOOSE Structural
Mechanics Tools for Modeling Core Bowing Phenomena in Fast Reactors*, Argonne,
2022. Common duct geometry (Table 2.1) and the VP1 gradient (Fig. 2-2).

| Common geometry (Table 2.1) | Value |
|---|---|
| Duct corner-to-corner / length / wall | 150 / 4000 / 3 mm |
| Duct-to-duct gap | 6 mm |
| ACLP pad thickness / ACLP gap | 2.75 / 0.5 mm |
| TLP pad thickness | N/A (TLP gap = 6 mm, as the duct wall) |
| Duct pitch (VP3A) | 138.9 mm |
| Core bottom / top / ACLP / TLP | 1500 / 2500 / 3000 / 4000 mm |

The restrained follow-ups VP3A and VP3B and their ANL/IAEA results are in
[`../vp3/README.md`](../vp3/README.md).

### Cited but not in `0.biblio`

1. J. J. Grudzinski, T. Moran, C. Grandy, *Supplement to the NUBOW-3D Manual*,
   ANL/NE-15/9 Rev. 1, 2017.

## 5. Files

| File | Content |
|---|---|
| `geometry.yaml` | VP1 duct dimensions, `orientation: vertical`, mesh density |
| `boundary_conditions.yaml` | corner-to-corner T(y, z) on inner/outer faces, clamp at z = 0 |
| `input.yaml` | steady thermo-mechanical solve, linear solvers |
| `steel_iaea.yaml` | material card (IWGFR/75, T_ref = 400 °C), also used by `flats/` |
| `mesh.geo` | hexagonal shell, parameters from the variant's `geometry.yaml` |
| `reference/anl_nse-21-30_vp1.csv` | corner temperatures and tip bow from ANL/NSE-21/30 |
| `Allrun` | gmsh, z3st, `plots_vp1.py`, for this folder or a variant |
| `plots_vp1.py` | profiles against the beam solution and the references, 3D figure |
| `view_vp1.py` | interactive pyvista viewer, `.vtu` export |
| `animate_vp1.py` | GIF of the bow as the gradient is applied and removed |
| `show_mesh.py` | display `mesh.msh` (or a variant's, as argument) with pyvista, faces coloured by tag |
| `figures/` | figures taken from the literature |
| `flats/`, `flats_1515Ti/` | variants, each with its own `README.md` |

## Next steps

- VP2 (`../vp2`): the same duct, non-linear gradient (hot half, cold half).
- VP3 (`../vp3`): the same duct and gradient, restrained by neighbouring
  ducts at the ACLP and TLP.
