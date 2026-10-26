# VP2: free thermal bow of a single hexagonal duct, non-linear gradient

IAEA Verification Problem 2 = **Example 2** of IWGFR/75 Part I (1990): the
VP1 duct ([`../vp1`](../vp1/README.md)), clamped at the bottom and free at the
top, with a non-linear temperature distribution across the duct: one half
hot, the other half cold. The bow has a closed-form beam solution; the
non-linear part of T only adds self-equilibrated stresses, which do not bend
the duct. The tip bow is compared with the beam and with the eleven
participants of the IAEA benchmark.

```bash
conda activate z3st
./Allrun                    # mesh, solve, plots_vp2.py
python show_mesh.py         # mesh, faces coloured by tag
```

| Tip bow, z = 4 m | u_y (mm) |
|---|---|
| z3st | **−16.17** |
| IWGFR/75 linear theory | −16.333 (z3st −1.0 %) |
| IWGFR/75, 11 codes | −16.020 … −16.444 |

## 1. Data (IWGFR/75 Part I, `0.biblio/22012674.pdf`)

| Quantity | Value | Source |
|---|---|---|
| Specification | "As for EXAMPLE 1 but with the non-linear temperature distribution shown for EXAMPLE 2 of Fig. 1" | App. 2, p. 108 |
| Duct | sharp-cornered hexagon, **150 mm between mid-wall corners**, wall 3 mm, length 4000 mm, at 20 °C | App. 2, Ex. 1, p. 107 |
| Core / LRP / URP | z = 1500–2500 / 3000 / 4000 mm | App. 2, Ex. 1, p. 107 |
| T below the core | 400 °C | Fig. 1, p. 12 |
| T rise | linear over the core, constant above | Fig. 1 |
| T at the top of the core | corners 3, 4, 5: 550 °C; corners 6, 1, 2: 500 °C; linear between corners | Fig. 1 |
| α, E, ν | 18·10⁻⁶ /°C, 1.7·10⁵ MPa, 0.3 | App. 2, Ex. 1, p. 108 |
| Bow, linear theory | core top 1.333, LRP 4.333, **URP 16.333 mm** | DT1BOW, p. 51 |
| Bow, 11 codes at the URP | 16.020 (Japan, Mitsubishi) … 16.444 mm (Japan, NAIG/Toshiba), spread 1.6 % | DT1BOW, p. 51; §22, p. 13 |

Same values in `reference/iwgfr75_ex2.csv`.

### Temperature field

Corners 1 and 4 on the y axis (−y cold, +y hot), corners 2, 6 and 3, 5 at
y = ∓R/2 with R = 75 mm the mid-wall corner radius. T is constant on the
four flats between corners of the same half and linear in y on the two side
flats 2–3 and 5–6:

```
T(y, z) = 400 + r(z) · (125 + 25 · clip(y / 0.0375, −1, 1))     [°C]
r(z)    = clip((z − 1.5) / 1.0, 0, 1)
```

Check against the worked example of Fig. 1: at y = 20 mm, z = 2.5 m,
T = 525 + 25 · 20 / 37.5 = 538.33 °C.

### Beam solution

T is not linear in y, so the curvature is the section average (plane
sections):

```
κ = α ∫ T y dA / ∫ y² dA
```

On the thin-walled hexagon this gives **4/3 of Example 1** (both on the
contour and lumped on the six corners), so all bows are 4/3 of VP1's: 1.333,
4.333, 16.333 mm. `plots_vp2.py` integrates κ over the meshed ring:
1.3332 × Example 1, tip −16.331 mm.

## 2. Implementation in z3st

As VP1 (same `mesh.geo`, `input.yaml`, `steel_iaea.yaml`, clamp), except:

1. **Geometry** (`geometry.yaml`). `D = 0.1534641` m outer corner to corner,
   i.e. 150 mm between mid-wall corners as the spec (D = 150 + 2 t/√3 mm).
   VP1 uses 150 mm outer, as ANL; for VP1's linear T the bow does not depend
   on the section size, for VP2 it does: with VP1's geometry the beam gives
   1.361 × Example 1 (−16.67 mm, +2.1 %), because the knees of T at
   y = ±37.5 mm then sit outside the mid-wall corners (±36.6 mm).
2. **Temperature** (`boundary_conditions.yaml`). The field above, in K, on
   `inner` and `outer`.
3. **Post-processing** (`plots_vp2.py`). T across the section at the top of
   the core against Fig. 1; section-mean u_y along z against the beam and
   DT1BOW: `output/vp2_profiles.png`, `output/vp2_3d.png`.

Mesh: 10 per flat × 2 through the wall × 160 along z (25 mm), 19 200
hexahedra, as VP1.

## 3. Results

| Location | z (m) | z3st section mean | Beam (meshed section) | IWGFR/75 linear theory | z3st vs. IWGFR/75 |
|---|---|---|---|---|---|
| core top | 2.5 | −1.304 mm | −1.333 mm | −1.333 mm | −2.2 % |
| LRP | 3.0 | −4.277 mm | −4.333 mm | −4.333 mm | −1.3 % |
| URP (top) | 4.0 | **−16.166 mm** | −16.331 mm | **−16.333 mm** | **−1.0 %** |

The offsets against the beam are those of VP1 (−2.1, −1.2, −0.9 %), not
specific to the non-linear profile. z3st lies within the range of the IAEA
codes (−16.020 … −16.444 mm).

![profiles](output/vp2_profiles.png)
![3d](output/vp2_3d.png)

## 4. Files

| File | Content |
|---|---|
| `geometry.yaml` | duct, 150 mm between mid-wall corners, `orientation: vertical`, mesh density |
| `boundary_conditions.yaml` | non-linear T(y, z) on inner/outer faces, clamp at z = 0 |
| `input.yaml` | steady thermo-mechanical solve, linear solvers (as VP1) |
| `steel_iaea.yaml` | material card (IWGFR/75, T_ref = 400 °C, as VP1) |
| `mesh.geo` | hexagonal shell (as VP1), parameters from `geometry.yaml` |
| `reference/iwgfr75_ex2.csv` | corner temperatures, linear theory and the codes' range |
| `Allrun` | gmsh, z3st, `plots_vp2.py` |
| `plots_vp2.py` | T and deflection against the beam and DT1BOW, 3D figure |
| `show_mesh.py` | display `mesh.msh` with pyvista, faces coloured by tag |
