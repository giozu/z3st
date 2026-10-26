# VP1 gradient across the flats: free bow of VP3B duct 1

The VP1 duct and inputs of `../` with the gradient turned to point at a flat,
as on duct 1 of IAEA VP3B. Without neighbours this is the free bow that the
VP3B restraint reduces, so it is the upper bound for `../../vp3`.

```bash
conda activate z3st
../Allrun flats             # mesh, solve, ../plots_vp1.py flats
python ../view_vp1.py --case flats
python ../animate_vp1.py --case flats
```

## 1. Data

### From the literature

ANL/NSE-22/38 (`sandbox/U_hex_bowing/0.biblio/ANL.pdf`), §2.2, p. 13: in VP3B
"the center duct has a thermal gradient similar to VP3A but the orientation
is rotated such that the duct bows across the flats. All the other ducts
remain at 400 °C." The VP3A gradient is the VP1 gradient (§2.1, Fig. 2-2).
The values across the flats are not stated; §2.4 reuses "the temperature
gradient equivalent to Duct 1 from VP3B" without giving them either.

### Interpretation

The same field as VP1, 50 K / 0.150 m = 333.3 K/m, along the flat normal y:

```
T(y, z) = 400 + r(z) · (125 + 333.3 · y)     [°C],   r(z) = clip((z − 1.5) / 1.0, 0, 1)
```

At the top of the core this gives 525 °C on the axis, 546.1 / 503.9 °C at
mid-wall of the +y / −y flats and 546.6 / 503.4 °C on the outer flat faces,
i.e. 43.3 K across the flats instead of the 50 K across the corners. This is
variant A ("the benchmark as specified") of the VP3B design
(`docs/superpowers/specs/2026-10-13-vp3b-design.md`), which notes that the
maximum of ANL Fig. 2-14, about 548 °C, matches the pad face of this field.
The alternative, 550 / 500 °C on the two flats (394 K/m at mid-wall), is the
gradient of `../flats_1515Ti`.

Material: `../steel_iaea.yaml` (α = 1.8·10⁻⁵ 1/K, IWGFR/75 Example 1;
T_ref = 400 °C). See `../README.md`, §1.

## 2. Implementation in z3st

As `../README.md`, §2, with:

- `geometry.yaml`: `orientation: flat` (flats at ±y; hot flat +y, cold −y);
- `boundary_conditions.yaml`: the field above on `inner` and `outer`,
  clamp at `bottom`;
- `input.yaml`: material `../steel_iaea.yaml`, output
  `output/fields_vp1_flats.{xdmf,h5}`.

`../Allrun flats` meshes `../mesh.geo` with this `geometry.yaml` into
`mesh.msh` here and runs z3st from this folder.

## 3. Results

| Location | z (m) | z3st section mean | z3st corner mean | Beam | Diff. |
|---|---|---|---|---|---|
| core top | 2.5 | −0.979 mm | −0.979 mm | −1.000 mm | −2.1 % |
| ACLP | 3.0 | −3.211 mm | −3.211 mm | −3.250 mm | −1.2 % |
| TLP (top) | 4.0 | **−12.134 mm** | −12.136 mm | **−12.250 mm** | −0.9 % |

The result is the VP1 result of `../` to four digits: a regular hexagon has
the same second moment of area about every axis, so the same gradient bows
the duct by the same amount across the flats as across the corners.

Against VP3B (ANL/NSE-22/38, Table 2.4), duct 1:

| z | Free, z3st (here) | Restrained, IAEA | Restrained, MOOSE |
|---|---|---|---|
| ACLP (3 m) | 3.21 mm | 0.61 mm | 0.68 mm |
| TLP (4 m) | 12.13 mm | 8.10 mm | 8.09 mm |

The neighbours take about a third of the free tip bow.

![profiles](output/vp1_profiles.png)
![3d](output/vp1_3d.png)

## 4. Files

| File | Content |
|---|---|
| `geometry.yaml` | VP1 dimensions, `orientation: flat` |
| `boundary_conditions.yaml` | VP1 gradient along the flat normal, clamp at z = 0 |
| `input.yaml` | as `../input.yaml`, material `../steel_iaea.yaml` |
| `mesh.msh`, `log_*.md`, `output/` | written by `../Allrun flats` |
