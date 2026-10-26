# VP1-style free bow, 15-15Ti, 50 °C across the flats

The first version of the VP1 case: the VP1 duct with 550 / 500 °C on the two
flats and the 15-15Ti material card. It is not the IAEA problem (see `../`
for that) and has no published reference; it checks z3st against the beam
solution with a real cladding steel. It is variant B of the VP3B design
(`docs/superpowers/specs/2026-10-13-vp3b-design.md`).

```bash
conda activate z3st
../Allrun flats_1515Ti      # mesh, solve, ../plots_vp1.py flats_1515Ti
python ../view_vp1.py --case flats_1515Ti
python ../animate_vp1.py --case flats_1515Ti
```

## 1. Data

| Quantity | Value | Source |
|---|---|---|
| Duct geometry, core elevations | as `../` | ANL/NSE-22/38, Table 2.1 |
| T profile | 400 °C below the core, linear rise over the core, constant above, 50 °C hot − cold | ANL/NSE-22/38, Fig. 2-2 |
| Hot / cold side | +y / −y **flats**, 550 / 500 °C at mid-wall | this case |
| Material | `15_15Ti_local.yaml`: E = 170 GPa, ν = 0.3, α = 1.2·10⁻⁵ 1/K, T_ref = 300 K, k(T) from `k_1515.py` | local 15-15Ti card |

In ANL Fig. 2-2 the 550 / 500 °C pair sits on two opposite **corners**; here
it is applied to the **flats**. The gradient is then 50 °C over
2a = 126.9 mm (a = 63.45 mm, mid-wall apothem), 394 K/m, against 333 K/m in
VP1.

```
T(y, z) = 673.15 + r(z) · (125 + 25 · y / a)     [K]
r(z)    = clip((z − 1.5) / 1.0, 0, 1)            axial ramp over the core
a       = 0.0634519 m                            mid-wall apothem
```

| Region | z (m) | T, +y flat | T, −y flat |
|---|---|---|---|
| below core | 0 – 1.5 | 400 °C | 400 °C |
| core | 1.5 – 2.5 | linear 400 → 550 °C | linear 400 → 500 °C |
| above core | 2.5 – 4.0 | 550 °C | 500 °C |

## 2. Implementation in z3st

As `../README.md`, §2, with:

- `geometry.yaml`: `orientation: flat` (flats at ±y);
- `boundary_conditions.yaml`: the field above on `inner` and `outer`, clamp
  at `bottom`;
- `15_15Ti_local.yaml`, `k_1515.py`: k(T) is evaluated but does not change
  T, which is imposed on both faces. It makes the thermal solve nonlinear:
  the staggered loop takes 3 iterations on the first step, 1 afterwards;
- `input.yaml`: output `output/fields_vp1_1515Ti.{xdmf,h5}`.

**Animation.** With T_ref = 300 K the duct also expands uniformly from 300 K
to 400 °C, which is not part of the bow. `../animate_vp1.py` removes it with
the solution at uniform 400 °C: `iso/` is this case with
`temperature: "673.15 + 0.0 * z"` on both faces, run by `../Allrun
flats_1515Ti` after the main solve (with `PYTHONPATH=..` so that z3st finds
`k_1515.py`), output `iso/output/fields_vp1_1515Ti_iso.h5`.

## 3. Results

### Temperature field

| z (m) | T +y flat | T −y flat | ΔT across flats (mid-wall) | T max | T min |
|---|---|---|---|---|---|
| 0 – 1.5 | 400.0 °C | 400.0 °C | 0.0 °C | 400.0 | 400.0 |
| 2.0 (mid core) | 475.0 | 450.0 | 25.0 | 475.3 | 449.7 |
| 2.5 – 4.0 | 550.0 | 500.0 | **50.0** | 550.6 | 499.4 |

The extremes 550.6 / 499.4 °C are on the outer surface, y = ±64.95 mm.

### Deflection

```
κ(z) = α · dT/dy = α · 25 · r(z) / a        (κ = 4.73·10⁻³ 1/m above the core)
u_y(z) = −∫₀ᶻ κ(s) (z − s) ds               (clamped at z = 0)
```

| Location | z (m) | z3st section mean | z3st corner mean | Beam | Diff. |
|---|---|---|---|---|---|
| core top | 2.5 | −0.772 mm | −0.772 mm | −0.788 mm | −2.0 % |
| ACLP | 3.0 | −2.530 mm | −2.530 mm | −2.561 mm | −1.2 % |
| TLP (top) | 4.0 | **−9.562 mm** | −9.563 mm | **−9.653 mm** | **−0.9 %** |

The top rises by u_z = 21.0 mm (expansion from T_ref = 300 K).

![profiles](output/vp1_profiles.png)
![3d](output/vp1_3d.png)

### Stress

The median von Mises stress is 0.17 MPa. The first cell row at the clamp
reaches 669 MPa (z = 12.5 mm): the clamp holds the nozzle at its 300 K size
while the steel there is at 400 °C. This is an artefact of T_ref = 300 K with
a fully clamped bottom; `../` uses T_ref = 400 °C and has no such stress.

## 4. Files

| File | Content |
|---|---|
| `geometry.yaml` | VP1 dimensions, `orientation: flat` |
| `boundary_conditions.yaml` | 550 / 500 °C across the flats, clamp at z = 0 |
| `input.yaml` | steady thermo-mechanical solve, output `fields_vp1_1515Ti` |
| `15_15Ti_local.yaml`, `k_1515.py` | material card and k(T) |
| `iso/` | the same case at uniform 400 °C (`input.yaml`, `boundary_conditions.yaml`), for the animation |
| `mesh.msh`, `log_*.md`, `output/` | written by `../Allrun flats_1515Ti` |
