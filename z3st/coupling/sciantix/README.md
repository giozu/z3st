# Z3ST ↔ SCIANTIX coupling (prototype)

Drives SCIANTIX (mesoscale fission-gas behaviour) from Z3ST to compute gaseous
swelling (→ the eigenstrain bus) and fission gas release per fuel point.
The binding wraps the C-linkage coupling entry of SCIANTIX. No SCIANTIX
physics is reimplemented in Python.

## 1. Build SCIANTIX as a shared library

This is the only supported recipe. Run it from the SCIANTIX repository root:

```bash
cd <sciantix>
mkdir -p build
g++ -O2 -std=c++17 -DCOUPLING_TU -fPIC -shared $(find include -type d | sed 's/^/-I/') \
    $(find src -name '*.C') -o build/libsciantix_tu.so
export SCIANTIX_LIB=$PWD/build/libsciantix_tu.so   # add this line to ~/.bashrc
```

### `Allmake.sh` and CMake do not produce this library

The SCIANTIX `CMakeLists.txt` builds an executable by default and, under `COUPLING_TU`, a
static library (`add_library(sciantix STATIC ${SOURCES})` with
`CMAKE_STATIC_LIBRARY_SUFFIX ".a"`). The binding loads a `.so` through `ctypes`
and cannot use a `.a`. Use the `g++` line above.

### `-DCOUPLING_TU` sets the physics

The `extern "C"` entry points live in `src/coupling/TUSrcCoupling.C` and are
exported unconditionally. A library built without the macro loads and every
call returns, with different physics. With the macro:
- `Simulation::execute()` skips `Burnup()`, `EffectiveBurnup()` and
  `Densification()` (Z3ST owns those).
- `SetVariables.C` takes burnup from `Sciantix_history[7]/[8]`, i.e. from Z3ST.

Without it, SCIANTIX computes its own burnup and ignores the value Z3ST passes in.
This command confirms the exported symbols (it does not detect the macro):
```bash
nm -D --defined-only build/libsciantix_tu.so | grep -E 'callSciantix|getSciantixOptions'
```

### Library path

`SCIANTIX_LIB` is read by `sciantix_binding.py` and nothing else resolves it.
A library under `/tmp` is removed on reboot and the next run fails with
`OSError: cannot open shared object file`. Put the library on a persistent path
and the `export` in `~/.bashrc`.

## 2. Validate the binding (standalone)

```bash
python3 smoke_test.py     # ramps T at fixed fission rate, prints swelling/bu/FGR
```
Compare the output with a SCIANTIX standalone run using the same
`input_history.txt` (same T, fission rate, dt). The two must match.

## 3. Array layout (verified against the SCIANTIX v2.2.1 source)

`include/MainVariables.h`: `options[40]`, `history[20]`, `variables[300]`,
`scaling_factors[20]`, `diffusion_modes[720]` (= 18 mode blocks × 40 modes).

| host writes: `history[]` | idx | reads: `variables[]` | idx |
|---|---|---|---|
| Temperature old/new (K) | 0,1 | Xe produced (at/m³) | 1 |
| Fission rate old/new (fiss/m³s) | 2,3 | Xe released (at/m³) | 6 |
| Hydrostatic stress old/new (MPa) | 4,5 | intragranular gas swelling (/) | 24 |
| time step Δt (s) | 6 | intergranular gas swelling (/) | 36 |
| steam pressure old/new (atm) | 9,10 | Burnup (MWd/kgUO₂) | 38 |

Sources: `src/operations/SetVariablesFunctions.C` (history and variable slots),
`src/operations/SetVariables.C:47` (`history[6]` → `physics_variable["Time step"]`,
seconds). The models integrate on the time step at `history[6]`.

### Burnup ownership (`history[7]`/`[8]`)

In a plain build these two slots hold time (h) and step number and are output only.
In a `-DCOUPLING_TU` build SCIANTIX skips its own `Burnup()`, `EffectiveBurnup()` and
`Densification()` (`Simulation.C:43`) and reads burnup from `history[7]` (old) and
`history[8]` (new) (`SetVariables.C:74`). Z3ST computes burnup with its RADAR
model and feeds it in. `advance(..., burnup_old=, burnup_new=)` writes those slots.
`spine.update_state` passes the per-dof burnup pair, converted from MWd/kgU to MWd/kgUO₂.

## 4. Z3ST integration (the eigenstrain bus)

SCIANTIX gaseous swelling is a numerical, stateful per-point field, not a UFL
expression. It is carried on the state bus, like burnup and creep. Default off.

1. `SciantixField` (in `sciantix_binding.py`) holds one SCIANTIX point per `V_t`
   dof of the fissile region. The library and model settings are read once and shared.
2. `spine.initialize_fields` builds the field (when `models.fission_gas.enabled`)
   and a `gas_swelling` Function on `V_t`. `spine.update_state(dt)` calls
   `field.step(dt, T, fission_rate, burnup_old=..., burnup_new=...)` with `T` from the
   temperature field, `fission_rate = q''' / E_fission`, and the host burnup pair
   (§3). The hydrostatic stress argument of `step` is not passed. The returned ΔV/V
   is written into `gas_swelling`.
3. `materials/sciantix_swelling.py::gaseous_swelling` is the eigenstrain callable.
   It returns `(gas_swelling/3)·I` from that field. A fuel card opts in with
   `eigenstrain: materials.sciantix_swelling.gaseous_swelling`.

The field has `snapshot()`/`restore()` for adaptive-timestep rollback, called from
`spine.snapshot_state`/`restore_state`. Config:

```yaml
models:
  fission_gas:
    enabled: true
    lib: /path/to/libsciantix_tu.so      # else $SCIANTIX_LIB ; build with -DCOUPLING_TU
    initial_conditions: input_initial_conditions.txt
    energy_per_fission: 3.2e-11          # J/fission (≈ 200 MeV)
```
The run directory needs `input_settings.txt` and `input_initial_conditions.txt` (the
files a SCIANTIX standalone run uses).

## 5. Initial conditions in coupling mode

Both handled in `load_initial_conditions`:
1. In coupling mode SCIANTIX does not read `input_initial_conditions.txt`
   (standalone only, `file_manager/InputReading.C`). The host seeds `variables[]`.
2. The standalone one-time `Initialization()` (`file_manager/Initialization.C`),
   skipped by the coupling entry, sets grain-boundary defaults absent from the
   IC file (`variables[25]`=2e13, `[35]`=0.5, `[37]`=1.0) and converts U% → at/m³.
   Without the grain-boundary defaults the intergranular model returns `nan` and
   releases nothing.

## 6. Validation against the SCIANTIX Baker gold

Binding written against SCIANTIX 2.2.1. Build the library as in §1, then run:
```bash
cd <sciantix>/regression/baker/test_Baker1977__1273K
SCIANTIX_LIB=<sciantix>/build/libsciantix_tu.so \
  PYTHONPATH=<z3st>/z3st/coupling/sciantix python3 \
  <z3st>/z3st/coupling/sciantix/validate_baker.py
```
The four engineering outputs match the standalone gold to about 1e-7 relative error:
FGR 0.132097, intragranular swelling 3.07e-4, intergranular swelling 0.0417,
burnup 6.719.

`validate_baker.py` feeds the gold burnup trajectory through `history[7]/[8]`. It passes
against both a plain and a `-DCOUPLING_TU` build. With the coupling build the gas
outputs match to about 1e-7 and burnup matches exactly.

## 7. Case in the suite

`z3st/cases/regression/fg_test_2D` is the PWR rod run with `models.fission_gas.enabled`.
It needs `SCIANTIX_LIB` built as in §1. It carries a gold (`output/non-regression_gold.json`)
and is in the local suite. There are no unit tests of `SciantixField`.

## Limits

- Fresh fuel only: no diffusion-mode projection for a pre-irradiated restart.
- One SCIANTIX point per dof.
- Effective burnup (the SCIANTIX high-burnup-structure input, Khvostov,
  temperature-gated) is computed inside the skipped `EffectiveBurnup()` and is not
  transferred by the coupling. In a `-DCOUPLING_TU` build it stays at 0, so the
  high-burnup-structure model (`iHighBurnupStructureFormation=1`) is not supported.
- `Irradiation time` and `FIMA` are also computed inside the skipped `Burnup()` and
  stay uncomputed in a coupling build. No Z3ST model reads them.
