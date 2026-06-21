# Z3ST ↔ SCIANTIX coupling (prototype)

Drive SCIANTIX (mesoscale fission-gas behaviour) from Z3ST to compute **gaseous
swelling** (→ the eigenstrain bus) and **fission gas release** per fuel point.
The binding wraps SCIANTIX's existing C-linkage coupling entry — **no SCIANTIX
physics is reimplemented in Python**.

## 1. Build SCIANTIX as a shared library

The binding `ctypes`-loads a `.so`, but SCIANTIX currently builds a *static*
lib / executable (`add_library(sciantix STATIC ...)`). Add a shared target that
includes the coupling shim (`src/coupling/TUSrcCoupling.C`, which defines the
`extern "C"` `callSciantix` / `getSciantixOptions`):

```cmake
# in SCIANTIX CMakeLists.txt, alongside the existing targets
add_library(sciantix_shared SHARED ${SOURCES})
set_target_properties(sciantix_shared PROPERTIES OUTPUT_NAME sciantix)
```

Build, then point the binding at it:
```bash
cmake -B build -DCMAKE_POSITION_INDEPENDENT_CODE=ON && cmake --build build
export SCIANTIX_LIB=$PWD/build/libsciantix.so
```
(`callSciantix` is already `extern "C"`, so the symbol is undecorated and
`ctypes`-callable; `-fPIC` is the only build flag the shared lib needs.)

**Build with `-DCOUPLING_TU`.** Build to a persistent path (not `/tmp`, which is
wiped on reboot) and put the export in your shell profile (`~/.bashrc`) so the
binding finds it across sessions:
```bash
g++ -O2 -std=c++17 -DCOUPLING_TU -fPIC -shared $(find include -type d | sed 's/^/-I/') \
    $(find src -name '*.C') -o build/libsciantix_tu.so
export SCIANTIX_LIB=$PWD/build/libsciantix_tu.so   # add this line to ~/.bashrc
```

## 2. Validate the binding (standalone)

```bash
python3 smoke_test.py     # ramps T at fixed fission rate, prints swelling/bu/FGR
```
Diff the output against a SCIANTIX **standalone** run with the same
`input_history.txt` (same T, fission rate, dt) — they must match. That closes the
binding correctness gap before any Z3ST integration.

## 3. Array layout (verified against SCIANTIX v2.2.1, 2026-06-18)

`include/MainVariables.h`: `options[40]`, `history[20]`, `variables[300]`,
`scaling_factors[20]`, `diffusion_modes[720]` (= 18 mode blocks × 40 modes).

| host writes — `history[]` | idx | reads — `variables[]` | idx |
|---|---|---|---|
| Temperature old/new (K) | 0,1 | Xe produced (at/m³) | 1 |
| Fission rate old/new (fiss/m³s) | 2,3 | Xe released (at/m³) | 6 |
| Hydrostatic stress old/new (MPa) | 4,5 | intragranular gas swelling (/) | 24 |
| time step Δt (s) | 6 | **intergranular gas swelling (/)** | **36** |
| steam pressure old/new (atm) | 9,10 | Burnup (MWd/kgUO₂) | 38 |

Sources: `src/operations/SetVariablesFunctions.C` (history + variable slots),
`src/operations/SetVariables.C:47` (`history[6]` → `physics_variable["Time step"]`,
seconds). The time step at `history[6]` is what the models integrate on.

**Burnup ownership (`history[7]`/`[8]`).** In a plain build these two slots are
Time (h) / step-number and are output-only. In a **`-DCOUPLING_TU`** build SCIANTIX
skips its own `Burnup()`/`EffectiveBurnup()`/`Densification()` (`Simulation.C:43`)
and instead **reads burnup from `history[7]` (old) / `history[8]` (new)**
(`SetVariables.C:74`). So the host owns burnup — Z3ST computes it with its RADAR
model and feeds it in. `advance(..., burnup_old=, burnup_new=)` writes those slots;
Z3ST's `spine.update_state` passes the per-dof burnup pair automatically.

## 4. Z3ST integration — IMPLEMENTED (the eigenstrain bus)

SCIANTIX gaseous swelling is a **numerical, stateful per-point field**, not a UFL
expression — so it rides the **state bus**, exactly like burnup/creep. Wired in
2026-06-19 (punch-list **CODE-FEATURE-4**), default OFF:

1. `SciantixField` (in `sciantix_binding.py`) holds one SCIANTIX point per `V_t`
   dof of the fissile region; the library + model settings are read once and shared.
2. `spine.initialize_fields` builds the field (when `models.fission_gas.enabled`)
   and a `gas_swelling` Function on `V_t`; `spine.update_state(dt)` calls
   `field.step(dt, T, fission_rate, burnup_old, burnup_new)` with `T` from the
   temperature field, `fission_rate = q''' / E_fission`, and the host burnup pair
   (Z3ST's RADAR model owns burnup; SCIANTIX consumes it — §3). The returned ΔV/V
   is written into `gas_swelling`.
3. `materials/sciantix_swelling.py::gaseous_swelling` is the eigenstrain callable —
   returns `(gas_swelling/3)·I` reading that field (same pattern as the burnup-fed
   swelling laws). A fuel card opts in with
   `eigenstrain: materials.sciantix_swelling.gaseous_swelling`.

No solver change was needed (the eigenstrain/state buses already existed). The
field also has `snapshot()`/`restore()` for adaptive-timestep rollback, hooked into
`spine.snapshot_state`/`restore_state`. Config:

```yaml
models:
  fission_gas:
    enabled: true
    lib: /path/to/libsciantix_tu.so      # else $SCIANTIX_LIB ; build with -DCOUPLING_TU
    initial_conditions: input_initial_conditions.txt
    energy_per_fission: 3.2e-11          # J/fission (≈ 200 MeV)
```
The run directory needs `input_settings.txt` + `input_initial_conditions.txt` (same
files a SCIANTIX standalone run uses). Current scope: fresh fuel, one point per dof.

## 5. Effective burnup for HBS — optional SCIANTIX patch

A `-DCOUPLING_TU` build skips `Burnup()`, `EffectiveBurnup()` and `Densification()`
(`Simulation.C:43`). Z3ST owns total burnup (fed in, §3) and densification (it is a
mechanical eigenstrain), so those two are correctly SCIANTIX-off. But **effective
burnup** is a SCIANTIX-internal HBS input (Khvostov, temperature-gated) with no Z3ST
equivalent — and the coupling does not transfer it, so in a stock `-DCOUPLING_TU`
build it stays frozen at 0 and HBS would be silently wrong.

Fix (SCIANTIX-side, in `effective_burnup_coupling.patch`): take `EffectiveBurnup()`
out of the guard and have it accumulate the temperature-gated **burnup increment**
(`Burnup.getIncrement()`) instead of `Specific power / 86400`. This needs no Specific
power (which the coupling build does not compute), keeps `dBu_eff ≤ dBu`, and is
numerically identical in a standalone build. Apply with:
```bash
cd <sciantix> && patch -p1 < <z3st>/z3st/coupling/sciantix/effective_burnup_coupling.patch
```
Verified: applies cleanly; builds standalone + `-DCOUPLING_TU`; the standalone HBS
regression `test_UO2HBS` reproduces effective burnup / restructured fraction / HBS
porosity exactly (rel err 0). Only needed if HBS (`iHighBurnupStructureFormation=1`)
is used; not required for the Baker validation. Note: `Irradiation time` and `FIMA`
remain uncomputed in a coupling build (also inside the skipped `Burnup()`); they have
no Z3ST consumer yet, but a model needing them would face the same gap.

## Status
Binding written against SCIANTIX 2.2.1; array map verified against source
2026-06-18 and **validated end-to-end the same day**.

Build the shared lib (no SCIANTIX repo edit needed):
```bash
cd <sciantix>
g++ -O2 -std=c++17 -fPIC -shared $(find include -type d | sed 's/^/-I/') \
    $(find src -name '*.C') -o /tmp/libsciantix.so
```
Then validate against SCIANTIX's own Baker regression gold:
```bash
cd <sciantix>/regression/baker/test_Baker1977__1273K
SCIANTIX_LIB=/tmp/libsciantix.so \
  PYTHONPATH=<z3st>/z3st/coupling/sciantix python3 \
  <z3st>/z3st/coupling/sciantix/validate_baker.py
```
Result (2026-06-18): all four engineering outputs match the standalone gold to
~1e-7 relative error — FGR 0.132097, intragranular swelling 3.07e-4, intergranular
swelling 0.0417, burnup 6.719. **[VALIDATION] PASS.**

Two non-obvious facts this surfaced (both handled in `load_initial_conditions`):
1. In coupling mode SCIANTIX does **not** read `input_initial_conditions.txt`
   (standalone-only, `file_manager/InputReading.C`) — the host seeds `variables[]`.
2. The standalone's one-time `Initialization()` (`file_manager/Initialization.C`),
   also skipped by the coupling entry, sets grain-boundary defaults absent from the
   IC file (`variables[25]`=2e13, `[35]`=0.5, `[37]`=1.0) and converts U% → at/m³.
   Without the grain-boundary defaults the intergranular model returns `nan` and
   releases nothing.

Coupling-build check (2026-06-19): with a `-DCOUPLING_TU` lib SCIANTIX skips its own
burnup and consumes Z3ST's (fed via `history[7]/[8]`); driving the gold burnup
trajectory reproduces the same gas outputs (~1e-7) and burnup matches exactly — so
the host-owns-burnup design is verified. `validate_baker.py` feeds the burnup
trajectory and passes against both a plain and a `-DCOUPLING_TU` build.

Z3ST integration (§4) is DONE and exercised by the `SciantixField` unit checks.
Remaining: fresh-fuel only (no diffusion-mode projection for a pre-irradiated
restart); a full pwr-rod case run with `models.fission_gas.enabled` + a
`-DCOUPLING_TU` lib is the next end-to-end step.
