cdcdcdc# Z3ST ↔ SCIANTIX coupling (prototype)

Drive SCIANTIX (mesoscale fission-gas behaviour) from Z3ST to compute **gaseous
swelling** (→ the eigenstrain bus) and **fission gas release** per fuel point.
The binding wraps SCIANTIX's existing C-linkage coupling entry — **no SCIANTIX
physics is reimplemented in Python**; the C++ stays the single source of truth.

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

## 2. Validate the binding (standalone)

```bash
python3 smoke_test.py     # ramps T at fixed fission rate, prints swelling/bu/FGR
```
Diff the output against a SCIANTIX **standalone** run with the same
`input_history.txt` (same T, fission rate, dt) — they must match. That closes the
binding correctness gap before any Z3ST integration.

## 3. Array layout (read from SCIANTIX v2.2.1)

`include/MainVariables.h`: `options[40]`, `history[20]`, `variables[300]`.

| host writes — `history[]` | idx | reads — `variables[]` | idx |
|---|---|---|---|
| Temperature old/new (K) | 0,1 | Xe produced (at/m³) | 1 |
| Fission rate old/new (fiss/m³s) | 2,3 | Xe released (at/m³) | 6 |
| Hydrostatic stress old/new (MPa) | 4,5 | intragranular gas swelling (/) | 24 |
| time step Δt | 6 *(confirm)* | **intergranular gas swelling (/)** | **36** |
| steam pressure old/new (atm) | 9,10 | Burnup (MWd/kgUO₂) | 38 |

Source: `src/operations/SetVariablesFunctions.C`. Indices marked *(confirm)*
should be checked by the SCIANTIX author against the current source.

## 4. Z3ST integration (next step — the eigenstrain bus)

SCIANTIX gaseous swelling is a **numerical, stateful per-point field**, not a UFL
expression — so it rides the **state bus**, exactly like burnup/creep:

1. one `SciantixSolver` per fuel integration region (radial ring first), kept on
   the material;
2. in `spine.update_state(dt)`, for each fuel region call `solver.advance(dt, T,
   T, fission_rate, ...)` with `T` from the temperature field and
   `fission_rate ∝ q'''`; write the returned `gaseous_swelling` into a DG0
   field `model.gas_swelling`;
3. a `materials/sciantix_swelling.py` eigenstrain callable returns
   `(gas_swelling/3)·I` reading that field — same pattern as
   `fuel_swelling.solid_swelling` reading `model.burnup`. FGR feeds a future
   plenum-pressure model.

No solver change is needed (the eigenstrain/state buses already exist); the work
is the per-region state bookkeeping. This is punch-list **CODE-FEATURE-4**.

## Status
Draft. Binding written against SCIANTIX 2.2.1; **not yet run end-to-end** (needs
the shared lib). The single ctypes call and the array map are the only things to
verify; once the smoke test matches a standalone run, the Z3ST integration is
mechanical.
