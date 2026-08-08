#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: validation run for the develop merge
#
# Completes the verification left unfinished on the laptop:
#   1. environment check      - refuses to run on the wrong interpreter
#   2. full local suite       - 49 cases, regenerates non-regression_summary.txt
#   3. pwr_rod_2D             - gold regression, the 1800-day PWR rod case
#   4. fg_test_2D serial      - reference for the MPI comparison
#   5. fg_test_2D on N ranks  - compared column-by-column against the serial run
#   6. thermal_gradient_3D    - unbuffered, to settle whether it converges
#
# Phase 1 is mandatory: the laptop run was invalidated by a PATH where
# /usr/bin preceded the conda env, so gmsh resolved correctly while python3
# did not, and half the work ran on dolfinx 0.10.0 / numpy 1.26 instead of
# 0.11.0 / 2.4. This script fails fast rather than produce results that look
# fine and are not comparable.
#
# Usage:
#   ./z3st_validate.sh                  # all phases
#   ./z3st_validate.sh 3 5              # only phases 3 and 5
#   RANKS=8 ./z3st_validate.sh 5        # MPI phase on 8 ranks
#   nohup setsid ./z3st_validate.sh > run.log 2>&1 &    # survives a disconnect
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

set -uo pipefail

# --. Configuration: override from the environment if the paths differ --..
Z3ST_ROOT="${Z3ST_ROOT:-$HOME/z3st}"
CONDA_ENV="${CONDA_ENV:-z3st}"
RANKS="${RANKS:-4}"
PULL="${PULL:-1}"          # PULL=0 to skip the git pull
RESULTS="${RESULTS:-$Z3ST_ROOT/validation_$(date +%Y%m%d_%H%M%S)}"

CASES="$Z3ST_ROOT/z3st/cases"
PHASES=("$@")

log()  { printf '\n\033[1m=== %s ===\033[0m\n' "$*"; }
ok()   { printf '  \033[32mOK\033[0m    %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; }
info() { printf '  ..    %s\n' "$*"; }

want() {   # phase selection: no arguments means every phase
    [ ${#PHASES[@]} -eq 0 ] && return 0
    for p in "${PHASES[@]}"; do [ "$p" = "$1" ] && return 0; done
    return 1
}

mkdir -p "$RESULTS"
echo "Z3ST validation - $(date)"
echo "repo    : $Z3ST_ROOT"
echo "env     : $CONDA_ENV"
echo "ranks   : $RANKS"
echo "results : $RESULTS"

# --.. Phase 1: environment. Always runs, never skipped. ..--
log "1. Environment"

if [ ! -d "$Z3ST_ROOT/z3st" ]; then
    bad "no z3st package under $Z3ST_ROOT - set Z3ST_ROOT"; exit 1
fi

# conda activate, so PATH ordering is the shell's problem and not ours
CONDA_BASE="$(conda info --base 2>/dev/null)"
if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" || { bad "cannot activate env '$CONDA_ENV'"; exit 1; }
    ok "conda env '$CONDA_ENV' activated"
else
    info "conda not found, falling back to PATH prepend"
    export PATH="$HOME/miniconda3/envs/$CONDA_ENV/bin:$PATH"
fi

info "python3 -> $(command -v python3)"
info "mpirun  -> $(command -v mpirun)"
info "gmsh    -> $(command -v gmsh)"

# The interpreter must be the env's, hold numpy >= 2, and import z3st.
python3 - <<'PY' || exit 1
import sys
fail = []

try:
    import numpy
    major = int(numpy.__version__.split(".")[0])
    print(f"  ..    numpy   {numpy.__version__}")
    if major < 2:
        fail.append(f"numpy {numpy.__version__} < 2: np.trapezoid is missing and "
                    "most verification cases will fail")
except ImportError as e:
    fail.append(f"numpy not importable: {e}")

try:
    import dolfinx
    print(f"  ..    dolfinx {dolfinx.__version__}")
    if not dolfinx.__version__.startswith("0.11"):
        fail.append(f"dolfinx {dolfinx.__version__}: Z3ST targets 0.11.x, "
                    "gold references are not comparable across versions")
except ImportError as e:
    fail.append(f"dolfinx not importable: {e}")

try:
    import z3st
    print(f"  ..    z3st    {z3st.__file__}")
except ImportError as e:
    fail.append(f"z3st not importable: {e} - pip install -e . in the repo root")

if fail:
    print("\n  environment is not usable for this validation:")
    for f in fail:
        print(f"    - {f}")
    sys.exit(1)
PY
ok "interpreter usable"

# fg_test_2D drives SCIANTIX through a -DCOUPLING_TU shared library.
if [ -n "${SCIANTIX_LIB:-}" ]; then
    if [ -f "$SCIANTIX_LIB" ]; then ok "SCIANTIX_LIB $SCIANTIX_LIB"
    else bad "SCIANTIX_LIB set but missing: $SCIANTIX_LIB"; exit 1; fi
else
    info "SCIANTIX_LIB unset: phases 4 and 5 rely on the loader finding libsciantix.so"
fi

info "cores: $(nproc)"

if [ "$PULL" = "1" ]; then
    git -C "$Z3ST_ROOT" pull --ff-only && ok "pulled: $(git -C "$Z3ST_ROOT" log --oneline -1)"
fi

STATUS=0

# --.. Phase 2: the full local suite ..--
if want 2; then
    log "2. Local suite (49 cases)"
    ( cd "$CASES" && bash non-regression_local.sh ) > "$RESULTS/suite.log" 2>&1
    rc=$?
    cp "$CASES/non-regression_summary.txt" "$RESULTS/" 2>/dev/null
    n_ok=$(grep -c '\[OK\]'   "$CASES/non-regression_summary.txt" 2>/dev/null || echo 0)
    n_fail=$(grep -c '\[FAIL\]' "$CASES/non-regression_summary.txt" 2>/dev/null || echo 0)
    if [ "$rc" -eq 0 ]; then ok "suite passed: $n_ok OK"
    else
        bad "suite: $n_ok OK, $n_fail FAIL (exit $rc)"
        grep '^Case:.*FAIL' "$CASES/non-regression_summary.txt" 2>/dev/null | sed 's/^/        /'
        STATUS=1
    fi
    info "summary regenerated in the repo, ready to commit"
fi

# --.. Phase 3: pwr_rod_2D against its gold ..--
if want 3; then
    log "3. pwr_rod_2D (gold regression)"
    C="$CASES/regression/pwr_rod_2D"
    ( cd "$C" && ./Allclean >/dev/null 2>&1; ./Allrun ) > "$RESULTS/pwr_rod_2D.log" 2>&1
    cp "$C/output/non-regression.json" "$RESULTS/pwr_rod_2D_nonreg.json" 2>/dev/null
    if grep -q 'PASS Regression within tolerance' "$RESULTS/pwr_rod_2D.log"; then
        ok "gold regression PASS"
    else
        bad "gold regression did not pass - see $RESULTS/pwr_rod_2D.log"
        grep -A8 'Regression check vs GOLD' "$RESULTS/pwr_rod_2D.log" | sed 's/^/        /'
        STATUS=1
    fi
fi

# --.. Phase 4: fg_test_2D serial, the MPI reference ..--
if want 4; then
    log "4. fg_test_2D serial (MPI reference)"
    C="$CASES/regression/fg_test_2D"
    ( cd "$C" && ./Allclean >/dev/null 2>&1; ./Allrun ) > "$RESULTS/fg_serial.log" 2>&1
    if [ -s "$C/output/history.csv" ]; then
        cp "$C/output/history.csv" "$RESULTS/fg_serial_history.csv"
        ok "serial reference: $(( $(wc -l < "$RESULTS/fg_serial_history.csv") - 1 )) steps"
    else
        bad "no history.csv - see $RESULTS/fg_serial.log"; STATUS=1
    fi
fi

# --.. Phase 5: fg_test_2D in parallel, compared against the serial run.
# This is what validates the collective reductions added to the case's
# diagnostics.py: burnup and the fission-gas averages are reduced across
# ranks over owned dofs, so they must match the serial values. ..--
if want 5; then
    log "5. fg_test_2D on $RANKS ranks (vs serial)"
    C="$CASES/regression/fg_test_2D"
    REF="$RESULTS/fg_serial_history.csv"
    if [ ! -s "$REF" ]; then
        bad "no serial reference: run phase 4 first"; STATUS=1
    else
        rm -f "$C/output/history.csv"
        ( cd "$C" && Z3ST_PLAIN_LOG=1 mpirun -n "$RANKS" python3 -m z3st ) \
            > "$RESULTS/fg_mpi${RANKS}.log" 2>&1
        if [ -s "$C/output/history.csv" ]; then
            cp "$C/output/history.csv" "$RESULTS/fg_mpi${RANKS}_history.csv"
            n_hdr=$(grep -c '^step,' "$C/output/history.csv")
            if [ "$n_hdr" -eq 1 ]; then ok "single header row (the rank-0 write guard holds)"
            else bad "$n_hdr header rows: ranks are interleaving writes"; STATUS=1; fi

            python3 - "$REF" "$C/output/history.csv" <<'PY' || STATUS=1
import csv, sys

s = list(csv.DictReader(open(sys.argv[1])))
p = list(csv.DictReader(open(sys.argv[2])))

if len(s) != len(p):
    print(f"  FAIL  step count differs: serial {len(s)}, parallel {len(p)}")
    sys.exit(1)
if s[0].keys() != p[0].keys():
    print("  FAIL  column names differ")
    sys.exit(1)

# Quantities diagnostics.py reduces itself: any ghost double-counting or
# rank-local value shows up here as an O(1) error, not a rounding difference.
REDUCED = ("burnup_avg_MWdkgU", "burnup_max_MWdkgU", "fg_produced_avg",
           "fg_ingrain_avg", "fg_gb_avg", "fg_released_avg", "fgr_frac")
# Read verbatim from solver-side global scalars: differences here come from
# the parallel iteration path, not from the diagnostics module.
SOLVER = ("gap_um", "contact_pressure_MPa", "T_max_K", "T_min_K")

worst = {}
for rs, rp in zip(s, p):
    if rs["step"] != rp["step"]:
        print(f"  FAIL  steps misaligned at {rs['step']}")
        sys.exit(1)
    for c in s[0].keys():
        a, b = float(rs[c]), float(rp[c])
        rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
        if rel > worst.get(c, (-1,))[0]:
            worst[c] = (rel, rs["step"], a, b)

print(f"\n  {'column':<24}{'max rel diff':>14}  step   serial            parallel")
print("  " + "-" * 86)
for c in s[0].keys():
    rel, st, a, b = worst[c]
    kind = "reduced" if c in REDUCED else ("solver" if c in SOLVER else "")
    print(f"  {c:<24}{rel:14.3e}  {st:>4}   {a:<17.8g} {b:<17.8g} {kind}")
print("  " + "-" * 86)

# The reduced quantities are the acceptance criterion. Burnup follows the power
# history and is independent of the mechanical solution, so it must be exact;
# the fission-gas averages inherit the solution's parallel drift, hence a
# looser bound.
TOL = {"burnup_avg_MWdkgU": 1e-12, "burnup_max_MWdkgU": 1e-12}
bad = []
for c in REDUCED:
    if c not in worst:
        continue
    tol = TOL.get(c, 1e-4)
    if worst[c][0] > tol:
        bad.append(f"{c}: {worst[c][0]:.3e} > {tol:.0e} (step {worst[c][1]})")

if bad:
    print("\n  FAIL  reduced quantities disagree between serial and parallel:")
    for b in bad:
        print(f"          {b}")
    print("        an O(1) error means ghost double-counting or a missing")
    print("        allreduce; a small one means the solution itself drifted.")
    sys.exit(1)

print("\n  OK    reduced quantities agree: collective reductions are correct")
mx = max(worst[c][0] for c in SOLVER if c in worst)
print(f"  ..    solver-side scalars differ by up to {mx:.3e} (parallel iteration path)")
PY
        else
            bad "no history.csv from the parallel run"; STATUS=1
        fi
    fi
fi

# --.. Phase 6: thermal_gradient_3D, unbuffered.
# On the laptop this case sat 33 minutes on its first solve at 923% CPU with a
# buffered log, so slow and non-converging were indistinguishable. Unbuffered
# output makes the iterations visible. ..--
if want 6; then
    log "6. thermal_gradient_3D (unbuffered)"
    C="$CASES/verification/mechanics/thermal_gradient_3D"
    ( cd "$C" && ./Allclean >/dev/null 2>&1
      gmsh mesh.geo -3 > log_mesh.md 2>&1 || gmsh mesh.geo -2 > log_mesh.md 2>&1
      Z3ST_PLAIN_LOG=1 PYTHONUNBUFFERED=1 timeout 3600 python3 -m z3st ) \
        > "$RESULTS/thermal_gradient_3D.log" 2>&1
    rc=$?
    iters=$(grep -c 'Iteration' "$RESULTS/thermal_gradient_3D.log" 2>/dev/null || echo 0)
    if [ "$rc" -eq 124 ]; then
        bad "timed out after 1 h at $iters iterations: not converging"
        STATUS=1
    elif [ "$rc" -eq 0 ]; then
        ok "completed, $iters iterations: it was slow, not stuck"
        grep -E 'Simulation completed|Total time steps' "$RESULTS/thermal_gradient_3D.log" | sed 's/^/        /'
    else
        bad "exit $rc - see $RESULTS/thermal_gradient_3D.log"; STATUS=1
    fi
fi

log "Summary"
if [ "$STATUS" -eq 0 ]; then ok "every phase run passed"
else bad "at least one phase failed - see the logs"; fi
echo "  logs and artefacts: $RESULTS"
exit "$STATUS"
