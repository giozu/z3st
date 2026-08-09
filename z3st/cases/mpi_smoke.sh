#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: MPI smoke test — 10 s, run by non-regression_github.sh after the cases.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Guards the one property no case verifies, because it is a property of the
# runtime and not of the physics: z3st stays rank-symmetric and gates its
# output. Both failure modes below were real, and both were invisible in serial.
#
#   1. Exit deadlock. dolfinx's PETSc wrappers do collective work in __del__.
#      A rank holding a different set of live Python objects has the collector
#      order those destructors differently, and the run hangs at exit with every
#      result already written. Three designs hung this way before the current
#      one; a `timeout` is therefore part of the check, not a convenience.
#
#   2. Output multiplied by the rank count. Every rank runs the same code, so an
#      ungated print appears once per process. The gate that fixes this is a
#      per-write flag inside the stdout wrapper -- and an earlier version was
#      installed only when markdown was active, which made it silently INERT
#      under mpirun, because there stdout is a tty even when the shell redirects
#      to a file.
#
# Not a new case: a case directory would duplicate a mesh, an input and a gold
# to test something none of them describes. This reuses the cheapest existing
# case and asserts on its log.

set -u
CASE="verification/cluster/mass_conservation_1D"
NP=${NP:-2}
LOG=$(mktemp)
trap 'rm -f "$LOG"' EXIT

cd "${0%/*}/$CASE" || { echo "MPI smoke: case not found: $CASE"; exit 1; }

echo "--.. MPI smoke test: $CASE on $NP ranks"

# 300 s is ~30x the serial cost: generous for a slow runner, still bounded, so a
# deadlock fails the job instead of hanging it until the CI timeout.
timeout 300 mpirun -n "$NP" python3 -m z3st > "$LOG" 2>&1
status=$?

if [[ $status -eq 124 ]]; then
    echo "MPI smoke FAILED: timed out on $NP ranks — the classic rank-asymmetry"
    echo "  deadlock in __del__. Results may be complete; the exit is what hangs."
    tail -n 20 "$LOG"
    exit 1
fi
if [[ $status -ne 0 ]]; then
    echo "MPI smoke FAILED: exit $status on $NP ranks"
    tail -n 30 "$LOG"
    exit 1
fi

# Every step banner must appear exactly once. On N ranks without the gate it
# appears N times, which is the whole failure mode -- and counting duplicates is
# what makes this check able to fail: a plain "did it run" test passes just as
# happily with the output multiplied.
dupes=$(grep -oE '^\[STEP [0-9]+/[0-9]+\]|^## Step [0-9]+' "$LOG" | sort | uniq -d)
if [[ -n "$dupes" ]]; then
    echo "MPI smoke FAILED: step banners duplicated — output is not rank-gated:"
    echo "$dupes" | head -5
    exit 1
fi

steps=$(grep -cE '^\[STEP [0-9]+/[0-9]+\]|^## Step [0-9]+' "$LOG")
if [[ $steps -eq 0 ]]; then
    echo "MPI smoke FAILED: no step banner in the log at all. Either the run"
    echo "  produced nothing, or rank 0 was gated off along with the others --"
    echo "  a gate that silences everyone passes a duplicate check trivially."
    tail -n 20 "$LOG"
    exit 1
fi

echo "MPI smoke PASSED: $NP ranks, clean exit, $steps step banners, none duplicated"
