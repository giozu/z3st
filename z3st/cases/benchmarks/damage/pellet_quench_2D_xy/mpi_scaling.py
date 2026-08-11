#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: MPI sweep for the thermal-shock fracture benchmark
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Re-runs this case at several MPI rank counts and checks that the result does not
depend on the decomposition.

Why this exists
---------------
The non-regression suite runs serially throughout, so a result that depends on
the domain decomposition would pass every case and reach a user unnoticed. This
script is the check that would catch it: each rank count must reproduce the
serial verdict exactly.

On timings
----------
Wall times are recorded in the JSON because they are useful while developing,
but they are deliberately not plotted and are not reported in the paper. On a
shared desktop the wall time cannot be separated from whatever else the machine
is doing -- the same serial run of this case has been measured at 220 s and at
288 s -- and with a threaded BLAS left unpinned an earlier sweep measured thread
contention rather than MPI. A speed-up obtained that way says more about the
machine than about the code. Use a quiet machine with enough cores, and a mesh
large enough to be worth decomposing, before quoting any of these numbers.

Ordering note: the sweep runs the highest rank count first and finishes at one
rank, so ``output/`` is left in the serial state the gold expects.

Usage
-----
    python3 mpi_scaling.py                 # 1, 2, 4 ranks
    python3 mpi_scaling.py 1 2 4 8         # explicit list; 8 oversubscribes
                                           # a machine with fewer cores
"""

import json
import os
import re
import subprocess
import sys


CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
RECORD = os.path.join(OUT, "mpi_scaling.json")

TIME_RE = re.compile(r"completed in ([0-9.]+)\s*s")
CELLS_RE = re.compile(r"Num cells:\s*([0-9]+)")


def physical_cores():
    """Physical cores, so the plot can mark oversubscribed points honestly."""
    try:
        out = subprocess.run(["lscpu", "-p=Core,Socket"], capture_output=True,
                             text=True, check=True).stdout
        pairs = {ln for ln in out.splitlines() if not ln.startswith("#")}
        return len(pairs)
    except Exception:
        return os.cpu_count() or 1


def ensure_mesh():
    """Generate mesh.msh if absent or older than mesh.geo.

    Allclean removes *.msh, and unlike Allrun this script does not rebuild it,
    so a sweep launched after a clean used to die on the first rank with a bare
    FileNotFoundError from every rank at once.
    """
    geo = os.path.join(CASE_DIR, "mesh.geo")
    msh = os.path.join(CASE_DIR, "mesh.msh")
    if os.path.exists(msh) and os.path.getmtime(msh) >= os.path.getmtime(geo):
        return
    print("[mpi_scaling] mesh.msh missing or stale, running gmsh", flush=True)
    with open(os.path.join(CASE_DIR, "log_mesh.md"), "w") as fh:
        r = subprocess.run(["gmsh", "-2", "mesh.geo", "-format", "msh2"],
                           cwd=CASE_DIR, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0 or not os.path.exists(msh):
        raise SystemExit("[mpi_scaling] gmsh failed; see log_mesh.md")


def run(nranks):
    """One run at nranks. Returns (wall_seconds, n_cells, results_dict)."""
    log = os.path.join(CASE_DIR, f"log_z3st_np{nranks}.md")
    cmd = ["python3", "-m", "z3st"] if nranks == 1 else \
          ["mpirun", "-n", str(nranks), "--oversubscribe", "python3", "-m", "z3st"]

    # One thread per rank, always, including the serial baseline. OpenBLAS
    # otherwise opens a thread pool inside every rank -- measured at ~25 threads
    # and 250 % CPU per rank on this case -- so four ranks ask for ten cores on a
    # six-core machine and the sweep measures thread thrashing instead of MPI.
    # Pinning also makes the baseline honest: without it the serial run gets a
    # threaded BLAS and the parallel runs are compared against a faster reference
    # than the one the speed-up implies.
    env = dict(os.environ,
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

    print(f"[mpi_scaling] {nranks} rank(s), 1 thread each: {' '.join(cmd)}",
          flush=True)
    with open(log, "w") as fh:
        r = subprocess.run(cmd, cwd=CASE_DIR, stdout=fh,
                           stderr=subprocess.STDOUT, env=env)
    text = open(log).read()
    if r.returncode != 0:
        tail = "\n".join(text.splitlines()[-15:])
        raise SystemExit(f"[mpi_scaling] run at {nranks} rank(s) failed:\n{tail}")
    t = TIME_RE.search(text)
    c = CELLS_RE.search(text)
    if not t:
        raise SystemExit(f"[mpi_scaling] no wall time in {log}")

    # The verdict is what establishes rank independence, so it is part of the
    # measurement rather than a separate step. It must therefore be fatal: an
    # earlier version warned and carried on, which left the stale verdict from
    # the previous rank count in place and made the drift comparison read zero
    # by comparing a file with itself. A check that cannot fail is not a check.
    jpath = os.path.join(OUT, "non-regression.json")
    if os.path.exists(jpath):
        os.remove(jpath)
    nr = subprocess.run(["python3", "non-regression.py"], cwd=CASE_DIR,
                        capture_output=True, text=True, env=env)
    if nr.returncode != 0 or not os.path.exists(jpath):
        tail = "\n".join((nr.stdout + nr.stderr).splitlines()[-15:])
        raise SystemExit(
            f"[mpi_scaling] non-regression failed after the {nranks}-rank run "
            f"(exit {nr.returncode}); rank independence cannot be established:"
            f"\n{tail}")
    res = {k: v["numerical"]
           for k, v in json.load(open(jpath))["results"].items()}

    # Keep the per-rank verdict, so the comparison is against archived files
    # rather than against whatever happens to be in output/ at the end.
    os.makedirs(os.path.join(OUT, "mpi"), exist_ok=True)
    with open(os.path.join(OUT, "mpi", f"non-regression_np{nranks}.json"), "w") as fh:
        json.dump(res, fh, indent=4)

    return float(t.group(1)), int(c.group(1)) if c else None, res


def main():
    ranks = [int(a) for a in sys.argv[1:]] or [1, 2, 4]
    ranks = sorted(set(ranks))
    cores = physical_cores()
    ensure_mesh()

    # Descending, so the final run is serial and output/ is left canonical.
    data = {}
    cells = None
    for n in sorted(ranks, reverse=True):
        t, c, res = run(n)
        cells = cells or c
        data[n] = {"wall_s": t, "results": res}
        print(f"[mpi_scaling]   {n} rank(s): {t:.1f} s", flush=True)

    # Rank independence: every rank count must reproduce the serial result.
    ref = data[min(ranks)]["results"]
    drift = {}
    for n in ranks:
        for k, v in data[n]["results"].items():
            if k in ref and ref[k]:
                d = abs(v - ref[k]) / abs(ref[k])
                drift[k] = max(drift.get(k, 0.0), d)
    print("\n[mpi_scaling] rank independence (max relative drift vs serial):")
    for k, d in sorted(drift.items()):
        print(f"  {k:28s} {d:.2e}{'   <-- NOT rank independent' if d > 1e-9 else ''}")
    rank_independent = all(d <= 1e-9 for d in drift.values())

    t1 = data[min(ranks)]["wall_s"]
    record = {"cells": cells, "physical_cores": cores,
              "threads_per_rank": 1,
              "rank_independent": rank_independent,
              "runs": {str(n): data[n]["wall_s"] for n in ranks},
              "speedup": {str(n): t1 / data[n]["wall_s"] for n in ranks}}
    with open(RECORD, "w") as fh:
        json.dump(record, fh, indent=4)
    print(f"\n[mpi_scaling] record written: {RECORD}")

    print("\n[mpi_scaling] summary")
    print(f"  cells             {cells}")
    print(f"  physical cores    {cores}, one thread per rank")
    for n in sorted(ranks):
        print(f"  {n} rank(s)        {data[n]['wall_s']:.1f} s"
              f"   (indicative only, see the module docstring)")
    print(f"  rank independent  {rank_independent}")
    if not rank_independent:
        raise SystemExit("[mpi_scaling] result depends on the rank count")


if __name__ == "__main__":
    main()
