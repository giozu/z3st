#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Mesh refinement series of the spherical_shell case.

Each level doubles the cells along every patch edge (Na) and across the
thickness (Nr), and takes the square root of the radial growth ratio, so the
radial grading is the same at every level. For each level the case is copied
next to itself, meshed, run and compared with the exact solution by
non-regression.py; the errors are collected in convergence_data.txt, which
plots.py turns into output/mesh_convergence.png.

Not part of Allrun: the finest level needs about 150 s and 2.5 GB. A fourth
level (32 x 32 x 80, 245760 cells) needs about 14 GB and is left out.

Run:  python3 mesh_convergence.py
"""

import json
import os
import re
import shutil
import subprocess
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "convergence_data.txt")
METRICS = ["T_max", "L2_error_T", "Linf_error_T",
           "L2_error_sigma_rr", "L2_error_sigma_tt", "L2_error_u_r"]

# (Na, Nr, q): Na cells along each patch edge, Nr radial cells, growth ratio q
LEVELS = [(4, 10, 1.3), (8, 20, 1.3 ** 0.5), (16, 40, 1.3 ** 0.25)]


def run_level(na, nr, q):
    work = os.path.join(os.path.dirname(CASE_DIR), f"_spherical_shell_conv_{na}")
    shutil.rmtree(work, ignore_errors=True)
    shutil.copytree(CASE_DIR, work, ignore=shutil.ignore_patterns("output", "__pycache__"))
    os.makedirs(os.path.join(work, "output"))
    try:
        geo = os.path.join(work, "mesh.geo")
        txt = open(geo).read()
        for name, val in (("Na", na), ("Nr", nr), ("q ", q)):
            txt = re.sub(rf"^({name.strip()}\s*= DefineNumber\[ )[0-9.]+", rf"\g<1>{val}", txt, flags=re.M)
        open(geo, "w").write(txt)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        subprocess.run(["gmsh", "mesh.geo", "-3"], cwd=work, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with open(os.path.join(work, "log_z3st.md"), "w") as log:
            subprocess.run([sys.executable, "-m", "z3st"], cwd=work, env=env, check=True,
                           stdout=log, stderr=subprocess.STDOUT)
        subprocess.run([sys.executable, "non-regression.py"], cwd=work, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cells = int(re.search(r"Num cells:\s*(\d+)", open(os.path.join(work, "log_z3st.md")).read()).group(1))
        res = json.load(open(os.path.join(work, "output", "non-regression.json")))["results"]
        return cells, [res[m]["rel_error"] for m in METRICS]
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    rows = []
    for na, nr, q in LEVELS:
        cells, errs = run_level(na, nr, q)
        rows.append([na, nr, q, cells] + errs)
        print(f"Na={na:3d} Nr={nr:3d} cells={cells:7d}  " + "  ".join(f"{e:.3e}" for e in errs))
    header = "Na Nr q cells " + " ".join(METRICS)
    with open(OUT, "w") as f:
        f.write("# spherical_shell mesh refinement series, written by mesh_convergence.py\n")
        f.write("# errors relative to the exact solution, as defined in non-regression.py\n")
        f.write("# " + header + "\n")
        for r in rows:
            f.write(f"{r[0]} {r[1]} {r[2]:.6f} {r[3]} " + " ".join(f"{e:.6e}" for e in r[4:]) + "\n")
    print(f"[INFO] wrote {OUT}")


if __name__ == "__main__":
    main()
