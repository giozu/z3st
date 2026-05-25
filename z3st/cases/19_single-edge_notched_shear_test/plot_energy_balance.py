#!/usr/bin/env python3
"""Standalone energy-balance plot for the SENS shear test.

Reads ``energies.txt`` and writes ``energy_balance.png``: three traces
(E_el, E_frac, E_tot vs Step) with two horizontal references — the
notch baseline ``Gc * Dn`` (step 0) and the Ambati Fig. 12d arrest target
``Gc * (Dn + 0.55 mm)``.

This is a post-hoc copy of the energy-balance block inside
``non-regression.py`` (which runs as part of ``Allrun`` and writes the
plot into ``./output/`` against the live run). Use this script instead
when you want to regenerate the plot for a backed-up output directory
(e.g. ``output_starconvex_g00/``, ``output_backup/``) without disturbing
the live state.

Usage:
    python3 plot_energy_balance.py                       # default: live run
                                                          # reads CASE_DIR/energies.txt
                                                          # writes CASE_DIR/output/energy_balance.png

    python3 plot_energy_balance.py output_starconvex_g1  # post-hoc on a backup
                                                          # reads output_starconvex_g1/energies.txt
                                                          #   (or CASE_DIR/energies.txt if absent)
                                                          # writes output_starconvex_g1/energy_balance.png

Relative paths are resolved against the script's own directory.
"""

import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml


CASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Optional positional argument: the output directory (default = CASE_DIR/output).
# Absolute paths are honoured; relative paths are resolved against CASE_DIR.
if len(sys.argv) > 1:
    arg = sys.argv[1]
    OUTPUT_DIR = arg if os.path.isabs(arg) else os.path.join(CASE_DIR, arg)
else:
    OUTPUT_DIR = os.path.join(CASE_DIR, "output")

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"[INFO] Using OUTPUT_DIR = {OUTPUT_DIR}")

# Resolve which energies.txt to read. Preference order:
#   1. <OUTPUT_DIR>/energies.txt   (snapshotted alongside a backup)
#   2. <CASE_DIR>/energies.txt     (the live-run file at the case root)
_candidates = [
    os.path.join(OUTPUT_DIR, "energies.txt"),
    os.path.join(CASE_DIR,   "energies.txt"),
]
energy_file = next((p for p in _candidates if os.path.exists(p)), None)
if energy_file is None:
    raise FileNotFoundError(
        "No energies.txt found in either "
        f"{OUTPUT_DIR} or {CASE_DIR}."
    )
print(f"[INFO] Using energy file = {energy_file}")


# ----- Read material + geometry constants for the reference lines -----------
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    cfg = yaml.safe_load(f)
mat_rel = list(cfg["materials"].values())[0]
MATERIAL_FILE = os.path.normpath(os.path.join(CASE_DIR, mat_rel))

with open(MATERIAL_FILE) as f:
    mat = yaml.safe_load(f)
Gc = float(mat["Gc"])

with open(os.path.join(CASE_DIR, "mesh.geo")) as f:
    Dn = float(re.search(r"Dn\s*=\s*([\d\.eE+-]+)", f.read()).group(1))


# ----- Read energies and plot -----------------------------------------------
data = np.genfromtxt(energy_file, names=True, skip_header=0)
if data.ndim == 0:
    # genfromtxt returns 0-d when there's a single row; promote to 1-d.
    data = data.reshape(1)

# Ambati Fig. 12d arrest target: notch length + ~0.55 mm of additional propagation.
ambati_arc    = 0.55e-3
E_frac_notch  = Gc * Dn
E_frac_target = Gc * (Dn + ambati_arc)

plt.figure(figsize=(8, 5))
plt.plot(data["Step"], data["E_el"],   "b-o", markersize=3, label=r"Elastic $E_{el}$")
plt.plot(data["Step"], data["E_frac"], "r-s", markersize=3, label=r"Fracture $E_{frac}$")
plt.plot(data["Step"], data["E_tot"],  "k--", lw=1.5,        label=r"Total $E_{tot}$")
plt.axhline(E_frac_notch, color="gray", ls="--", alpha=0.5,
            label=rf"$G_c \cdot D_n = {E_frac_notch:.2f}$ J (notch baseline, step 0)")
plt.axhline(E_frac_target, color="gray", ls=":", alpha=0.7,
            label=rf"$G_c \cdot (D_n + 0.55\,\mathrm{{mm}}) = {E_frac_target:.2f}$ J (Ambati Fig. 12d total)")
plt.xlabel("Step")
plt.ylabel("Energy (J)")
plt.title("Global energy balance")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="upper left", fontsize="small")
plt.tight_layout()

png_path = os.path.join(OUTPUT_DIR, "energy_balance.png")
plt.savefig(png_path, dpi=200)
plt.close()
print(f"[INFO] energy_balance.png saved -> {png_path}")
print(f"[INFO]   final  E_el = {float(data['E_el'][-1]):.3f} J, "
      f"E_frac = {float(data['E_frac'][-1]):.3f} J  "
      f"(step {int(data['Step'][-1])})")
