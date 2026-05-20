#!/usr/bin/env python3
"""Force-displacement curve for SENS.

Two-stage workflow:
  1. If `output/fields_*.vtu` files are present, extract `(step, u_x, F_x)`
     per step and write `output/force_displacement.csv` (overwriting any
     previous version).
  2. Always read `output/force_displacement.csv` and plot `F_x` vs `u_x`
     into `output/force_displacement.png`.

If no VTU files exist (e.g. after `./Allclean`), stage 1 is skipped and
the script simply re-renders the plot from the saved CSV.
"""

import os
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

CASE_DIR = os.path.dirname(os.path.abspath(__file__))

if len(sys.argv) > 1:
    arg = sys.argv[1]
    OUTPUT_DIR = arg if os.path.isabs(arg) else os.path.join(CASE_DIR, arg)
else:
    OUTPUT_DIR = os.path.join(CASE_DIR, "output")

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"[INFO] Using OUTPUT_DIR = {OUTPUT_DIR}")
VTU_FILES  = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
CSV_FILE   = os.path.join(OUTPUT_DIR, "force_displacement.csv")
PNG_FILE   = os.path.join(OUTPUT_DIR, "force_displacement.png")

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])


# ----- Stage 1: extract per-step (u, F) from VTUs and write CSV -------------
if VTU_FILES:
    print(f"[INFO] Extracting (u, F) per step from {len(VTU_FILES)} VTU files...")
    top_tol = Ly / 200.0
    rows = []
    for step, vtufile in enumerate(VTU_FILES):
        try:
            m = pv.read(vtufile)
        except Exception as e:
            print(f"[WARN] skipping {os.path.basename(vtufile)}: {e}")
            continue
        pts = m.points
        x_pts, y_pts = pts[:, 0], pts[:, 1]
        top_mask = np.abs(y_pts - Ly) < top_tol
        if not np.any(top_mask):
            continue

        # u_x on top edge.
        if "Displacement" in m.point_data:
            u_x_top = float(np.max(np.asarray(m.point_data["Displacement"])[top_mask, 0]))
        else:
            u_x_top = float("nan")

        # F_x = integral of sigma_xy along top edge (per unit out-of-plane depth).
        if "Stress_steel (points)" in m.point_data:
            S_top = np.asarray(m.point_data["Stress_steel (points)"])[top_mask, 1]
            x_top = x_pts[top_mask]
            order = np.argsort(x_top)
            F_per_depth = float(np.trapezoid(S_top[order], x_top[order]))
            F_kN = F_per_depth * 1e-3 * 1e-3   # × Lz_ambati × kN_conversion
        else:
            F_kN = float("nan")

        rows.append((step, u_x_top * 1e3, F_kN))

    arr = np.array(rows)
    np.savetxt(CSV_FILE, arr,
               header="Step,u_x_top_mm,F_x_kN",
               delimiter=",", fmt=["%d", "%.6e", "%.6e"], comments="")
    print(f"[INFO] CSV saved: {CSV_FILE}  ({len(rows)} rows)")
else:
    print(f"[INFO] No VTU files in {OUTPUT_DIR}; skipping extraction, using existing CSV.")


# ----- Stage 2: plot from CSV -----------------------------------------------
if not os.path.exists(CSV_FILE):
    raise FileNotFoundError(
        f"No VTU files and no existing CSV ({CSV_FILE}); nothing to plot."
    )

data = np.genfromtxt(CSV_FILE, delimiter=",", skip_header=1)
if data.ndim == 1:
    data = data.reshape(1, -1)
u_mm = data[:, 1]
F_kN = data[:, 2]
finite = np.isfinite(u_mm) & np.isfinite(F_kN)
u_mm, F_kN = u_mm[finite], F_kN[finite]

plt.figure(figsize=(7, 5))
plt.plot(u_mm, F_kN, "C0-o", markersize=3, label="Numerical")
if F_kN.size > 0:
    peak = int(np.nanargmax(F_kN))
    plt.plot(u_mm[peak], F_kN[peak], "r*", markersize=12,
             label=f"Peak: F = {F_kN[peak]:.3f} kN at u = {u_mm[peak]:.4f} mm")
plt.xlabel(r"Top-edge displacement $u_x$ (mm)")
plt.ylabel(r"Force $F_x$ (kN, per 1 mm depth)")
plt.title("Force-displacement curve (cf. Ambati Fig. 13)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(PNG_FILE, dpi=200)
plt.close()
print(f"[INFO] force_displacement.png saved -> {PNG_FILE}")
