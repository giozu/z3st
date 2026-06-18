# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Post-processing: render the damage field at the final time step of the
# plate thermal-shock run to a PNG (the nucleated edge-crack array).
# Reads the XDMF/.h5 written by the solver; no dolfinx import needed, so it
# runs even with the file lock held (HDF5_USE_FILE_LOCKING=FALSE).

import glob
import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# Plate dimensions (must match geometry.yaml).
L, H = 0.025, 0.0098

h5_files = sorted(glob.glob("output/*.h5"))
if not h5_files:
    raise SystemExit("[plot_damage] no output/*.h5 found — run the solver first.")

with h5py.File(h5_files[0], "r") as h:
    X = h["Mesh/mesh/geometry"][:]      # (N, 3) node coordinates
    topo = h["Mesh/mesh/topology"][:]   # (M, 3) triangle connectivity
    steps = sorted(h["Function/Damage"].keys())
    D = h["Function/Damage/" + steps[-1]][:].ravel()

# Localisation diagnostics (printed to the run log).
d_edge = np.minimum.reduce([X[:, 1], H - X[:, 1], X[:, 0]])
print(f"[plot_damage] final step '{steps[-1]}': "
      f"max D = {D.max():.3f}, crack-core fraction (D>0.9) = {(D > 0.9).mean():.3f}, "
      f"interior max D = {D[d_edge > 1.5e-3].max():.2e}")

triang = mtri.Triangulation(X[:, 0] * 1e3, X[:, 1] * 1e3, topo)

fig, ax = plt.subplots(figsize=(13, 5.6))
tc = ax.tripcolor(triang, D, shading="gouraud", cmap="inferno", vmin=0.0, vmax=1.0)
ax.set_aspect("equal")
ax.set_xlim(0, L * 1e3)
ax.set_ylim(0, H * 1e3)
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_title("Plate thermal shock — damage D (final step)")
fig.colorbar(tc, ax=ax, fraction=0.025, pad=0.01, label="D")
fig.tight_layout()

out = "output/damage_field.png"
fig.savefig(out, dpi=130)
print(f"[plot_damage] saved {out}")
