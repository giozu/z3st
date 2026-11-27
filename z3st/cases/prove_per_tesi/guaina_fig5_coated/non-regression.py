#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
guaina_fig5_coated: stresses at mid-height in Zircaloy and Cr.

  output/stress_midheight.png          radial profiles at the last step
  output/stress_midheight_history.png  mean over the wall vs time
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np

Lz = 77.0e-3     # m  height
R_int = 4.75e-3  # m  Zircaloy / Cr interface

# Read the stress of every step
with h5py.File("output/simulation.h5", "r") as f:
    nodes = f["Mesh/mesh/geometry"][:]  # node coordinates, (n_nodes, 3)
    cells = f["Mesh/mesh/topology"][:]  # node indices of each cell, (n_cells, 4)
    # One dataset per saved step, named by its time with "_" for the decimal
    # point (e.g. "87_5"); sort them numerically, not alphabetically
    steps = sorted(f["Function/Stress"], key=lambda s: float(s.replace("_", ".")))
    times = [float(s.replace("_", ".")) for s in steps]
    # S[step, cell, component]: one value per cell; the 9 components are the
    # 3x3 tensor in (r, θ, z) flattened by rows (rr = 0, θθ = 4, zz = 8)
    S = np.array([f["Function/Stress"][s][:] for s in steps]) * 1e-6  # MPa

# Cell centres = mean of the cell's node coordinates (x = r, y = z)
r, z = nodes[cells].mean(axis=1)[:, :2].T
# Mid-height row: the cell centres closest to Lz/2 (Lz/2 itself may fall on a cell edge)
z_mid = z[np.argmin(abs(z - Lz / 2))]
row = np.flatnonzero(np.isclose(z, z_mid, atol=1e-9))  # indices of the cells in that row
row = row[np.argsort(r[row])]  # sorted from the inner to the outer radius
r, S = r[row], S[:, row, :] # select only the cells at mid-height

# Stress components: tensor ordered (r, θ, z), so rr = 0, θθ = 4, zz = 8
s_rr, s_tt, s_zz = S[:, :, 0], S[:, :, 4], S[:, :, 8]
zr = r < R_int  # Zircaloy cells
cr = r > R_int  # Chromium cells

# Radial profiles at the last step
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(r[zr] * 1e3, s_rr[-1, zr], "o-", label=r"$\sigma_{rr}$")
ax1.plot(r[zr] * 1e3, s_tt[-1, zr], "o-", label=r"$\sigma_{\theta\theta}$")
ax1.plot(r[zr] * 1e3, s_zz[-1, zr], "o-", label=r"$\sigma_{zz}$")
ax1.set_title("Zircaloy")
ax1.set_xlabel("Radius (mm)")
ax1.set_ylabel("Stress (MPa)")
ax1.grid(alpha=0.5)
ax1.legend()

ax2.plot(r[cr] * 1e3, s_rr[-1, cr], "o-", label=r"$\sigma_{rr}$")
ax2.plot(r[cr] * 1e3, s_tt[-1, cr], "o-", label=r"$\sigma_{\theta\theta}$")
ax2.plot(r[cr] * 1e3, s_zz[-1, cr], "o-", label=r"$\sigma_{zz}$")
ax2.set_title("Chromium")
ax2.set_xlabel("Radius (mm)")
ax2.set_ylabel("Stress (MPa)")
ax2.grid(alpha=0.5)
ax2.legend()

fig.tight_layout()
fig.savefig("output/stress_midheight.png", dpi=200)

# Mean over the wall (weighted by r, cylindrical) vs time
zr_tt = np.average(s_tt[:, zr], axis=1, weights=r[zr])
zr_zz = np.average(s_zz[:, zr], axis=1, weights=r[zr])
cr_tt = np.average(s_tt[:, cr], axis=1, weights=r[cr])
cr_zz = np.average(s_zz[:, cr], axis=1, weights=r[cr])

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(times, zr_tt, "-", color="C1", label=r"Zircaloy $\sigma_{\theta\theta}$")
ax.plot(times, zr_zz, "-", color="C2", label=r"Zircaloy $\sigma_{zz}$")
ax.plot(times, cr_tt, "--", color="C1", label=r"Chromium $\sigma_{\theta\theta}$")
ax.plot(times, cr_zz, "--", color="C2", label=r"Chromium $\sigma_{zz}$")
ax.set(xlabel="Time (s)", ylabel="Mean stress at mid-height (MPa)")
ax.grid(alpha=0.5)
ax.legend()
fig.tight_layout()
fig.savefig("output/stress_midheight_history.png", dpi=200)

print(f"z = {z_mid * 1e3:.1f} mm, t = {times[-1]:.1f} s")
print(f"Zircaloy: hoop = {zr_tt[-1]:.2f} MPa, axial = {zr_zz[-1]:.2f} MPa")
print(f"Chromium: hoop = {cr_tt[-1]:.2f} MPa, axial = {cr_zz[-1]:.2f} MPa")
