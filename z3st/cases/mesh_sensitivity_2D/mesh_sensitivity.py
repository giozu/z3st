#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST script --.. ..- .-.. .-.. ---
"""
Mesh sensitivity script
-----------------------

"""

import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.utils_extract_vtu import extract_field

# Configuration
GEO_FILE = "mesh.geo"
VTU_FILE = os.path.join("output", "fields.vtu")
nx_values = [11, 21, 41, 81]
results = []

# Parameters
Lx, Ly = 0.1, 1.0
k, E, nu, alpha = 48.1, 1.77e11, 0.3, 1.7e-5
Ti, To = 490.0, 480.0
q0, mu = 2.0e6, 24

y_target = Ly / 2


def update_geo(nx):
    """To update nx in mesh.geo"""
    with open(GEO_FILE, "r") as f:
        content = f.read()
    content = re.sub(r"nx\s*=\s*\d+;", f"nx = {nx};", content)
    with open(GEO_FILE, "w") as f:
        f.write(content)


def analytic_T(x):
    """Analytical temperature profile, exponential source, slab"""
    x_rel = x
    term1 = Ti + (To - Ti) * (x_rel / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x_rel / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x_rel) - 1))
    return term1 + term2


def sigma_th_ana(x, T_num, c=1.0):
    """Analytical thermal stress"""
    T_mean = np.trapezoid(T_num, x) / (x.max() - x.min())
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


# --- Sensitivity loop ---
for nx in nx_values:
    current_h = Lx / (nx - 1)
    mask_tol = 0.05

    print(f"Testing mesh: {nx-1} elements | h={current_h:.4e}")

    update_geo(nx)

    subprocess.run(
        ["gmsh", GEO_FILE, "-2", "-o", "mesh.msh", "-format", "msh2"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(["python3", "-m", "z3st"], check=True, stdout=subprocess.DEVNULL)

    x_n, y_n, _, T_all = extract_field(VTU_FILE, field_name="Temperature")
    mask_n = np.abs(y_n - y_target) < mask_tol
    idx_n = np.argsort(x_n[mask_n])
    xn_p, Tn_p = x_n[mask_n][idx_n], T_all[mask_n][idx_n]

    x_c, y_c, _, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
    mask_c = np.abs(y_c - y_target) < mask_tol
    idx_c = np.argsort(x_c[mask_c])
    xc_p, Sn_p = x_c[mask_c][idx_c], S_all[mask_c, 4][idx_c]

    T_ref = analytic_T(xn_p)
    l2_err_T = np.sqrt(np.mean((Tn_p - T_ref) ** 2)) / (Ti - To)

    S_ref = sigma_th_ana(xc_p, analytic_T(xc_p))
    l2_err_S = np.sqrt(np.mean((Sn_p - S_ref) ** 2)) / np.max(np.abs(S_ref))

    print(f"   -> Err T: {l2_err_T:.2e} | Err S: {l2_err_S:.2e}")
    results.append([current_h, l2_err_T, l2_err_S])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
    ax1.plot(xn_p, Tn_p, "ro", label="Numerical T")
    ax1.plot(xn_p, T_ref, "k--", label="Analytical T")
    ax1.set_ylabel("Temperature (K)")
    ax1.legend()

    ax2.plot(xc_p, Sn_p / 1e6, "bo", label="Numerical sigma_yy")
    ax2.plot(xc_p, S_ref / 1e6, "k--", label="Analytical sigma_yy")
    ax2.set_ylabel("Stress (MPa)")
    ax2.set_xlabel("x (m)")
    ax2.legend()

    plt.savefig(f"output/profile_nx_{nx-1}.png")
    plt.close()

res = np.array(results)
h, eT, eS = res[:, 0], res[:, 1], res[:, 2]

plt.figure(figsize=(10, 8))
plt.loglog(h, eT, "bo-", label=r"Error $L_2$ Temperature", linewidth=2)
plt.loglog(h, eS, "ro-", label=r"Error $L_2$ Stress ($\sigma_{yy}$)", linewidth=2)

plt.loglog(h, h**2 * (eT[0] / h[0] ** 2), "k--", alpha=0.5, label="Slope 2 (Theoretical T)")
plt.loglog(h, h * (eS[0] / h[0]), "k:", alpha=0.5, label=r"Slope 1 (Theoretical $\sigma$)")

plt.grid(True, which="both", ls="-", alpha=0.5)
plt.xlabel("h (mesh size)")
plt.ylabel("Relative L2 error")
plt.legend()
plt.title("Convergence study: 2D slab")
plt.savefig("output/convergence_study.png", dpi=300)
plt.show()
