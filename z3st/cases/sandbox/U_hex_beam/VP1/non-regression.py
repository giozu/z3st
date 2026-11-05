#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- beam1d non-regression script --.. ..- .-.. .-.. ---
"""
beam1d -- IAEA VP1 (z-gradient) as an Euler-Bernoulli beam, axial and
bending, clamped at the base, tip free.

References: Oechsner, "Elasto-Plasticity of Frame Structure Elements"
(Springer, 2014), Sects. 5.2-5.3; IWGFR/75 linear theory (DT1BOW); thesis
log, Eq. eigencurv. Layout and checks as in z3st teaching/01_1D.

What the FE actually solves
---------------------------
Section: the undeformed hexagonal wall of geometry.yaml (beam1d.section), the
same D, t, orientation the Z3ST mesh is built from; checked here against the
closed forms and against the 'bottom' facets of the Z3ST mesh.

Two Hermite bending planes (dofs (u, u') per node) and a linear rod for u_z
(beam1d.element). The loads are the free strains of the section, from the T
expression of boundary_conditions.yaml (the string of the Z3ST run):
    eps0*(z)  = 1 / A_w int alpha (T - T_ref) dA
    kappa*(z) = 1 / I_w int alpha (T - T_ref) (x, y) dA .
VP1: T linear in y at every height, slope -25 K / 0.075 m at the core top, so
    kappa*_y = -kappa_0,   kappa_0 = alpha 25 / 0.075 = 6.0e-3 1/m
    eps0*    = alpha 125 K
above the core, both 0 below it and linear across it.

Clamped base and free tip: the beam is statically determinate and nothing
restrains the free strains, so N = 0, M = 0 and the bow towards the cold
corner (+y) is w = u_y with
    w'' = kappa_0 ramp(z),   w(0) = w'(0) = 0 .
Closed form, core from z_b to z_t, L_c = z_t - z_b:
    w = 0                                                       z <= z_b
    w = kappa_0 (z - z_b)^3 / (6 L_c)                           z_b < z <= z_t
    w = kappa_0 [L_c^2 / 6 + L_c (z - z_t) / 2 + (z - z_t)^2 / 2]   z > z_t
    w(core top, LRP, URP) = 1.000, 3.250, 12.250 mm  (IWGFR/75 linear theory)
and u_z = alpha 125 K int_0^z ramp: 1.125, 2.250, 4.500 mm. EI and EA do not
enter w and u_z, only K, M, N and the reactions. T is linear in the section,
so the fibre stress sigma_zz is 0 everywhere.

The exact w is cubic in the core and parabolic above, C1 at the kinks: it
lies in the Hermite space once z_b and z_t are nodes. The FE is then exact
everywhere, already with 3 elements (one per zone), as the bar of
teaching/01_1D is exact with linear elements. The linear rod is exact at the
nodes only: u_z is quadratic in the core, linear between the nodes in the FE
(error h^2 u_z'' / 8 at mid-element, 1.6e-4 of u_z(URP) with h = 0.05 m), and
N = EA (u_z' - eps0*) is exact at the element midpoints.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from z3st.utils.non_regression import error_metric, finish, metric  # noqa: E402

from beam1d import case as case_mod  # noqa: E402
from beam1d import element, solver, thermal  # noqa: E402
from beam1d.mesh3d import section_from_msh  # noqa: E402
from beam1d.section import HexSection  # noqa: E402

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT_DIR, "non-regression.json")

c = case_mod.load_case(CASE_DIR)
geo, mat, sec = c.geo, c.mat, c.section
CASE_3D = os.path.normpath(os.path.join(CASE_DIR, c.inp["case_3d"]))
alpha, E = float(mat["alpha"]), float(mat["E"])
Lz = c.H
z_b, z_t = c.stations["core_bottom"], c.stations["core_top"]
L_c = z_t - z_b
STATIONS = {"core top": z_t, "LRP": c.stations["LRP"], "URP": Lz}

d = np.load(os.path.join(OUT_DIR, "fields.npz"))
z, EI, EA = d["z"], float(d["EI"]), float(d["EA"])
z_s, w_s, v_s, uz_s = d["z_s"], d["uy_s"], d["ux_s"], d["uz_s"]

TOLERANCE = 1e-8

# --.. ..- .-.. .-.. --- analytical reference --.. ..- .-.. .-.. ---
kappa_0 = alpha * 25.0 / 0.075
dT_mean = 125.0


def ramp(zq):
    return np.clip((np.asarray(zq) - z_b) / L_c, 0.0, 1.0)


def w_cf(zq):
    zq = np.asarray(zq, dtype=float)
    core = kappa_0 * (zq - z_b) ** 3 / (6.0 * L_c)
    above = kappa_0 * (L_c**2 / 6.0 + L_c * (zq - z_t) / 2.0 + (zq - z_t) ** 2 / 2.0)
    return np.where(zq <= z_b, 0.0, np.where(zq <= z_t, core, above))


def uz_cf(zq):
    zq = np.asarray(zq, dtype=float)
    core = (zq - z_b) ** 2 / (2.0 * L_c)
    above = L_c / 2.0 + (zq - z_t)
    return alpha * dT_mean * np.where(zq <= z_b, 0.0, np.where(zq <= z_t, core, above))


print(f"[INFO] Geometry  : H = {Lz} m, core {z_b}-{z_t} m, D = {float(geo['D']) * 1e3:.6f} mm (outer, corners)")
print(f"[INFO] Section   : A_w = {sec.A:.6e} m^2, I_w = {sec.I:.6e} m^4, EI = {EI:.4e} N m^2, EA = {EA:.4e} N")
print(f"[INFO] Reference : kappa_0 = alpha 25 / 0.075 = {kappa_0:.6e} 1/m, eps0* = alpha 125 K = {alpha * dT_mean:.6e}")
for name, zz in STATIONS.items():
    print(f"[INFO]             w({name}, z = {zz} m) = {float(w_cf(zz)) * 1e3:.4f} mm, "
          f"u_z = {float(uz_cf(zz)) * 1e3:.4f} mm")

# --.. ..- .-.. .-.. --- FE fields --.. ..- .-.. .-.. ---
w_n = d["U_y"][0::2]
w_at = {name: float(solver.evaluate(z, d["U_y"], zz)[0][0]) for name, zz in STATIONS.items()}
uz_at = {name: float(solver.evaluate_rod(z, d["U_z"], zz)[0][0]) for name, zz in STATIONS.items()}
k_s, M_s = d["ky_s"], d["My_s"]
for name in STATIONS:
    print(f"[INFO] FE        : w({name}) = {w_at[name] * 1e3:.6f} mm, u_z({name}) = {uz_at[name] * 1e3:.6f} mm")

# --.. ..- .-.. .-.. --- plots --.. ..- .-.. .-.. ---
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(19, 5.5), sharey=True)

ax1.plot(k_s * 1e3, z_s, "-", color="#0072B2", label=r"FE $-\kappa^*_y$ (section integral)")
ax1.plot(kappa_0 * ramp(z) * 1e3, z, "o", color="#D55E00", mfc="none", label=r"$\kappa_0$ ramp at the nodes")
ax1.set_xlabel(r"free curvature of the $w$ plane (1/km)", fontsize=12)
ax1.set_ylabel("z (m)", fontsize=12)
ax1.set_title("Free curvature of the section", fontsize=12)

ax2.plot(w_s * 1e3, z_s, "-", color="#0072B2", label=r"FE $w = u_y$ (Hermite)")
ax2.plot(w_n * 1e3, z, "o", color="#0072B2", ms=4)
z_d = np.linspace(0.0, Lz, 400)
ax2.plot(w_cf(z_d) * 1e3, z_d, "k--", lw=1.2, label="closed form")
ax2.set_xlabel(r"$w$ towards the cold corner (mm)", fontsize=12)
ax2.set_title("Bow along the axis", fontsize=12)

ax3.plot((w_s - w_cf(z_s)) / float(w_cf(Lz)), z_s, "-", color="#0072B2",
         label=r"$(w - w_{cf}) / w_{cf}(URP)$")
ax3.plot(M_s / (EI * kappa_0), z_s, "-", color="#D55E00", label=r"$M / (EI\,\kappa_0)$")
ax3.set_xlabel("residual (-)", fontsize=12)
ax3.set_title("FE - closed form, bending moment", fontsize=12)

ax4.plot(uz_s * 1e3, z_s, "-", color="#009E73", label=r"FE $u_z$ (linear rod)")
ax4.plot(d["U_z"] * 1e3, z, "o", color="#009E73", ms=4)
ax4.plot(uz_cf(z_d) * 1e3, z_d, "k--", lw=1.2, label=r"closed form $\alpha\,125\,$K$\int$ramp")
ax4.set_xlabel(r"$u_z$ (mm)", fontsize=12)
ax4.set_title("Axial displacement", fontsize=12)

for ax in (ax1, ax2, ax3, ax4):
    for zz in (z_b, z_t):
        ax.axhline(zz, color="grey", ls=":", lw=0.8)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="best", fontsize=9)
fig.suptitle(rf"IAEA VP1 as an Euler-Bernoulli beam  |  clamped base, {len(z) - 1} elements,  "
             rf"$\kappa_0$ = {kappa_0 * 1e3:.2f} 1/km,  $w$(URP) = {w_at['URP'] * 1e3:.3f} mm,  "
             rf"$u_z$(URP) = {uz_at['URP'] * 1e3:.3f} mm", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.94))
plot_path = os.path.join(OUT_DIR, "bow.png")
plt.savefig(plot_path, dpi=200)
plt.close(fig)
print(f"[INFO] Plot saved to: {plot_path}")

# --.. ..- .-.. .-.. --- section: closed forms and the Z3ST mesh --.. ..- .-.. .-.. ---
# The beam section is the undeformed cross-section of the 3D model: the same
# integrals over the 'bottom' facets of the Z3ST mesh (prismatic wall).
print("\n" + "=" * 72)
print("Section properties: quadrature, closed forms, Z3ST mesh")
print("=" * 72)
msh = section_from_msh(os.path.join(CASE_3D, "mesh.msh"), "bottom")
a_o = sec.a_o
print(f"  {'':>6} {'beam':>14} {'closed form':>14} {'Z3ST mesh':>14}")
for key, val, cf in (("A", sec.A, d["A_w_cf"]), ("I_yy", sec.I_yy, d["I_w_cf"]), ("I_xx", sec.I_xx, d["I_w_cf"]),
                     ("I_xy", sec.I_xy, 0.0), ("S_x", sec.S_x, 0.0), ("S_y", sec.S_y, 0.0)):
    print(f"  {key:>6} {val:14.6e} {float(cf):14.6e} {msh[key]:14.6e}")
flat = HexSection(float(geo["D"]), float(geo["t"]), "flat", c.inp["section"]["n_gauss"])
iso_err = max(abs(sec.I_xx - sec.I_yy), abs(flat.I_xx - sec.I), abs(flat.I_yy - sec.I), abs(sec.I_xy),
              abs(flat.I_xy)) / sec.I
S_err = max(abs(sec.S_x), abs(sec.S_y)) / (sec.A * a_o)
B1_err = float(np.max(np.abs(d["B1"]))) / (float(d["B0"]) * a_o)
mesh_err = max(abs(msh["A"] - sec.A) / sec.A,
               max(abs(msh[k] - getattr(sec, k)) for k in ("I_xx", "I_yy", "I_xy")) / sec.I,
               max(abs(msh["S_x"]), abs(msh["S_y"])) / (sec.A * a_o))
print(f"  isotropy (I_xx = I_yy, I_xy = 0, flat = vertical): {iso_err:.3e}; S / (A a_o) = {S_err:.3e}; "
      f"B1 / (B0 a_o) = {B1_err:.3e}")
print(f"  beam vs Z3ST mesh (max relative): {mesh_err:.3e}")
print("=" * 72)

# --.. ..- .-.. .-.. --- stiffness matrix extraction --.. ..- .-.. .-.. ---
# K_e by Gauss integration of EI N'' N''^T (what beam1d assembles) against the
# closed form of Oechsner Eq. (5.236), on the first element of the mesh; the
# rod against Eq. (5.79). The unconstrained global K of a bending plane has
# exactly two zero eigenvalues (rigid translation and rotation), the rod one.
print("\n" + "=" * 72)
print("Stiffness matrix extraction (Oechsner Eqs. 5.236 and 5.79)")
print("=" * 72)
h = float(z[1] - z[0])
K_e = element.stiffness(EI, h)
K_e_cf = element.stiffness_closed(EI, h)
K_e_err = float(np.max(np.abs(K_e - K_e_cf)))
print(f"Element length h = {h:.4f} m, unit EI / h^3 = {EI / h**3:.4e} N/m")
print("K_e * h^3 / EI (dofs u_1, phi_1, u_2, phi_2; closed form: 12, 6h, 4h^2, 2h^2):")
with np.printoptions(precision=4, linewidth=120, suppress=True):
    print(K_e * h**3 / EI)
print(f"Max |K_e - K_e(5.236)| = {K_e_err:.3e} (max entry {np.max(np.abs(K_e_cf)):.3e})")

xi_r, w_r = element.gauss(2)
dN_r = element.rod_shape(xi_r, h)[1]
K_rod_q = EA * h * (dN_r * w_r) @ dN_r.T
K_rod_err = float(np.max(np.abs(K_rod_q - element.rod_stiffness(EA, h)))) / (EA / h)
print(f"Rod: max |int EA N' N'^T - K_e(5.79)| / (EA / h) = {K_rod_err:.3e}")

K = d["K"]
eig = np.sort(np.linalg.eigvalsh(K))
n_zero = int(np.sum(np.abs(eig) < 1e-10 * eig[-1]))
eig_r = np.sort(np.linalg.eigvalsh(d["K_rod"]))
n_zero_rod = int(np.sum(np.abs(eig_r) < 1e-10 * eig_r[-1]))
print(f"\nGlobal K (bending plane): {K.shape[0]} x {K.shape[1]}, unconstrained")
print(f"  smallest 3 eigenvalues: {eig[:3]}")
print(f"  zero eigenvalues: {n_zero} (rigid-body translation and rotation); rod: {n_zero_rod}")
print("=" * 72)

# --.. ..- .-.. .-.. --- one element, tip load --.. ..- .-.. .-.. ---
# A cantilever of one element with a force P at the tip: the Hermite cubic is
# the exact solution, w(H) = P H^3 / 3EI and w'(H) = P H^2 / 2EI.
P = 1.0e3
U1, _ = solver.solve(element.stiffness(EI, Lz), np.array([0.0, 0.0, P, 0.0]), [0, 1])
tip_err = max(abs(U1[2] - P * Lz**3 / (3 * EI)) / (P * Lz**3 / (3 * EI)),
              abs(U1[3] - P * Lz**2 / (2 * EI)) / (P * Lz**2 / (2 * EI)))
print(f"\n[INFO] 1 element, tip force {P:.0f} N: w(H) = {U1[2] * 1e3:.6f} mm "
      f"(P H^3 / 3EI = {P * Lz**3 / (3 * EI) * 1e3:.6f} mm), max rel. error {tip_err:.3e}")

# --.. ..- .-.. .-.. --- shape functions --.. ..- .-.. .-.. ---
# Hermite functions on one element: N_iu interpolate u, N_iphi the rotation.
# Kronecker checks at the two nodes: N_iu(z_j) = delta_ij, N_iu'(z_j) = 0,
# N_iphi(z_j) = 0, N_iphi'(z_j) = delta_ij.
xi_d = np.linspace(0.0, 1.0, 200)
N, dN, _ = element.shape(xi_d, h)
fig, ax = plt.subplots(figsize=(10, 4.5))
labels = [r"$N_{1u}$", r"$N_{1\varphi}/h$", r"$N_{2u}$", r"$N_{2\varphi}/h$"]
scale = [1.0, 1.0 / h, 1.0, 1.0 / h]
for i, col in enumerate(plt.cm.viridis(np.linspace(0.1, 0.9, 4))):
    ax.plot(xi_d, N[i] * scale[i], "-", color=col, lw=2.0, label=labels[i])
ax.axhline(0.0, color="grey", ls=":", lw=0.6)
ax.axhline(1.0, color="grey", ls=":", lw=0.6)
ax.set_xlabel(r"$\xi = (z - z_1) / h$", fontsize=12)
ax.set_ylabel(r"$N(\xi)$", fontsize=12)
ax.set_title(f"Cubic Hermite shape functions (Oechsner Eqs. 5.222-5.225), h = {h:.3f} m", fontsize=12)
ax.legend(loc="upper center", ncol=4, fontsize=10)
ax.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
sf_path = os.path.join(OUT_DIR, "shape_functions.png")
plt.savefig(sf_path, dpi=200)
plt.close(fig)
print(f"\n[INFO] Shape functions plotted at: {sf_path}")

N_n, dN_n, _ = element.shape(np.array([0.0, 1.0]), h)
kron_err = float(max(np.max(np.abs(N_n[[0, 2]] - np.eye(2))), np.max(np.abs(dN_n[[0, 2]])),
                     np.max(np.abs(N_n[[1, 3]])), np.max(np.abs(dN_n[[1, 3]] - np.eye(2)))))
print(f"       Kronecker check at the nodes: max deviation = {kron_err:.3e}")

# --.. ..- .-.. .-.. --- three elements, one per zone --.. ..- .-.. .-.. ---
# Nodes at 0, z_b, z_t, Lz only, same section integral and T: the nodal values
# are exact, and so is w in between (cubic in the core, parabola above).
_, eps = thermal.build(c, sec)
z3 = solver.mesh(Lz, 1, [z_b, z_t])
z3_g, xi, wg = solver.gauss_points(z3, c.inp["element"]["n_gauss"])
k3 = -sec.curvature(eps, z3_g.ravel())[:, 1].reshape(z3_g.shape)
K3, f3 = solver.assemble(z3, EI, k3, xi, wg)
U3, _ = solver.solve(K3, f3, [0, 1])
w3_s = solver.evaluate(z3, U3, z_d)[0]
w3_err = float(np.max(np.abs(w3_s - w_cf(z_d))) / float(w_cf(Lz)))
print(f"\n[INFO] 3 elements (nodes {z3.tolist()} m): w(URP) = {U3[-2] * 1e3:.6f} mm, "
      f"max |w - w_cf| / w_cf(URP) = {w3_err:.3e}")

# Galerkin reconstruction on the 3 elements: w(z) = sum_i (u_i N_iu + phi_i N_iphi)
fig, ax = plt.subplots(figsize=(10, 5))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(U3)))
for j in range(len(U3)):
    if U3[j] == 0.0:
        continue
    Uj = np.zeros_like(U3)
    Uj[j] = U3[j]
    kind = "u" if j % 2 == 0 else r"\varphi"
    ax.plot(z_d, solver.evaluate(z3, Uj, z_d)[0] * 1e3, ":", color=colors[j], lw=1.4,
            label=rf"${kind}_{{{j // 2}}}\,N_{{{j // 2}{kind}}}$")
ax.plot(z_d, w3_s * 1e3, "k-", lw=2.5, label=r"$w = \sum_i (u_i N_{iu} + \varphi_i N_{i\varphi})$")
ax.plot(z_d, w_cf(z_d) * 1e3, "--", color="#D55E00", lw=1.5, label="closed form")
ax.plot(z3, U3[0::2] * 1e3, "ko", mfc="white", ms=7, mew=1.5)
ax.set_xlabel("z (m)", fontsize=12)
ax.set_ylabel(r"$w$ (mm)", fontsize=12)
ax.set_title("Galerkin reconstruction on 3 Hermite elements (one per zone)", fontsize=12)
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
decomp_path = os.path.join(OUT_DIR, "displacement_decomposition.png")
plt.savefig(decomp_path, dpi=200)
plt.close(fig)
print(f"[INFO] Displacement decomposition plotted at: {decomp_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
i3 = int(np.argmin(np.abs(z_s - 3.0)))
k_above, e0_above = float(k_s[i3]), float(d["e0_s"][i3])
R = np.concatenate([d["R_y"], d["R_x"]])
R_scale = np.array([EI * kappa_0 / Lz, EI * kappa_0] * (len(R) // 2))   # force, moment
N_scale = EA * alpha * dT_mean
sig_scale = E * alpha * 50.0
errors = {
    "A_w": metric(d["A_w"], d["A_w_cf"]),
    "I_w": metric(d["I_w"], d["I_w_cf"]),
    "section_isotropy": error_metric(iso_err),
    "section_centroid": error_metric(S_err),
    "B1": error_metric(B1_err),
    "section_vs_3d_mesh": error_metric(mesh_err),
    "kappa_above_core": metric(k_above, kappa_0),
    "kappa_perp": error_metric(np.max(np.abs(d["kx_s"])), rel=np.max(np.abs(d["kx_s"])) / kappa_0),
    "eps0_above_core": metric(e0_above, alpha * dT_mean),
    **{f"w_{name.replace(' ', '_')}": metric(w_at[name], w_cf(zz)) for name, zz in STATIONS.items()},
    "w_along_z": error_metric(np.max(np.abs(w_s - w_cf(z_s))), rel=np.max(np.abs(w_s - w_cf(z_s))) / w_cf(Lz)),
    "v_URP": error_metric(abs(v_s[-1]), rel=abs(v_s[-1]) / w_cf(Lz)),
    **{f"uz_{name.replace(' ', '_')}": metric(uz_at[name], uz_cf(zz)) for name, zz in STATIONS.items()},
    # linear rod: exact at the nodes, linear in between (u_z is quadratic in the core)
    "uz_nodes": error_metric(np.max(np.abs(d["U_z"] - uz_cf(z))),
                             rel=np.max(np.abs(d["U_z"] - uz_cf(z))) / uz_cf(Lz)),
    "M_max": error_metric(np.max(np.abs(M_s)), rel=np.max(np.abs(M_s)) / (EI * kappa_0)),
    "N_max": error_metric(np.max(np.abs(d["N_mid"])), rel=np.max(np.abs(d["N_mid"])) / N_scale),
    "clamp_reactions": error_metric(np.max(np.abs(R)), rel=np.max(np.abs(R) / R_scale)),
    "clamp_reaction_axial": error_metric(abs(float(d["R_z"][0])), rel=abs(float(d["R_z"][0])) / N_scale),
    "sigma_zz_max": error_metric(np.max(np.abs(d["sigma_mid"])), rel=np.max(np.abs(d["sigma_mid"])) / sig_scale),
    "sigma_resultants": error_metric(max(np.max(np.abs(d["N_sig_mid"])) / N_scale,
                                         np.max(np.abs(d["Mx_sig_mid"])) / (EI * kappa_0),
                                         np.max(np.abs(d["My_sig_mid"])) / (EI * kappa_0))),
    "K_e_vs_5.236": error_metric(K_e_err, rel=K_e_err / np.max(np.abs(K_e_cf))),
    "K_rod_vs_5.79": error_metric(K_rod_err),
    "rigid_modes": metric(n_zero, 2),
    "rigid_modes_rod": metric(n_zero_rod, 1),
    "tip_load_1_element": error_metric(tip_err),
    "hermite_kronecker": error_metric(kron_err),
    "w_3_elements": error_metric(w3_err),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
