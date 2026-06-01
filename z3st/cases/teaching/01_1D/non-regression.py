#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: teaching/01_1D  --  1D bar in uniaxial tension on a true 1D mesh.

Reference: Bower, "Applied Mechanics of Solids" (CRC, 2010),
ch. 7.2 + ch. 8.1.5.

What the FE actually solves
---------------------------
The mesh has tdim = 1 (line elements). The displacement function space is a
1-component vector field u(x) = u_x(x). The strain tensor is therefore
non-zero only in the axial direction:
    eps_xx = du_x/dx,   eps_yy = eps_zz = eps_xy = eps_xz = eps_yz = 0
There is no transverse motion modelled at all, so the bar is *implicitly
constrained* in y and z (think: a 1D fibre inside an infinitely rigid sleeve).

Applying the 3D isotropic Hooke's law sigma = lambda tr(eps) I + 2G eps to
this constrained strain gives:
    sigma_xx = (lambda + 2G) eps_xx
    sigma_yy = sigma_zz = lambda eps_xx      (Poisson reaction stresses)
    sigma_xy = sigma_yz = sigma_xz = 0

With sigma_xx = P (the applied traction):
    eps_xx     = P / (lambda + 2G)
    sigma_yy   = sigma_zz = nu/(1 - nu) * P
    u_x(x)     = eps_xx * x
    u_x(L)     = P * L / (lambda + 2G)

Equivalent in (E, nu): lambda + 2G = E (1 - nu) / [(1 + nu)(1 - 2nu)].
For E = 200 GPa, nu = 0.3:
    lambda + 2G = 269.23 GPa
    nu/(1 - nu) = 0.4286 -> sigma_yy = sigma_zz = 42.86 MPa
    u_x(L = 1 m, P = 100 MPa) = 0.3714 mm

This is *not* the engineering bar formula u_x = P L / E that Bower's ch. 7.2
hand-derives -- that one assumes free transverse contraction (sigma_yy =
sigma_zz = 0). The engineering bar is recovered in case 01_3D, where the
mesh is truly 3D and the transverse Poisson contraction is free.

The pedagogical pair (01_1D vs 01_3D) shows that a "1D" mesh in FE software
is not the same as the engineering 1D model -- it carries an implicit lateral
constraint that bumps the axial stiffness from E to (lambda + 2G).
"""

import os
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import dolfinx
import ufl
from dolfinx.fem.petsc import assemble_matrix

from z3st.core.mesh import load_mesh
from z3st.utils.utils_extract_vtu import (
    list_fields,
    extract_field,
    extract_displacement,
)
from z3st.utils.utils_verification import pass_fail_check, regression_check

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(GEOMETRY_FILE, "r") as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get("Lx"))

with open(MATERIAL_FILE, "r") as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get("E"))
nu = float(mat_data.get("nu"))

# Lame parameters derived the same way z3st does (see spine.py:97-98).
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
M_constrained = lmbda + 2.0 * G                # axial "constrained" modulus

with open(BC_FILE, "r") as f:
    bc_data = yaml.safe_load(f)
mech_list = bc_data.get("mechanical", {}).get("steel", [])
P = next(
    (float(bc["traction"]) for bc in mech_list
     if bc.get("type") == "Neumann" and bc.get("region") == "xmax"),
    None,
)

with open(MESH_GEO_FILE, "r") as f:
    content = f.read()
# In mesh.geo `nx` is the Transfinite Curve parameter, which is the number of
# *nodes* on the curve (endpoints included). The number of line elements is
# therefore (nx - 1).
nx = int(re.search(r"nx\s*=\s*(\d+);", content).group(1)) - 1   # line elements

print(f"[INFO] Geometry  : Lx = {Lx} m (true 1D mesh)")
print(f"[INFO] Material  : E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO]            lambda = {lmbda:.3e}, G = {G:.3e}, "
      f"lambda + 2G = {M_constrained:.3e} Pa")
print(f"[INFO] Load      : P = {P:.2e} Pa ({P * 1e-6:.1f} MPa)")
print(f"[INFO] Mesh      : nx = {nx} line elements")

TOLERANCE = 1e-4

# --.. ..- .-.. .-.. --- analytical reference (constrained bar) ----------------
sigma_xx_ref = P
sigma_yy_ref = nu / (1.0 - nu) * P             # = lambda * eps_xx
sigma_zz_ref = sigma_yy_ref
eps_xx_ref = P / M_constrained
u_xL_ref = eps_xx_ref * Lx

print(f"[INFO] Constrained-bar ref: eps_xx = {eps_xx_ref:.4e}, "
      f"u_x(L) = {u_xL_ref * 1e3:.4f} mm, "
      f"sigma_yy = sigma_zz = {sigma_yy_ref * 1e-6:.2f} MPa")

# --.. ..- .-.. .-.. --- inspect VTU --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- extract FE fields --.. ..- .-.. .-.. ---
# Stress is per-cell. On a 1D mesh every cell is on the "axis", so we just
# extract everything and sort by x.
x_S, _y_S, _z_S, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
sort_idx = np.argsort(x_S)
x_s = x_S[sort_idx]
sigma_xx = S_all[sort_idx, 0]
sigma_yy = S_all[sort_idx, 4]
sigma_zz = S_all[sort_idx, 8]

# Displacement is per-node. On a 1D mesh u_vec is a 1D array of shape
# (n_nodes,) -- the single axial component. On 2D/3D it would be 2D.
x_n, _y_n, _z_n, u_vec = extract_displacement(VTU_FILE)
u_x = u_vec if u_vec.ndim == 1 else u_vec[:, 0]
sort_n = np.argsort(x_n)
x_n_sorted = x_n[sort_n]
u_x_sorted = u_x[sort_n]

# Tip displacement (rightmost node).
u_xL_num = float(u_x_sorted[-1])

print(f"[INFO] FE tip disp: u_x(L) = {u_xL_num * 1e3:.6f} mm "
      f"(ref {u_xL_ref * 1e3:.6f} mm)")

# --.. ..- .-.. .-.. --- plots --.. ..- .-.. .-.. ---
Pa_to_MPa = 1e-6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(x_s, sigma_xx * Pa_to_MPa, "bo-", label=r"FE $\sigma_{xx}$",
         markersize=6, alpha=0.8)
ax1.plot(x_s, sigma_yy * Pa_to_MPa, "rs-", label=r"FE $\sigma_{yy}$",
         markersize=6, alpha=0.8)
ax1.plot(x_s, sigma_zz * Pa_to_MPa, "m^-", label=r"FE $\sigma_{zz}$",
         markersize=6, alpha=0.8)
ax1.axhline(P * Pa_to_MPa, color="k", linestyle="--", linewidth=1.2,
            label=rf"Analytic $\sigma_{{xx}} = P = {P*Pa_to_MPa:.0f}$ MPa")
ax1.axhline(sigma_yy_ref * Pa_to_MPa, color="r", linestyle=":", linewidth=1.0,
            label=rf"Analytic $\sigma_{{yy}}=\sigma_{{zz}}=\frac{{\nu}}{{1-\nu}}P"
                  rf"={sigma_yy_ref*Pa_to_MPa:.1f}$ MPa")
ax1.set_xlabel("x (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.set_title("Stress along the bar (1D mesh)", fontsize=13)
ax1.legend(loc="best", fontsize=9)
ax1.grid(True, linestyle="--", alpha=0.6)

ax2.plot(x_n_sorted, u_x_sorted * 1e3, "bo", label=r"FE $u_x$", markersize=7)
x_dense = np.linspace(0.0, Lx, 200)
ax2.plot(x_dense, eps_xx_ref * x_dense * 1e3, "k--", linewidth=1.5,
         label=r"Analytic $u_x = P\,x/(\lambda+2G)$")
ax2.set_xlabel("x (m)", fontsize=12)
ax2.set_ylabel(r"$u_x$ (mm)", fontsize=12)
ax2.set_title("Axial displacement along the bar", fontsize=13)
ax2.legend(loc="best")
ax2.grid(True, linestyle="--", alpha=0.6)

fig.suptitle(
    rf"1D bar in tension (true 1D mesh)  |  $P$ = {P*Pa_to_MPa:.0f} MPa,  "
    rf"$L$ = {Lx} m,  $E$ = {E*1e-9:.0f} GPa,  $\nu$ = {nu}",
    fontsize=13,
)
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "bar_tension.png")
plt.savefig(plot_path, dpi=200)
print(f"[INFO] Plot saved to: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
sigma_xx_err = float(np.max(np.abs(sigma_xx - sigma_xx_ref)) / P)
sigma_yy_err = float(np.max(np.abs(sigma_yy - sigma_yy_ref)) / P)
sigma_zz_err = float(np.max(np.abs(sigma_zz - sigma_zz_ref)) / P)
u_xL_err = float(abs(u_xL_num - u_xL_ref) / u_xL_ref)

print("\n[RESULT] Relative errors vs constrained-bar analytical:")
print(f"   sigma_xx ~ P                    : {sigma_xx_err:.3e}")
print(f"   sigma_yy ~ nu/(1-nu)*P          : {sigma_yy_err:.3e}")
print(f"   sigma_zz ~ nu/(1-nu)*P          : {sigma_zz_err:.3e}")
print(f"   u_x(L)  ~ P*L / (lambda + 2G)   : {u_xL_err:.3e}")

errors = {
    "sigma_xx_max_err": {
        "numerical": float(np.max(sigma_xx)),
        "reference": float(sigma_xx_ref),
        "abs_error": float(np.max(np.abs(sigma_xx - sigma_xx_ref))),
        "rel_error": sigma_xx_err,
    },
    "sigma_yy_max_err": {
        "numerical": float(np.max(sigma_yy)),
        "reference": float(sigma_yy_ref),
        "abs_error": float(np.max(np.abs(sigma_yy - sigma_yy_ref))),
        "rel_error": sigma_yy_err,
    },
    "sigma_zz_max_err": {
        "numerical": float(np.max(sigma_zz)),
        "reference": float(sigma_zz_ref),
        "abs_error": float(np.max(np.abs(sigma_zz - sigma_zz_ref))),
        "rel_error": sigma_zz_err,
    },
    "u_xL": {
        "numerical": u_xL_num,
        "reference": u_xL_ref,
        "abs_error": abs(u_xL_num - u_xL_ref),
        "rel_error": u_xL_err,
    },
}

# --.. ..- .-.. .-.. --- stiffness matrix extraction ------------------------
# Re-assemble the global stiffness matrix using exactly the same mesh, function
# space and bilinear form that z3st used during the run. With the default
# nx = 11 (10 line elements) you get an 11x11 tridiagonal matrix; with nx = 1
# (gmsh clamps to 2 nodes -> 1 line element) you get the canonical 2x2 K_e
# from Bower ch. 7.2 / 8.1.5:
#
#     K_e = (lambda + 2G) / L_e * [[+1, -1], [-1, +1]]
#
# This block is purely diagnostic / pedagogical -- it does not affect the
# pass/fail of the case, but it makes z3st's FE assembly visible and lets the
# student verify it by hand against the textbook formula.

print("\n" + "=" * 72)
print("Stiffness matrix extraction (Bower ch. 7.2 / 8.1.5)")
print("=" * 72)

MESH_FILE = os.path.join(CASE_DIR, "mesh.msh")
mesh, _ct, _ft = load_mesh(MESH_FILE, comm=MPI.COMM_WORLD)
print(f"Mesh  : tdim = {mesh.topology.dim}, "
      f"{mesh.topology.index_map(0).size_local} nodes, "
      f"{mesh.topology.index_map(mesh.topology.dim).size_local} cells")

# Same V_m z3st builds (see core/finite_element_setup.py:32):
V = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (mesh.topology.dim,)))
n_dofs = V.dofmap.index_map.size_local * V.dofmap.bs
print(f"V_m   : {V.ufl_element()}")
print(f"#dofs : {n_dofs}")

# Same bilinear form z3st builds on a 1D mesh: see our tdim==1 patch in
# models/mechanical_model.py::epsilon (3x3 padded strain with only eps_xx),
# combined with the default lame branch (lambda * tr(eps) * I + 2G * eps).
u_tr = ufl.TrialFunction(V)
v_te = ufl.TestFunction(V)

def _epsilon_1d(w):
    return ufl.as_tensor([[w[0].dx(0), 0.0, 0.0],
                          [0.0,        0.0, 0.0],
                          [0.0,        0.0, 0.0]])

eps_u = _epsilon_1d(u_tr)
sigma_u = lmbda * ufl.tr(eps_u) * ufl.Identity(3) + 2.0 * G * eps_u
a_form = dolfinx.fem.form(ufl.inner(sigma_u, _epsilon_1d(v_te)) * ufl.dx(domain=mesh))

K_petsc = assemble_matrix(a_form)
K_petsc.assemble()

# Convert sparse PETSc matrix to dense numpy for printing.
indptr, indices, data = K_petsc.getValuesCSR()
K = np.zeros((n_dofs, n_dofs))
for i in range(n_dofs):
    for j_idx in range(indptr[i], indptr[i + 1]):
        K[i, indices[j_idx]] = data[j_idx]

# Spatial coords of each dof, so the user can see what each row/col refers to.
dof_coords = V.tabulate_dof_coordinates()
order = np.argsort(dof_coords[:, 0])   # left-to-right node ordering
print("Dof -> coordinate (sorted by x):")
for i in order:
    print(f"   dof {i:2d}  ->  x = {dof_coords[i, 0]:.4f} m")

# Pretty-print K in element-stiffness units (M/L_e_uniform).
# For a uniform mesh, h = Lx / n_elements, so the "natural unit" of K is M/h.
n_elements = mesh.topology.index_map(mesh.topology.dim).size_local
h = Lx / n_elements
unit = M_constrained / h
print()
print(f"Element length h        = Lx / N = {Lx} / {n_elements} = {h:.4f} m")
print(f"Element stiffness unit  = (lambda + 2G) / h = {unit:.4e} N/m^2")

print(f"\nGlobal stiffness K (raw, N/m^2):")
with np.printoptions(precision=4, linewidth=120, suppress=True):
    if n_dofs <= 12:
        print(K)
    else:
        print(f"  [{n_dofs}x{n_dofs} matrix omitted; showing K * h / M instead]")

print(f"\nNormalised K * h / (lambda + 2G)   (analytical: 1 elem -> [[+1,-1],[-1,+1]];")
print(f"                                                 N elems -> tridiag with +2 on inner diag, -1 off):")
K_norm = K / unit
with np.printoptions(precision=4, linewidth=120, suppress=True):
    print(K_norm)

# Eigenvalues: there must be exactly one zero eigenvalue (rigid-body
# translation along x); all others are strictly positive.
eigvals = np.sort(np.linalg.eigvalsh(K))
print(f"\nEigenvalues of K (sorted):")
print(f"  smallest 3: {eigvals[:3]}")
print(f"  largest  3: {eigvals[-3:]}")
print(f"  (one zero = rigid-body translation; rank deficiency = {np.sum(np.abs(eigvals) < unit * 1e-10)})")

# For the single-element case, verify against the closed-form K_e by hand.
if n_dofs == 2:
    K_analytic = unit * np.array([[+1.0, -1.0],
                                  [-1.0, +1.0]])
    err = float(np.max(np.abs(K - K_analytic)))
    print(f"\n[1 element] Closed-form K_e = (lambda+2G)/L_e * [[+1,-1],[-1,+1]]:")
    with np.printoptions(precision=4, linewidth=120, suppress=True):
        print(K_analytic)
    print(f"   Max |K_FE - K_analytic| = {err:.3e}")

print("=" * 72)

# --.. ..- .-.. .-.. --- shape functions visualisation -----------------------
# The Lagrange P1 shape functions N_i(x) are the "hat functions" that the FE
# uses as a basis. They are never built explicitly in z3st/dolfinx -- basix
# evaluates them on the fly during assembly. But we can recover N_i by setting
# the i-th dof to 1 and all others to 0: the resulting Function *is* N_i.

from dolfinx.geometry import (
    bb_tree, compute_collisions_points, compute_colliding_cells,
)

# Sample many points along the bar.
xs_sample = np.linspace(0.0, Lx, 400)
points = np.column_stack([xs_sample,
                          np.zeros_like(xs_sample),
                          np.zeros_like(xs_sample)])

# Find which cell contains each sample point (geometric search).
tree = bb_tree(mesh, mesh.topology.dim)
candidates = compute_collisions_points(tree, points)
collisions = compute_colliding_cells(mesh, candidates, points)
cells_for_pt = np.array([
    collisions.links(i)[0] if len(collisions.links(i)) > 0 else 0
    for i in range(len(points))
], dtype=np.int32)

# Build N_i for each dof and sample it.
N_func = dolfinx.fem.Function(V)
N_samples = np.zeros((n_dofs, len(xs_sample)))
for i in range(n_dofs):
    N_func.x.array[:] = 0.0
    N_func.x.array[i] = 1.0
    vals = N_func.eval(points, cells_for_pt)
    # vals shape is (n_points, value_size) -- value_size = 1 on a 1D mesh.
    N_samples[i, :] = vals[:, 0] if vals.ndim == 2 else vals

# Plot.
fig, ax = plt.subplots(figsize=(10, 4.5))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_dofs))
for i in range(n_dofs):
    # Walk through dofs in left-to-right order for the label.
    rank_i = int(np.where(order == i)[0][0])
    ax.plot(xs_sample, N_samples[i, :], color=colors[rank_i], linewidth=2.0,
            label=rf"$N_{{{rank_i}}}(x)$  (node at $x={dof_coords[i,0]:.3f}$ m)")
    ax.axvline(dof_coords[i, 0], color=colors[rank_i], linestyle=":",
               linewidth=0.8, alpha=0.5)

# Mark element boundaries.
for x_node in np.sort(dof_coords[:, 0]):
    ax.plot(x_node, 1.0, "k|", markersize=15)
    ax.plot(x_node, 0.0, "ko", markersize=5)

ax.set_xlabel("x (m)", fontsize=12)
ax.set_ylabel(r"$N_i(x)$", fontsize=12)
ax.set_title(f"Lagrange P1 (CG1) shape functions on the current mesh "
             f"({n_dofs} dofs, {n_elements} element"
             f"{'s' if n_elements > 1 else ''})", fontsize=12)
ax.axhline(1.0, color="grey", linestyle=":", linewidth=0.6)
ax.axhline(0.0, color="grey", linestyle=":", linewidth=0.6)
ax.set_ylim(-0.1, 1.15)
ax.legend(loc="upper center", ncol=min(n_dofs, 4), fontsize=9, frameon=True)
ax.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

sf_path = os.path.join(CASE_DIR, "output", "shape_functions.png")
plt.savefig(sf_path, dpi=200)
print(f"\n[INFO] Shape functions plotted at: {sf_path}")

# Kronecker delta check, evaluated at the EXACT node positions (no sampling
# bias): build a (n_dofs x n_dofs) table N_ij = N_i(x_j). Should be Identity.
node_points = np.column_stack([dof_coords[:, 0],
                               np.zeros(n_dofs),
                               np.zeros(n_dofs)])
candidates_n = compute_collisions_points(tree, node_points)
collisions_n = compute_colliding_cells(mesh, candidates_n, node_points)
cells_node = np.array([collisions_n.links(j)[0] for j in range(n_dofs)], dtype=np.int32)
N_at_nodes = np.zeros((n_dofs, n_dofs))
for i in range(n_dofs):
    N_func.x.array[:] = 0.0
    N_func.x.array[i] = 1.0
    vals = N_func.eval(node_points, cells_node)
    N_at_nodes[i, :] = vals[:, 0] if vals.ndim == 2 else vals
kron_err = float(np.max(np.abs(N_at_nodes - np.eye(n_dofs))))
print(f"       Kronecker delta check at exact node coords: "
      f"max |N_i(x_j) - delta_ij| = {kron_err:.3e}")

# --.. ..- .-.. .-.. --- displacement reconstruction via N_i u_i -------------
# Show that the FE displacement u_x(x) is literally a weighted sum of the
# shape functions, with weights = nodal displacements u_i. This is the
# Galerkin trial-function representation u(x) = sum_i u_i N_i(x).

# Pick the nodal displacements at the (already-sorted) dof positions.
u_i = np.array([u_x_sorted[np.argmin(np.abs(x_n_sorted - dof_coords[i, 0]))]
                for i in range(n_dofs)])

print(f"\nNodal displacements u_i (the 'weights'):")
for i in order:
    print(f"   u_{int(np.where(order == i)[0][0])}  (at x = {dof_coords[i, 0]:.3f} m)  "
          f"=  {u_i[i] * 1e3:.4f} mm")

# Reconstruct u(x) = sum_i u_i N_i(x) and verify it matches the FE field.
u_reconstructed = np.sum(u_i[:, None] * N_samples, axis=0)

# Sample the FE u_x at the same xs_sample for a direct overlay.
u_func = dolfinx.fem.Function(V)
u_func.x.array[:] = 0.0
for i in range(n_dofs):
    u_func.x.array[i] = u_i[i]
u_fe_at_xs = u_func.eval(points, cells_for_pt)
u_fe_at_xs = u_fe_at_xs[:, 0] if u_fe_at_xs.ndim == 2 else u_fe_at_xs

recon_err = float(np.max(np.abs(u_reconstructed - u_fe_at_xs)))
print(f"\nReconstruction check: max |sum(u_i N_i) - u_FE| = {recon_err:.3e} m")

fig, ax = plt.subplots(figsize=(10, 5))

# Each scaled shape function u_i * N_i(x) (dotted, light)
for i in range(n_dofs):
    rank_i = int(np.where(order == i)[0][0])
    ax.plot(xs_sample, u_i[i] * N_samples[i, :] * 1e3,
            color=colors[rank_i], linestyle=":", linewidth=1.3, alpha=0.8,
            label=rf"$u_{{{rank_i}}}\,N_{{{rank_i}}}(x)$"
                  rf"  ($u_{{{rank_i}}}={u_i[i]*1e3:.3f}$ mm)")

# Sum = reconstructed FE displacement (thick black)
ax.plot(xs_sample, u_reconstructed * 1e3, "k-", linewidth=2.5,
        label=r"$u_x(x) = \sum_i u_i\, N_i(x)$  (sum)")

# Analytical reference (constrained-bar PL/M * x)
ax.plot(xs_sample, eps_xx_ref * xs_sample * 1e3, "r--", linewidth=1.5,
        alpha=0.8, label=r"Analytical  $u_x = P\,x/(\lambda+2G)$")

# Mark nodes
for i in range(n_dofs):
    ax.plot(dof_coords[i, 0], u_i[i] * 1e3, "ko", markersize=7,
            markerfacecolor="white", markeredgewidth=1.5)

ax.set_xlabel("x (m)", fontsize=12)
ax.set_ylabel(r"$u_x$ (mm)", fontsize=12)
ax.set_title(f"Galerkin reconstruction: $u_x(x) = \\sum_i u_i\\,N_i(x)$ "
             f"({n_dofs} dofs, {n_elements} element"
             f"{'s' if n_elements > 1 else ''})", fontsize=12)
ax.axhline(0.0, color="grey", linestyle=":", linewidth=0.6)
ax.legend(loc="upper left", fontsize=9, frameon=True)
ax.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

decomp_path = os.path.join(CASE_DIR, "output", "displacement_decomposition.png")
plt.savefig(decomp_path, dpi=200)
print(f"[INFO] Displacement decomposition plotted at: {decomp_path}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
