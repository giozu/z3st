# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Fragmented-pellet relocation solver (Phase 1: mechanics-only, prescribed T).
#
# An explicitly cracked (fragmented) UO2 pellet inside a zircaloy cladding, 2D
# plane-strain r-theta section. Each fuel fragment is a free body held only by
# multi-body contact against the cladding (dolfinx_contact, frictionless
# Nitsche). A prescribed parabolic T(r) drives differential thermal expansion;
# the fragments relocate outward and close the fuel-clad gap. This is the
# discrete reference the smeared Barani softening model is later compared with.
#
# Phase 1 scope: fragment-to-cladding contact only (crack faces are free
# surfaces); prescribed T(r); no thermal solve. Fragment-to-fragment
# self-contact and the coupled thermal solve arrive in later phases.
#
# The fragments carry no Dirichlet constraint; global and per-fragment rigid
# body motion is handled by a per-subdomain near-null space in the CG/GAMG
# solve. The thermal load on each free fragment is self-equilibrated, so the
# singular tangent is consistent.
#
#   python relocation_solve.py

from pathlib import Path

import numpy as np
import yaml
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem, mesh as dmesh, la
from dolfinx.fem import (extract_function_spaces, form, functionspace, Function)
from dolfinx.fem.petsc import (apply_lifting, assemble_matrix, assemble_vector,
                               create_vector, set_bc)
from dolfinx.io import gmsh as dgmsh, VTXWriter
from dolfinx.graph import adjacencylist

import dolfinx_contact.cpp
from dolfinx_contact.helpers import (lame_parameters, epsilon, sigma_func,
                                     rigid_motions_nullspace_subdomains)
from dolfinx_contact.general_contact.contact_problem import (ContactProblem,
                                                             FrictionLaw)
from dolfinx_contact.newton_solver import NewtonSolver

from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent

# --. physical-group tags (from make_mesh.py) --..
TAG_FUEL, TAG_CLAD = 1, 2
TAG_OUTER = 10        # fuel fragment outer arcs (candidate contact surface)
TAG_CLAD_INNER = 60   # cladding inner surface

# --. prescribed thermal field (parabolic, centreline-hot) --..
R_PELLET = 5.0e-3     # fuel outer radius (m)
GAP0 = 1.0e-4         # as-fabricated fuel-clad gap (m) = RCI - R_out
T_CENTRE = 1400.0     # K, pellet centreline
T_SURF = 700.0        # K, pellet rim and cladding
T_REF = 600.0         # K, stress-free reference

# --. Nitsche / solver knobs --..
ORDER = 1
GAMMA = 10.0
THETA = 1.0           # symmetric Nitsche
Q_DEGREE = 3
N_LOAD_STEPS = 6
# Weak elastic foundation to regularise the rigid-body modes of fragments that
# are not yet in contact (otherwise the tangent is singular). Small relative to
# the elastic stiffness E/R^2, so it perturbs the relocation negligibly. A
# Phase-1 expedient; later phases replace it with proper null-space handling.
K_FOUND_FACTOR = 1.0e-5


def load_material(fname):
    with open(HERE / fname) as f:
        m = yaml.safe_load(f)
    return float(m["E"]), float(m["nu"]), float(m["alpha"])


def split_fragments(msh, cell_tags):
    """Cell MeshTags with one unique id per connected body: each fuel fragment
    gets id 0..N-1, the cladding gets id N. Fragments are separated by the open
    crack gaps, so connected components on the fuel dual graph isolates them."""
    tdim = msh.topology.dim
    clad_cells = cell_tags.indices[cell_tags.values == TAG_CLAD]
    fuel_cells = cell_tags.indices[cell_tags.values == TAG_FUEL]
    fuel_set = set(int(c) for c in fuel_cells)
    local = {int(c): i for i, c in enumerate(fuel_cells)}

    msh.topology.create_connectivity(tdim, tdim - 1)
    c2f = msh.topology.connectivity(tdim, tdim - 1)
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    rows, cols = [], []
    for c in fuel_cells:
        for fac in c2f.links(c):
            for nb in f2c.links(fac):
                if int(nb) != int(c) and int(nb) in fuel_set:
                    rows.append(local[int(c)])
                    cols.append(local[int(nb)])
    adj = coo_matrix((np.ones(len(rows)), (rows, cols)),
                     shape=(fuel_cells.size, fuel_cells.size))
    n_frag, labels = connected_components(adj, directed=False)

    idx = np.concatenate([fuel_cells, clad_cells]).astype(np.int32)
    val = np.concatenate([labels.astype(np.int32),
                          np.full(clad_cells.size, n_frag, dtype=np.int32)])
    o = np.argsort(idx)
    return dmesh.meshtags(msh, tdim, idx[o], val[o]), n_frag


def temperature(x, scale=1.0):
    """Parabolic profile in the fuel, clamped to T_SURF in clad/gap; ramped."""
    r = np.sqrt(x[0] ** 2 + x[1] ** 2)
    full = T_SURF + (T_CENTRE - T_SURF) * np.maximum(0.0, 1.0 - (r / R_PELLET) ** 2)
    return T_REF + scale * (full - T_REF)


def main():
    comm = MPI.COMM_WORLD
    md = dgmsh.read_from_msh(str(HERE / "mesh.msh"), comm, rank=0, gdim=2)
    msh, cell_tags, facet_tags = md.mesh, md.cell_tags, md.facet_tags
    tdim = msh.topology.dim

    msh.topology.create_connectivity(tdim - 1, 0)
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_connectivity(tdim, tdim)

    domain_marker, n_frag = split_fragments(msh, cell_tags)
    print(f"[relocation] {n_frag} fuel fragments + 1 cladding, "
          f"{domain_marker.values.size} cells, {msh.geometry.x.shape[0]} nodes")

    E_f, nu_f, alpha_f = load_material("fuel.yaml")
    E_c, nu_c, alpha_c = load_material("clad.yaml")

    # --. spaces --..
    V = functionspace(msh, ("Lagrange", ORDER, (tdim,)))
    u = Function(V, name="u")          # accumulated displacement
    du = Function(V, name="du")        # Newton increment
    v = ufl.TestFunction(V)
    w = ufl.TrialFunction(V)

    # --. per-body material fields (DG-0) --..
    V0 = functionspace(msh, ("DG", 0))
    mu0, lmbda0, beta0, fric = (Function(V0) for _ in range(4))
    mu_fn, lmbda_fn = lame_parameters(plane_strain=True)
    for cells, E, nu, al in ((cell_tags.indices[cell_tags.values == TAG_FUEL],
                              E_f, nu_f, alpha_f),
                             (cell_tags.indices[cell_tags.values == TAG_CLAD],
                              E_c, nu_c, alpha_c)):
        mu, lm = mu_fn(E, nu), lmbda_fn(E, nu)
        mu0.x.array[cells] = mu
        lmbda0.x.array[cells] = lm
        beta0.x.array[cells] = (3.0 * lm + 2.0 * mu) * al    # thermoelastic modulus
    fric.x.array[:] = 0.0

    # --. prescribed temperature field --..
    VT = functionspace(msh, ("Lagrange", 1))
    T = Function(VT, name="T")
    T.interpolate(lambda x: temperature(x, 1.0))

    # --. bulk thermoelastic residual (increment form) --..
    sigma = sigma_func(mu0, lmbda0)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    k_found = K_FOUND_FACTOR * E_f / R_PELLET ** 2
    F = ufl.inner(sigma(u), epsilon(v)) * dx
    F -= beta0 * (T - T_REF) * ufl.div(v) * dx
    F += k_found * ufl.inner(u, v) * dx          # weak foundation (regularisation)
    F = ufl.replace(F, {u: u + du})
    J = ufl.derivative(F, du, w)
    F_form, J_form = form(F), form(J)

    bcs = []   # no Dirichlet: rigid modes handled by the per-subdomain nullspace

    # --. contact: fuel outer arcs (10) <-> cladding inner (60), unbiased --..
    surfaces = adjacencylist(np.array([TAG_OUTER, TAG_CLAD_INNER], dtype=np.int32),
                             np.array([0, 2], dtype=np.int32))
    contact_pairs = [(0, 1), (1, 0)]
    search = [dolfinx_contact.cpp.ContactMode.ClosestPoint for _ in contact_pairs]
    cp = ContactProblem([facet_tags], surfaces, contact_pairs, msh, Q_DEGREE, search)
    cp.generate_contact_data(FrictionLaw.Frictionless, V,
                             {"u": u, "du": du, "mu": mu0, "lambda": lmbda0,
                              "fric": fric},
                             E_f * GAMMA * ORDER ** 2, THETA)

    # --. Newton callbacks --..
    def compute_coefficients(x, coeffs):
        nloc = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
        du.x.array[:nloc] = x.array_r[:nloc]
        du.x.scatter_forward()
        cp.update_contact_data(du)

    def compute_residual(x, b, coeffs):
        b.zeroEntries()
        b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        cp.assemble_vector(b, V)
        assemble_vector(b, F_form)
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    def compute_jacobian(x, A, coeffs):
        A.zeroEntries()
        cp.assemble_matrix(A, V)
        assemble_matrix(A, J_form, bcs=bcs)
        A.assemble()

    A = cp.create_matrix(J_form)
    b = create_vector(extract_function_spaces(F_form))

    newton = NewtonSolver(msh.comm, A, b, cp.coeffs)
    newton.set_residual(compute_residual)
    newton.set_jacobian(compute_jacobian)
    newton.set_coefficients(compute_coefficients)

    tags = np.unique(domain_marker.values)
    null_space = rigid_motions_nullspace_subdomains(V, domain_marker, tags)
    newton.A.setNearNullSpace(null_space)
    # contact linearisation floors the relative residual near 1e-4 (the absolute
    # residual still drops ~15 orders); a tight rtol is unreachable for contact.
    newton.set_newton_options({"atol": 1e-7, "rtol": 1e-4, "max_it": 30,
                               "relaxation_parameter": 1.0,
                               "convergence_criterion": "residual"})
    # direct solve: the foundation makes the tangent SPD and small enough for LU
    newton.set_krylov_options({"ksp_type": "preonly", "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps"})

    vtx = VTXWriter(msh.comm, str(HERE / "output" / "relocation_u.bp"), [u],
                    engine="BP4")

    # --. load stepping: ramp the thermal field --..
    for i in range(N_LOAD_STEPS):
        scale = (i + 1) / N_LOAD_STEPS
        T.interpolate(lambda x, s=scale: temperature(x, s))
        n_it, converged = newton.solve(du, write_solution=True)
        du.x.scatter_forward()
        u.x.array[:] += du.x.array[:]
        cp.update_contact_detection(u)
        A = cp.create_matrix(J_form)
        A.setNearNullSpace(null_space)
        newton.set_petsc_matrix(A)
        du.x.array[:] = 0.1 * du.x.array[:]
        cp.update_contact_data(du)
        vtx.write(float(i + 1))
        print(f"[relocation] step {i+1}/{N_LOAD_STEPS}  T_scale={scale:.2f}  "
              f"Newton_its={n_it}  converged={converged}")
    vtx.close()

    report_relocation(u, msh, facet_tags)
    plot_relocation(u, msh)


def report_relocation(u, msh, facet_tags):
    """Mean outward radial displacement of the fuel rim = relocation, and the
    fraction of the as-fabricated gap it closes."""
    tdim = msh.topology.dim
    V = u.function_space
    outer_facets = facet_tags.indices[facet_tags.values == TAG_OUTER]
    dofs = fem.locate_dofs_topological(V, tdim - 1, outer_facets)  # blocked dofs
    coords = V.tabulate_dof_coordinates()[dofs]
    uvals = u.x.array.reshape(-1, tdim)[dofs]
    r = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    rhat = coords[:, :2] / r[:, None]
    ur = np.sum(uvals[:, :2] * rhat, axis=1)
    print(f"[relocation] fuel rim: mean u_r = {ur.mean()*1e6:.2f} um, "
          f"max u_r = {ur.max()*1e6:.2f} um  (gap0 = {GAP0*1e6:.0f} um, "
          f"mean closure = {100*ur.mean()/GAP0:.1f}%)")


def plot_relocation(u, msh, scale=25.0):
    """Deformed fragments (displacement x scale) coloured by |u|, over the
    undeformed cladding, saved to relocation.png."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.tri import Triangulation
        import dolfinx.plot
    except Exception as e:
        print(f"[relocation] plot skipped ({e})")
        return
    tdim = msh.topology.dim
    V = u.function_space
    cells, _types, x = dolfinx.plot.vtk_mesh(V)
    tris = cells.reshape(-1, 4)[:, 1:]
    uu = u.x.array.reshape(-1, tdim)
    umag = np.sqrt(uu[:, 0] ** 2 + uu[:, 1] ** 2) * 1e6            # um
    xd = (x[:, 0] + scale * uu[:, 0]) * 1e3                         # mm
    yd = (x[:, 1] + scale * uu[:, 1]) * 1e3
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    tpc = ax.tripcolor(Triangulation(xd, yd, tris), umag,
                       shading="gouraud", cmap="magma")
    ax.triplot(Triangulation(x[:, 0] * 1e3, x[:, 1] * 1e3, tris),
               color="0.7", lw=0.15, alpha=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(f"Fragment relocation (displacement x{scale:.0f})")
    fig.colorbar(tpc, ax=ax, label="|u| (µm)", shrink=0.85)
    fig.tight_layout()
    fig.savefig(HERE / "relocation.png", dpi=140)
    print(f"[relocation] wrote {HERE / 'relocation.png'}")


if __name__ == "__main__":
    main()
