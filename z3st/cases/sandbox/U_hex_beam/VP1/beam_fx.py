#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam_fx: IAEA VP1 Euler-Bernoulli beam in weak form (FEniCSx prototype)
# Version: 0.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
The beam of beam1d (same case files, section, T and eps*) solved with dolfinx
on the 1D gmsh mesh mesh.msh (mesh.geo): nz elements on [0, H] along x, the
stations of geometry.yaml as named physical points (the regions of
boundary_conditions.yaml). Writes output/fields_fx.npz and compares it with
output/fields.npz of the Hermite solver.

Bending plane a = x, y: Lagrange elements are C0, so u'' is not available and
the plane is solved in mixed form, unknowns u (m) and the bending moment
M = EI (u'' - k) (N m), k = -kappa*_a (1/m), with M'' = q (N/m):

    int M tau / EI + int u' tau' + int k tau = sum_ends theta tau n
    int M' v' + int q v + sum_i F_i v(z_i)  = 0

u Dirichlet in V (v = 0 there). A prescribed rotation theta at an end is
natural (right-hand side, n = -1 at z = 0, +1 at z = H); at an end with free
rotation M = 0 is imposed in Q (tau = 0 there). Point forces F_i (N) act at
the stations. With q = F = 0 and a free end, M = 0 and u'' = k: the free bow.

Axial rod, unknown u_z (m), eps0* (-), EA (N):

    int EA (u_z' - eps0*) w' = int q_z w + sum_i F_z,i w(z_i)

Not implemented here: pads (Contact), moments (Neumann moment), a rotation
prescribed away from the ends.

Run from VP1 with U_hex_beam on PYTHONPATH (not part of Allrun):
    gmsh $(python -m z3st.utils.geo_args) mesh.geo -1 -o mesh.msh > log_mesh.md
    python beam_fx.py [degree]
"""

import os
import sys

import basix.ufl
import numpy as np
import ufl
from dolfinx import fem, geometry
from dolfinx.fem import petsc as fpetsc
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
from petsc4py import PETSc

from beam1d import case as case_mod
from beam1d import thermal

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")


def load_mesh(c, path):
    """Mesh and vertex tags of mesh.msh, checked against nz and the stations of geometry.yaml."""
    d = gmshio.read_from_msh(path, MPI.COMM_SELF, rank=0, gdim=1)
    msh, ft = d.mesh, d.facet_tags
    n_el = msh.topology.index_map(1).size_local
    if n_el != int(c.geo["nz"]):
        raise ValueError(f"{path}: {n_el} elements, nz = {c.geo['nz']} in geometry.yaml (re-run gmsh)")
    tags = {name: g.tag for name, g in d.physical_groups.items() if g.dim == 0}
    x = msh.geometry.x[:, 0]
    for name, zc in c.stations.items():
        v = ft.find(tags[name]) if name in tags else []
        if len(v) != 1 or not np.isclose(x[v[0]], zc):
            raise ValueError(f"{path}: station '{name}' at z = {zc} m is not a tagged vertex (check mesh.geo)")
    msh.topology.create_connectivity(0, 1)
    return msh, ft, tags


def evaluate(expr, V, zq):
    """expr (UFL) interpolated into V and evaluated at the heights zq."""
    f = fem.Function(V)
    f.interpolate(fem.Expression(expr, V.element.interpolation_points))
    msh = V.mesh
    pts = np.column_stack([zq, np.zeros_like(zq), np.zeros_like(zq)])
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cells = geometry.compute_colliding_cells(msh, geometry.compute_collisions_points(tree, pts), pts)
    return f.eval(pts, np.array([cells.links(i)[0] for i in range(len(zq))])).ravel()


def solve_linear(a, L, bcs, loads):
    """Solve a = L with the Dirichlet bcs and point loads [(dof, value)] added to the vector (MUMPS)."""
    a, L = fem.form(a), fem.form(L)
    A = fpetsc.assemble_matrix(a, bcs)
    A.assemble()
    b = fpetsc.assemble_vector(L)
    fpetsc.apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    for dof, F in loads:
        b.setValue(dof, F, addv=PETSc.InsertMode.ADD_VALUES)
    b.assemble()
    fpetsc.set_bc(b, bcs)
    ksp = PETSc.KSP().create(MPI.COMM_SELF)
    ksp.setOperators(A)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.getPC().setFactorSolverType("mumps")
    x = A.createVecRight()
    ksp.solve(b, x)
    if ksp.getConvergedReason() < 0:
        raise RuntimeError(f"MUMPS failed (reason {ksp.getConvergedReason()}): a mechanism in the supports?")
    return x


ENDS = ("bottom", "top")


def parse_bcs(c):
    """Per component, keyed by region: fixed displacements, end rotations, point forces; line load q."""
    fix = {p: {} for p in "xyz"}
    rot = {p: {} for p in "xy"}
    force = {p: {} for p in "xyz"}
    for bc in case_mod.mechanical(c, "Dirichlet"):
        r = bc["region"]
        case_mod.station(c, r)
        for p, val in zip("xyz", bc.get("displacement") or []):
            if val is not None:
                fix[p][r] = float(val)
        for p, val in zip("xy", bc.get("rotation") or []):
            if val is not None:
                if r not in ENDS:
                    raise NotImplementedError(f"rotation at '{r}': only at the ends in the mixed form")
                rot[p][r] = float(val)
    for bc in case_mod.mechanical(c, "Neumann"):
        if bc.get("moment"):
            raise NotImplementedError("Neumann moment: not in beam_fx")
        r = bc["region"]
        case_mod.station(c, r)
        for p, F in zip("xyz", bc.get("force") or []):
            force[p][r] = force[p].get(r, 0.0) + float(F)
    if case_mod.mechanical(c, "Contact"):
        raise NotImplementedError("Contact pads: not in beam_fx")
    g = np.zeros(3)
    for bc in case_mod.mechanical(c, "Body"):
        g += np.asarray(bc.get("gravity") or [0.0, 0.0, 0.0], dtype=float)
    q = float(c.mat.get("rho", 0.0)) * c.section.A * g
    return fix, rot, force, q


def bending_plane(msh, ft, tags, EI, k, q, fix, rot, force, degree):
    """u, M of one bending plane in mixed form; returns the solution and its space."""
    el = basix.ufl.element("Lagrange", "interval", degree)
    W = fem.functionspace(msh, basix.ufl.mixed_element([el, el]))
    (u, M), (v, tau) = ufl.TrialFunctions(W), ufl.TestFunctions(W)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    n = ufl.FacetNormal(msh)[0]

    a = (M * tau / EI + u.dx(0) * tau.dx(0) + M.dx(0) * v.dx(0)) * ufl.dx
    L = -(k * tau + fem.Constant(msh, q) * v) * ufl.dx
    for r, val in rot.items():
        L += fem.Constant(msh, val) * tau * n * ds(tags[r])

    Wu, Wm = W.sub(0), W.sub(1)
    dofs = lambda S, r: fem.locate_dofs_topological(S, 0, ft.find(tags[r]))   # noqa: E731
    bcs = [fem.dirichletbc(PETSc.ScalarType(val), dofs(Wu, r), Wu) for r, val in fix.items()]
    bcs += [fem.dirichletbc(PETSc.ScalarType(0.0), dofs(Wm, r), Wm) for r in ENDS if r not in rot]
    loads = [(dofs(Wu, r)[0], -F) for r, F in force.items()]

    w = fem.Function(W)
    w.x.array[:] = solve_linear(a, L, bcs, loads).array
    return w, W


def axial_rod(msh, ft, tags, EA, e0, q, fix, force, degree):
    """u_z of the axial rod."""
    V = fem.functionspace(msh, ("Lagrange", degree))
    u, w = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = EA * u.dx(0) * w.dx(0) * ufl.dx
    L = (EA * e0 * w.dx(0) + fem.Constant(msh, q) * w) * ufl.dx
    dofs = lambda r: fem.locate_dofs_topological(V, 0, ft.find(tags[r]))   # noqa: E731
    bcs = [fem.dirichletbc(PETSc.ScalarType(val), dofs(r), V) for r, val in fix.items()]
    loads = [(dofs(r)[0], F) for r, F in force.items()]
    uz = fem.Function(V)
    uz.x.array[:] = solve_linear(a, L, bcs, loads).array
    return uz


def run(case_dir=CASE_DIR, degree=3):
    """Solve the case on case_dir/mesh.msh; returns the case and the fields at the heights z_s."""
    c = case_mod.load_case(case_dir)
    sec = c.section
    sec.modulus_weighted(c.mat["E"])   # stops on an E varying over the section
    E = float(c.mat["E"])
    EI, EA = E * sec.I, E * sec.A
    _, eps = thermal.build(c, sec)
    fix, rot, force, q = parse_bcs(c)

    msh, ft, tags = load_mesh(c, os.path.join(c.dir, "mesh.msh"))
    z = np.sort(msh.geometry.x[:, 0])
    zs = np.unique(np.concatenate([np.linspace(0.0, c.H, 801), z]))
    out = {"z": z, "z_s": zs, "EI": EI, "EA": EA, "degree": degree}

    # free strains on DG of the element degree: exact when kappa*, eps0* are polynomials
    # of degree <= p per element (T kinks on nodes, constant alpha, as VP1)
    D = fem.functionspace(msh, ("DG", degree))
    for a, ia in (("x", 0), ("y", 1)):
        k = fem.Function(D)
        k.interpolate(lambda x, ia=ia: -sec.curvature(eps, x[0])[:, ia])
        w, W = bending_plane(msh, ft, tags, EI, k, q[ia], fix[a], rot[a], force[a], degree)
        u, M = ufl.split(w)
        Vu, VM = W.sub(0).collapse()[0], W.sub(1).collapse()[0]
        Vd = fem.functionspace(msh, ("DG", degree - 1))
        out.update({f"u{a}_s": evaluate(u, Vu, zs), f"du{a}_s": evaluate(u.dx(0), Vd, zs),
                    f"M{a}_s": evaluate(M, VM, zs)})

    e0 = fem.Function(D)
    e0.interpolate(lambda x: sec.mean(eps, x[0]))
    uz = axial_rod(msh, ft, tags, EA, e0, q[2], fix["z"], force["z"], degree)
    out.update({"uz_s": evaluate(uz, uz.function_space, zs),
                "N_s": evaluate(EA * (uz.dx(0) - e0), D, zs)})
    return c, out


def main(degree=3):
    c, out = run(CASE_DIR, degree)
    zs = out["z_s"]
    print(f"[INFO] Case      : {c.dir}")
    print(f"[INFO] Mesh      : {len(out['z']) - 1} elements, P{degree}, stations {c.stations}")
    for name in ("core_top", "LRP", "top"):
        i = int(np.argmin(np.abs(zs - c.stations[name])))
        print(f"[INFO] {name:9s} : u_y = {out['uy_s'][i] * 1e3:.6f} mm, u_x = {out['ux_s'][i] * 1e3:.2e} mm, "
              f"u_z = {out['uz_s'][i] * 1e3:.6f} mm")
    print(f"[INFO] Tip       : u_y'(H) = {out['duy_s'][-1] * 1e3:.6f} mrad")
    print(f"[INFO] Resultants: max |M_x|, |M_y| = {np.abs(out['Mx_s']).max():.2e}, "
          f"{np.abs(out['My_s']).max():.2e} N m, max |N| = {np.abs(out['N_s']).max():.2e} N")

    herm = os.path.join(OUT_DIR, "fields.npz")
    if os.path.exists(herm):
        d = np.load(herm)
        for key in ("ux_s", "uy_s", "uz_s", "duy_s"):
            ref = np.interp(zs, d["z_s"], d[key])
            print(f"[INFO] vs beam1d : max |{key} fx - Hermite| = {np.abs(out[key] - ref).max():.3e}")

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "fields_fx.npz")
    np.savez(path, **out)
    print(f"[INFO] Fields    : {path}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
