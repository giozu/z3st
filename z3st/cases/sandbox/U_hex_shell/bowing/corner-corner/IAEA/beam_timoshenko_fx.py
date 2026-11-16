#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 / VP2 as a 1D beam solved by FEM (dolfinx), against the closed form
of hex_post (VP1, VP1Tdep) / vp2_post (VP2) and z3st 3D: the section-mean
profiles of the fields file of post.yaml when it exists (hex_post.analyse),
and the stations of the non-regression gold.

Inputs are those of the case folder: post.yaml (material card, closed_form
module), the temperature module (T1, T2, T_vpN, PHI_HOT, D) and its
geometry.yaml.

Section sums on the mid-wall line (thin wall, dA = t dp), at every z, with
c the coordinate towards the cold corner from the section centre:

    EA = int E dA            S  = int E c dA          N* = int E eps_th dA
    EI = int E c^2 dA        M* = - int E eps_th c dA
    GA = K_SHEAR int E / (2 (1 + nu)) dA

The reference axis is the section centre, as the section-mean displacement
of z3st and the closed forms; with E(T) the E-weighted centroid moves off it
and S couples stretching and bending (S = 0 for a constant card).

eps_th = alpha (T - T_ref) for a constant card, the card's thermal_strain
hook otherwise, E constant or the card's E hook.

Timoshenko beam with stretching, mixed (u, w, theta), P2 / P2 / P2, with
eps = u' - c theta' on the section:

    int EA u' du' - S (theta' du' + u' dtheta') + EI theta' dtheta'
        + GA (w' - theta)(dw' - dtheta) dz = int N* du' + M* dtheta' dz

Clamped at z = 0 (u = w = theta = 0), free at the top. Without lateral
loads the shear force is zero and the solution is the Euler-Bernoulli one,
the closed form. Lateral loads, pad springs
and other supports enter as extra terms of the same forms.

--element hermite: Euler-Bernoulli with stretching, C1, (u, w) in
P2 / cubic Hermite (Basix), theta = w':

    int EA u' du' - S (w'' du' + u' dw'') + EI w'' dw'' dz = int N* du' + M* dw'' dz

Clamped by the value and the slope dofs of w at z = 0. dolfinx does not scale
the Hermite slope dofs by the cell size: C1 only on a uniform mesh (checked).

Usage (z3st conda env)::

    python beam_timoshenko_fx.py [VP1 VP1/1515Ti_Tdependent VP2 ...] [--n 400] [--element mixed|hermite]

Outputs in ``U_hex_shell/results/``: <case>_beam_timoshenko_fx[_hermite].png and _summary.json.
"""

import argparse
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 3))
sys.path.insert(0, U_HEX)
sys.path.insert(0, os.path.join(HERE, "VP2"))
import hex_post as hp  # noqa: E402
import vp2_post  # noqa: E402

K_SHEAR = 0.5          # shear coefficient of a thin-walled closed tube
N_P = 6 * 400          # points on the mid-wall line (400 a flat)
N_ZK = 4001            # z levels of the section sums
GOLD_STATIONS = {"uy_core_top": 2.5, "uy_ACLP": 3.0, "uy_duct_top": 4.0}


def _resolve(case_dir):
    """post.yaml, material card, temperature module and closed form of a case folder."""
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    mat = hp.load_material(os.path.join(case_dir, cfg["material"]))
    kind, mod_path = str(cfg["closed_form"]).split(":", 1)
    mod_dir, module = os.path.split(os.path.normpath(os.path.join(case_dir, mod_path)))
    m = hp._import(module, mod_dir)
    fr = geometry_frame(mod_dir)
    cf = vp2_post.closed_form({**cfg, "closed_form": f"{kind}:{module}"}, mat, mod_dir, fr)
    return cfg, mat, m, getattr(m, f"T_{kind}"), fr, cf


def geometry_frame(mod_dir):
    """The hex_post frame entries used by perimeter_point and the closed forms, from geometry.yaml."""
    with open(os.path.join(mod_dir, "geometry.yaml")) as fh:
        g = yaml.safe_load(fh)
    D, t, H = float(g["D"]), float(g["t"]), float(g["H"])
    a_mid = 0.5 * D * np.cos(np.pi / 6) - 0.5 * t
    a0 = np.pi / 6 if g.get("orientation", "flat") == "vertical" else 0.0   # corners on +-y / +-x
    return {"a0": a0, "a_mid": a_mid, "L_mid": 2 * a_mid * np.tan(np.pi / 6), "t": t,
            "D": D, "z0": 0.0, "H": H}


def section_sums(mat, m, T_field, fr):
    """EA, EI, GA, N*, M* on N_ZK levels along z."""
    p = (np.arange(N_P) + 0.5) * 6 * fr["L_mid"] / N_P
    x, y, _ = hp.perimeter_point(fr, p, 0.0)
    dA = fr["t"] * 6 * fr["L_mid"] / N_P
    e = -np.array([np.cos(m.PHI_HOT), np.sin(m.PHI_HOT)])   # towards the cold corner
    c = x * e[0] + y * e[1]
    zk = np.linspace(fr["z0"], fr["z0"] + fr["H"], N_ZK)
    T = np.stack([T_field(np.array([x, y, np.full_like(x, z)])) for z in zk])   # (N_ZK, N_P)
    if "alpha" in mat:
        E = np.full_like(T, mat["E"])
        eps = mat["alpha"] * (T - mat["T_ref"])
    else:
        E = hp._hook(mat["E"])(T)
        eps = hp._hook(mat["eigenstrain"].rsplit(".", 1)[0] + ".thermal_strain")(
            T, mat["T_ref"], mat["expansion_law"] or "Sobolev")
    EA = E.sum(1) * dA
    return {"z": zk, "EA": EA, "S": (E * c).sum(1) * dA, "EI": (E * c**2).sum(1) * dA,
            "GA": K_SHEAR * EA / (2 * (1 + mat["nu"])),
            "N": (E * eps).sum(1) * dA, "M": -(E * eps * c).sum(1) * dA}


def solve(sec, H, n_el, element="mixed"):
    """1D FEM: u_z, w, theta sorted by z, at the P2 nodes (mixed) or at 2 n_el + 1
    uniform points (hermite)."""
    from mpi4py import MPI
    import basix
    import basix.ufl
    import ufl
    from dolfinx import fem, geometry, mesh
    from dolfinx.fem.petsc import LinearProblem
    from petsc4py import PETSc

    dom = mesh.create_interval(MPI.COMM_SELF, n_el, [0.0, H])
    Q1 = fem.functionspace(dom, ("Lagrange", 1))
    zq = Q1.tabulate_dof_coordinates()[:, 0]

    def coef(key):
        f = fem.Function(Q1)
        f.x.array[:] = np.interp(zq, sec["z"], sec[key])
        return f

    EA, S, EI, GA, Ns, Ms = (coef(k) for k in ("EA", "S", "EI", "GA", "N", "M"))
    clamp = mesh.locate_entities_boundary(dom, 0, lambda x: np.isclose(x[0], 0.0))
    opts = {"ksp_type": "preonly", "pc_type": "lu"}
    P2 = basix.ufl.element("Lagrange", "interval", 2)

    if element == "hermite":
        h = np.diff(np.sort(dom.geometry.x[:, 0]))
        if np.ptp(h) > 1e-12 * h.mean():
            raise ValueError("Hermite slope dofs are not scaled by dolfinx: uniform mesh only")
        HM = basix.ufl.element(basix.ElementFamily.Hermite, basix.CellType.interval, 3)
        W = fem.functionspace(dom, basix.ufl.mixed_element([P2, HM]))
        u, w = ufl.TrialFunctions(W)
        du, dw = ufl.TestFunctions(W)
        k, dk = w.dx(0).dx(0), dw.dx(0).dx(0)
        a = (EA * u.dx(0) * du.dx(0) - S * (k * du.dx(0) + u.dx(0) * dk) + EI * k * dk) * ufl.dx
        L = (Ns * du.dx(0) + Ms * dk) * ufl.dx
        bcs = [fem.dirichletbc(PETSc.ScalarType(0.0), fem.locate_dofs_topological(W.sub(i), 0, clamp), W.sub(i))
               for i in (0, 1)]   # u; value and slope of w
        sol = LinearProblem(a, L, bcs=bcs, petsc_options_prefix="beam_c1_", petsc_options=opts).solve()
        DG2 = fem.functionspace(dom, ("DG", 2))
        slope = fem.Function(DG2)
        slope.interpolate(fem.Expression(sol.sub(1).dx(0), DG2.element.interpolation_points))
        z = np.linspace(0.0, H, 2 * n_el + 1)
        p = np.stack([z, 0 * z, 0 * z], axis=1)
        tree = geometry.bb_tree(dom, 1)
        cells = geometry.compute_colliding_cells(dom, geometry.compute_collisions_points(tree, p), p)
        cells = np.array([cells.links(j)[0] for j in range(len(z))], dtype=np.int32)
        ev = lambda f: np.ravel(f.eval(p, cells)).copy()
        return {"u_z": (z, ev(sol.sub(0))), "w": (z, ev(sol.sub(1))), "theta": (z, ev(slope))}

    W = fem.functionspace(dom, basix.ufl.mixed_element([P2, P2, P2]))
    u, w, th = ufl.TrialFunctions(W)
    du, dw, dth = ufl.TestFunctions(W)
    a = (EA * u.dx(0) * du.dx(0) - S * (th.dx(0) * du.dx(0) + u.dx(0) * dth.dx(0))
         + EI * th.dx(0) * dth.dx(0) + GA * (w.dx(0) - th) * (dw.dx(0) - dth)) * ufl.dx
    L = (Ns * du.dx(0) + Ms * dth.dx(0)) * ufl.dx
    bcs = [fem.dirichletbc(PETSc.ScalarType(0.0), fem.locate_dofs_topological(W.sub(i), 0, clamp), W.sub(i))
           for i in (0, 1, 2)]
    sol = LinearProblem(a, L, bcs=bcs, petsc_options_prefix="beam_", petsc_options=opts).solve()

    out = {}
    for name, i in (("u_z", 0), ("w", 1), ("theta", 2)):
        f = sol.sub(i).collapse()
        z = f.function_space.tabulate_dof_coordinates()[:, 0]
        i = np.argsort(z)
        out[name] = (z[i], f.x.array[i].copy())
    return out


def gold_3d(case_dir, mat):
    """z3st 3D bow at the stations from the non-regression gold: output/, or
    <expansion_law>/output/ for a T-dependent case (its Allrun links output/ there)."""
    paths = [os.path.join(case_dir, d, "non-regression_gold.json")
             for d in ("output", os.path.join(str(mat.get("expansion_law")), "output"))]
    path = next((p for p in paths if os.path.exists(p)), None)
    if path is None:
        return {}
    with open(path) as fh:
        res = json.load(fh)["results"]
    return {GOLD_STATIONS[k]: v["numerical"] for k, v in res.items() if k in GOLD_STATIONS}


def fields_3d(case_dir, cfg, mat, cf):
    """Section-mean bow, u_z and curvature of the z3st 3D run (hex_post.analyse), None without the fields file."""
    path = os.path.join(case_dir, cfg["fields"])
    if not os.path.exists(path):
        return None
    d = hp.load_z3st(path, temperature=cfg["temperature"], case_dir=case_dir)
    res = hp.analyse(d, mat, cfg.get("support", "clamped"), cf)
    return {"z": res["z_nodes"], "w": res["w"], "u_z": res["u_bar"][:, 2],
            "z_k": res["z_sec"], "kappa": res["kappa"]}


def run(case_dir, n_el, out_dir, element="mixed"):
    case = os.path.relpath(os.path.abspath(case_dir), HERE)
    tag = case.replace(os.sep, "_")
    cfg, mat, m, T_field, fr, cf = _resolve(case_dir)
    sec = section_sums(mat, m, T_field, fr)
    fem_ = solve(sec, fr["H"], n_el, element)
    label = {"mixed": "Timoshenko + stretching, C0 P2/P2/P2",
             "hermite": "Euler-Bernoulli + stretching, C1 Hermite / P2"}[element]
    suffix = "" if element == "mixed" else f"_{element}"
    zw, w = fem_["w"]
    zu, uz = fem_["u_z"]
    zt, th = fem_["theta"]
    kap = np.gradient(th, zt)
    w_cf, uz_cf = cf.w(zw), cf.u_z(zu)
    f3 = fields_3d(case_dir, cfg, mat, cf)
    g3d = gold_3d(case_dir, mat)

    at = lambda zq, z, v: 1e3 * float(np.interp(zq, z, v))
    stations = sorted(set(cfg.get("stations_z", [])) | set(g3d))
    rows = []
    for z in stations:
        r = {"z_m": z, "w_1d_mm": at(z, zw, w), "w_cf_mm": 1e3 * float(cf.w(np.array([z]))[0]),
             "w_3d_mm": at(z, f3["z"], f3["w"]) if f3 else (1e3 * g3d[z] if z in g3d else None),
             "w_3d_gold_mm": 1e3 * g3d[z] if z in g3d else None,
             "u_z_1d_mm": at(z, zu, uz), "u_z_cf_mm": 1e3 * float(cf.u_z(np.array([z]))[0]),
             "u_z_3d_mm": at(z, f3["z"], f3["u_z"]) if f3 else None}
        rows.append(r)
    err = {"w_max_abs_mm": 1e3 * np.max(np.abs(w - w_cf)), "u_z_max_abs_mm": 1e3 * np.max(np.abs(uz - uz_cf)),
           "w_tip_rel": abs(w[-1] / w_cf[-1] - 1)}
    if f3:
        w1_on3 = np.interp(f3["z"], zw, w)
        err.update({"w_1d_vs_3d_max_abs_mm": 1e3 * np.max(np.abs(w1_on3 - f3["w"])),
                    "w_1d_vs_3d_tip_rel": abs(f3["w"][-1] / w1_on3[-1] - 1),
                    "u_z_1d_vs_3d_max_abs_mm": 1e3 * np.max(np.abs(np.interp(f3["z"], zu, uz) - f3["u_z"]))})

    src3 = "fields" if f3 else ("gold" if g3d else None)
    print(f"\n{case}: 1D FEM ({label}, {n_el} elements) vs closed form ({type(cf).__name__})"
          f"{f' and z3st 3D ({src3})' if src3 else ''}")
    print(f"  EI = {sec['EI'][0]:.4e} N m^2 at z = 0, {sec['EI'][-1]:.4e} at the top;"
          f"  max |S| = {np.abs(sec['S']).max():.3e} N;  max M* = {sec['M'].max():.4e} N m")
    dash = f"{'-':>10}"
    print(f"  {'z (m)':>6} {'w 1D':>10} {'w closed':>10} {'w 3D':>10} {'u_z 1D':>10} {'u_z closed':>11} {'u_z 3D':>10}   (mm)")
    for r in rows:
        w3 = f"{r['w_3d_mm']:10.4f}" if r["w_3d_mm"] is not None else dash
        u3 = f"{r['u_z_3d_mm']:10.4f}" if r["u_z_3d_mm"] is not None else dash
        print(f"  {r['z_m']:6.2f} {r['w_1d_mm']:10.4f} {r['w_cf_mm']:10.4f} {w3} {r['u_z_1d_mm']:10.4f} {r['u_z_cf_mm']:11.4f} {u3}")
    print(f"  max |w 1D - closed| = {err['w_max_abs_mm']:.2e} mm, max |u_z 1D - closed| = {err['u_z_max_abs_mm']:.2e} mm")
    if f3:
        print(f"  max |w 1D - 3D| = {err['w_1d_vs_3d_max_abs_mm']:.4f} mm (tip {100 * err['w_1d_vs_3d_tip_rel']:.2f} %),"
              f" max |u_z 1D - 3D| = {err['u_z_1d_vs_3d_max_abs_mm']:.4f} mm")

    fig, axs = plt.subplots(1, 4, figsize=(18, 6.5), sharey=True)
    ax = axs[0]
    ax.plot(1e3 * w, zw, "-", color="tab:blue", lw=2, label="1D FEM")
    ax.plot(1e3 * w_cf, zw, "--", color="k", lw=1, label=f"closed form ({type(cf).__name__})")
    if f3:
        ax.plot(1e3 * f3["w"], f3["z"], ":", color="tab:red", lw=2, label="z3st 3D, section mean")
    if g3d:
        ax.plot([1e3 * v for v in g3d.values()], list(g3d), "o", mfc="none", mec="tab:red", ms=8,
                label="z3st 3D, non-regression gold")
    ax.set(xlabel="bow w towards the cold corner (mm)", ylabel="z (m)", title="Bow")
    ax = axs[1]
    ax.plot(kap, zt, "-", color="tab:blue", lw=2, label="1D FEM, theta'")
    ax.plot(cf.kappa(zt), zt, "--", color="k", lw=1, label="closed form")
    if f3:
        ax.plot(f3["kappa"], f3["z_k"], ":", color="tab:red", lw=2, label="z3st 3D, plane-section fit")
    ax.set(xlabel="curvature (1/m)", title="Curvature")
    ax = axs[2]
    if f3:
        ax.plot(1e3 * (np.interp(f3["z"], zw, w) - f3["w"]), f3["z"], "-", color="tab:red", label="w: 1D - 3D")
        ax.plot(1e3 * (np.interp(f3["z"], zu, uz) - f3["u_z"]), f3["z"], "-", color="tab:orange", label="u_z: 1D - 3D")
    ax.plot(1e3 * (w - w_cf), zw, "-", color="tab:blue", label="w: 1D - closed")
    ax.plot(1e3 * (uz - uz_cf), zu, "-", color="tab:green", label="u_z: 1D - closed")
    ax.axvline(0, color="k", lw=0.8)
    ax.set(xlabel="difference (mm)", title="Differences")
    ax = axs[3]
    ax.plot(1e3 * uz, zu, "-", color="tab:green", lw=2, label="1D FEM")
    ax.plot(1e3 * uz_cf, zu, "--", color="k", lw=1, label="closed form")
    if f3:
        ax.plot(1e3 * f3["u_z"], f3["z"], ":", color="tab:red", lw=2, label="z3st 3D, section mean")
    ax.set(xlabel=f"axial u_z (mm), T_ref {mat['T_ref'] - 273.15:.0f} C", title="Axial growth")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(f"{case}: 1D beam FEM ({label}, {n_el} el.) vs closed form and z3st 3D")
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"{tag}_beam_timoshenko_fx{suffix}.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)
    with open(os.path.join(out_dir, f"{tag}_beam_timoshenko_fx{suffix}_summary.json"), "w") as fh:
        json.dump({"case": case, "element": element, "n_el": n_el, "closed_form": type(cf).__name__, "K_SHEAR": K_SHEAR,
                   "z3st_3d": src3, "stations": rows, "errors": err}, fh, indent=2, default=float)
    print(f"[INFO] {png}")
    return err


def main():
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cases", nargs="*", default=["VP1", os.path.join("VP1", "1515Ti_Tdependent"), "VP2"],
                    help="case folders, relative to this file")
    ap.add_argument("--n", type=int, default=400, help="axial elements (multiple of 8: nodes at 1.5 and 2.5 m)")
    ap.add_argument("--element", choices=("mixed", "hermite"), default="mixed",
                    help="mixed: C0 Timoshenko (any mesh); hermite: C1 Euler-Bernoulli (uniform mesh)")
    ap.add_argument("--out", default=hp.RESULTS)
    args = ap.parse_args()
    for c in args.cases:
        run(os.path.join(HERE, c), args.n, args.out, args.element)


if __name__ == "__main__":
    main()
