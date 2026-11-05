# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
python -m beam1d [case_dir]: solve the case in case_dir (default: the current
directory) and write output/fields.npz. The case files are described in
beam1d.case; the section is the undeformed hexagonal wall (hex_shell) of geometry.yaml
(beam1d.section), the one z3st meshes from the same numbers.

Unknowns: the centroid displacements u_x(z), u_y(z) (two Hermite bending
planes) and u_z(z) (linear rod). Fibre strain eps = u_z' - x u_x'' - y u_y'',
stress sigma = E (eps - eps*). With E uniform over the section (B1 = 0) the
three problems decouple (Ma et al. 2019, Eq. 17):

    (EA (u_z' - eps0*))' + q_z = 0
    (EI (u_a'' - k_a))''     = q_a ,     k_a = -kappa*_a ,  a = x, y

eps0* and kappa* are the section integrals of eps*; supports, loads and pads
come from boundary_conditions.yaml, the pads are solved by an active set
(beam1d.solver).
"""

import os
import sys

import numpy as np

from beam1d import case as case_mod
from beam1d import solver, thermal

PLANES = (("x", 0), ("y", 1))


def _bcs(c, z):
    """Fixed dofs and values per problem, point loads, pads, distributed load per unit length."""
    fix = {p: ([], []) for p in ("x", "y", "z")}
    for bc in case_mod.mechanical(c, "Dirichlet"):
        i = solver.node_at(z, case_mod.station(c, bc["region"]))
        for p, val in zip("xyz", bc.get("displacement") or []):
            if val is not None:
                fix[p][0].append(i if p == "z" else 2 * i)
                fix[p][1].append(float(val))
        for p, val in zip("xy", bc.get("rotation") or []):
            if val is not None:
                fix[p][0].append(2 * i + 1)
                fix[p][1].append(float(val))
    loads = {p: [] for p in ("x", "y", "z")}
    for bc in case_mod.mechanical(c, "Neumann"):
        i = solver.node_at(z, case_mod.station(c, bc["region"]))
        for p, F in zip("xyz", bc.get("force") or []):
            loads[p].append((i if p == "z" else 2 * i, float(F)))
        for p, m in zip("xy", bc.get("moment") or []):
            loads[p].append((2 * i + 1, float(m)))
    pads = {"x": [], "y": []}
    for bc in case_mod.mechanical(c, "Contact"):
        zc = case_mod.station(c, bc["region"])
        pads[bc["direction"]].append({"name": bc["region"], "z": zc, "dof": 2 * solver.node_at(z, zc),
                                      "gap": float(bc["gap"]), "side": int(bc["side"])})
    g = np.zeros(3)
    for bc in case_mod.mechanical(c, "Body"):
        g += np.asarray(bc.get("gravity") or [0.0, 0.0, 0.0], dtype=float)
    if np.any(g) and "rho" not in c.mat:
        raise KeyError(f"gravity needs rho in the material card {c.mat_path}")
    q = float(c.mat.get("rho", 0.0)) * c.section.A * g
    return fix, loads, pads, q


def run(case_dir=".", section_model=None, contacts=True):
    """Solve the case; returns the case and the fields written to output/fields.npz.

    section_model overrides the section model of a T table (contour_linear |
    ma12); contacts=False drops the pads (the free bow of a restrained case).
    """
    c = case_mod.load_case(case_dir)
    if not contacts:
        mech = c.bcs.get("mechanical", {})
        mech[c.mat_name] = [bc for bc in mech.get(c.mat_name) or [] if bc.get("type") != "Contact"]
    sec = c.section
    B = sec.modulus_weighted(c.mat["E"])   # stops on an E varying over the section
    E = float(c.mat["E"])
    EI, EA = E * sec.I, E * sec.A
    models = c.inp.get("models", {"bending": True, "axial": True})
    T, eps = thermal.build(c, sec, section_model)

    z = solver.mesh(c.H, int(c.geo["nz"]), list(c.stations.values()))
    z_g, xi, wg = solver.gauss_points(z, c.inp.get("element", {}).get("n_gauss", 3))
    kap_g = sec.curvature(eps, z_g.ravel()).reshape(*z_g.shape, 2)
    e0_g = sec.mean(eps, z_g.ravel()).reshape(z_g.shape)
    fix, loads, pads, q = _bcs(c, z)
    max_iters = c.inp.get("contact", {}).get("max_iters", 20)

    zs = np.unique(np.concatenate([np.linspace(0.0, c.H, 801), z]))
    z_mid = 0.5 * (z[:-1] + z[1:])
    kap_s, e0_s = sec.curvature(eps, zs), sec.mean(eps, zs)
    cf = sec.closed_form()
    out = {"z": z, "z_s": zs, "z_mid": z_mid, "E": E, "EI": EI, "EA": EA,
           "A_w": sec.A, "I_w": sec.I, "S_x": sec.S_x, "S_y": sec.S_y, "I_xx": sec.I_xx, "I_yy": sec.I_yy,
           "I_xy": sec.I_xy, "A_w_cf": cf["A"], "I_w_cf": cf["I"], "r": sec.radius_of_gyration(),
           "B0": B["B0"], "B1": B["B1"], "B2": B["B2"], "e0_s": e0_s, "sec_xy": sec.xy, "sec_dA": sec.dA}

    d2u_mid = {}
    pad_rows = []
    for a, ia in PLANES:
        k_s = -kap_s[:, ia]
        if not models.get("bending", True):
            d2u_mid[a] = np.zeros_like(z_mid)
            continue
        K, f = solver.assemble(z, EI, -kap_g[..., ia], xi, wg, q[ia])
        for dof, F in loads[a]:
            f[dof] += F
        U, R, F_pad, active, it = solver.solve_contact(K, f, *fix[a], pads[a], max_iters)
        u, du, d2u = solver.evaluate(z, U, zs)
        d2u_mid[a] = solver.evaluate(z, U, z_mid)[2]
        out.update({f"U_{a}": U, f"R_{a}": R, f"fixed_{a}": np.asarray(fix[a][0], dtype=int),
                    f"u{a}_s": u, f"du{a}_s": du, f"d2u{a}_s": d2u, f"k{a}_s": k_s,
                    f"M{a}_s": EI * (d2u - k_s), f"contact_iters_{a}": it})
        out["K"] = K   # unconstrained, the same for both planes
        for p, F, on in zip(pads[a], F_pad, active):
            pad_rows.append((p["name"], a, p["z"], p["gap"], p["side"], F, on, U[p["dof"]]))

    duz_mid = np.zeros_like(z_mid)
    if models.get("axial", True):
        K_r, f_r = solver.assemble_rod(z, EA, e0_g, xi, wg, q[2])
        for dof, F in loads["z"]:
            f_r[dof] += F
        U_z, R_z = solver.solve(K_r, f_r, *fix["z"])
        uz, duz = solver.evaluate_rod(z, U_z, zs)
        duz_mid = solver.evaluate_rod(z, U_z, z_mid)[1]
        out.update({"U_z": U_z, "R_z": R_z, "fixed_z": np.asarray(fix["z"][0], dtype=int),
                    "uz_s": uz, "duz_s": duz, "K_rod": K_r,
                    "N_mid": EA * (duz_mid - sec.mean(eps, z_mid))})

    # fibre stress at the element midpoints and its resultants
    sigma = sec.fibre_stress(E, sec.field(eps, z_mid), duz_mid, d2u_mid["x"], d2u_mid["y"])
    N_sig, M_sig = sec.resultants(sigma)
    out.update({"sigma_mid": sigma, "N_sig_mid": N_sig, "Mx_sig_mid": M_sig[0], "My_sig_mid": M_sig[1],
                "q": q})
    if pad_rows:
        cols = list(zip(*pad_rows))
        out.update({"pad_name": np.array(cols[0]), "pad_dir": np.array(cols[1]), "pad_z": np.array(cols[2]),
                    "pad_gap": np.array(cols[3]), "pad_side": np.array(cols[4]), "pad_force": np.array(cols[5]),
                    "pad_active": np.array(cols[6]), "pad_u": np.array(cols[7])})
    return c, out


def main(case_dir="."):
    c, out = run(case_dir)
    sec = c.section
    print(f"[INFO] Case      : {c.dir}")
    print(f"[INFO] Section   : a_o = {sec.a_o * 1e3:.4f} mm, a_i = {sec.a_i * 1e3:.4f} mm, "
          f"orientation = {sec.orientation}")
    print(f"[INFO]             A_w = {sec.A:.6e} m^2 (closed form {out['A_w_cf']:.6e}), "
          f"I_w = {sec.I:.6e} m^4 (closed form {out['I_w_cf']:.6e}), I_xy = {sec.I_xy:.1e}, "
          f"S = ({sec.S_x:.1e}, {sec.S_y:.1e}) m^3, r = {out['r'] * 1e3:.3f} mm")
    print(f"[INFO] Material  : {c.mat_path}")
    print(f"[INFO]             E = {out['E']:.3e} Pa, alpha = {c.mat['alpha']}, T_ref = {float(c.mat['T_ref']):.2f} K, "
          f"EI = {out['EI']:.4e} N m^2, EA = {out['EA']:.4e} N")
    z = out["z"]
    print(f"[INFO] Mesh      : {len(z) - 1} elements on [0, {c.H}] m, h_max = {np.diff(z).max():.4f} m, "
          f"stations {c.stations}")
    for a in ("x", "y", "z"):
        if f"U_{a}" not in out:
            continue
        U = out[f"U_{a}"]
        tip = f"u_{a}(H) = {U[-2 if a != 'z' else -1] * 1e3:.6f} mm"
        if a != "z":
            tip += f", u_{a}'(H) = {U[-1] * 1e3:.6f} mrad"
        print(f"[INFO] Problem {a} : {tip}, fixed dofs {out[f'fixed_{a}'].tolist()}, "
              f"reactions {np.array2string(out[f'R_{a}'], precision=4)} (N, N m)")
    for row in zip(*(out[k] for k in ("pad_name", "pad_dir", "pad_z", "pad_gap", "pad_side", "pad_force",
                                       "pad_active", "pad_u"))) if "pad_name" in out else ():
        name, a, zc, gap, side, F, on, u = row
        print(f"[INFO] Pad       : {name} at z = {zc} m, u_{a} = {u * 1e3:+.6f} mm, gap {side:+d} x {gap * 1e3:.3f} mm, "
              f"{'closed' if on else 'open'}, force {F:+.4f} N")
    os.makedirs(os.path.join(c.dir, "output"), exist_ok=True)
    path = os.path.join(c.dir, "output", "fields.npz")
    np.savez(path, **out)
    print(f"[INFO] Fields    : {path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
