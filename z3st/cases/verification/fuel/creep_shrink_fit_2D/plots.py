#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import glob
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_extract_xdmf import extract_field_xdmf

pv.OFF_SCREEN = True
# start_xvfb is deprecated in recent PyVista but is the working headless path
# here; the deprecation notice is silenced.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        pv.start_xvfb()
    except Exception:
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
HISTORY_CSV = os.path.join(OUT, "history.csv")

# parameters
# Geometry, material and solver constants come from this case's YAML files.
# See case_params.py.
from case_params import (
    R_GAS, R_PELLET, R_CLAD_I, R_CLAD_O, E_EL, NU, E_FUEL, NU_FUEL,
    CREEP_A0, CREEP_N, CREEP_Q, elastic_factor, report,
)






TOL = 2.0e-5






# ----------------------------------------------------------------------
# 1. MESH PLOT
# ----------------------------------------------------------------------
def plot_mesh():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    grid = pv.read(files[0])
    # region by cell-centre radius: pellet vs cladding (gap is unmeshed)
    cc = grid.cell_centers().points[:, 0]
    region = np.where(cc <= R_PELLET + 1e-6, 0, 1)  # 0 pellet, 1 clad
    grid.cell_data["region"] = region

    p = pv.Plotter(off_screen=True, window_size=(700, 1100))
    p.add_mesh(grid, scalars="region", show_edges=True, line_width=1,
               cmap=["#4C72B0", "#DD8452"], show_scalar_bar=False)
    p.add_text("UO2 pellet | gap 65 um | Zircaloy clad",
               font_size=12, position=(0.04, 0.965), viewport=True, color="black")
    p.view_xy()
    p.camera.zoom(0.95)
    p.add_axes(xlabel="r", ylabel="z")
    p.screenshot(os.path.join(OUT, "mesh.png"))
    p.close()
    print("  wrote output/mesh.png")








# ----------------------------------------------------------------------
# 2b. RADIAL PROFILE at mid-height
# ----------------------------------------------------------------------
_HISTORY = None


def history():
    """history.csv as a named record array, read once per process."""
    global _HISTORY
    if _HISTORY is None:
        _HISTORY = np.genfromtxt(HISTORY_CSV, delimiter=",", names=True, dtype=float)
    return _HISTORY


def _step_days(n):
    """Solve-step times in days, read from history.csv (falls back to indices)."""
    if os.path.exists(HISTORY_CSV):
        s = np.atleast_1d(history()["time_s"]) / 86400.0
        if len(s) >= n:
            return s[:n]
    return np.arange(n, dtype=float)


def plot_radial_profile():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    z_mid = 0.005          # mid-height (m)
    dz = 3.0e-4            # band half-width to pick a horizontal cut

    days = _step_days(len(files))      # solve-step times, for the legend
    n = len(files)
    label_every = max(1, n // 6)       # 6 labelled curves, rest unlabelled
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    cmap = plt.cm.viridis
    for i, f in enumerate(files):
        g = pv.read(f)
        r = g.points[:, 0]
        z = g.points[:, 1]
        ur = g.point_data["Displacement"][:, 0] * 1e6  # um
        band = np.abs(z - z_mid) < dz
        rb, ub = r[band], ur[band]
        color = cmap(i / (len(files) - 1))
        for lo, hi in [(0.0, R_PELLET + 1e-9), (R_CLAD_I - 1e-9, R_CLAD_O + 1e-9)]:
            sel = (rb >= lo) & (rb <= hi)
            o = np.argsort(rb[sel])
            ax.plot(rb[sel][o] * 1e3, ub[sel][o], "-", color=color, lw=1.8)
        if i % label_every == 0 or i == n - 1:
            ax.plot([], [], "-", color=color, lw=1.8, label=f"{days[i]:.0f} d")

    # shade the gap region
    ax.axvspan(R_PELLET * 1e3, R_CLAD_I * 1e3, color="0.85", alpha=0.7)
    ax.text((R_PELLET + R_CLAD_I) / 2 * 1e3, ax.get_ylim()[1] * 0.05,
            "gap", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("radial displacement u_r (um)")
    ax.set_title("Radial displacement profile (mid-height)\n"
                 "clad creeps under the shrink-fit interference, contact pressure relaxes")
    ax.legend(title="time", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.text(0.97, 0.06, "fresh (0 MWd/kgU) curve is non-zero: thermal\n"
            "expansion from T_ref (293 K) to 580 K coolant",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right", color="#444",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "radial_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/radial_profile.png")


# ----------------------------------------------------------------------
# 2c. STRESS PROFILE (radial + hoop) at mid-height, last step
# ----------------------------------------------------------------------
def _stress_profile_data():
    """(r, sigma_rr, sigma_theta) in MPa on a mid-height band, last step.

    Prefers nodal VTU stress; falls back to the cell-wise XDMF field
    (piecewise constant, sampled at cell centres). The stress tensor is
    stored in cylindrical order (r, theta, z), so in the flattened
    9-component array column 0 is sigma_rr and column 4 is sigma_theta.
    """
    z_mid, dz = 0.005, 3.0e-4
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    if files:
        g = pv.read(files[-1])
        if "Stress (points)" not in g.point_data:
            return None
        r, z = g.points[:, 0], g.points[:, 1]
        s = np.asarray(g.point_data["Stress (points)"]).reshape(-1, 9)
    else:
        # the solver keeps fields.h5 locked while running; read-only access to
        # flushed steps is safe, so disable HDF5 locking before h5py loads
        os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        xdmf = os.path.join(OUT, "fields.xdmf")
        if not os.path.exists(xdmf):
            return None
        r, z, _, s = extract_field_xdmf(xdmf, "Stress", step_index=-1)
        s = np.asarray(s).reshape(len(r), 9)
    band = np.abs(z - z_mid) < dz
    return r[band], s[band, 0] / 1e6, s[band, 4] / 1e6


def plot_stress_profile():
    data = _stress_profile_data()
    if data is None:
        print("  [plots] no stress field in output: skipping stress profile")
        return
    r, s_rr, s_tt = data
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    first = True
    for lo, hi in [(0.0, R_PELLET + 1e-9), (R_CLAD_I - 1e-9, R_CLAD_O + 1e-9)]:
        sel = (r >= lo) & (r <= hi)
        o = np.argsort(r[sel])
        ax.plot(r[sel][o] * 1e3, s_rr[sel][o], "o-", color="#4C72B0", ms=3,
                lw=1.8, label=r"$\sigma_{rr}$ (radial)" if first else None)
        ax.plot(r[sel][o] * 1e3, s_tt[sel][o], "s-", color="#C44E52", ms=3,
                lw=1.8, label=r"$\sigma_{\theta\theta}$ (hoop)" if first else None)
        first = False

    ax.axvspan(R_PELLET * 1e3, R_CLAD_I * 1e3, color="0.85", alpha=0.7)
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.text((R_PELLET + R_CLAD_I) / 2 * 1e3, ax.get_ylim()[1] * 0.05,
            "gap", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title("Radial and hoop stress at mid-height (last step)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/stress_profile.png")

    # console summary at the four radii of interest
    def at(rt):
        i = np.argmin(np.abs(r - rt))
        return s_rr[i], s_tt[i]
    print("    mid-height stresses (MPa):")
    for name, rt in [("fuel centre", 0.0), ("fuel surface", R_PELLET),
                     ("clad inner", R_CLAD_I), ("clad outer", R_CLAD_O)]:
        srr, stt = at(rt)
        print(f"      {name:12s} r={rt * 1e3:6.3f} mm   "
              f"sigma_rr={srr:9.1f}   sigma_theta={stt:9.1f}")






# ----------------------------------------------------------------------
# 3. FINAL-STEP FIELDS: temperature, radial displacement, clad von Mises
# ----------------------------------------------------------------------
def _bar_args(title):
    return dict(title=title, vertical=True, position_x=0.80, position_y=0.20,
                height=0.6, width=0.10, title_font_size=18, label_font_size=15,
                n_labels=4, fmt="%.0f", color="black")


def plot_fields():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    grid = pv.read(files[-1])

    # 3D displacement vector for warping (u_r, u_z, 0)
    u2 = grid.point_data["Displacement"]
    u3 = np.zeros((u2.shape[0], 3))
    u3[:, 0] = u2[:, 0]
    u3[:, 1] = u2[:, 1]
    grid.point_data["u3"] = u3
    grid.point_data["u_r_um"] = u2[:, 0] * 1e6
    grid.point_data["T"] = grid.point_data["Temperature"]

    # cell-centre radius to split pellet / clad robustly
    cc = grid.cell_centers().points[:, 0]
    clad = grid.extract_cells(np.where(cc > R_CLAD_I - 1e-6)[0])

    # ---- (a) temperature ----
    p = pv.Plotter(off_screen=True, window_size=(560, 1050))
    p.add_mesh(grid, scalars="T", cmap="inferno", show_scalar_bar=True,
               scalar_bar_args=_bar_args("T (K)"))
    p.add_text("Temperature (centreline 1706 K -> coolant 580 K)",
               font_size=11, position=(0.03, 0.965), viewport=True, color="black")
    p.view_xy(); p.camera.zoom(0.92)
    p.screenshot(os.path.join(OUT, "field_temperature.png")); p.close()
    print("  wrote output/field_temperature.png")

    # ---- (b) radial displacement field u_r on the (undeformed) mesh ----
    # A warped mesh is misleading for contact: each body warps by its own
    # displacement, so the pellet (large u_r) visually overruns the clad.
    p = pv.Plotter(off_screen=True, window_size=(560, 1050))
    p.add_mesh(grid, scalars="u_r_um", cmap="viridis", show_edges=False,
               scalar_bar_args=_bar_args("u_r (um)"))
    p.add_text("Radial displacement (pellet expands into the clad)",
               font_size=11, position=(0.03, 0.965), viewport=True, color="black")
    p.view_xy(); p.camera.zoom(0.92)
    p.screenshot(os.path.join(OUT, "field_radial_disp.png")); p.close()
    print("  wrote output/field_radial_disp.png")

    # ---- (c) cladding von Mises ----
    vm_key = "VonMises (points)" if "VonMises (points)" in grid.point_data else None
    if vm_key is not None:
        vm = np.asarray(clad.point_data[vm_key]) / 1e6  # MPa
        clad.point_data["vM_MPa"] = vm
        finite = vm[np.isfinite(vm)]
        clim = [float(finite.min()), float(finite.max())]
        p = pv.Plotter(off_screen=True, window_size=(560, 1050))
        p.add_mesh(clad, scalars="vM_MPa", cmap="plasma", clim=clim,
                   show_edges=True, line_width=0.4,
                   scalar_bar_args=_bar_args("vM (MPa)"))
        p.add_text(f"Cladding von Mises ({clim[0]:.0f}-{clim[1]:.0f} MPa)",
                   font_size=11, position=(0.03, 0.965), viewport=True, color="black")
        p.view_xy(); p.camera.zoom(0.92)
        p.screenshot(os.path.join(OUT, "field_clad_vonmises.png")); p.close()
        print("  wrote output/field_clad_vonmises.png")


def load_history_interference():
    """Load gap and interference from history.csv."""
    if not os.path.exists(HISTORY_CSV):
        raise FileNotFoundError(f"{HISTORY_CSV} not found. Run the case first.")

    data = history()
    time_s = data["time_s"]
    time_h = time_s / 3600.0
    time_days = time_s / 86400.0
    gap_um = data["gap_um"]
    pressure = data["contact_pressure_MPa"]
    interference_um = np.maximum(0.0, -gap_um)  # negative gap -> interference
    # No truncation to a nominal step count: with time_adaptivity enabled the
    # number of steps actually written is a runtime outcome, not sum(n_steps).
    # The history is aligned against the XDMF step count below instead.
    return time_h, time_days, gap_um, interference_um, pressure

def contact_pressure(t, Delta, a, b, c, E1, nu1, E2, nu2, A, n, T):
    """
    Contact pressure Pk(t) of a shrink-fit joint under creep.

    Parameters
    ----------
    t : array_like
        Times [h].
    Delta : float or array_like
        Radial interference [mm]. Scalar or array.
    a, b, c : float
        Inner, common and outer radii [mm].
    E1, nu1 : float
        Young's modulus and Poisson's ratio of the shaft.
    E2, nu2 : float
        Young's modulus and Poisson's ratio of the hub.
    A, n : float
        Norton-law parameters (hub).

    Returns
    -------
    Pk : ndarray
        Contact pressure [MPa] at each time in t.
    f, k1, k2, phi : float
        Intermediate quantities.
    """
    t = np.asarray(t, dtype=float)
    Delta = np.asarray(Delta, dtype=float)
    # Ensure Young moduli are in MPa for the formula (E may be defined in Pa)
    E1 = float(E1) 
    E2 = float(E2) 
    Delta_arr = np.full_like(t, float(Delta[0]), dtype=float) if Delta.ndim == 0 else np.broadcast_to(Delta, t.shape).astype(float)
    Delta_safe = np.maximum(Delta_arr, 0)

    # --- Step 1: elastic factor f (eq. 4-5) ---
    # Pk = Delta_u_el / [ b/E1*((a^2+b^2)/(b^2-a^2)+nu1) + b/E2*((c^2+b^2)/(c^2-b^2)-nu2) ]
    f = elastic_factor(a, b, c, E1, nu1, E2, nu2)  # [MPa/mm]

    # --- Step 2: k1, shaft creep contribution (eq. 16) ---
    # k1 = 0 for a solid shaft (a = 0): near-hydrostatic compression
    if a == 0:
        k1 = 0.0
    else:
        k1 = (np.sqrt(3) / 2) * (
            (np.sqrt(3) / n) * (a**(2/n) * b**(2/n)) / (b**(2/n) - a**(2/n))
        ) ** n

    # --- Step 3: k2, hub creep contribution (eq. 20) ---
    k2=(np.sqrt(3) / 2) * (
            (np.sqrt(3) / n) * (c**(2/n) * b**(2/n)) / (c**(2/n) - b**(2/n))
        ) ** n

    # --- Step 4: phi, relaxation rate (eq. 21-22) ---
    phi = -(CREEP_A0 * np.exp(-CREEP_Q / (R_GAS * T)) / b) * f**n * (k2 - k1)
    contact_mask = Delta_arr > 1e-9          # seuil physique, pas 1e-12
    Pk = np.zeros_like(t)

    if not np.any(contact_mask):
        # No interference at any recorded time, so there is no assembled shrink
        # fit to relax and the analytical pressure is identically zero. idx0
        # and t_rel below are only meaningful once contact exists.
        return Pk / 1e6, f, k1, k2, phi

    idx0 = np.argmax(contact_mask)        # premier indice de contact
    t0 = t[idx0]                          # time at which creep starts
    t_rel = t[idx0:] - t0


    # --- Step 5: elastic interference Delta_u_el(t) (eq. 21) ---
    if n == 1:
        # n = 1 is the special case: no relaxation (see the paper's Discussion)
        Delta_u_el = Delta_arr.copy()
    else:
        # Avoid dividing by zero when Delta_safe == 0 (no interference)
        Delta_u_el = np.zeros_like(Delta_arr, dtype=float)


        numer = 1 + phi[idx0:] * (1 - n) * t_rel / (Delta_safe[idx0:] ** (1 - n))
        Delta_u_el[idx0:] = Delta_arr[idx0:] * (numer) ** (1 / (1 - n))

    # --- Step 6: contact pressure Pk(t) = Delta_u_el(t) * f (eq. 5/21) ---
    Pk = Delta_u_el * f
    return Pk/1e6, f, k1, k2, phi



def plot_contact_pressure_evolution():
    """Esposito eq. (21) against the Z3ST contact pressure over the relaxation."""
    _problems = report()
    if _problems:
        warnings.warn(
            "case is inconsistent with the Esposito comparison; see the warnings "
            "above. The analytical curve below is not a valid verification.",
            RuntimeWarning,
        )

    # --- Joint geometry [mm] ---
    a = 0.0      # shaft inner radius (a = 0 for a solid shaft)
    b = R_CLAD_I    # rayon commun (interface de contact)
    c = R_CLAD_O     # hub outer radius

    # --- Elastic properties ---
    E1, nu1 = E_FUEL, NU_FUEL   # shaft  [MPa, -]
    E2, nu2 = E_EL, NU   # hub  [MPa, -]

    # --- Creep properties (Norton law: eps_eq_dot = A * sigma_eq^n) ---
    # Hub properties (A2, n2). For a solid shaft (a = 0) k1 = 0 and the shaft
    # properties do not enter.
    n = CREEP_N                  # exponent used in the law (n1 = n2 assumed)
    A = CREEP_A0               # Norton parameter used in phi (hub)

    # --- Time / interference table read from history.csv ---
    # t_array : time [days]
    # Delta   : radial interference [mm] (scalar or array)


    t_array, time_days, gap_um, interference_um, pressure = load_history_interference()
    Delta = np.maximum(interference_um * 1e-6, 0.0)  # [mm] : 1 µm = 1e-6 m

    t_array=t_array*3600
    # =====================================================================
    # 2. CONTACT-PRESSURE CLOSED FORM
    # =====================================================================





    # =====================================================================
    # 3. EVALUATION AND OUTPUT
    # =====================================================================

    # Clad temperature per step, from the same CSV as everything else: the
    # closed form uses it only through exp(-Q/(R*T)), and every series below is
    # then aligned by construction.
    Temp = history()["T_max_K"]

    Pk, f, k1, k2, phi = contact_pressure(
        t_array, Delta, a, b, c, E1, nu1, E2, nu2, A, n, Temp
    )



    print(f"Elastic factor f     = {f:.5f} MPa/m")
    print(f"k1                   = {k1:.5e}")
    print(f"k2                   = {k2:.5e}")
    print()
    print("Table read from history.csv:")
    print(f"{'t [days]':>10} | {'gap [um]':>10} | {'interference [um]':>18}")
    print("-" * 45)
    step = max(1, len(t_array) // 20)
    for ti, gap_i, inter_i in zip(time_days[::step], gap_um[::step], interference_um[::step]):
        print(f"{ti:10.1f} | {gap_i:10.3f} | {inter_i:18.3f}")
    print()
    print(f"{'t [days]':>10} | {'Pk [MPa]':>10}")
    print("-" * 25)
    for ti, Pi in zip(time_days[::5], Pk[::5]):
        print(f"{ti:10.1f} | {Pi:10.3f}")

    # --- Figure: Pk = f(t) ---
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(time_days, Pk, "-", color="#4C72B0", lw=1.8,
            label=f"Esposito eq. (21) — Tresca hub, n = {n:g}")
    ax.plot(time_days, pressure, "o", ms=5, mfc="none", mec="#C44E52", mew=1.6,
            label="Z3ST — J2 (von Mises), penalty contact")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("contact pressure $P_k$ (MPa)")
    ax.set_title("Shrink-fit contact pressure relaxing under Norton creep\n"
                 "clad hub on a pellet shaft, isothermal at 580 K")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    output_path = os.path.join(OUT, "contact_pressure_evolution.png")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Figure written to {output_path}")


def plot_interference_pressure():
    """Contact pressure against the interference that produces it.

    p = f * Delta is linear while the joint stays elastic, so the slope of this
    curve is the elastic factor of Esposito eq. (4). Interference is the closed
    gap: gap < 0 in history.csv means the surfaces overlap.
    """
    d = history()
    pressure_mpa = np.asarray(d["contact_pressure_MPa"], dtype=float)
    gap_um = np.asarray(d["gap_um"], dtype=float)
    time_days = np.asarray(d["time_days"], dtype=float)
    interference = np.maximum(0.0, -gap_um * 1e-6)   # um -> m, positive part

    # Create plot
    fig, ax2 = plt.subplots(figsize=(6.4, 4.6))
    
    # Middle plot: Contact pressure vs interference
    ax2.plot(interference * 1e6, pressure_mpa, "b-s", lw=2, ms=4, label="p vs interference")
    ax2.set_xlabel("Interference (μm)")
    ax2.set_ylabel("Contact pressure (MPa)")
    ax2.set_title("Contact pressure as function of interference")
    ax2.grid(True, ls=":", alpha=0.6)
    ax2.legend()
    
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "interference_pressure.png"), dpi=150)
    plt.close(fig)
    print("[INFO] interference_pressure.png saved")
    
    # Print summary statistics
    contact_steps = pressure_mpa > 0
    if contact_steps.any():
        print(f"[INFO] Contact initiated at time: {time_days[contact_steps][0]:.1f} days")
        print(f"[INFO] Contact initiated at interference: {interference[contact_steps][0]*1e6:.3f} μm")
        print(f"[INFO] Final interference: {interference[-1]*1e6:.3f} μm")
        print(f"[INFO] Final contact pressure: {pressure_mpa[-1]:.2f} MPa")
    else:
        print("[INFO] No contact detected during simulation")


if __name__ == "__main__":
    have_vtu = bool(glob.glob(os.path.join(OUT, "*.vtu")))
    have_csv = os.path.exists(os.path.join(OUT, "history.csv"))
    if not have_vtu and not have_csv:
        raise SystemExit("No output found - run the case first (Allrun).")
    print("[plots] generating figures...")
    plot_stress_profile()        # VTU-preferred; falls back to the XDMF h5
    plot_interference_pressure()
    plot_contact_pressure_evolution()   # Esposito eq. (21) overlay
    if have_vtu:
        plot_mesh()
        plot_radial_profile()
        plot_fields()
    else:
        print("  [plots] XDMF-only run: skipping mesh / radial-profile / field "
              "renders (pyvista cannot read dolfinx XDMF). Open output/fields.xdmf "
              "in ParaView, or re-run with 'output: {format: vtu}' for those figures.")
    print("[plots] done.")
