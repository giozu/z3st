#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST parametric study --.. ..- .-.. .-.. ---
"""
GB-fracture parametric study for two_elliptical_cavities_2D.

Physics
-------
Fission-gas bubbles sitting on a grain boundary (GB) in oxide fuel do two
things: they concentrate stress at their tips and they reduce the load-bearing
GB ligament. This script predicts the critical bubble pressure ``p_crit`` to
fail the GB, as a function of the bubble areal coverage ``Fc`` and the GB
toughness ``Gc`` — the actual deliverable of this case (a GB-failure map),
rather than a single ``max_damage`` number.

First-approximation reference model
-----------------------------------
- AT2 critical stress (GB):          sigma_c = sqrt(27 * E * Gc / (256 * lc))
- elliptical-tip stress concentration: K_t   = 2*(ax/ay) - 1   (same factor the
  non-regression already prints)
- load-bearing ligament: the intact fraction (1 - Fc) carries the load, so the
  ligament stress is K_t * p / (1 - Fc); the GB fails when it reaches sigma_c:

      p_crit(Fc, Gc) = sigma_c(Gc) * (1 - Fc) / K_t

  -> p_crit = sigma_c/K_t at Fc=0 (isolated bubble, tip nucleation) and -> 0 as
     Fc -> 1 (the ligament vanishes). Higher Gc shifts the whole curve up.

Coverage <-> spacing (areal, GB-surface convention)
---------------------------------------------------
      Fc = pi * a^2 / lambda^2   ->   lambda = a * sqrt(pi / Fc)
This is the convention mesh.geo uses (Lx = 2*ax*sqrt(pi/Fc_target)); it is the
physically correct one for a *surface* (the GB), distinct from the linear
coverage 2a/lambda and from any volumetric porosity (which involves the lens
height ay, not the footprint ax).

Usage
-----
    python3 parametric_study.py            # analytical reference map only (fast)
    python3 parametric_study.py --fem      # + FEM sweep overlay (slow; see below)
"""
import os
import sys
import shutil
import subprocess
import yaml
import numpy as np
import h5py
import matplotlib.pyplot as plt

from z3st.materials.oxide import GC_GB  # (J/m² == pJ/µm²)

CASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE, "output")

# --. material + numerics (read, never hardcode) --..
E = float(yaml.safe_load(open(os.path.join(CASE, "../../../materials/oxide.yaml")))["E"])  # MPa
lc = float(yaml.safe_load(open(os.path.join(CASE, "input.yaml")))["damage"]["lc"])         # µm
_g = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
ax, ay = float(_g["ax"]), float(_g["ay"])             # µm
K_t = 2.0 * (ax / ay) - 1.0                           # elliptical-tip concentration


def sigma_c(Gc):
    """AT2 critical stress (MPa) for a toughness Gc (pJ/µm² = J/m²)."""
    return np.sqrt(27.0 * E * Gc / (256.0 * lc))


def p_crit_ref(Fc, Gc):
    """First-approximation critical bubble pressure (MPa) — see module docstring."""
    return sigma_c(Gc) * (1.0 - Fc) / K_t


def spacing_from_Fc(Fc):
    """Bubble centre-to-centre spacing lambda (µm) for an areal coverage Fc."""
    return ax * np.sqrt(np.pi / Fc)


def reference_map(gc_values=(GC_GB, 1.0, 5.0), fc=np.linspace(0.02, 0.9, 200), ax_plot=None):
    """Plot p_crit vs Fc for several GB toughnesses (the ligament reference)."""
    own = ax_plot is None
    if own:
        plt.figure(figsize=(9, 6)); ax_plot = plt.gca()
    for Gc in gc_values:
        ax_plot.plot(fc, p_crit_ref(fc, Gc),
                     label=f"Gc = {Gc:g} J/m²  (σc = {sigma_c(Gc):.0f} MPa)")
    ax_plot.set_xlabel("GB bubble coverage  Fc  (areal)")
    ax_plot.set_ylabel("critical bubble pressure  p_crit  (MPa)")
    ax_plot.set_title(f"GB fracture map — Kt={K_t:.2f}, lc={lc} µm, E={E/1e3:.0f} GPa")
    ax_plot.grid(True, ls=":"); ax_plot.legend()
    if own:
        os.makedirs(OUT, exist_ok=True)
        path = os.path.join(OUT, "pcrit_vs_Fc_reference.png")
        plt.tight_layout(); plt.savefig(path, dpi=200)
        print(f"[INFO] reference map -> {path}")
    # report the case's operating point
    Fc0 = 4.0 * np.pi * ax**2 / float(_g["Lx"])**2 if "Lx" in _g else 0.40
    print(f"[INFO] operating point Fc≈{Fc0:.2f}, GB Gc={GC_GB} J/m²: "
          f"p_crit ≈ {p_crit_ref(Fc0, GC_GB):.0f} MPa "
          f"(bubble-spacing λ ≈ {spacing_from_Fc(Fc0):.3f} µm)")
    return ax_plot


# ---------------------------------------------------------------------------
# FEM sweep (opt-in, --fem): regenerate the mesh at each Fc, ramp the bubble
# pressure to p_crit, run z3st, and read off the percolation pressure (the
# first ramp step at which damage forms across the GB ligament). The case's
# mesh/BC/input are backed up and restored, so it is left untouched.
# ---------------------------------------------------------------------------
def _percolation_pressure(h5_path, pressures, band, threshold):
    """First ramp pressure at which max Damage in the GB band |y|<band reaches
    `threshold` (a proxy for ligament percolation), else NaN."""
    with h5py.File(h5_path, "r") as f:
        y = np.array(f["Mesh/mesh/geometry"])[:, 1]
        D = f["Function/Damage"]
        steps = sorted(D.keys(), key=lambda s: float(s.replace("_", ".")))
        gb = np.abs(y) < band
        dmax = []
        for s in steps:
            dmax.append(float(np.array(D[s]).reshape(-1)[gb].max()))
        dmax = np.array(dmax)
        hit = np.where(dmax >= threshold)[0]
        p = float(pressures[hit[0]]) if hit.size else float("nan")
        return p, dmax


def run_fem_point(Fc, p_max, n_steps=15, band=None, threshold=0.9, verbose=True):
    """Regenerate the mesh at coverage Fc, ramp the bubble pressure 0->p_max in
    n_steps, run z3st, and return (p_crit, max_damage_history). Backs up and
    restores mesh.msh / boundary_conditions.yaml / input.yaml."""
    band = band if band is not None else 2.0 * lc
    files = ["mesh.msh", "boundary_conditions.yaml", "input.yaml"]
    bk = {f: f"/tmp/bk_te_{f}" for f in files}
    for f in files:
        p = os.path.join(CASE, f)
        if os.path.exists(p):
            shutil.copy(p, bk[f])
    pressures = np.linspace(0.0, p_max, n_steps)
    try:
        # 1) mesh at this coverage
        subprocess.run(["gmsh", "-setnumber", "Fc_target", f"{Fc}",
                        os.path.join(CASE, "mesh.geo"), "-2",
                        "-o", os.path.join(CASE, "mesh.msh")],
                       check=True, capture_output=True)
        # 2) pressure ramp + matching n_steps
        bc = yaml.safe_load(open(os.path.join(CASE, "boundary_conditions.yaml")))
        for e in bc["mechanical"]["solid"]:
            if e.get("type") == "Neumann" and e.get("region") == "cavity":
                e["traction"] = [float(-p) for p in pressures]
        yaml.safe_dump(bc, open(os.path.join(CASE, "boundary_conditions.yaml"), "w"), sort_keys=False)
        inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
        inp["n_steps"] = int(n_steps)
        yaml.safe_dump(inp, open(os.path.join(CASE, "input.yaml"), "w"), sort_keys=False)
        # 3) run
        if verbose:
            print(f"[fem] Fc={Fc}: ramp 0->{p_max:.0f} MPa in {n_steps} steps ...", flush=True)
        subprocess.run([sys.executable, "-m", "z3st"], cwd=CASE, check=True,
                       capture_output=True, env={**os.environ, "HDF5_USE_FILE_LOCKING": "FALSE"})
        # 4) percolation pressure
        p_crit, dmax = _percolation_pressure(os.path.join(OUT, "results.h5"),
                                             pressures, band, threshold)
        if verbose:
            print(f"[fem] Fc={Fc}: max-damage per step = {np.round(dmax, 3)}")
            print(f"[fem] Fc={Fc}: p_crit = {p_crit if p_crit == p_crit else 'not reached'} MPa")
        return p_crit, dmax
    finally:
        for f in files:
            if os.path.exists(bk[f]):
                shutil.copy(bk[f], os.path.join(CASE, f))


def plot_case(Fc, p_crit, out_png):
    """Representative figure for one swept Fc: the final Damage field (the GB
    crack) next to sigma_yy, read straight from this run's results.h5."""
    from matplotlib.tri import Triangulation
    with h5py.File(os.path.join(OUT, "results.h5"), "r") as f:
        g = np.array(f["Mesh/mesh/geometry"]); xn, yn = g[:, 0], g[:, 1]
        topo = np.array(f["Mesh/mesh/topology"])

        def last(name, comp=None):
            grp = f[f"Function/{name}"]
            k = sorted(grp.keys(), key=lambda s: float(s.replace("_", ".")))[-1]
            a = np.array(grp[k]).reshape(np.array(grp[k]).shape[0], -1)
            return a[:, comp] if comp is not None else a.reshape(-1)

        d = last("Damage")
        syy = last("Stress", 4)

    tri = Triangulation(xn, yn, topo)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    cf = ax[0].tricontourf(tri, d, levels=np.linspace(0, 1, 21), cmap="inferno")
    fig.colorbar(cf, ax=ax[0]); ax[0].set_title(f"Damage (max = {d.max():.2f})")
    pc = ax[1].tripcolor(tri, facecolors=syy, cmap="viridis")
    fig.colorbar(pc, ax=ax[1]); ax[1].set_title(r"$\sigma_{yy}$ (MPa)")
    for a in ax:
        a.set_aspect("equal"); a.set_xlabel(r"x ($\mu$m)"); a.set_ylabel(r"y ($\mu$m)")
    tag = f"{p_crit:.0f} MPa" if p_crit == p_crit else "no fracture"
    fig.suptitle(f"Fc = {Fc:.2f}    |    $p_{{crit}}$ = {tag}    "
                 f"(λ = {spacing_from_Fc(Fc):.3f} µm)", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fem_sweep(fc_values, p_max=3000.0, n_steps=20):
    """Sweep coverage. Each point fractures at the energy-regime load (~p_max),
    saves a per-Fc figure, and returns its percolation pressure."""
    pts = []
    for Fc in fc_values:
        try:
            p_crit, _ = run_fem_point(Fc, p_max=p_max, n_steps=n_steps)
        except subprocess.CalledProcessError as e:
            print(f"[fem] Fc={Fc}: z3st failed: {e.stderr.decode()[-300:] if e.stderr else e}")
            p_crit = float("nan")
        out_png = os.path.join(OUT, f"sweep_Fc{Fc:.2f}.png")
        try:
            plot_case(Fc, p_crit, out_png)
            print(f"[fem] Fc={Fc}: figure -> {out_png}")
        except Exception as ex:
            print(f"[fem] Fc={Fc}: per-case figure skipped: {ex}")
        pts.append((Fc, p_crit))
    return pts


if __name__ == "__main__":
    ax_plot = reference_map()
    if "--fem" in sys.argv:
        # values after --fem are the Fc points to sweep (default: 0.2/0.4/0.6)
        rest = [a for a in sys.argv[sys.argv.index("--fem") + 1:] if not a.startswith("-")]
        fc_values = [float(a) for a in rest] if rest else [0.2, 0.4, 0.6]
        pts = fem_sweep(fc_values)
        ok = [(f, p) for f, p in pts if p == p]  # drop NaN
        if ok:
            fc, pc = zip(*ok)
            ax_plot.plot(fc, pc, "ks-", ms=9, label="phase-field $p_{crit}$")
            ax_plot.set_ylim(0, 1.1 * max(pc))
            ax_plot.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(OUT, "pcrit_vs_Fc_sweep.png"), dpi=200)
            print(f"[INFO] overlaid {len(ok)} FEM point(s) -> pcrit_vs_Fc_sweep.png")
