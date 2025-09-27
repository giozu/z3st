#!/usr/bin/env python3
"""
Verification + non-regression checker for ZEST VTU outputs.

- Auto-detect .VTU file under output/
- Analytic verification for 1D Dirichlet-Dirichlet slab:
    T(x) = T0 + (T1 - T0) * (x - xmin)/L
  with T0, T1, L defined inside main().
- Optional T(x) plots: centerline and cross-section-averaged
- Baseline freeze/compare with tolerances
"""
import os
import json
import sys
import argparse
import pathlib
from glob import glob
from datetime import datetime

import numpy as np

# Headless-safe plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyvista as pv
pv.OFF_SCREEN = True  # safe on headless/HPC

OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ANSI colors
class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    RESET = "\033[0m"



# ---------------------- helpers ----------------------
def _detect_temp_field(grid):
    cand_names = ["Temperature", "temperature", "T", "temp"]
    for name in cand_names:
        if name in grid.point_data:
            return ("point", name)
        else:
            print("Temperature field not found in point_data.")

    for name in cand_names:
        if name in grid.cell_data:
            return ("cell", name)
        else:
            print("Temperature field not found in cell_data.")

    return (None, None)


def _mag(arr: np.ndarray) -> np.ndarray:
    """ Return magnitude of array, handling both 1D and multi-D cases. """
    arr = np.asarray(arr)
    if arr.ndim == 1:
        return arr
    return np.linalg.norm(arr.reshape(arr.shape[0], -1), axis=1)


# ---------------------- metrics ----------------------
def compute_metrics_from_vtu(vtu_path: str, fields=("Temperature", "Displacement")) -> dict:
    # Legge il file VTU con pyvista
    grid = pv.read(vtu_path)

    # Array di coordinate di tutti i punti della mesh
    pts = grid.points

    # Inizializza il dizionario dei risultati con il numero di punti
    out = {"n_points": int(len(pts))}

    # Recupera i dati associati ai punti della mesh (point_data)
    # → es. Temperature, Displacement, PhaseField
    avail = {k: np.asarray(v) for k, v in grid.point_data.items()}

    # Se non c'è nulla in point_data ma ci sono dati a livello di cella,
    # li converte in point_data (così possiamo trattarli in modo uniforme)
    if not avail and grid.cell_data:
        grid = grid.cell_data_to_point_data()
        avail = {k: np.asarray(v) for k, v in grid.point_data.items()}

    # Loop sui campi richiesti (Temperature, Displacement, …)
    for name in fields:
        # Cerca corrispondenze: esatto o contenuto nel nome
        matches = [k for k in avail if k == name] or [k for k in avail if name in k]

        # Loop sui campi trovati
        for key in matches:
            # Applica _mag → così se è vettoriale (es. Displacement),
            # calcola il modulo in ogni punto, se è scalare lo lascia invariato
            values = _mag(avail[key].astype(float))

            # Elimina eventuali valori NaN o infiniti
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue

            # Calcola e salva nel dizionario alcune metriche di sintesi:
            out[f"{key}:min"]  = float(values.min())   # minimo
            out[f"{key}:max"]  = float(values.max())   # massimo
            out[f"{key}:mean"] = float(values.mean())  # media
            out[f"{key}:L2"]   = float(np.sqrt(np.mean(values**2))) # norma L2 media

    return out

# ---------------------- verification ----------------------
def verification(vtu_path: str, axis: int = 0,
                       T0_known: float = 600.0, T1_known: float = 800.0,
                       L_domain: float = 1.0,
                       tol_len: float = 1e-12):
    """
    Verification vs. analytic T(x) = T0 + (T1 - T0) * (x - xmin)/L_domain.
    No rescaling: requires mesh length == L_domain, otherwise FAIL.
    """
    grid = pv.read(vtu_path)
    loc, field = _detect_temp_field(grid) # "point", "Temperature"

    if loc == "point":
        T = np.asarray(grid.point_data[field], dtype=float)

    X = grid.points[:, axis]
    xmin, xmax = float(X.min()), float(X.max())
    L_mesh = xmax - xmin

    # check mesh length
    if abs(L_mesh - L_domain) > tol_len:
        return False, {"error": f"Mesh length {L_mesh} != expected L={L_domain}"}

    # analytic solution
    T_ana = T0_known + (T1_known - T0_known) * (X - xmin) / L_domain
    diff = T - T_ana
    L2 = np.sqrt(np.mean(diff ** 2))
    Linf = np.max(np.abs(diff))
    
    return  {"L2": L2, "Linf": Linf,
             "T0": T0_known, "T1": T1_known,
             "L": L_domain, "L_mesh": L_mesh}


# ---------------------- plots ----------------------
def plot_T_vs_x_centerline(vtu_path: str, out_png: str, out_csv: str = None,
                           resolution: int = 400, axis: int = 0,
                           T0_known: float = 600.0, T1_known: float = 800.0,
                           L_domain: float = 1.0, tol_len: float = 1e-12) -> None:
    grid = pv.read(vtu_path)
    loc, field = _detect_temp_field(grid)

    # bounds
    mins = np.array(grid.bounds)[[0, 2, 4]] # xmin, ymin, zmin
    maxs = np.array(grid.bounds)[[1, 3, 5]] # xmax, ymax, zmax
    xmin, xmax = float(mins[axis]), float(maxs[axis])
    L_mesh = xmax - xmin
    if abs(L_mesh - L_domain) > tol_len:
        print(f"[PLOT] WARNING: Mesh length {L_mesh} != expected L={L_domain}.")

    # line through center
    mid = 0.5 * (mins + maxs)
    p0, p1 = mid.copy(), mid.copy()
    shrink = 1e-4 * L_mesh
    p0[axis] = xmin + shrink
    p1[axis] = xmax - shrink

    try:
        samp = grid.sample_over_line(tuple(p0), tuple(p1), resolution=resolution)
        x = samp.points[:, axis]
        if loc == "point" and field in samp.point_data:
            T = np.asarray(samp.point_data[field], dtype=float)
        else:
            arrs = samp.point_data if field in samp.point_data else samp.cell_data
            T = np.asarray(arrs[field], dtype=float)
    except TypeError:
        line = pv.Line(tuple(p0), tuple(p1), resolution=resolution)
        samp = grid.sample(line)
        x = samp.points[:, axis]
        arrs = samp.point_data if field in samp.point_data else samp.cell_data
        T = np.asarray(arrs[field], dtype=float)

    ix = np.argsort(x)
    x, T = x[ix], T[ix]

    # analytic
    T_ana = T0_known + (T1_known - T0_known) * (x - xmin) / L_domain
    L2 = float(np.sqrt(np.mean((T - T_ana) ** 2)) /
                   (np.sqrt(np.mean(T_ana ** 2)) + 1e-30))

    plt.figure()
    plt.plot(x, T, lw=2, label="numerical")
    plt.plot(x, T_ana, lw=2, ls="--", label=f"analytic")
    plt.xlabel("x [m]"); plt.ylabel("Temperature [K]")
    plt.title("Centerline temperature profile")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig(out_png, dpi=180); plt.close()
    print(f"[PLOT] Saved: {out_png}")

    if out_csv:
        np.savetxt(out_csv, np.column_stack([x, T, T_ana]),
                   delimiter=",", header="x,T,T_analytic", comments="")
        print(f"[PLOT] Saved: {out_csv}")


def plot_T_vs_x_avg(vtu_path: str, out_png: str, out_csv: str = None,
                    nbins: int = 100, axis: int = 0,
                    T0_known: float = 600.0, T1_known: float = 800.0,
                    L_domain: float = 1.0, tol_len: float = 1e-12) -> None:
    grid = pv.read(vtu_path)
    loc, field = _detect_temp_field(grid)
    if not field:
        print("[PLOT] Temperature field not found; skip averaged plot.")
        return

    # bounds
    mins = np.array(grid.bounds)[[0, 2, 4]]
    maxs = np.array(grid.bounds)[[1, 3, 5]]
    xmin, xmax = float(mins[axis]), float(maxs[axis])
    L_mesh = xmax - xmin
    if abs(L_mesh - L_domain) > tol_len:
        print(f"[PLOT] ERROR: Mesh length {L_mesh} != expected L={L_domain}. No plot generated.")
        return

    # values
    x = grid.points[:, axis]
    if loc == "point":
        T = np.asarray(grid.point_data[field], dtype=float)

    # binning
    edges = np.linspace(xmin, xmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    T_avg = np.empty(nbins)
    for i in range(nbins):
        mask = (x >= edges[i]) & (x < edges[i + 1])
        T_avg[i] = T[mask].mean() if mask.any() else np.nan
    valid = np.isfinite(T_avg)
    centers, T_avg = centers[valid], T_avg[valid]

    # analytic on bin centers
    T_ana = T0_known + (T1_known - T0_known) * (centers - xmin) / L_domain

    plt.figure()
    plt.plot(centers, T_avg, lw=2, label="numerical")
    plt.plot(centers, T_ana, lw=2, ls="--", label="analytic")
    plt.xlabel("x [m]"); plt.ylabel("⟨Temperature⟩ over y–z [K]")
    plt.title("Cross-section averaged temperature")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig(out_png, dpi=180); plt.close()
    print(f"[PLOT] Saved: {out_png}")

    if out_csv:
        np.savetxt(out_csv, np.column_stack([centers, T_avg, T_ana]),
                   delimiter=",", header="x,T_avg,T_analytic", comments="")
        print(f"[PLOT] Saved: {out_csv}")


# ---------------------- baseline ----------------------
def write_baseline(path: str, data: dict) -> None:
    data_with_date = dict(data)
    data_with_date["_baseline_date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "w") as f:
        json.dump(data_with_date, f, indent=2, sort_keys=True)
    print(f"[BASELINE] wrote {path} (date: {data_with_date['_baseline_date']})")


def read_json(path: str) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def compare_metrics(cur: dict, ref: dict, atol=1e-6):
    """
    Compare current vs baseline metrics with tolerances.
    Returns (ok_all, rows) where rows is a list of (key, cur, ref, diff, status).
    """
    rows = []
    keys = set(k for k, v in cur.items() if isinstance(v, (int, float))) | \
           set(k for k, v in ref.items() if isinstance(v, (int, float)))
    ok_all = True
    for k in sorted(keys):
        if k not in cur or k not in ref:
            rows.append((k, None, None, None, "MISSING"))
            ok_all = False
            continue
        cv, rv = float(cur[k]), float(ref[k])
        tol_ok = abs(cv - rv) <= atol
        status = "OK" if tol_ok else "FAIL"
        if not tol_ok:
            ok_all = False
        rows.append((k, cv, rv, cv - rv, status))
    return ok_all, rows



# ---------------------- main ----------------------
def main():
    # ------------------ user constants ------------------
    T0_known = 600.0   # left Dirichlet temperature [K]
    T1_known = 800.0   # right Dirichlet temperature [K]
    L_domain  = 0.3    # domain length [m]
    # ----------------------------------------------------

    ap = argparse.ArgumentParser(description="Verification + non-regression checker for ZEST VTU outputs.")
    ap.add_argument("--vtu", help="Path to VTU (e.grid., output/fields.vtu)")
    ap.add_argument("--write-baseline", action="store_true", help="Write baseline.json from current metrics.")
    ap.add_argument("--baseline", default="baseline.json", help="Path to baseline JSON.")
    ap.add_argument("--atol", type=float, default=1e-4, help="Absolute tolerance for regression.")
    ap.add_argument("--plot-centerline", action="store_true", help="Save centerline T(x) plot.")
    ap.add_argument("--plot-avg", action="store_true", help="Save cross-section-averaged T(x) plot.")
    ap.add_argument("--csv", action="store_true", help="Also write CSV files for generated plots.")
    ap.add_argument("--line-resolution", type=int, default=400, help="Samples along the centerline.")
    ap.add_argument("--avg-bins", type=int, default=100, help="Number of x-bins for the averaged plot.")
    ap.add_argument("--axis", type=int, default=0, choices=(0, 1, 2), help="Axis for 1D treatment (0=x, 1=y, 2=z).")
    args = ap.parse_args()

    # Auto-detect VTU if not provided
    vtu_path = args.vtu
    if not vtu_path:
        cands = sorted(glob(os.path.join(OUTPUT_DIR, "fields_step*.vtu")) +
                       glob(os.path.join(OUTPUT_DIR, "fields.vtu")))
        assert cands, "No VTU found in output/. Run the case first."
        vtu_path = cands[-1]
        print(f"[INFO] Auto-selected VTU: {vtu_path}")

    vtu = pathlib.Path(vtu_path)
    assert vtu.exists(), f"VTU not found: {vtu}"

    # 1) Base metrics
    metrics = compute_metrics_from_vtu(str(vtu))
    print("[METRICS]")
    for k in sorted(metrics.keys()):
        print(f"  {k}: {metrics[k]}")

    # 2) Analytic verification
    info = verification(
        str(vtu),
        axis=args.axis,
        T0_known=T0_known,
        T1_known=T1_known,
        L_domain=L_domain
    )
    print(f"[VERIFICATION]: L2={info['L2']:.3e}, Linf={info['Linf']:.3e}")
    metrics.update({f"verif:{k}": v for k, v in info.items() if isinstance(v, (int, float))})

    # 3) Plots
    if args.plot_centerline:
        plot_T_vs_x_centerline(
            str(vtu),
            out_png=os.path.join(OUTPUT_DIR, "T_vs_x_centerline.png"),
            out_csv=os.path.join(OUTPUT_DIR, "T_vs_x_centerline.csv") if args.csv else None,
            resolution=args.line_resolution,
            axis=args.axis,
            T0_known=T0_known,
            T1_known=T1_known,
            L_domain=L_domain
        )
    if args.plot_avg:
        plot_T_vs_x_avg(
            str(vtu),
            out_png=os.path.join(OUTPUT_DIR, "T_vs_x_avg.png"),
            out_csv=os.path.join(OUTPUT_DIR, "T_vs_x_avg.csv") if args.csv else None,
            nbins=args.avg_bins,
            axis=args.axis,
            T0_known=T0_known,
            T1_known=T1_known,
            L_domain=L_domain
        )

    # 4) Baseline logic
    if args.write_baseline:
        write_baseline(args.baseline, metrics)
        return 0

    if pathlib.Path(args.baseline).exists():
        baseline = read_json(args.baseline)
        ok, rows = compare_metrics(metrics, baseline, atol=args.atol)

        print("\n" + "─"*29 + " Regression test " + "─"*29)
        print(f"{'Metric':20s} {'Current':>12s} {'Baseline':>12s} {'Δ (cur-ref)':>12s} {'Status':>8s}")
        print("─"*70)
        for k, cv, rv, diff, status in rows:
            if cv is None or rv is None:
                print(f"{k:20s} {'-':>12s} {'-':>12s} {'-':>12s} {status:>8s}")
            else:
                color = Colors.GREEN if status == "OK" else Colors.RED
                print(f"{k:20s} {cv:12.4e} {rv:12.4e} {diff:12.2e} {color}{status:>8s}{Colors.RESET}")
        print("─"*70)
        overall_status = f"{Colors.GREEN}PASSED{Colors.RESET}" if ok else f"{Colors.RED}FAILED{Colors.RESET}"
        print(f"Overall: {overall_status}")

    else:
        print(f"[WARN] baseline not found: {args.baseline}. Use --write-baseline to create it.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
