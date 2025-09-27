#!/usr/bin/env python3
"""
Verification + non-regression checker for ZEST VTU outputs.

- Auto-detect .VTU file under output/
- Analytic verification for 1D Dirichlet-Dirichlet slab:
    T(x) = T0 + (T1 - T0) * (x / L)
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

# ---------------------- helpers: fields & plots ----------------------
def _detect_temp_field(grid):
    variable_names = ["Temperature", "temperature", "T", "temp"]
    for name in variable_names:
        if name in grid.point_data:
            return ("point", name)
    for name in variable_names:
        if name in grid.cell_data:
            return ("cell", name)
    return (None, None)

def plot_T_vs_x_centerline(vtu_path: str, out_png: str, out_csv: str = None,
                           resolution: int = 400, axis: int = 0,
                           T0_known: float = 600.0, T1_known: float = 800.0,
                           L_domain: float = 0.3) -> None:
    
    g = pv.read(vtu_path)
    loc, field = _detect_temp_field(g)
    if not field:
        print("[PLOT] Temperature field not found; skip centerline plot.")
        return

    # line along chosen axis
    mins = np.array(g.bounds)[[0, 2, 4]]
    maxs = np.array(g.bounds)[[1, 3, 5]]
    mid = 0.5 * (mins + maxs)

    p0 = mid.copy(); p1 = mid.copy()
    shrink = 1e-4 * (maxs[axis] - mins[axis])
    p0[axis] = mins[axis] + shrink
    p1[axis] = maxs[axis] - shrink

    try:
        samp = g.sample_over_line(tuple(p0), tuple(p1), resolution=resolution)
        x_raw = samp.points[:, axis]
        if loc == "point" and field in samp.point_data:
            T = np.asarray(samp.point_data[field], dtype=float)
        else:
            arrs = samp.point_data if field in samp.point_data else samp.cell_data
            T = np.asarray(arrs[field], dtype=float)
    except TypeError:
        line = pv.Line(tuple(p0), tuple(p1), resolution=resolution)
        samp = g.sample(line)
        x_raw = samp.points[:, axis]
        arrs = samp.point_data if field in samp.point_data else samp.cell_data
        T = np.asarray(arrs[field], dtype=float)

    # sort and rescale to [0, L_domain]
    ix = np.argsort(x_raw)
    x = (x_raw[ix] - x_raw.min()) / max(x_raw.max() - x_raw.min(), 1e-12) * L_domain
    T = T[ix]

    # analytic solution
    T_ana = T0_known + (T1_known - T0_known) * (x / L_domain)
    rel_L2 = float(np.sqrt(np.mean((T - T_ana) ** 2)) /
                   (np.sqrt(np.mean(T_ana ** 2)) + 1e-30))

    plt.figure()
    plt.plot(x, T, lw=2, label="numerical")
    plt.plot(x, T_ana, lw=2, ls="--", label=f"analytic (rel L2={rel_L2:.1e})")
    plt.xlabel("x [m]"); plt.ylabel("Temperature [K]")
    plt.title("Centerline temperature profile")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig(out_png, dpi=180); plt.close()
    print(f"[PLOT] Saved: {out_png}")

    if out_csv:
        np.savetxt(out_csv, np.column_stack([x, T, T_ana]),
                   delimiter=",", header="x,T,T_analytic", comments="")
        print(f"[PLOT] Saved: {out_csv}")


# ---------------------- metrics & verification ----------------------
def _mag(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 1:
        return arr
    return np.linalg.norm(arr.reshape(arr.shape[0], -1), axis=1)


def compute_metrics_from_vtu(vtu_path: str, fields=("Temperature", "Displacement")) -> dict:
    g = pv.read(vtu_path)
    pts = g.points
    out = {"n_points": int(len(pts))}
    avail = {k: np.asarray(v) for k, v in g.point_data.items()}
    if not avail and g.cell_data:
        g = g.cell_data_to_point_data()
        avail = {k: np.asarray(v) for k, v in g.point_data.items()}

    for name in fields:
        matches = [k for k in avail if k == name] or [k for k in avail if name in k]
        for key in matches:
            values = _mag(avail[key].astype(float))
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            out[f"{key}:min"] = float(values.min())
            out[f"{key}:max"] = float(values.max())
            out[f"{key}:mean"] = float(values.mean())
            out[f"{key}:L2"] = float(np.sqrt(np.mean(values**2)))
    return out

def verify_linear_slab(vtu_path: str, axis: int = 0,
                       T0_known: float = 600.0, T1_known: float = 800.0,
                       L_domain: float = 0.3,
                       tol_rel_L2: float = 1e-3, tol_Linf: float = 5e-3,
                       tol_len: float = 1e-12):
    g = pv.read(vtu_path)
    loc, field = _detect_temp_field(g)
    if not field:
        return False, {"error": "Temperature field not found in VTU."}

    if loc == "point":
        T = np.asarray(g.point_data[field], dtype=float)
    else:
        g = g.cell_data_to_point_data()
        T = np.asarray(g.point_data[field], dtype=float)

    # coordinate raw
    X = g.points[:, axis]
    xmin, xmax = float(X.min()), float(X.max())
    L_mesh = xmax - xmin

    # check consistency with L_domain
    if abs(L_mesh - L_domain) > tol_len:
        return False, {"error": f"Mesh length {L_mesh} != expected L={L_domain}"}

    # analytic line
    T_ana = T0_known + (T1_known - T0_known) * (X - xmin) / L_domain
    diff = T - T_ana
    rel_L2 = float(np.sqrt(np.mean(diff ** 2)) / (np.sqrt(np.mean(T_ana ** 2)) + 1e-30))
    Linf = float(np.max(np.abs(diff)))
    ok = (rel_L2 < tol_rel_L2) and (Linf < tol_Linf)
    return ok, {"rel_L2": rel_L2, "Linf": Linf, "T0": T0_known, "T1": T1_known,
                "L": L_domain, "L_mesh": L_mesh}



# ---------------------- baseline I/O & comparison ----------------------
def write_baseline(path: str, data: dict) -> None:
    # Add a creation date
    data_with_date = dict(data)
    data_with_date["_baseline_date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(path, "w") as f:
        json.dump(data_with_date, f, indent=2, sort_keys=True)
    print(f"[BASELINE] wrote {path} (date: {data_with_date['_baseline_date']})")



def read_json(path: str) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def compare_metrics(cur: dict, ref: dict, rtol=5e-4, atol=1e-6):
    diffs = {}
    keys = set(k for k, v in cur.items() if isinstance(v, (int, float))) | \
           set(k for k, v in ref.items() if isinstance(v, (int, float)))
    ok_all = True
    for k in sorted(keys):
        if k not in cur or k not in ref:
            ok_all = False
            diffs.setdefault(k, []).append(f"missing in {'current' if k not in cur else 'baseline'}")
            continue
        cv, rv = float(cur[k]), float(ref[k])
        if not (abs(cv - rv) <= atol + rtol * max(1.0, abs(cv), abs(rv))):
            ok_all = False
            diffs.setdefault(k, []).append(f"got {cv:.6e}, ref {rv:.6e}")
    return ok_all, diffs


# ---------------------- main ----------------------
def main():
    # ------------------ user constants ------------------
    T0_known = 600.0   # left boundary temperature [K]
    T1_known = 800.0   # right boundary temperature [K]
    L_domain  = 0.3    # domain length [m]
    # ----------------------------------------------------

    ap = argparse.ArgumentParser(description="Verification + non-regression checker for ZEST VTU outputs.")
    ap.add_argument("--vtu", help="Path to VTU (e.g., output/fields.vtu)")
    ap.add_argument("--verify-slab", action="store_true", help="Analytic check for 1D Dirichlet–Dirichlet slab.")
    ap.add_argument("--write-baseline", action="store_true", help="Write baseline.json from current metrics.")
    ap.add_argument("--baseline", default="baseline.json", help="Path to baseline JSON.")
    ap.add_argument("--rtol", type=float, default=5e-4, help="Relative tolerance for regression.")
    ap.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance for regression.")
    ap.add_argument("--plot-centerline", action="store_true", help="Save centerline T(x) plot.")
    ap.add_argument("--csv", action="store_true", help="Also write CSV files for generated plots.")
    ap.add_argument("--line-resolution", type=int, default=400, help="Samples along the centerline.")
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
    if args.verify_slab:
        ok_ver, info = verify_linear_slab(
            str(vtu),
            axis=args.axis,
            T0_known=T0_known,
            T1_known=T1_known,
            L_domain=L_domain
        )
        if ok_ver:
            print(f"[VERIF] linear slab: PASS  rel_L2={info['rel_L2']:.3e}, "
                  f"Linf={info['Linf']:.3e}  (T0={T0_known:.2f}, T1={T1_known:.2f}, L={L_domain:.2f})")
        else:
            print(f"[VERIF] linear slab: FAIL  details={info}")
            return 1
        metrics.update({f"verif:{k}": v for k, v in info.items() if isinstance(v, (int, float))})

    # 3) Plots
    if args.verify_slab or args.plot_centerline:
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

    # 4) Baseline logic
    if args.write_baseline:
        write_baseline(args.baseline, metrics)
        return 0

    if pathlib.Path(args.baseline).exists():
        baseline = read_json(args.baseline)
        ok, diffs = compare_metrics(metrics, baseline, rtol=args.rtol, atol=args.atol)
        if ok:
            print("[OK] regression test PASSED.")
            return 0
        print("[FAIL] regression test FAILED. Differences:")
        for key, issues in diffs.items():
            print(f"  - {key}:")
            for line in issues:
                print(f"      * {line}")
        return 1
    else:
        print(f"[WARN] baseline not found: {args.baseline}. Use --write-baseline to create it.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
