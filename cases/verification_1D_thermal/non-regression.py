#!/usr/bin/env python3
"""
Verification + non-regression checker for ZEST VTU outputs.

- Auto-detect VTU under output/
- Analytic verification for 1D Dirichlet–Dirichlet slab (linear T(x))
  using KNOWN boundary temperatures from boundary_conditions.yaml (or --T0/--T1)
- Optional T(x) plots: centerline and cross-section-averaged
- Baseline freeze/compare with tolerances
"""
import os
import json
import sys
import argparse
import pathlib
from glob import glob

import numpy as np
import yaml

# Headless-safe plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pyvista as pv
pv.OFF_SCREEN = True  # safe on headless/HPC

OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------- helpers: BC reading ----------------------
def resolve_bc_path(input_path: str = "input.yaml", bc_override: str = None) -> str:
    """Locate boundary_conditions.yaml, preferring explicit override, then input.yaml."""
    if bc_override:
        return bc_override
    try:
        with open(input_path, "r") as f:
            data = yaml.safe_load(f) or {}
        return data.get("boundary_conditions_path", "boundary_conditions.yaml")
    except Exception:
        return "boundary_conditions.yaml"


def read_T0_T1_from_bc(bc_path: str, left_region: str = "xmin", right_region: str = "xmax"):
    """
    Parse boundary_conditions.yaml for Dirichlet temperatures on given regions.
    Expected structure:
      thermal_bcs:
        <block_name>:
          - {type: Dirichlet, region: xmin, temperature: 600.0}
          - {type: Dirichlet, region: xmax, temperature: 800.0}
    Returns (T0, T1) as floats or (None, None) if not found.
    """
    try:
        with open(bc_path, "r") as f:
            bc = yaml.safe_load(f) or {}
    except Exception:
        return None, None

    tbcs = bc.get("thermal_bcs", {})
    T0 = T1 = None
    if isinstance(tbcs, dict):
        for _blk, lst in tbcs.items():
            if not isinstance(lst, list):
                continue
            for entry in lst:
                if not isinstance(entry, dict):
                    continue
                if str(entry.get("type", "")).lower() != "dirichlet":
                    continue
                reg = entry.get("region")
                T = entry.get("temperature", None)
                if T is None:
                    continue
                try:
                    T = float(T)
                except Exception:
                    continue
                if reg == left_region:
                    T0 = T
                elif reg == right_region:
                    T1 = T
    return T0, T1


# ---------------------- helpers: fields & plots ----------------------
def _detect_temp_field(grid):
    # Prefer point data called "Temperature"; fallbacks allowed
    if "Temperature" in grid.point_data:
        return ("point", "Temperature")
    if "Temperature" in grid.cell_data:
        return ("cell", "Temperature")
    for k in grid.point_data.keys():
        if "Temp" in k:
            return ("point", k)
    for k in grid.cell_data.keys():
        if "Temp" in k:
            return ("cell", k)
    return (None, None)


def plot_T_vs_x_centerline(vtu_path: str, out_png: str, out_csv: str = None,
                           resolution: int = 400, axis: int = 0,
                           T0_known: float = None, T1_known: float = None) -> None:
    g = pv.read(vtu_path)
    loc, field = _detect_temp_field(g)
    if not field:
        print("[PLOT] Temperature field not found; skip centerline plot.")
        return

    xmin, xmax, ymin, ymax, zmin, zmax = g.bounds
    if axis != 0:
        print("[PLOT] Centerline currently implemented for axis=0 (x). Using axis=0.")

    # shrink the sampled segment so endpoints are inside the domain
    L = xmax - xmin
    shrink = 1e-4 * max(L, 1.0)
    p0 = (xmin + shrink, 0.5*(ymin + ymax), 0.5*(zmin + zmax))
    p1 = (xmax - shrink, 0.5*(ymin + ymax), 0.5*(zmin + zmax))

    # sample along the line (PyVista API compatibility)
    try:
        samp = g.sample_over_line(p0, p1, resolution=resolution)
        x = samp.points[:, 0]
        if loc == "point" and field in samp.point_data:
            T = np.asarray(samp.point_data[field], dtype=float)
        elif loc == "cell" and field in samp.cell_data:
            T = np.asarray(samp.cell_data[field], dtype=float)
        else:
            arrs = samp.point_data if samp.point_data else samp.cell_data
            T = np.asarray(arrs[field], dtype=float)
    except TypeError:
        line = pv.Line(p0, p1, resolution=resolution)
        samp = g.sample(line)
        x = samp.points[:, 0]
        arrs = samp.point_data if field in samp.point_data else samp.cell_data
        T = np.asarray(arrs[field], dtype=float)

    # sort
    ix = np.argsort(x)
    x = x[ix]
    T = T[ix]

    # Use known boundary temps if provided; otherwise fall back to line averages near ends
    xs0, xs1 = float(x.min()), float(x.max())
    Ls = max(xs1 - xs0, 1e-12)
    if T0_known is None or T1_known is None:
        # robust fallback: average small zones near ends of the SAMPLED line
        eps = 0.01 * Ls
        left_mask = x <= (xs0 + eps)
        right_mask = x >= (xs1 - eps)
        T0_est = float(np.nanmean(T[left_mask])) if left_mask.any() else float(T[0])
        T1_est = float(np.nanmean(T[right_mask])) if right_mask.any() else float(T[-1])
        T0 = T0_known if T0_known is not None else T0_est
        T1 = T1_known if T1_known is not None else T1_est
    else:
        T0, T1 = float(T0_known), float(T1_known)

    # analytic line defined over the sampled segment
    T_ana = T0 + (T1 - T0) * (x - xs0) / Ls
    denom = np.sqrt(np.mean(T_ana**2) + 1e-30)
    rel_L2 = float(np.sqrt(np.mean((T - T_ana)**2)) / denom)

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


def plot_T_vs_x_avg(vtu_path: str, out_png: str, out_csv: str = None,
                    nbins: int = 100, axis: int = 0) -> None:
    g = pv.read(vtu_path)
    loc, field = _detect_temp_field(g)
    if not field:
        print("[PLOT] Temperature field not found; skip averaged plot.")
        return

    x = g.points[:, axis]
    # choose the correct data container
    if loc == "point":
        T = np.asarray(g.point_data[field], dtype=float)
    else:
        # for cell data, expand to points via interpolation (simple fallback)
        T = np.asarray(g.point_data.get(field, g.cell_data[field]), dtype=float)
    xmin, xmax = float(x.min()), float(x.max())
    edges = np.linspace(xmin, xmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    T_avg = np.empty(nbins)

    for i in range(nbins):
        right_edge = edges[i + 1] if i < nbins - 1 else edges[i + 1] + 1e-15
        mask = (x >= edges[i]) & (x < right_edge)
        T_avg[i] = T[mask].mean() if mask.any() else np.nan

    valid = np.isfinite(T_avg)
    centers, T_avg = centers[valid], T_avg[valid]

    plt.figure()
    plt.plot(centers, T_avg, lw=2)
    plt.xlabel("x [m]")
    plt.ylabel("⟨Temperature⟩ over y–z [K]")
    plt.title("Cross-section averaged temperature")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"[PLOT] Saved: {out_png}")

    if out_csv:
        np.savetxt(out_csv, np.column_stack([centers, T_avg]),
                   delimiter=",", header="x,T_avg", comments="")
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
                       tol_rel_L2: float = 1e-3, tol_Linf: float = 5e-3,
                       T0_known: float = None, T1_known: float = None):
    """
    Verification vs. analytic T(x) = T0 + (T1 - T0) * (x - xmin) / (xmax - xmin)
    using KNOWN boundary temps when provided. Otherwise falls back to robust inference over the mesh.
    """
    g = pv.read(vtu_path)
    loc, field = _detect_temp_field(g)
    if not field:
        return False, {"error": "Temperature field not found in VTU."}

    # use point data if available; otherwise interpolate cell data to points via PyVista (simple)
    if loc == "point":
        T = np.asarray(g.point_data[field]).astype(float)
    else:
        # cell -> point interpolation (approx)
        g_point = g.point_data.copy()
        if field in g_point:
            T = np.asarray(g_point[field]).astype(float)
        else:
            # fallback: use cell data replicated to points by nearest cell center (rough)
            T = np.repeat(np.asarray(g.cell_data[field]).astype(float), 1)[:len(g.points)]

    X = g.points[:, axis]
    xmin, xmax = float(X.min()), float(X.max())
    L = xmax - xmin
    if L <= 0:
        return False, {"error": "Degenerate geometry (L<=0) on chosen axis."}

    # Known boundary temperatures preferred
    if T0_known is None or T1_known is None:
        left_mask = np.isclose(X, xmin, atol=1e-8)
        right_mask = np.isclose(X, xmax, atol=1e-8)
        if left_mask.any():
            T0_est = float(T[left_mask].mean())
        else:
            k = max(5, int(0.5 * (len(X) ** (1 / 3))))
            idx_left = np.argpartition(np.abs(X - xmin), k)[:k]
            T0_est = float(T[idx_left].mean())
        if right_mask.any():
            T1_est = float(T[right_mask].mean())
        else:
            k = max(5, int(0.5 * (len(X) ** (1 / 3))))
            idx_right = np.argpartition(np.abs(X - xmax), k)[:k]
            T1_est = float(T[idx_right].mean())

        T0 = T0_known if T0_known is not None else T0_est
        T1 = T1_known if T1_known is not None else T1_est
    else:
        T0, T1 = float(T0_known), float(T1_known)

    eps = 1e-12
    T_ana = T0 + (T1 - T0) * (X - xmin) / max(L, eps)
    diff = T - T_ana
    denom = np.sqrt(np.mean(T_ana ** 2) + 1e-30)
    rel_L2 = float(np.sqrt(np.mean(diff ** 2)) / denom)
    Linf = float(np.max(np.abs(diff)))
    ok = (rel_L2 < tol_rel_L2) and (Linf < tol_Linf)
    return ok, {"rel_L2": rel_L2, "Linf": Linf, "T0": T0, "T1": T1, "xmin": xmin, "xmax": xmax}


# ---------------------- baseline I/O & comparison ----------------------
def write_baseline(path: str, data: dict) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
    print(f"[BASELINE] wrote {path}")


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
    ap = argparse.ArgumentParser(description="Verification + non-regression checker for ZEST VTU outputs.")
    ap.add_argument("--vtu", help="Path to VTU (e.g., output/fields.vtu)")
    ap.add_argument("--verify-slab", action="store_true", help="Analytic check for 1D Dirichlet–Dirichlet slab.")
    ap.add_argument("--write-baseline", action="store_true", help="Write baseline.json from current metrics.")
    ap.add_argument("--baseline", default="baseline.json", help="Path to baseline JSON.")
    ap.add_argument("--rtol", type=float, default=5e-4, help="Relative tolerance for regression.")
    ap.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance for regression.")
    ap.add_argument("--plot-centerline", action="store_true", help="Save centerline T(x) plot.")
    ap.add_argument("--plot-avg", action="store_true", help="Save cross-section-averaged T(x) plot.")
    ap.add_argument("--csv", action="store_true", help="Also write CSV files for generated plots.")
    ap.add_argument("--line-resolution", type=int, default=400, help="Samples along the centerline.")
    ap.add_argument("--avg-bins", type=int, default=100, help="Number of x-bins for the averaged plot.")
    ap.add_argument("--axis", type=int, default=0, choices=(0, 1, 2), help="Axis for 1D treatment (0=x, 1=y, 2=z).")
    # NEW: known BC values and paths
    ap.add_argument("--input", default="input.yaml", help="Path to input.yaml (to locate BCs).")
    ap.add_argument("--bc-file", help="Path to boundary_conditions.yaml (overrides --input).")
    ap.add_argument("--left-region", default="xmin", help="Region name for left Dirichlet boundary.")
    ap.add_argument("--right-region", default="xmax", help="Region name for right Dirichlet boundary.")
    ap.add_argument("--T0", type=float, help="Left boundary temperature [K] (overrides BC file).")
    ap.add_argument("--T1", type=float, help="Right boundary temperature [K] (overrides BC file).")
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

    # Resolve known boundary temperatures T0,T1
    T0_known = args.T0
    T1_known = args.T1
    if T0_known is None or T1_known is None:
        bc_path = resolve_bc_path(args.input, args.bc_file)
        T0_bc, T1_bc = read_T0_T1_from_bc(bc_path, args.left_region, args.right_region)
        if T0_known is None:
            T0_known = T0_bc
        if T1_known is None:
            T1_known = T1_bc
        if T0_known is not None and T1_known is not None:
            print(f"[INFO] Using BC temperatures: T0={T0_known} K, T1={T1_known} K from {bc_path}")
        else:
            print("[INFO] Could not read both T0/T1 from BC file; falling back to inference where needed.")

    # If baseline has verif:* keys, auto-enable verification
    if pathlib.Path(args.baseline).exists():
        baseline_meta = read_json(args.baseline)
        if any(k.startswith("verif:") for k in baseline_meta) and not args.verify_slab:
            print("[INFO] Baseline includes verification keys; enabling --verify-slab automatically.")
            args.verify_slab = True

    # 1) Base metrics
    metrics = compute_metrics_from_vtu(str(vtu))
    print("[METRICS]")
    for k in sorted(metrics.keys()):
        print(f"  {k}: {metrics[k]}")

    # 2) Analytic verification (1D slab) with KNOWN BCs
    if args.verify_slab:
        ok_ver, info = verify_linear_slab(
            str(vtu),
            axis=args.axis,
            T0_known=T0_known,
            T1_known=T1_known
        )
        if ok_ver:
            print(f"[VERIF] linear slab: PASS  rel_L2={info['rel_L2']:.3e}, "
                  f"Linf={info['Linf']:.3e}  (T0={info['T0']:.2f}, T1={info['T1']:.2f})")
        else:
            print(f"[VERIF] linear slab: FAIL  details={info}")
            return 1
        metrics.update({f"verif:{k}": v for k, v in info.items() if isinstance(v, (int, float))})

    # 3) Plots (centerline uses known T0/T1 if available)
    if args.verify_slab or args.plot_centerline:
        plot_T_vs_x_centerline(
            str(vtu),
            out_png=os.path.join(OUTPUT_DIR, "T_vs_x_centerline.png"),
            out_csv=os.path.join(OUTPUT_DIR, "T_vs_x_centerline.csv") if args.csv else None,
            resolution=args.line_resolution,
            axis=args.axis,
            T0_known=T0_known,
            T1_known=T1_known
        )
    if args.plot_avg:
        plot_T_vs_x_avg(
            str(vtu),
            out_png=os.path.join(OUTPUT_DIR, "T_vs_x_avg.png"),
            out_csv=os.path.join(OUTPUT_DIR, "T_vs_x_avg.csv") if args.csv else None,
            nbins=args.avg_bins,
            axis=args.axis,
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
