#!/usr/bin/env python3
import json, sys, argparse, pathlib
import numpy as np
import pyvista as pv

# Default list of field names to look for in the VTU file
DEFAULT_FIELDS = [
    "Temperature", "T", 
    "Displacement", "u", "U",
    "Von Mises (points)", "Von Mises (cells)",
    "Stress (points)", "Stress (cells)",
    "Strain (points)", "Strain (cells)",
]

def _mag(arr):
    """
    Compute magnitude if the array is vector/tensor data (2D array: N x dim).
    If array is scalar data (1D), return it unchanged.
    """
    if arr.ndim == 2:
        return np.linalg.norm(arr, axis=1)
    return arr

def compute_metrics(vtu_path: str, fields_wanted=None):
    """
    Reads a VTU file and computes simple metrics (min, max, mean, L2 norm, count)
    for the specified fields.

    Parameters:
        vtu_path (str): Path to the VTU file.
        fields_wanted (list of str): Optional override list of field names to search for.

    Returns:
        dict: mapping "location::field_name" → metrics dictionary.
    """
    grid = pv.read(vtu_path)
    fields_wanted = fields_wanted or DEFAULT_FIELDS

    # Collect all available point and cell data into a single dictionary
    avail = {}
    for k, v in grid.point_data.items():
        avail[f"point::{k}"] = v
    for k, v in grid.cell_data.items():
        avail[f"cell::{k}"] = v

    metrics = {}
    for name in fields_wanted:
        # Try to find exact matches for "::name"
        matches = [k for k in avail.keys() if k.endswith(f"::{name}")]
        # If no exact match, look for partial matches
        if not matches:
            matches = [k for k in avail.keys() if name in k]

        # For each matching field, compute metrics
        for key in matches:
            data = avail[key]
            arr = np.asarray(data)
            values = _mag(arr.astype(float))     # magnitude if vector/tensor
            values = values[np.isfinite(values)] # ignore NaN/inf

            if values.size == 0:
                continue

            m = {
                "min": float(values.min()),
                "max": float(values.max()),
                "mean": float(values.mean()),
                "l2norm": float(np.sqrt(np.dot(values, values))),
                "count": int(values.size),
            }
            metrics[key] = m
    return metrics

def compare_metrics(current, baseline, rtol=1e-3, atol=1e-8):
    """
    Compare current metrics against baseline metrics.

    Returns:
        ok_all (bool): True if all comparisons are within tolerances.
        diffs (dict): key → list of differences for failed checks.
    """
    diffs = {}
    ok_all = True
    for key, bvals in baseline.items():
        if key not in current:
            ok_all = False
            diffs[key] = [f"missing field in current output"]
            continue

        cvals = current[key]
        viol = []
        for k in ["min", "max", "mean", "l2norm"]:
            bv = bvals[k]
            cv = cvals[k]
            # Check if values are close within tolerances
            if not np.isclose(cv, bv, rtol=rtol, atol=atol):
                rel = (cv - bv) / (abs(bv) + 1e-30)
                viol.append(f"{k}: current={cv:.6e}, baseline={bv:.6e}, rel={rel:.3e}")
        if viol:
            ok_all = False
            diffs[key] = viol
    return ok_all, diffs

def main():
    """
    Entry point for the regression test script.
    Can either:
        - Write/update the baseline metrics from a VTU file
        python regression.py --vtu output/fields.vtu --write-baseline

        - Compare current VTU metrics with an existing baseline
        python regression.py --vtu output/fields.vtu

    """

    ap = argparse.ArgumentParser(description="Regression test (VTU).")
    ap.add_argument("--vtu", default="output/fields.vtu", help="Path to the VTU file.")
    ap.add_argument("--baseline", default="baseline.json",
                    help="Path to the baseline JSON file.")
    ap.add_argument("--write-baseline", action="store_true",
                    help="Write/update the baseline with current metrics.")
    ap.add_argument("--rtol", type=float, default=1e-3, help="Relative tolerance.")
    ap.add_argument("--atol", type=float, default=1e-3, help="Absolute tolerance.")
    ap.add_argument("--fields", nargs="*", default=None,
                    help="Optional list of fields to consider (overrides default).")
    args = ap.parse_args()

    # Compute metrics from the provided VTU file
    metrics = compute_metrics(args.vtu, args.fields)

    base_path = pathlib.Path(args.baseline)

    # Write baseline mode
    if args.write_baseline:
        base_path.parent.mkdir(parents=True, exist_ok=True)
        with open(base_path, "w") as f:
            json.dump(metrics, f, indent=2)
        print(f"[baseline]: {base_path}")
        return 0

    # Baseline must exist for comparison mode
    if not base_path.exists():
        print(f"[ERROR] baseline missing: {base_path}\n"
              f"→ generate baseline with: python {pathlib.Path(__file__).name} --write-baseline")
        return 2

    # Load baseline and compare
    with open(base_path) as f:
        baseline = json.load(f)

    ok, diffs = compare_metrics(metrics, baseline, rtol=args.rtol, atol=args.atol)
    if ok:
        print("[OK] regression test PASSED.")
        return 0
    else:
        print("[FAIL] regression test FAILED. Differences:")
        for key, issues in diffs.items():
            print(f"  - {key}:")
            for line in issues:
                print(f"      * {line}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
