# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import json
import os
import sys
import collections.abc
import numpy as np

USE_COLOR = sys.stdout.isatty()

GREEN = "\033[92m" if USE_COLOR else ""
RED = "\033[91m" if USE_COLOR else ""
BOLD = "\033[1m" if USE_COLOR else ""
END = "\033[0m" if USE_COLOR else ""


def pass_fail_check(errors, tolerance, out_json, case_dir):
    """
    Evaluate relative errors in 'errors' dict and mark PASS/FAIL.

    Parameters
    ----------
    errors : dict
        Contains entries with 'rel_error' values.
    tolerance : float
        Relative tolerance for pass/fail check.
    out_json : str
        Path to output JSON file.
    case_dir : str
        Case directory (for logs).
    """
    print(f"\nPass/Fail check (tolerance = {tolerance:.1e})")

    all_pass = True
    for key, val in errors.items():
        rel_err = val["rel_error"]

        if isinstance(rel_err, (list, tuple)):
            rel_arr = np.array(rel_err, dtype=float)
            passed = np.all(rel_arr < tolerance)
            rel_err_display = float(np.max(rel_arr))
        else:
            passed = rel_err < tolerance
            rel_err_display = float(rel_err)

        color = GREEN if passed else RED
        status = "PASS" if passed else "FAIL"
        all_pass &= passed

        print(f"  {key:18s} → rel err = {rel_err_display:.2e}  → {color}{status}{END}")
        val["status"] = status

    summary = (
        f"{GREEN}PASS{END} All checks within tolerance"
        if all_pass
        else f"{RED}FAIL{END} Some checks exceeded tolerance"
    )
    print(f"\n[SUMMARY] {BOLD}{summary}{END}")

    # Save results
    with open(out_json, "w") as f:
        json.dump(
            {"results": errors, "tolerance": tolerance, "summary": "PASS" if all_pass else "FAIL"},
            f,
            indent=4,
        )

    print(f"[INFO] non-regression results written to: {out_json}")
    return all_pass


def _write_regression_verdict(case_dir, verdict):
    """Persist the gold-regression verdict ('PASS'/'FAIL') into the run's
    non-regression.json so the suite drivers (non-regression*.sh) can surface
    it — a regression is otherwise invisible outside stdout."""
    out_json = os.path.join(case_dir, "output", "non-regression.json")
    try:
        with open(out_json, "r") as f:
            data = json.load(f)
        data["regression"] = verdict
        with open(out_json, "w") as f:
            json.dump(data, f, indent=4)
    except (FileNotFoundError, json.JSONDecodeError):
        pass


def regression_check(errors, case_dir, regression_tol=1e-3, near_zero_factor=5e-2):
    """
    This function evaluates the qualitative performance of the current run compared to the 
    GOLD benchmark. It normalizes error data (handling both scalars and sequences) and 
    determines if the solver's accuracy has improved or degraded, providing immediate 
    visual feedback on the impact of recent code changes or mesh refinements.

    Parameters
    ----------
    errors : dict
        Current non-regression results.
    case_dir : str
        Directory of the case (expects output/non-regression_gold.json inside).
    regression_tol : float
        Relative tolerance for regression comparison.
    """
    gold_file = os.path.join(case_dir, "output", "non-regression_gold.json")

    if not os.path.exists(gold_file):
        print(f"\n[INFO] No GOLD file found → skipping regression check.")
        return None

    print(f"\nRegression check vs GOLD reference:")

    try:
        with open(gold_file, "r") as f:
            gold_data = json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        # A corrupt / unparseable gold must FAIL loudly, never be silently
        # skipped: otherwise a broken gold disables the regression check and the
        # suite reports a false pass.
        print(f"  {RED}[ERROR] GOLD file unreadable ({exc}); marking regression FAIL.{END}")
        _write_regression_verdict(case_dir, "FAIL")
        return False

    gold_results = gold_data.get("results", gold_data)
    reg_pass = True

    for key in sorted(set(errors) | set(gold_results)):
        if key not in errors:
            print(f"  {key:18s} → present in GOLD but MISSING from run → REGRESSION")
            reg_pass = False
            continue

        num_now = errors[key]["numerical"]
        num_gold = gold_results.get(key, {}).get("numerical", None)

        rel_err_now = errors[key].get("rel_error", None)
        rel_err_gold = gold_results.get(key, {}).get("rel_error", None)

        if num_gold is None:
            print(f"  {key:18s} → missing in GOLD file, skipped")
            continue

        num_now_arr = (
            np.array(num_now, dtype=float) if isinstance(num_now, list) else np.array([num_now])
        )
        num_gold_arr = (
            np.array(num_gold, dtype=float) if isinstance(num_gold, list) else np.array([num_gold])
        )

        if num_now_arr.shape != num_gold_arr.shape:
            print(f"  {key:18s} → shape {num_now_arr.shape} vs GOLD {num_gold_arr.shape} → REGRESSION")
            reg_pass = False
            continue

        abs_diff_arr = np.abs(num_now_arr - num_gold_arr)
        rel_diff_arr = abs_diff_arr / np.maximum(np.abs(num_gold_arr), 1e-12)
        rel_diff_display = float(np.max(rel_diff_arr))

        # Fields whose analytical reference is exactly zero are 'should-be-zero'
        # residuals: their numerical value is floating-point noise (build- and
        # BLAS-dependent), so a RELATIVE comparison to ~0 is ill-defined -- a few
        # parts in 1e8 of a ~1e-6 residual spuriously flips the gold. Compare these
        # by an ABSOLUTE tolerance scaled to the residual's own (near-zero)
        # magnitude instead. This can only rescue brittle near-zero fields: a field
        # already inside the relative band stays inside this looser absolute band.
        ref_val = gold_results.get(key, {}).get("reference", errors[key].get("reference"))
        near_zero_ref = ref_val is not None and np.all(
            np.abs(np.atleast_1d(np.asarray(ref_val, dtype=float))) <= 0.0
        )
        if near_zero_ref:
            atol_eff = max(1e-12, near_zero_factor * float(np.max(np.abs(num_gold_arr))))
            passed = bool(np.all(abs_diff_arr <= atol_eff))
        else:
            passed = bool(np.all(rel_diff_arr < regression_tol))


        # --- Accuracy trend analysis ---

        # Normalize relative errors: if the error is a sequence (e.g., an error field over a mesh), 
        # extract the maximum value to represent the 'worst-case' scenario.
        if rel_err_now is not None and isinstance(rel_err_now, collections.abc.Sequence): # if rel_err_now is a list
            err_now_val = float(np.max(rel_err_now))
        else: # if rel_err_now is a scalar
            err_now_val = rel_err_now

        if rel_err_gold is not None and isinstance(rel_err_gold, collections.abc.Sequence): # if rel_err_gold is a list
            err_gold_val = float(np.max(rel_err_gold))
        else: # if rel_err_gold is a scalar
            err_gold_val = rel_err_gold

        # Determine if the current accuracy has improved (lower error) 
        # or worsened (higher error) relative to the GOLD reference.
        # This logic triggers visual indicators ('better'/'worse') in the terminal output.
        improvement = (
            err_now_val is not None and err_gold_val is not None and err_now_val < err_gold_val
        )

        worsening = (
            err_now_val is not None and err_gold_val is not None and err_now_val > err_gold_val
        )

        color = GREEN if passed else RED
        status = "PASS" if passed else "FAIL"

        if improvement:
            extra = f" ({GREEN}better{END})"
        elif worsening:
            extra = f" ({RED}worse{END})"
        else:
            extra = ""

        print(f"  {key:18s} → rel diff = {rel_diff_display:.2e}  → {color}{status}{END}{extra}")
        reg_pass &= passed

    summary_reg = (
        f"{GREEN}PASS{END} Regression within tolerance"
        if reg_pass
        else f"{RED}FAIL{END} Regression mismatch"
    )
    print(f"\n[SUMMARY] {BOLD}{summary_reg}{END}")

    _write_regression_verdict(case_dir, "PASS" if reg_pass else "FAIL")
    return reg_pass
