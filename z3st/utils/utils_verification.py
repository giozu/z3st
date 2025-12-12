# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import json
import os
import sys

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


def regression_check(errors, case_dir, regression_tol=1e-3):
    """
    Compare numerical results with a non-regression_gold.json reference file.

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

    with open(gold_file, "r") as f:
        gold_data = json.load(f)

    gold_results = gold_data.get("results", gold_data)
    reg_pass = True

    for key in errors:
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

        rel_diff_arr = np.abs(num_now_arr - num_gold_arr) / np.maximum(np.abs(num_gold_arr), 1e-12)

        passed = np.all(rel_diff_arr < regression_tol)
        rel_diff_display = float(np.max(rel_diff_arr))

        improvement = (
            rel_err_now is not None and rel_err_gold is not None and rel_err_now < rel_err_gold
        )

        worsening = (
            rel_err_now is not None and rel_err_gold is not None and rel_err_now > rel_err_gold
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

    return reg_pass
