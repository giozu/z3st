#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST Non-regression test suite for github
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date: 2025-12-12
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

set -e
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# CI case-inclusion policy (see CONTEXT.md §6 / punch_list.md CASES-P1-5):
# This list is a tight subset of the local non-regression.sh suite, chosen
# for fast CI turnaround under the dolfinx/dolfinx:stable container budget.
# Each case here:
#   - has a stable output/non-regression_gold.json,
#   - runs in well under one minute,
#   - exercises one physics path that's not covered by another listed case
#     (linear elasticity, axisymmetric thermal, 3D thermal, volumetric
#     heating, multi-material with Robin pair gap, 2D damage, J2 plasticity,
#     custom-plasticity hook).
# Active WIP cases are intentionally excluded, they are too long for CI 
# and their golds are still being blessed. Run them locally via cases/non-regression.sh.
CASES=(
    "1_thin_slab_2D"
    "2_thin_cylindrical_shell_2D"
    "3_thick_slab_adiabatic_3D"
    "7_box_heated"
    "14_full_cylinder"
    "16_coaxial_cylinders_3D"
    "18_box_crack_2D"
    "20_plasticity_2D"
    "demo_CP_single_grain"
    "V_coaxial_contact_verification"
    "V_swelling_verification"
    "V_burnup_verification"
    "V_fuel_swelling_verification"
    "V_axial_power_verification"
    "V_axial_table_verification"
    "V_creep_verification"
    "V_creep_relaxation_verification"
)
SUMMARY_FILE="${ROOT_DIR}/cases/non-regression_GH_summary.txt"

echo "Running minimal non-regression suite..."
echo "---------------------------------------"
echo "Root directory: ${ROOT_DIR}"
echo "Selected cases:"
printf '  - %s\n' "${CASES[@]}"
echo ""

echo "Z3ST non-regression summary - $(date)" > "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

global_status=0

for case_name in "${CASES[@]}"; do
    case_dir="${ROOT_DIR}/cases/${case_name}"

    echo "-------------------------------------------------"
    echo "Running case: ${case_name}"
    echo "-------------------------------------------------"

    if [[ ! -d "$case_dir" ]]; then
        echo "(!) Skipping ${case_name}: directory not found."
        echo "Case: $case_name (NOT FOUND)" >> "$SUMMARY_FILE"
        global_status=1
        continue
    fi

    cd "$case_dir"
    chmod +x Allrun || true
    chmod +x Allclean || true

    ./Allclean || true
    ./Allrun > "${case_name}_log.txt" 2>&1
    exit_code=$?

    # Allrun's exit code only catches crashes: the per-case non-regression.py
    # always exits 0. Gate on the verdicts it writes into non-regression.json
    # ("summary" = analytic tolerance check, "regression" = vs GOLD reference)
    # so a numerical regression actually fails CI.
    fail_reason=""
    if [[ $exit_code -eq 0 && -f "output/non-regression.json" ]]; then
        summary_status=$(grep -o '"summary": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
        regression_status=$(grep -o '"regression": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
        if [[ "$summary_status" == "FAIL" ]]; then
            exit_code=1
            fail_reason=" (analytic tolerance)"
        elif [[ "$regression_status" == "FAIL" ]]; then
            exit_code=1
            fail_reason=" (regression vs gold)"
        fi
    fi

    if [[ $exit_code -ne 0 ]]; then
        echo "Case ${case_name} FAILED (exit code ${exit_code})${fail_reason}"
        echo "Case: $case_name -> FAIL${fail_reason}" >> "$SUMMARY_FILE"
        global_status=1
    else
        echo "Case ${case_name} PASSED"
        echo "Case: $case_name -> OK" >> "$SUMMARY_FILE"
    fi

    cd "$ROOT_DIR"
done

echo ""
echo "---------------------------------------"
echo "Non-regression summary written to: ${SUMMARY_FILE}"
echo "---------------------------------------"
cat "$SUMMARY_FILE"

exit $global_status
