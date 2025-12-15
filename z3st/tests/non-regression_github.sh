#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST Non-Regression test suite for github
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date: 2025-12-12
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

set -e
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
CASES=(
    "box_heated"
    "coaxial_cylinders"
    "cylindrical_shell_thick_GPS"
    "plate"
    "plate_non_linear"
    "thick_cylindrical_thermal_shield"
    "thin_thermal_slab_adiabatic"
    "thin_thermal_slab_with_neumann"
)
SUMMARY_FILE="${ROOT_DIR}/tests/non-regression_summary.txt"

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
    case_dir="${ROOT_DIR}/tests/${case_name}"

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

    if [[ $exit_code -ne 0 ]]; then
        echo "Case ${case_name} FAILED (exit code ${exit_code})"
        echo "Case: $case_name -> FAIL" >> "$SUMMARY_FILE"
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
