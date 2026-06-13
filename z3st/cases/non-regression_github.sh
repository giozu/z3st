#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST Non-regression test suite for GitHub
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date: 2026-06-13
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

set -e
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# The CI case list lives in cases/cases_ci.txt
mapfile -t CASES < <(grep -vE '^[[:space:]]*(#|$)' "${ROOT_DIR}/cases/cases_ci.txt" | awk '{print $1}')
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
    ./Allrun > "${case_name//\//_}_log.txt" 2>&1
    exit_code=$?

    # A numerical regression fails CI.
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
