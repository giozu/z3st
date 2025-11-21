#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Performing Z3ST non-regression
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date  : 2025-10-13
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SUMMARY_FILE="${ROOT_DIR}/non-regression_summary.txt"

CASES=(
    "1_thin_thermal_slab"
    "2_thin_cylindrical_thermal_shield"
    "3_thin_thermal_slab_adiabatic"
    "4_thin_cylindrical_thermal_shield_adiabatic"
    "5_thin_thermal_slab_non_adiabatic"
    "6_thin_cylindrical_thermal_shield_non_adiabatic"
    "7_box_heated"
    "8_cylindrical_shell_thick_plane_strain"
    "9_cylindrical_shell_thick_GPS"
    "10_cylindrical_shell_thick_plane_stress"
    "11_cylindrical_shell_thin"
    "12_thick_cylindrical_thermal_shield"
    "13_thick_linear_thermal_shield_adiabatic"
    "14_thick_cylindrical_thermal_shield_adiabatic"
    "15_thick_linear_thermal_shield_non_adiabatic"
    "16_thick_cylindrical_thermal_shield_non_adiabatic"
    "17_annular_cylinder"
    "18_full_cylinder"
    "19_plate"
    "20_coaxial_cylinders"
)

echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo "Z3ST - Automated non-regression suite"
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo "Root directory : ${ROOT_DIR}"
echo "Non-regression cases:"
printf '  - %s\n' "${CASES[@]}"
echo ""

echo "Z3ST non-regression summary - $(date)" > "$SUMMARY_FILE"
echo "--.. ..- .-.. .-.. -----.. ..- .-.. .-.. -----.. ..- .-.. .-.. ---" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Loop through selected non-regression cases
for case_name in "${CASES[@]}"; do
    case_dir="${ROOT_DIR}/${case_name}"

    if [[ ! -d "$case_dir" ]]; then
        echo "(!)  Skipping ${case_name}: directory not found."
        echo "Case: $case_name (NOT FOUND)" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
        continue
    fi

    if [[ ! -f "${case_dir}/Allrun" ]]; then
        echo "(!)  Skipping ${case_name}: missing Allrun script."
        echo "Case: $case_name (NO Allrun)" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
        continue
    fi

    echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
    echo "Running case: $case_name"
    echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
    echo "Directory: $case_dir"

    cd "$case_dir" || continue

    ./Allclean

    start_time=$(date +%s)

    chmod +x Allrun
    ./Allrun

    exit_code=${PIPESTATUS[0]}
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    if [[ $exit_code -eq 0 ]]; then
        status="OK"
    else
        status="FAIL"
    fi

    {
        echo "Case: $case_name"
        printf "  %-15s : %s\n" "Run" "$status"

        if [[ -f "output/non-regression.json" ]]; then
            summary_status=$(grep -o '"summary": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
            printf "  %-15s : %s\n" "Non-regression" "$summary_status"
        else
            printf "  %-15s : %s\n" "Non-regression" "(no non-regression.json found)"
        fi

        printf "  %-15s : %ss\n" "Time" "$elapsed"
        echo ""
    } >> "$SUMMARY_FILE"

    echo "Case: $case_name"
    printf "  %-15s : %s\n" "Run" "$status"
    printf "  %-15s : %s\n" "Non-regression" "$summary_status"
    printf "  %-15s : %ss\n" "Time" "$elapsed"
    echo ""

    cd "$ROOT_DIR" || exit
done

echo ""
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo " All selected cases completed."
echo " Summary written to: ${SUMMARY_FILE}"
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo ""
