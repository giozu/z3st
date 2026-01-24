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
    "00_example"
    "1_thin_slab_2D"
    "1_thin_slab_neumann_3D"
    "2_thin_cylindrical_shell_2D"
    "3_thick_slab_adiabatic_3D"
    "3_thin_slab_adiabatic_2D"
    "4_thin_cylindrical_shell_adiabatic_2D"
    "4_thick_cylindrical_shell_adiabatic_2D"
    "5_thick_slab_non_adiabatic_3D"
    "5_thin_slab_non_adiabatic_2D"
    "6_thin_cylindrical_shell_non_adiabatic_2D"
    "6_thick_cylindrical_shell_non_adiabatic_2D"
    "7_box_heated"
    "8_thick_cylindrical_shell_plane_strain_2D"
    "9_thick_cylindrical_shell_GPS_2D"
    "9_thick_cylindrical_shell_GPS_3D"
    "11_thin_cylindrical_shell_Mariotte"
    "12_cylindrical_shell_thermal_gradient_2D"
    "12_cylindrical_shell_thermal_gradient_3D"
    "13_annular_cylinder"
    "14_full_cylinder"
    "20_coaxial_cylinders_3D"
    "21_plate_non_linear"
    "mesh_sensitivity_2D"
    "stress_strain_curve_stress"
    "stress_strain_curve_displacement"
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
