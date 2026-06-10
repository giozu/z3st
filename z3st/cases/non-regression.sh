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
    "1_thin_slab_neumann_2D"
    "1_thin_slab_neumann_3D"
    "1_thin_slab_non_linear"
    "2_thin_cylindrical_shell_2D"
    "3_thick_slab_adiabatic_3D"
    "3_thin_slab_adiabatic_2D"
    "4_thin_cylindrical_shell_adiabatic_2D"
    "4_thick_cylindrical_shell_adiabatic_2D"
    # "5_thick_slab_non_adiabatic_3D"
    "5_thin_slab_non_adiabatic_2D"
    "6_thin_cylindrical_shell_non_adiabatic_2D"
    "6_thick_cylindrical_shell_non_adiabatic_2D"
    "7_box_heated"
    "8_thick_cylindrical_shell_plane_strain_2D"
    "9_thick_cylindrical_shell_GPS_2D"
    # "9_thick_cylindrical_shell_GPS_3D"
    "11_thin_cylindrical_shell_Mariotte"
    "12_cylindrical_shell_thermal_gradient_2D"
    # "12_cylindrical_shell_thermal_gradient_3D"
    "13_annular_cylinder"
    "14_full_cylinder"
    "15_spherical_pressurized_cavity"
    "16_coaxial_cylinders_3D"
    # The following are present on disk but currently lack a
    # non-regression_gold.json (CASES-P1-6 / CASES-FOLLOWUP-5).
    # Uncomment after blessing a gold for each.
    # "15_single_elliptical_cavity_2D"
    # "15_two_elliptical_cavities_2D"
    # "17_stress_strain_curve_double_crack"
    "17_stress_strain_curve_displacement"
    # "17_stress_strain_curve_knotch"
    "17_stress_strain_curve_stress"
    # "18_box_knotch_2D"
    "18_box_crack_2D"
    # "19_single-edge_notched_tension_test"
    # "19_single-edge_notched_shear_test"
    "20_plasticity_2D"
    "demo_CP_single_grain"
    "V_swelling_verification"
    "V_fuel_swelling_verification"
    "V_burnup_verification"
    "V_coaxial_contact_verification"
    "U_pwr_rod_2D"
    "I_mesh_sensitivity_2D"
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
    # conda run -n z3st ./Allrun
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
            regression_status=$(grep -o '"regression": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
            printf "  %-15s : %s\n" "Non-regression" "$summary_status"
            printf "  %-15s : %s\n" "Gold regression" "${regression_status:-(no gold)}"
        else
            summary_status="(no non-regression.json — NOT VALIDATED)"
            regression_status=""
            printf "  %-15s : %s\n" "Non-regression" "$summary_status"
        fi

        printf "  %-15s : %ss\n" "Time" "$elapsed"
        echo ""
    } >> "$SUMMARY_FILE"

    echo "Case: $case_name"
    printf "  %-15s : %s\n" "Run" "$status"
    printf "  %-15s : %s\n" "Non-regression" "$summary_status"
    printf "  %-15s : %s\n" "Gold regression" "${regression_status:-(no gold)}"
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