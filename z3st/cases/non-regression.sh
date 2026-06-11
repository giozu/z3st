#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Performing Z3ST non-regression
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date  : 2025-10-13
#
# SUPERSEDED (2026-06-11) by non-regression_local.sh, which discovers
# its case set (every dir with Allrun + output/non-regression_gold.json,
# minus suite_exclude.txt) instead of this hand-maintained list, and
# exits non-zero on failures. Kept temporarily for comparison; delete
# once non-regression_local.sh has earned trust.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SUMMARY_FILE="${ROOT_DIR}/non-regression_summary.txt"

CASES=(
    "verification/mechanics/uniaxial_tension"
    "verification/thermal/thin_slab_dirichlet_2D"
    "verification/thermal/thin_slab_neumann_2D"
    "verification/thermal/thin_slab_neumann_3D"
    "verification/mechanics/uniaxial_tension_nonlinear"
    "verification/thermal/thin_cylindrical_shell_dirichlet_2D"
    "verification/thermal/thick_slab_adiabatic_3D"
    "verification/thermal/thin_slab_adiabatic_2D"
    "verification/thermal/thin_cylindrical_shell_adiabatic_2D"
    "verification/thermal/thick_cylindrical_shell_adiabatic_2D"
    # "verification/thermal/thick_slab_non_adiabatic_3D"
    "verification/thermal/thin_slab_non_adiabatic_2D"
    "verification/thermal/thin_cylindrical_shell_non_adiabatic_2D"
    "verification/thermal/thick_cylindrical_shell_non_adiabatic_2D"
    "verification/thermal/box_heated"
    "verification/mechanics/lame_plane_strain_2D"
    "verification/mechanics/lame_gps_2D"
    # "verification/mechanics/lame_gps_3D"
    "verification/mechanics/mariotte_thin_shell"
    "verification/mechanics/thermal_gradient_2D"
    # "verification/mechanics/thermal_gradient_3D"
    "verification/mechanics/annular_cylinder"
    "verification/mechanics/full_cylinder"
    "verification/mechanics/spherical_cavity"
    "regression/coaxial_gap_3D"
    # The following are present on disk but currently lack a
    # non-regression_gold.json (CASES-P1-6 / CASES-FOLLOWUP-5).
    # Uncomment after blessing a gold for each.
    # "regression/elliptical_cavity_2D"
    # "regression/two_elliptical_cavities_2D"
    # "benchmarks/double_crack_2D"
    "verification/plasticity/stress_strain_displacement"
    # "benchmarks/notched_plate_2D"
    "verification/plasticity/stress_strain_stress"
    # "regression/box_notch_2D"
    "regression/box_crack_2D"
    # "benchmarks/sen_tension"
    # "benchmarks/sen_shear"
    "verification/plasticity/j2_hardening_2D"
    "verification/plasticity/crystal_single_grain"
    "verification/fuel/swelling"
    "verification/fuel/fuel_swelling"
    "verification/fuel/burnup"
    "verification/fuel/axial_power"
    "verification/fuel/axial_table"
    "verification/fuel/creep"
    "verification/fuel/creep_relaxation"
    "verification/fuel/coaxial_contact"
    "verification/fuel/creep_law_discovery"
    "verification/thermal/spherical_shell"
    # "regression/pwr_rod_2D"
    "studies/mesh_sensitivity_2D"
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