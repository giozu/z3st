#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST local non-regression suite
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Author: Giovanni Zullo
# Date  : 2026-06-11
#
# Membership is discovered, not hand-maintained: every directory under
# cases/ that has BOTH an Allrun and a blessed output/non-regression_gold.json
# is part of the local suite. sandbox/ is never scanned (explicitly
# unprotected work lives there). Exceptions go in suite_exclude.txt,
# one case per line, with a reason as a trailing comment.
#
# Usage:
#   ./non-regression_local.sh            run the discovered suite
#   ./non-regression_local.sh --list     print the discovered case set and exit
#   ./non-regression_local.sh CASE...    run only the named cases (paths
#                                        relative to cases/, discovery rules
#                                        and exclusions still apply)
#
# A case fails if Allrun exits non-zero, if output/non-regression.json is
# missing afterwards, or if its "summary" or "regression" verdict is FAIL.
# The script exits non-zero if any case failed.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
EXCLUDE_FILE="${ROOT_DIR}/suite_exclude.txt"
SUMMARY_FILE="${ROOT_DIR}/non-regression_summary.txt"

# --.. ..- .-.. .-.. --- discovery --.. ..- .-.. .-.. ---
mapfile -t DISCOVERED < <(
    find "$ROOT_DIR" -type d -name sandbox -prune -o -type f -name Allrun -printf '%h\n' \
    | while read -r d; do
          [[ -f "$d/output/non-regression_gold.json" ]] && echo "${d#"$ROOT_DIR"/}"
      done | sort
)

# exclusions (strip blank lines and comments; first whitespace-separated field)
declare -A EXCLUDED
if [[ -f "$EXCLUDE_FILE" ]]; then
    while read -r entry _; do
        [[ -z "$entry" || "$entry" == \#* ]] && continue
        EXCLUDED["$entry"]=1
    done < "$EXCLUDE_FILE"
fi

CASES=()
SKIPPED=()
for c in "${DISCOVERED[@]}"; do
    if [[ -n "${EXCLUDED[$c]:-}" ]]; then SKIPPED+=("$c"); else CASES+=("$c"); fi
done

# optional positional filter
if [[ $# -gt 0 && "$1" != "--list" ]]; then
    FILTERED=()
    for want in "$@"; do
        hit=""
        for c in "${CASES[@]}"; do [[ "$c" == "$want" ]] && { FILTERED+=("$c"); hit=1; }; done
        [[ -z "$hit" ]] && echo "(!) '$want' is not in the discovered suite (no Allrun+gold, excluded, or misspelled)."
    done
    CASES=("${FILTERED[@]}")
fi

if [[ "${1:-}" == "--list" ]]; then
    echo "Discovered suite (${#CASES[@]} cases):"
    printf '  %s\n' "${CASES[@]}"
    if [[ ${#SKIPPED[@]} -gt 0 ]]; then
        echo "Excluded by ${EXCLUDE_FILE##*/} (${#SKIPPED[@]}):"
        printf '  %s\n' "${SKIPPED[@]}"
    fi
    exit 0
fi

echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo "Z3ST - Automated non-regression suite (discovery-based)"
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
echo "Root directory : ${ROOT_DIR}"
echo "Cases (${#CASES[@]}):"
printf '  - %s\n' "${CASES[@]}"
if [[ ${#SKIPPED[@]} -gt 0 ]]; then
    echo "Excluded (${#SKIPPED[@]}):"
    printf '  - %s\n' "${SKIPPED[@]}"
fi
echo ""

echo "Z3ST non-regression summary - $(date)" > "$SUMMARY_FILE"
echo "--.. ..- .-.. .-.. -----.. ..- .-.. .-.. -----.. ..- .-.. .-.. ---" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

global_status=0

for case_name in "${CASES[@]}"; do
    case_dir="${ROOT_DIR}/${case_name}"

    echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
    echo "Running case: $case_name"
    echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"

    cd "$case_dir" || { global_status=1; continue; }

    chmod +x Allrun Allclean 2>/dev/null
    [[ -x ./Allclean ]] && ./Allclean

    start_time=$(date +%s)
    ./Allrun
    exit_code=$?
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    run_status="OK"
    [[ $exit_code -ne 0 ]] && run_status="FAIL (exit ${exit_code})"

    if [[ -f "output/non-regression.json" ]]; then
        summary_status=$(grep -o '"summary": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
        regression_status=$(grep -o '"regression": *"[^"]*"' output/non-regression.json | head -1 | cut -d'"' -f4)
    else
        summary_status="MISSING non-regression.json"
        regression_status=""
    fi

    case_status="OK"
    if [[ $exit_code -ne 0 || "$summary_status" != "PASS" || "$regression_status" == "FAIL" ]]; then
        case_status="FAIL"
        global_status=1
    fi

    {
        echo "Case: $case_name [${case_status}]"
        printf "  %-15s : %s\n" "Run" "$run_status"
        printf "  %-15s : %s\n" "Non-regression" "$summary_status"
        printf "  %-15s : %s\n" "Gold regression" "${regression_status:-(no verdict)}"
        printf "  %-15s : %ss\n" "Time" "$elapsed"
        echo ""
    } | tee -a "$SUMMARY_FILE"

    cd "$ROOT_DIR" || exit 1
done

echo ""
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"
if [[ $global_status -eq 0 ]]; then
    echo " All ${#CASES[@]} cases PASSED."
else
    echo " FAILURES detected - see ${SUMMARY_FILE}"
fi
echo "--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---"

exit $global_status
