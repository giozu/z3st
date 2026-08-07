# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: shared case cleaner, sourced by every case Allclean.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# A case Allclean is two lines:
#
#     #!/bin/bash
#     . "$(python3 -c 'import z3st,os;print(os.path.dirname(z3st.__file__))')/utils/allclean.sh"
#
# Never "rm -rf output": the blessed output/non-regression_gold.json lives
# there, the suite runs Allclean before every case, and wiping the directory
# would disable that case's gold-regression check. Deleting by extension keeps
# the gold by construction rather than by everyone remembering to.

cd "${0%/*}" || exit 1

echo "Cleaning current directory..."

rm -f output/non-regression.json output/case_summary.json
rm -f output/conductivity_impact.json
rm -f output/*.vtu output/*.pvd output/*.xdmf output/*.h5
rm -f output/*.png output/*.csv
rm -f *.msh *.png log_*.md energies.txt input_check.txt
rm -rf __pycache__

echo "Done. Output and temporary files deleted (gold record preserved)."
