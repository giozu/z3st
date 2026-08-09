# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: shared case runner, sourced by every case Allrun.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# A vanilla case Allrun is three lines:
#
#     #!/bin/bash
#     DIM=2
#     . "$(python3 -c 'import z3st,os;print(os.path.dirname(z3st.__file__))')/utils/allrun.sh"
#
# Knobs, all optional:
#   DIM   gmsh dimension (1, 2 or 3)          default 2
#   POST  extra post-processing scripts, e.g. POST="plots.py diagnostics.py"
#
# $0 inside a sourced script is still the *caller* (the case Allrun), so the
# cd below lands in the case directory no matter where the suite invokes it.

set -e
cd "${0%/*}" || exit 1

gmsh mesh.geo -"${DIM:-2}" > log_mesh.md
python3 -m z3st > log_z3st.md
python3 non-regression.py

for script in ${POST:-}; do
    python3 "$script"
done

# Never fails the case: main() swallows its own errors, and a case with no
# parsable residuals simply gets no plot.
python3 -m z3st.utils.plot_convergence log_z3st.md
