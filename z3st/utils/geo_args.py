# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: print gmsh -setnumber / -setstring flags from geometry.yaml.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
#     gmsh $(python3 -m z3st.utils.geo_args) mesh.geo -2
#
# Only top-level scalars are passed; nested blocks (labels, ...) are skipped.
# The .geo picks them up through `If (!Exists(X)) X = default; EndIf`.

import shlex
import sys

import yaml


def geo_args(path="geometry.yaml"):
    with open(path) as f:
        geometry = yaml.safe_load(f) or {}
    args = []
    for key, value in geometry.items():
        if isinstance(value, bool) or key == "name":
            continue
        if isinstance(value, (int, float)):
            args += ["-setnumber", key, repr(value)]
        elif isinstance(value, str):
            args += ["-setstring", key, value]
    return args


if __name__ == "__main__":
    print(shlex.join(geo_args(*sys.argv[1:])))
