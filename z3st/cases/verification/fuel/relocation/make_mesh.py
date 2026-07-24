# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Regenerate the fragmented pellet + cladding mesh for the relocation case.
# A single 2D (plane-strain r-theta) section in SI metres, written to mesh.msh.
# Fuel fragments carry physical group 'fuel' (tag 1); the cladding ring is
# 'cladding' (tag 2). Boundary curves: 'outer' (10, fuel fragment outer arcs),
# 'cracks' (50, fragment-fragment faces), 'clad_inner' (60), 'clad_outer' (61).
#
#   python make_mesh.py        # writes mesh.msh in this directory

import sys
import shutil
from pathlib import Path

# the parametric generator lives in the shared utils
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "utils" / "geo_files"))
import pellet_fragmented as pf  # noqa: E402

HERE = Path(__file__).resolve().parent

# --. single-case parameters --..
# print_scale converts the generator's millimetre defaults to metres:
#   pellet radius 5 mm, fuel-clad gap 0.1 mm, crack aperture 0.1 mm.
PARAMS = dict(pf.DEFAULTS)
PARAMS.update({
    "SEED": 24,
    "N_main": 4,          # main radial cracks (branching adds a fragment -> 5)
    "circ_radii": [],     # no circumferential cracks for the first case
    "curl_amp": 0.22,
    "branch_prob": 0.35,
    "hub_offset": 0.15,
    "angle_jitter": 0.35,
    "print_scale": 1e-3,  # mm -> m
})


def main():
    build = HERE / "_mesh_build"
    stats = pf.generate_pellet(PARAMS, build)
    shutil.copy(build / "pellet.msh", HERE / "mesh.msh")
    shutil.copy(build / "top_view.png", HERE / "top_view.png")
    shutil.rmtree(build, ignore_errors=True)
    print(f"wrote {HERE / 'mesh.msh'}  ({stats['n_fragments']} fragments)")


if __name__ == "__main__":
    main()
