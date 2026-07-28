#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Single source of truth for the post-processing scripts of this case.

Every geometric, material and solver constant used by ``plots.py`` and
``pressure_creep_analysis.py`` is read here from the case's own YAML files.
Nothing is restated as a literal: an analytical overlay computed with
parameters that differ from the ones the solver actually ran is not a
verification, and hand-copied constants silently drift the moment the case is
retuned.

The only literal below is the universal gas constant.
"""

import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    with open(os.path.join(HERE, name), "r") as fh:
        return yaml.safe_load(fh)


_input = _load("input.yaml")
_geom = _load("geometry.yaml")
_clad = _load(_input["materials"]["clad"])
_fuel = _load(_input["materials"]["fuel"])

# --- Physical constant ------------------------------------------------------
R_GAS = 8.3145              # (J/(mol.K))

# --- Geometry (m) -----------------------------------------------------------
R_PELLET = float(_geom["outer_radius_1"])
R_CLAD_I = float(_geom["inner_radius_2"])
R_CLAD_O = float(_geom["outer_radius_2"])
G0 = R_CLAD_I - R_PELLET    # cold radial gap

# --- Contact ----------------------------------------------------------------
_contact = _input["models"]["contact"]
K_PEN = float(_contact["penalty_stiffness"])        # (Pa/m)
INITIAL_GAP = float(_contact["initial_gap"])        # (m)

# --- Time history -----------------------------------------------------------
_time = [float(t) for t in _input["time"]]
_n_steps = _input["n_steps"]
if not isinstance(_n_steps, list):
    _n_steps = [_n_steps]
T_TOTAL = _time[-1]                                 # (s)
NB_STEPS = int(sum(int(n) for n in _n_steps))

# --- Cladding (the creeping body: Esposito's hub) ---------------------------
E_CLAD = float(_clad["E"])                          # (Pa)
NU_CLAD = float(_clad["nu"])                        # (-)
ALPHA_CLAD = float(_clad["alpha"])                  # (1/K)
G = E_CLAD / (2.0 * (1.0 + NU_CLAD))                # shear modulus (Pa)

CREEP_A0 = float(_clad["creep_A0"])                 # (Pa^-n s^-1)
CREEP_N = float(_clad["creep_n"])                   # (-)
CREEP_Q = float(_clad["creep_Q"])                   # (J/mol)

# Irradiation creep must stay off for the Esposito comparison (see clad.yaml).
CREEP_IRR_B = float(_clad.get("creep_irr_B", 0.0))
FAST_FLUX = float(_clad.get("fast_flux", 0.0))

# --- Pellet (Esposito's shaft) ----------------------------------------------
E_FUEL = float(_fuel["E"])                          # (Pa)
NU_FUEL = float(_fuel["nu"])                        # (-)

# --- Backwards-compatible aliases used by the scripts -----------------------
E_EL = E_CLAD
NU = NU_CLAD
ALPHA = ALPHA_CLAD

# Reference stress for the uniaxial relaxation cross-check: an imposed total
# strain of 5e-5 m over the clad thickness, held constant.
SIGMA0 = E_CLAD * 5.0e-5 / (R_CLAD_I - R_CLAD_O)


def check_consistency():
    """Fail loudly on the inconsistencies that silently invalidate the overlay.

    Returns the list of problems found (empty when the case is coherent) so a
    caller can warn rather than abort if that is more useful.
    """
    problems = []

    if abs(G0 - INITIAL_GAP) > 1e-12:
        problems.append(
            f"cold gap from geometry.yaml ({G0:.3e} m) disagrees with "
            f"models.contact.initial_gap in input.yaml ({INITIAL_GAP:.3e} m)"
        )

    if CREEP_IRR_B * FAST_FLUX != 0.0:
        problems.append(
            "irradiation creep is active (creep_irr_B * fast_flux != 0); "
            "Esposito eq. (21) is a pure Norton law and the analytical "
            "overlay is not comparable with an in-pile term switched on"
        )

    if str(_clad.get("creep", "")).lower() != "norton":
        problems.append(f"clad creep law is '{_clad.get('creep')}', expected 'norton'")

    return problems


def report():
    """Print the resolved parameter set, then the consistency verdict."""
    print("[case_params] resolved from the case YAML files:")
    for key in (
        "R_PELLET", "R_CLAD_I", "R_CLAD_O", "G0", "K_PEN", "INITIAL_GAP",
        "T_TOTAL", "NB_STEPS", "E_CLAD", "NU_CLAD", "E_FUEL", "NU_FUEL",
        "CREEP_A0", "CREEP_N", "CREEP_Q", "CREEP_IRR_B", "FAST_FLUX",
    ):
        print(f"  {key:<12} = {globals()[key]}")

    problems = check_consistency()
    if problems:
        print("[case_params] INCONSISTENT:")
        for p in problems:
            print(f"  [WARNING] {p}")
    else:
        print("[case_params] consistency check: OK")
    return problems


if __name__ == "__main__":
    report()
