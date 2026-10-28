# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA CRP on LMFBR core mechanics codes, examples 1 (VP1) and 2 (VP2): the
results of the participants and the values declared by the IAEA, shared by
VP1/vp1_benchmark.py, VP2/vp2_benchmark.py and iwgfr_table.py.

Sources:
  * IWGFR/75 (Part I): Data Table DT1BOW (Y of every participant at core
    top, LRP and URP), Analysis Table AT1BOW (URP average, spread, theory),
    AT1LEN / DT1LEN (axial expansion, theory from the average temperature);
  * IWGFR/76 (Part II): the participant papers (page in ``page``), used
    where DT1BOW is not legible and for the sensitivity studies.

IAEA definitions (IWGFR/75, paras 17-19): spread = (max - min) / average and
accuracy = (code - average) / average, both in %. AVERAGE(10) leaves out
RAINBOW (Mitsubishi), 1.8 % below the average of all 11 (para 19).

Y and dimensions in mm; Z in mm at T_ref = 20 C (stress-free temperature of
the codes, DT1LEN); stations core top 2500, LRP 3000, URP 4000 mm.
"""

import numpy as np

STATIONS = ("Core top", "LRP", "URP")

# DT1BOW: Y (core top, LRP, URP) per participant, example 1 and 2
CODES = {
    "France (SOLO)":          {"page": "76 p. 25-27", "ex1": (1.000, 3.250, 12.250), "ex2": (1.333, 4.338, 16.334)},
    "Italy (SISCO)":          {"page": "76 p. 102",   "ex1": (1.010, 3.280, 12.310), "ex2": (1.350, 4.370, 16.410)},
    "UK (CRAMP)":             {"page": "76 p. 212",   "ex1": (0.989, 3.234, 12.213), "ex2": (1.318, 4.312, 16.284)},
    "CSSR (SANDRA)":          {"page": "76 p. 21",    "ex1": (1.000, 3.250, 12.250), "ex2": (1.333, 4.333, 16.333)},
    "India":                  {"page": "76 p. 66-67", "ex1": (1.002, 3.238, 12.180), "ex2": (1.332, 4.303, 16.190)},
    "Japan PNC (HIBEACON)":   {"page": "76 p. 192-193", "ex1": (1.005, 3.253, 12.249), "ex2": (1.339, 4.333, 16.332)},
    "Japan Mitsubishi (RAINBOW)": {"page": "76 p. 121-122", "ex1": (0.988, 3.195, 12.015), "ex2": (1.318, 4.260, 16.020)},
    "Belgium (STRAW)":        {"page": "76 p. 12-16", "ex1": (1.000, 3.249, 12.245), "ex2": (1.333, 4.332, 16.328)},
    "Germany (DDT)":          {"page": "76 p. 53-54", "ex1": (1.004, 3.265, 12.310), "ex2": (1.339, 4.354, 16.415)},
    "USSR (ACME)":            {"page": "76 p. 203",   "ex1": (1.000, 3.250, 12.250), "ex2": (1.330, 4.330, 16.330)},
    "Japan NAIG (ARKAS)":     {"page": "76 p. 135-136", "ex1": (1.011, 3.276, 12.333), "ex2": (1.348, 4.367, 16.444)},
}
REJECTED = "Japan Mitsubishi (RAINBOW)"   # out of AVERAGE(10)

# Sensitivity studies of the same problem and solutions of another T field (not in the IAEA averages)
SENSITIVITY = {
    "Japan NAIG (ARKAS), shell model": {"page": "76 p. 135-136", "same_case": True,
                                        "ex1": (1.017, 3.287, 12.165), "ex2": (1.356, 4.383, 16.219),
                                        "note": "shell model; mean of the 6 mid-wall corners"},
    "Japan PNC (HIBEACON), standard": {"page": "76 p. 158, 192", "same_case": False,
                                       "ex1": (0.89, 2.91, 11.0), "ex2": (1.20, 3.89, 14.7),
                                       "note": "face-averaged T: other T field"},
}

# Linear theory (DT1BOW / AT1BOW, small strain, 20 C geometry)
THEORY = {"ex1": (1.000, 3.250, 12.250), "ex2": (1.333, 4.333, 16.333)}

# AT1BOW at the URP: (average, spread %) declared by the IAEA
DECLARED_URP = {
    "ex1": {"AVERAGE(10)": (12.259, 1.241), "AVERAGE(ALL=11)": (12.237, 2.599)},
    "ex2": {"AVERAGE(10)": (16.340, 1.554)},
}

# Second-order solutions (DT1BOW, URP): not comparable with a small-strain run
SECOND_ORDER_URP = {"ex1": {"axial + lateral": 12.36, "axial only": 12.47, "lateral only": 12.14},
                    "ex2": {"axial + lateral": 16.480, "axial only": 16.627, "lateral only": 16.187}}

# Axial expansion Z at the URP, T_ref 20 C (AT1LEN, example 1); Belgium is at T_ref 0 C
AXIAL_URP = {"France (SOLO)": 31.97, "Italy (SISCO)": 31.87, "India": 31.95, "CSSR (SANDRA)": 31.86,
             "Germany (DDT)": 31.86, "Japan NAIG (ARKAS)": 31.89}
AXIAL_DECLARED = {"AVERAGE(6)": (31.90, 0.35)}
# DT1LEN theory from the average temperature inlet -> station (462.5 / 441.6 / 425 C), both examples
AXIAL_THEORY = (18.225, 22.77, 31.86)

# Across-flats outer dimension increase above the core, T_ref 20 C (mean T 525 C)
LATERAL = {"ex1": {"Italy (SISCO)": 1.21, "France (SOLO)": 1.208}, "ex2": {"Italy (SISCO)": 1.21}}


def accepted(ex):
    """Y of the codes in AVERAGE(10), (n, 3)."""
    return np.array([c[ex] for name, c in CODES.items() if name != REJECTED])


def band(ex):
    """Per station: average, spread (%), min, max of AVERAGE(10); at the
    URP the average and spread are the IAEA ones, elsewhere computed here."""
    Y = accepted(ex)
    out = {}
    for i, st in enumerate(STATIONS):
        avg = Y[:, i].mean()
        out[st] = {"average": avg, "spread": 100 * np.ptp(Y[:, i]) / avg, "min": Y[:, i].min(),
                   "max": Y[:, i].max(), "declared": False}
    avg, spread = DECLARED_URP[ex]["AVERAGE(10)"]
    out["URP"].update(average=avg, spread=spread, declared=True)
    return out


def outer_flats(d, z):
    """Outer across-flats distance at height z (m), undeformed and deformed:
    outermost nodes on the flats axis (x for corners on +-y)."""
    X, u = d["X"], d["u"]
    zn = X[np.argmin(np.abs(X[:, 2] - z)), 2]
    on = np.abs(X[:, 2] - zn) < 1e-6
    Xs, us = X[on], u[on]
    ic = int(np.argmax(np.abs(Xs[:, :2]).max(axis=0)))   # corner axis
    f = 1 - ic
    line = np.abs(Xs[:, ic]) < 1e-9 * np.abs(Xs[:, f]).max()
    hi, lo = np.argmax(np.where(line, Xs[:, f], -np.inf)), np.argmin(np.where(line, Xs[:, f], np.inf))
    return zn, Xs[hi, f] - Xs[lo, f], (Xs[hi, f] + us[hi, f]) - (Xs[lo, f] + us[lo, f])
