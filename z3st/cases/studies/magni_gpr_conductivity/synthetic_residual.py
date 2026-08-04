#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: analytic residual used to verify the GPR conductivity machinery
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Known-truth log residual for the synthetic GPR verification.

The GPR conductivity hook assimilates a correction on top of the Magni
correlation as a Gaussian process over the log residual

    r = log(k_data / k_magni).

To verify the machinery -- kernel evaluation, de-standardisation, the Newton
tangent dk/dT, the external-operator plumbing -- we do NOT need experimental
data. We need a residual whose value AND derivative are known in closed form,
exactly as `train_knet.py` fits the neural-network hook to the analytic law
k(T) = 1/(a + b*T) rather than to measurements.

This module defines that residual. `make_synthetic_gpr.py` samples it to build
a checkpoint; the non-regression check then asks whether the fitted GPR
reproduces it. Nothing here is experimental data and nothing here is a physical
claim: the correction is a deliberately smooth, bounded fiction whose only
purpose is to have an answer to compare against.

    r(T, Pu, p) = A * sin(pi * (T - T_LO) / (T_HI - T_LO))
                + B * (Pu - PU_REF)
                + C * (p  - P_REF)

Smooth and bounded by construction, so exp(r) stays a ~10% correction and the
conduction operator stays well posed over the whole sampling window.
"""

import numpy as np

# Sampling window. Chosen to cover the Dirichlet spans of the GPR cases with
# margin, so the cases never evaluate the fit outside the fitted region.
T_LO, T_HI = 500.0, 2000.0

# Residual shape. A is the peak log-correction (~10% in k), B and C tilt it
# with composition so the fit has to resolve more than a 1-D function of T.
A = 0.10
B = 0.25
C = 0.40
PU_REF = 0.20
P_REF = 0.05


def residual(T, Pu=0.0, Am=0.0, x=0.0, p=0.0):
    """Analytic log residual r = log(k_true / k_magni)."""
    T = np.asarray(T, dtype=float)
    w = np.pi / (T_HI - T_LO)
    return (A * np.sin(w * (T - T_LO))
            + B * (np.asarray(Pu, dtype=float) - PU_REF)
            + C * (np.asarray(p, dtype=float) - P_REF))


def dresidual_dT(T, Pu=0.0, Am=0.0, x=0.0, p=0.0):
    """dr/dT in closed form, for the tangent check."""
    T = np.asarray(T, dtype=float)
    w = np.pi / (T_HI - T_LO)
    return A * w * np.cos(w * (T - T_LO))


def k_true(T, Pu=0.0, Am=0.0, x=0.0, p=0.0):
    """The conductivity the fitted GPR is expected to reproduce."""
    from z3st.materials.magni_mox_thermal import k_numpy
    return k_numpy(T, Pu=Pu, Am=Am, x=x, p=p) * np.exp(
        residual(T, Pu=Pu, Am=Am, x=x, p=p))


def dk_true_dT(T, Pu=0.0, Am=0.0, x=0.0, p=0.0):
    """d(k_true)/dT in closed form, for the tangent check."""
    from z3st.materials.magni_mox_thermal import dk_dT_numpy, k_numpy
    base = k_numpy(T, Pu=Pu, Am=Am, x=x, p=p)
    dbase = dk_dT_numpy(T, Pu=Pu, Am=Am, x=x, p=p)
    r = residual(T, Pu=Pu, Am=Am, x=x, p=p)
    dr = dresidual_dT(T, Pu=Pu, Am=Am, x=x, p=p)
    return np.exp(r) * (dbase + base * dr)


__all__ = ["residual", "dresidual_dT", "k_true", "dk_true_dT",
           "T_LO", "T_HI", "A", "B", "C", "PU_REF", "P_REF"]
