# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


def k(T):
    """Thermal conductivity (W/m·K), T in K. Homework 2024-2025, NDT.

    Evaluated on the previous step's temperature field.
    """
    return 13.95 + 0.01163 * (T - 273.15)
