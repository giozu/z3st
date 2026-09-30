# ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---
# Z3ST: 15-15Ti thermal conductivity, case-local
# Author: Bianca Funaro
# ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. --- ..-. ..- -. .- .-. ---


def k(T):
    """Thermal conductivity (W/m·K), T in K. Homework 2024-2025, NDT.

    Evaluated on the previous step's temperature field.
    """
    return 13.95 + 0.01163 * (T - 273.15)
