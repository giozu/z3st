"""
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""

from importlib import metadata

try:
    __version__ = metadata.version("z3st")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

from .core import solver
from .models import mechanical_model, thermal_model, gap_model
from .utils import export_vtu, plot_convergence

__all__ = [
    "solver",
    "thermal_model",
    "mechanical_model",
    "gap_model",
    "export_vtu",
    "plot_convergence",
]

Solver = solver.Solver
Thermal = thermal_model.ThermalModel
Mechanical = mechanical_model.MechanicalModel
