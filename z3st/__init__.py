"""
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.2.0 (2026)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""

from importlib import metadata

try:
    __version__ = metadata.version("z3st")
except metadata.PackageNotFoundError:
    __version__ = "0.2.0"

def __getattr__(name):
    """Lazy import of heavy modules (dolfinx-dependent) to allow
    lightweight utilities (e.g., utils_extract_vtu) to be imported
    without requiring MPI / dolfinx at import time."""

    _lazy_map = {
        "solver": (".core", "solver"),
        "Spine": (".core.spine", "Spine"),
        "Solver": (".core.solver", "Solver"),
        "thermal_model": (".models", "thermal_model"),
        "Thermal": (".models.thermal_model", "ThermalModel"),
        "mechanical_model": (".models", "mechanical_model"),
        "Mechanical": (".models.mechanical_model", "MechanicalModel"),
        "gap_model": (".models", "gap_model"),
        "cluster_dynamic_model": (".models", "cluster_dynamic_model"),
        "Cluster": (".models.cluster_dynamic_model", "ClusterDynamicsModel"),
        "plot_convergence": (".utils", "plot_convergence"),
    }

    if name in _lazy_map:
        module_path, attr = _lazy_map[name]
        import importlib
        mod = importlib.import_module(module_path, __name__)
        obj = getattr(mod, attr)
        globals()[name] = obj
        return obj

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "solver",
    "Spine",
    "thermal_model",
    "mechanical_model",
    "gap_model",
    "cluster_dynamic_model",
    "plot_convergence",
]
