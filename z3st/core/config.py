# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


class Config:
    """
    Configuration manager for Z3ST simulations.

    This class parses the user input YAML file and initializes the
    global configuration used by all other modules (mesh, solvers, models).

    It loads:
      * active physical models (thermal, mechanical, gap conductance)
      * solver settings (linear/non-linear, tolerances, coupling scheme)
      * paths for geometry, mesh, and boundary conditions
      * number of time steps

    Attributes
    ----------
    input_file : dict
        Dictionary parsed from the YAML input file.
    on : dict
        Dictionary of active physical models (`thermal`, `mechanical`).
    gap_model : str
        Selected gap conductance model (e.g., "Fixed" or "Gas").
    h_gap_value : float
        Gap conductance coefficient (W/m²·K).
    geometry_path : str
        Path to the geometry YAML file.
    mesh_path : str
        Path to the mesh file (.msh or .py).
    boundary_conditions_path : str
        Path to the YAML file containing boundary condition definitions.
    n_steps : int
        Number of time steps.
    coupling : str
        Coupling strategy ("monolithic" or "staggered").
    mech_solver_type : str
        Type of mechanical solver ("linear" or "non-linear").
    mech_linear_solver : str
        Linear solver for the mechanical problem (e.g., "iterative_amg").
    rtol_mech : float
        Relative tolerance for mechanical solver.
    stag_tol_mech : float
        Staggered tolerance for the mechanical solver.
    mech_regime : str
        Mechanical regime ("3D" or "plane_stress").
    mech_convergence : str
        Convergence criterion ("norm" or "rel_norm").
    mech_debug : bool
        Debug flag for the mechanical model.
    th_solver_type : str
        Type of thermal solver ("linear" or "non-linear").
    th_linear_solver : str
        Linear solver for the thermal problem (e.g., "iterative_amg").
    rtol_th : float
        Relative tolerance for thermal solver.
    stag_tol_th : float
        Staggered tolerance for the thermal solver.
    th_convergence : str
        Convergence criterion for thermal solver.

    Notes
    -----
    This class is initialized at the start of each simulation.
    It centralizes configuration parsing so that all modules
    (thermal, mechanical, solver, gap model) share consistent parameters.
    """

    def __init__(self, input_file):
        """
        Initialize the Config class from a YAML input dictionary.

        Parameters
        ----------
        input_file : dict
            Parsed YAML configuration file containing all simulation parameters.
        """
        print("__Config initializer__")

        self.input_file = input_file

        # --. Active models --..
        models = self.input_file.get("models", {})
        self.on = {
            "thermal": models.get("thermal", False),
            "mechanical": models.get("mechanical", False),
        }

        # --. Gap conductance --..
        gap_config = self.input_file.get("models", {}).get("gap_conductance", {})
        self.gap_model = gap_config.get("type", None)
        self.h_gap_value = gap_config.get("value", 0.0)

        # --. Paths --..
        self.geometry_path = self.input_file.get("geometry_path", None)
        self.mesh_path = self.input_file.get("mesh_path", None)
        self.boundary_conditions_path = self.input_file.get("boundary_conditions_path", None)
        self.n_steps = self.input_file.get("n_steps", 10)

        # --. Debug print --..
        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")

        print(f"  → Gap conductance     : {self.gap_model} (value = {self.h_gap_value})")
        print("\n")
