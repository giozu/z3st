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
            "damage": models.get("damage", False),
        }
        gap_config = self.input_file.get("models", {}).get("gap_conductance", {})
        self.gap_model = gap_config.get("type", None)
        self.h_gap_value = gap_config.get("value", 0.0)

        # --. Paths --..
        self.geometry_path = self.input_file.get("geometry_path", None)
        self.mesh_path = self.input_file.get("mesh_path", None)
        self.boundary_conditions_path = self.input_file.get("boundary_conditions_path", None)
        self.n_steps = self.input_file.get("n_steps", 10)

        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")
        print(f"  → Gap conductance     : {self.gap_model} (value = {self.h_gap_value})")
        print("\n")
