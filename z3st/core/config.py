# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


class Config:
    """
    Configuration manager for Z3ST simulations.

    Parses the user input YAML and initializes the global configuration used by
    the other modules. Loads:
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
            "cluster": models.get("cluster", False),
            "plasticity": models.get("plasticity", False),
            "contact": bool(models.get("contact", False)),
        }

        # --. Fission-gas behaviour via SCIANTIX coupling (default OFF) --..
        # ``models.fission_gas`` may be a bool or a block:
        #   fission_gas:
        #     enabled: true
        #     lib: /path/to/libsciantix.so       # else $SCIANTIX_LIB
        #     initial_conditions: input_initial_conditions.txt
        #     energy_per_fission: 3.2e-11        # J/fission (≈ 200 MeV)
        # SCIANTIX also reads input_settings.txt from the run directory.
        fg = models.get("fission_gas", False)
        if fg is True:
            fg = {"enabled": True}            # bare "fission_gas: true" shorthand
        elif not isinstance(fg, dict):
            fg = {"enabled": bool(fg)}         # any other scalar (False/0/None)
        self.on["fission_gas"] = bool(fg.get("enabled", False))
        self.sciantix_lib = fg.get("lib", None)
        self.sciantix_ic = fg.get("initial_conditions", "input_initial_conditions.txt")
        self.sciantix_energy_per_fission = float(fg.get("energy_per_fission", 3.2e-11))

        gap_config = self.input_file.get("models", {}).get("gap_conductance", {})
        self.gap_model = gap_config.get("type", None)
        self.h_gap_value = gap_config.get("value", 0.0)

        # Contact-coupled gap conductance (Todreas & Kazimi, Nuclear Systems I,
        # 3rd ed., Eq. 8.141/8.142): on gap closure a contact term proportional
        # to the pellet-clad contact pressure is added to the open-gap value.
        cc = gap_config.get("contact_coupling", {})
        self.gap_contact_coupling = bool(cc.get("enabled", False))
        self.gap_meyer_hardness = float(cc.get("meyer_hardness", 9.65e8))   # Pa (Zircaloy ~ 14e4 psi)
        self.gap_contact_thickness = float(cc.get("gas_thickness", 4.0e-6))  # m (roughness-based gas space on contact)

        # Under-relaxation of the gap conductance between staggered iterations
        # (h ← ω·h_new + (1−ω)·h_prev). The contact-pressure → conductance →
        # temperature → expansion → pressure feedback is the loop that chatters
        # on gap closure; damping h attacks it at the source. 1.0 = off.
        self.gap_relax = float(gap_config.get("relax", 1.0))
        
        # --. Paths --..
        self.geometry_path = self.input_file.get("geometry_path", None)
        self.mesh_path = self.input_file.get("mesh_path", None)
        self.boundary_conditions_path = self.input_file.get("boundary_conditions_path", None)
        self.n_steps = self.input_file.get("n_steps", 10)
        # Normalised to lowercase ("2D" → "2d"); downstream regime branches
        # assume one of these five values, so reject anything else up front.
        self.regime = self.input_file.get("regime", "2d").lower()
        valid_regimes = {"1d", "2d", "3d", "axisymmetric", "plane_stress"}
        if self.regime not in valid_regimes:
            raise ValueError(
                f"Invalid regime '{self.regime}'. Must be one of {sorted(valid_regimes)}."
            )

        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Regime              : {self.regime}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")
        print(f"  → Gap conductance     : {self.gap_model} (value = {self.h_gap_value})")
        print("\n")
