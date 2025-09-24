class Config:
    def __init__(self, input_file):
        print("__Config initializer__")

        self.input_file = input_file

        # --- Active models ---
        models = self.input_file.get("models", {})
        self.on = {
            "thermal": models.get("thermal", False),
            "mechanical": models.get("mechanical", False),
            "fracture": models.get("fracture", False)
        }

        # --- Gap conductance ---
        gap_config = self.input_file.get("models", {}).get("gap_conductance", {})
        self.gap_model = gap_config.get("type", None)
        self.h_gap_value = gap_config.get("value", 0.0)

        # --- Paths ---
        self.geometry_path = self.input_file.get("geometry_path", None)
        self.mesh_path = self.input_file.get("mesh_path", None)
        self.boundary_conditions_path = self.input_file.get("boundary_conditions_path", None)
        self.n_steps = self.input_file.get("n_steps", 10)

        # --- Coupling regime ---
        self.coupling = self.input_file.get("solver_settings", {}).get("coupling", "staggered")

        # --- Mechanical solver options ---
        mech_cfg = self.input_file.get("mechanical", {})
        self.mech_solver_type = mech_cfg.get("solver", "linear")
        self.mech_linear_solver = mech_cfg.get("linear_solver", "iterative_amg")
        self.rtol_mech = float(mech_cfg.get("rtol", 1e-6))
        self.stag_tol_mech = float(mech_cfg.get("stag_tol", 1e-8))
        self.mech_regime = mech_cfg.get("mechanical_regime", "3D")
        self.mech_convergence = mech_cfg.get("convergence", "norm")

        # --- Thermal solver options ---
        th_cfg = self.input_file.get("thermal", {})
        self.th_solver_type = th_cfg.get("solver", "linear")
        self.th_linear_solver = th_cfg.get("linear_solver", "iterative_amg")
        self.rtol_th = float(th_cfg.get("rtol", 1e-6))
        self.stag_tol_th = float(th_cfg.get("stag_tol", 1e-8))
        self.th_convergence = th_cfg.get("convergence", "norm")

        # --- Phase-field solver options (if fracture is ON) ---
        pf_cfg = self.input_file.get("fracture", {})
        self.tol_pf = float(pf_cfg.get("tol", 1e-3))  # default se non definito

        # --- Debug print ---
        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")

        print(f"  → Gap conductance     : {self.gap_model} (value = {self.h_gap_value})")
        print(f"  → Mechanical    : {self.mech_regime}")
        print(f"  → Mechanical solver   : {self.mech_solver_type} ({self.mech_linear_solver})")
        print(f"     tolerances         : rtol={self.rtol_mech}, stag_tol={self.stag_tol_mech}")
        print(f"     convergence        : {self.mech_convergence}")
        print(f"  → Thermal solver      : {self.th_solver_type} ({self.th_linear_solver})")
        print(f"     tolerances         : rtol={self.rtol_th}, stag_tol={self.stag_tol_th}")
        print(f"     convergence        : {self.th_convergence}")
        print("\n")
