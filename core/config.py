class Config:
    def __init__(self, input_file):
        print("__Config initializer__")

        self.input_file = input_file
        
        models = self.input_file.get("models", {})
        self.on = {
            "thermal": models.get("thermal", False),
            "mechanical": models.get("mechanical", False),
            "fracture": models.get("fracture", False)
        }
        gap_config = self.input_file.get("models", {}).get("gap_conductance", {})
        self.gap_model = gap_config.get("type", None)
        self.h_gap_value = gap_config.get("value", 0.0)

        self.geometry_path = self.input_file.get("geometry_path", None)
        self.mesh_path = self.input_file.get("mesh_path", None)
        self.boundary_conditions_path = self.input_file.get("boundary_conditions_path", None)
        self.n_steps = self.input_file.get("n_steps", 10)
        self.coupling = self.input_file.get("solver_settings", {}).get("coupling", "staggered")
        self.regime = self.input_file.get("mechanical_regime", "3D")

        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")

        print(f"  → Gap conductance     : {self.gap_model}")
        if self.gap_model == "Fixed":
            print(f"    Value = {self.h_gap_value} W/m2-K")
        else:
            print(f"    Value = {self.h_gap_value}")
        
        print(f"  → Mechanical regime   : {self.regime}")
        print("\n")
