class Config:
    def __init__(self, input_file):
        print("__Config initializer__")
        self.config = input_file
        
        models = self.config.get("models", {})
        self.on = {
            "thermal": models.get("thermal", True),
            "mechanical": models.get("mechanical", True),
            "creep": models.get("creep", False)
        }
        
        self.geometry_path = self.config.get("geometry_path", None)
        self.mesh_path = self.config.get("mesh_path", None)
        self.boundary_conditions_path = self.config.get("boundary_conditions_path", None)
        self.n_steps = self.config.get("n_steps", 10)
        self.coupling = self.config.get("solver_settings", {}).get("coupling", "staggered")

        print(f"  → Geometry            : {self.geometry_path}")
        print(f"  → Mesh                : {self.mesh_path}")
        print(f"  → Boundary conditions : {self.boundary_conditions_path}")
        print(f"  → Time steps          : {self.n_steps}")
        print(f"  → Models active       :")
        for model, active in self.on.items():
            print(f"      {model:<10} → {'ON' if active else 'OFF'}")