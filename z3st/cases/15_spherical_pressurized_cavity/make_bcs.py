import yaml
import numpy as np

n_steps = 50
max_disp = 1.0e-6  # Max macroscopic displacement (1 micrometer)

# Create a linear ramp of displacements for the 50 steps
disp_history = np.linspace(0.0, max_disp, n_steps).tolist()

bcs = {
    "mechanical": {
        "solid": [
            {"type": "Clamp_x", "region": "xmin"}, # X plane symmetry
            {"type": "Clamp_y", "region": "ymin"}, # Y plane symmetry
            {"type": "Clamp_z", "region": "zmin"}, # Z plane symmetry (bottom clamp)
            {"type": "Dirichlet_z", "region": "zmax", "value": disp_history} # Top traction
        ]
    }
}

with open("boundary_conditions.yaml", "w") as f:
    yaml.dump(bcs, f, sort_keys=False)

print("boundary_conditions.yaml generated successfully!")