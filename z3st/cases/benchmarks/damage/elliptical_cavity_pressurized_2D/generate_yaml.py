import numpy as np
import yaml

# Ramp internal pressure (traction) from 0 up to a target value, over 401
# steps to match the reference load-stepping convergence config (case 15).
#
# TODO: pick a sensible P_MAX. As a first guess, size it so the expected
# critical cracking pressure (order of magnitude vs. sigma_c ~ 150 MPa and
# the bubble's stress concentration factor) falls comfortably inside the
# ramp -- e.g. start with P_MAX ~ 1.5e8 (150 MPa) and adjust once you see
# where damage actually starts growing.
P_MAX = 1.5e8  # Pa -- PLACEHOLDER, tune this
N_STEPS = 401 # number of load steps

steps = np.linspace(0, P_MAX, N_STEPS).tolist()

# Créez la structure YAML
config = {
    'mechanical': {
        'uo2': [
            {'type': 'Clamp_y', 'region': 'ymin'},
            {'type': 'Clamp_x', 'region': 'xmin'},
            {
                'type': 'Neumann',
                'region': 'cavity',
                'traction': steps
            }
        ]
    }
}

# Sauvegardez dans un fichier
with open('boundary_conditions.yaml', 'w') as f:
    yaml.dump(config, f, default_flow_style=False)

print(f"✅ Generated {len(steps)} pressure load steps")
print(f"First value: {steps[0]} Pa")
print(f"Last value: {steps[-1]} Pa")