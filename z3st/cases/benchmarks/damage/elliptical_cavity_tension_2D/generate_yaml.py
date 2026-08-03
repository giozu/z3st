import numpy as np
import yaml

# Générez 201 steps
steps = np.linspace(0, 1.2e-07, 401).tolist()

# Créez la structure YAML
config = {
    'mechanical': {
        'uo2': [
            {'type': 'Clamp_y', 'region': 'ymin'},
            {'type': 'Clamp_x', 'region': 'xmin'},
            {
                'type': 'Dirichlet_y',
                'region': 'ymax',
                'value': steps
            }
        ]
    }
}

# Sauvegardez dans un fichier
with open('generate.yaml', 'w') as f:
    yaml.dump(config, f, default_flow_style=False)

print(f"✅ Generated {len(steps)} load steps")
print(f"First value: {steps[0]}")
print(f"Last value: {steps[-1]}")