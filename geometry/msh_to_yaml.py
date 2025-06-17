import yaml

msh_file = "../mesh/thermal_shield_slab.msh"
yaml_file = "thermal_shield_slab.yaml"

def extract_physical_names(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    try:
        start = lines.index("$PhysicalNames\n") + 1
        end = lines.index("$EndPhysicalNames\n")
    except ValueError:
        raise RuntimeError("Section $PhysicalNames not found in the file.")

    section = lines[start:end]
    if section[0].strip().isdigit():
        section = section[1:]

    labels = {}
    for line in section:
        parts = line.strip().split()
        tag = int(parts[1])
        name = parts[2].strip('"')
        labels[name] = tag

    return labels

def write_yaml(labels, output_file):
    data = {
        "name": "thermal_shield_slab",
        "edge": "0.045",                    # meters
        "width": "0.1",                     # meters
        "height": "0.1",                    # meters
        "geometry_type": "rect",
        "labels": labels
    }
    with open(output_file, "w") as f:
        yaml.dump(data, f, sort_keys=False)
    print(f"YAML written to {output_file}")

labels = extract_physical_names(msh_file)
write_yaml(labels, yaml_file)
