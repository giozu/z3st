import yaml

msh_file = "../mesh/full_cylinder_2.msh"
yaml_file = "full_cylinder.yaml"

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
        "name": "full_cylinder",
        "radius": "0.004",              # meters
        "height": "0.01",               # meters
        "geometry_type": "cyl",
        "labels": labels
    }
    with open(output_file, "w") as f:
        yaml.dump(data, f, sort_keys=False)
    print(f"YAML written to {output_file}")

labels = extract_physical_names(msh_file)
write_yaml(labels, yaml_file)
