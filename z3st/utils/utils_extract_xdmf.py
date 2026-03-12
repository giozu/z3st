#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os
import h5py
import numpy as np

def list_fields_xdmf(xdmf_path):
    """
    List all fields available in the H5 file associated with the XDMF.
    """
    h5_path = xdmf_path.replace(".xdmf", ".h5")
    if not os.path.exists(h5_path):
        print(f"[ERROR] H5 file not found: {h5_path}")
        return []

    fields = []
    with h5py.File(h5_path, 'r') as f:
        if 'Function' in f:
            for field in f['Function'].keys():
                fields.append(field)
                
    print(f"[INFO] Available fields in {h5_path}:")
    for field in fields:
        print(f"  → {field}")
    return fields

def extract_field_xdmf(xdmf_path, field_name, step_index=-1, return_coords=True):
    """
    Extract a field and its coordinates from a Z3ST-generated XDMF/H5 file.
    
    Parameters
    ----------
    xdmf_path : str
        Path to the .xdmf file.
    field_name : str
        Name of the field to extract (e.g., 'Damage', 'Stress_solid').
    step_index : int, optional
        Index of the time step to extract (default -1, the last one).
    return_coords : bool, optional
        Whether to return coordinates (default True).
        
    Returns
    -------
    If return_coords is True:
        x, y, z, data : np.ndarray
    Else:
        data : np.ndarray
    """
    h5_path = xdmf_path.replace(".xdmf", ".h5")
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"H5 file not found: {h5_path}")

    with h5py.File(h5_path, 'r') as f:
        if 'Function' not in f or field_name not in f['Function']:
            available = list(f['Function'].keys()) if 'Function' in f else []
            raise ValueError(f"Field '{field_name}' not found in {h5_path}. Available: {available}")
            
        field_group = f[f"Function/{field_name}"]
        
        # Steps are names of datasets, usually strings representing time or step index
        # Z3ST seems to name them by time values or step indices (e.g., '0', '400', '800')
        # We sort them numerically if possible
        def try_float(s):
            try:
                return float(s.replace("_", "."))
            except ValueError:
                return s

        steps = sorted(field_group.keys(), key=try_float)
        target_step = steps[step_index]
        data = np.array(field_group[target_step])
        
        if not return_coords:
            return data
            
        if 'Mesh' not in f or 'mesh' not in f['Mesh']:
            raise ValueError(f"Mesh not found in {h5_path}")
            
        geometry = np.array(f["Mesh/mesh/geometry"])
        x, y, z = geometry[:, 0], geometry[:, 1], geometry[:, 2]
        
        # Check if data is nodal or cell-based
        num_nodes = geometry.shape[0]
        if data.shape[0] == num_nodes:
            # Nodal data (Displacement, Damage, Temperature)
            return x, y, z, data
        else:
            # Cell data (Stress, Strain) - compute cell centers
            if 'topology' not in f["Mesh/mesh"]:
                raise ValueError(f"Topology not found in {h5_path} for cell data extraction")
                
            topology = np.array(f["Mesh/mesh/topology"])
            # topology shape: (num_cells, nodes_per_cell)
            # cell_centers: mean of coordinates of nodes in each cell
            cell_centers = np.mean(geometry[topology], axis=1)
            
            xc = cell_centers[:, 0]
            yc = cell_centers[:, 1]
            zc = cell_centers[:, 2] if cell_centers.shape[1] > 2 else np.zeros_like(xc)
            
            return xc, yc, zc, data
