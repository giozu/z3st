#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST test script --.. ..- .-.. .-.. ---
"""
Z3ST case: cluster_dynamics_test

"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

# HDF5 file
h5_file = 'output/results.h5'

with h5py.File(h5_file, 'r') as f:
    # Mesh coordinates
    if 'Mesh' in f and 'mesh' in f['Mesh']:
        # Extracting the coordinates (x, y, z) of all the nodes
        coords = f['Mesh']['mesh']['geometry'][:]
        # Mesh is 1D, only first column
        x_coords_raw = coords[:, 0]
        print(f"Mesh loaded: n_points={len(x_coords_raw)}, range=[{x_coords_raw.min()}, {x_coords_raw.max()}]")
    else:
        raise KeyError("Mesh coordinates not found in HDF5 file!")

    # Cluster data
    if 'Function' in f and 'ClusterDensity' in f['Function']:
        cluster_group = f['Function']['ClusterDensity']
        
        # Helper to sort keys numerically by extracting the time value
        def get_time(key):
            try:
                return float(key.replace('_', '.'))
            except:
                return 0.0
        
        # Sort keys based on the numerical time value
        time_steps = sorted(cluster_group.keys(), key=get_time)
        print(f"Found {len(time_steps)} steps. Sorting numerically...")
        
        # Geometric reordering
        sort_idx = np.argsort(x_coords_raw)
        x_coords = x_coords_raw[sort_idx]
        
        # Matplotlib style improvements
        plt.style.use('seaborn-v0_8-muted') # or 'ggplot'
        fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
        
        frames = []
        all_max = 0
        for time_key in time_steps:
            all_max = max(all_max, np.max(cluster_group[time_key][:]))
        
        print("Rendering frames...")
        for time_key in time_steps:
            c_values = cluster_group[time_key][:][sort_idx]
            time_val = get_time(time_key)
            
            # Use a slightly more "premium" look
            line, = ax.plot(x_coords, c_values, color='#0077CC', lw=2.5, alpha=0.9)
            
            # Title as a centered artist
            title = ax.text(0.5, 1.05, f'Z3ST, cluster dynamics | Time: {time_val:.2f} s', 
                           transform=ax.transAxes, ha='center', fontsize=12, 
                           fontweight='bold', color='#333333')
            
            frames.append([line, title])

        ax.set_xlabel('Cluster size n', fontsize=11)
        ax.set_ylabel('Density c(n)', fontsize=11)
        ax.set_ylim([0, all_max * 1.1])
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_facecolor('#F8F9FA')
        
        # Create and save animation
        from matplotlib.animation import ArtistAnimation
        print("Creating high-quality GIF (output/cluster_evolution.gif)...")
        ani = ArtistAnimation(fig, frames, interval=250, blit=True)
        ani.save('output/cluster_evolution.gif', writer='pillow')
        plt.close()
        print("Done! GIF generated successfully.")