#!/usr/bin/env python3
"""
Plot cluster dynamics results from HDF5 file
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt

print("Reading HDF5 file...")

# Read HDF5 file directly
with h5py.File('output/results.h5', 'r') as f:
    print("Available groups:", list(f.keys()))
    
    # Navigate to ClusterDensity data
    if 'Function' in f and 'ClusterDensity' in f['Function']:
        cluster_group = f['Function']['ClusterDensity']
        
        # Get all time steps
        time_steps = sorted(cluster_group.keys())
        print(f"Found {len(time_steps)} time steps")
        
        # Read mesh coordinates from first dataset
        first_key = time_steps[0]
        first_data = cluster_group[first_key][:]
        n_points = len(first_data)
        
        # Assume uniform mesh from 0 to 100
        x_coords = np.linspace(0, 100, n_points)
        
        # Read data for selected time steps
        n_plots = min(6, len(time_steps))
        indices = np.linspace(0, len(time_steps)-1, n_plots, dtype=int)
        
        fig, ax = plt.subplots(figsize=(12, 7))
        colors = plt.cm.viridis(np.linspace(0, 1, n_plots))
        
        for i, idx in enumerate(indices):
            time_key = time_steps[idx]
            data = cluster_group[time_key][:].flatten()  # Flatten to 1D
            
            # Extract time from key (format: "0_1000000000000001" -> 0.1)
            try:
                time_val = float(time_key.replace('_', '.'))
            except:
                time_val = idx * 0.1
            
            # Relax mask and add diagnostics
            mask = data >= 0 # Plot everything to see if it's there
            if np.any(mask):
                ax.plot(x_coords, data, 
                        label=f't = {time_val:.2f} s', 
                        linewidth=2.5, color=colors[i], alpha=0.8)
                print(f"  Step {idx}: t={time_val:.2f}s, max={data.max():.2e}, avg={data.mean():.2e}")
        
        ax.set_xlabel('Cluster size n', fontsize=14, fontweight='bold')
        ax.set_ylabel('Cluster density c(n,t)', fontsize=14, fontweight='bold')
        ax.set_title('Cluster Dynamics Evolution', fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Only add legend if lines were added
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=11, loc='best', framealpha=0.9)
        ax.set_xlim([0, 100])
        
        # Use linear scale (data doesn't span many orders of magnitude after renormalization)
        ax.set_yscale('linear')
        
        plt.tight_layout()
        plt.savefig('output/cluster_evolution_plot.png', dpi=200, bbox_inches='tight')
        print('\n✓ Plot saved to output/cluster_evolution_plot.png')
        
    else:
        print("ERROR: ClusterDensity data not found in HDF5 file!")
        print("Available structure:")
        def print_structure(name, obj):
            print(f"  {name}")
        f.visititems(print_structure)
