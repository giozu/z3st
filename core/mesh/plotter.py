# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import numpy as np
import matplotlib.cm as cm
import pyvista
from dolfinx.plot import vtk_mesh
from core.diagnostic import log


class MeshPlotter:
    """Visualize mesh and boundary tags using PyVista."""

    def __init__(self, mesh, facet_tags, label_map: dict[str, int]):
        self.mesh = mesh
        self.facet_tags = facet_tags
        self.label_map = label_map

    def show(self):
        log.info("Rendering mesh with PyVista...")
        topology, cell_types, geometry = vtk_mesh(self.mesh, self.mesh.topology.dim - 1)
        surface = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        unique_tags = np.unique(self.facet_tags.values)
        cmap = cm.get_cmap("tab10", len(unique_tags))
        label_names = list(self.label_map.keys())

        plotter = pyvista.Plotter()
        for i, tag in enumerate(unique_tags):
            name = label_names[i] if i < len(label_names) else f"Tag {tag}"
            faces = surface.extract_cells(np.where(self.facet_tags.values == tag)[0])
            plotter.add_mesh(faces, color=cmap(i)[:3], show_edges=True, label=name)

        plotter.add_legend()
        plotter.show()
