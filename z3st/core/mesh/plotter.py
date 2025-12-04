# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import matplotlib.cm as cm
import numpy as np
import pyvista
from dolfinx.plot import vtk_mesh

from z3st.core.diagnostic import log


class MeshPlotter:
    """Visualize mesh and boundary tags using PyVista."""

    def __init__(self, mesh, facet_tags, label_map: dict[str, int]):
        self.mesh = mesh
        self.facet_tags = facet_tags
        self.label_map = label_map

    def show(self):
        log.info("Rendering mesh with PyVista...")

        topology, cell_types, _ = vtk_mesh(self.mesh, self.mesh.topology.dim - 1)
        surface = pyvista.UnstructuredGrid(topology, cell_types, self.mesh.geometry.x)

        unique_tags = np.unique(self.facet_tags.values)
        print("[INFO] Face labels present:", unique_tags)

        colors = cm.get_cmap("tab10", len(unique_tags))
        label_map = {}
        for i, tag in enumerate(unique_tags):
            label_name = list(self.label_map.keys())[list(self.label_map.values()).index(tag)]
            label_map[tag] = (label_name, colors(i)[:3])

        plotter = pyvista.Plotter()
        for tag, (name, color) in label_map.items():
            facets = self.facet_tags.find(tag)
            cells = surface.extract_cells(facets)
            plotter.add_mesh(cells, color=color, show_edges=True, label=name)

        plotter.add_legend()
        plotter.show()
