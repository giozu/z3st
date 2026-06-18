# z3st/models/contact_model.py
from dolfinx import fem, mesh
import numpy as np
import ufl
from mpi4py import MPI
from typing import Dict, Optional, List
import petsc4py.PETSc as PETSc

class ContactModel:
    """
    Modèle pour gérer le contact mécanique entre deux sous-domaines d'un même maillage.
    Utilise une méthode de pénalisation pour modéliser les interactions de contact.
    """

    def __init__(
        self,
        mesh: mesh.Mesh,
        subdomains: Dict[str, int],
        V_m: Dict[str, fem.FunctionSpace],
        facet_tags: mesh.MeshTags,
        geometry_config: Dict,
        contact_config: Optional[Dict] = None,
    ):
        """
        Args:
            mesh: Maillage unique contenant les deux sous-domaines.
            subdomains: Dictionnaire associant les noms des sous-domaines à leurs tags.
            V_m: Dictionnaire des espaces fonctionnels (un par sous-domaine).
            facet_tags: Tags des facettes du maillage.
            geometry_config: Configuration géométrique (ex: R_mid).
            contact_config: Configuration du contact (ex: penalty, gap).
        """
        self.mesh = mesh
        self.subdomains = subdomains
        self.V_m = V_m
        self.facet_tags = facet_tags
        self.geometry_config = geometry_config
        self.contact_config = {
            "penalty": 1e5,            # Reduced penalty for better stability (Pa/m)
            "gap": 0.001,          # Jeu initial (m)
            "gap_tolerance": 0.000001,  # Tolérance pour la détection de contact
        }
        if contact_config:
            self.contact_config.update(contact_config)

        # Defensive cast to float to handle YAML scientific notation quirk (e.g. '1e8' as string)
        self.contact_config["penalty"] = float(self.contact_config["penalty"])
        self.contact_config["gap"] = float(self.contact_config["gap"])
        self.contact_config["gap_tolerance"] = float(self.contact_config["gap_tolerance"])

        # Trouver les facettes de l'interface
        self.interface_facets = self._find_interface_facets()

        # Identifier les sommets de l'interface pour le calcul des forces
        self.interface_vertices = self._find_interface_vertices()

        # Fonctions auxiliaires pour stocker le déplacement "cible" de l'autre côté
        self.u_target = {
            name: fem.Function(V_m[name]) for name in self.subdomains
        }

    def _find_interface_facets(self) -> List[int]:
        """
        Trouve les facettes de l'interface entre les deux sous-domaines.
        """
        labels = self.geometry_config.get("labels", {})
        
        # On cherche d'abord le label générique 'interface'
        if "interface" in labels:
            return self.facet_tags.find(labels["interface"])
        
        # Sinon, on cherche les labels spécifiques définis dans geometry.yaml
        # (ex: steel_in_right et steel_o_left)
        found_facets = []
        for contact_label in ["steel_in_right", "steel_o_left"]:
            if contact_label in labels:
                found_facets.append(self.facet_tags.find(labels[contact_label]))
        
        if not found_facets:
            available = list(labels.keys())
            raise KeyError(f"Label 'interface' non trouvé. Labels disponibles: {available}")
            
        return np.unique(np.concatenate(found_facets))

    def _find_interface_vertices(self) -> np.ndarray:
        """
        Identifie les indices des sommets (vertices) rattachés aux facettes de l'interface.
        """
        self.mesh.topology.create_connectivity(self.mesh.topology.dim - 1, 0)
        facet_to_vertex = self.mesh.topology.connectivity(self.mesh.topology.dim - 1, 0)
        
        vertices = []
        for facet in self.interface_facets:
            vertices.extend(facet_to_vertex.links(facet))
        return np.unique(np.array(vertices, dtype=np.int32))

    def compute_contact_forces(
        self,
        u_new: Dict[str, fem.Function],
    ) -> Dict[str, fem.Function]:
        """
        Calcule les forces de contact par pénalisation entre les deux sous-domaines.
        """
        penalty = self.contact_config["penalty"]
        gap_tolerance = self.contact_config["gap_tolerance"]

        # Récupérer les coordonnées des sommets (restreint à la dimension du déplacement)
        tdim = self.mesh.topology.dim
        interface_coords = self.mesh.geometry.x[self.interface_vertices, :tdim]

        # Récupérer les déplacements des sommets de l'interface
        # Pour P1, le reshape par tdim aligne les composantes par nœud
        u_in_all = u_new["steel_in"].x.array.reshape((-1, tdim))
        u_o_all = u_new["steel_o"].x.array.reshape((-1, tdim))
        
        u_in = u_in_all[self.interface_vertices]
        u_o = u_o_all[self.interface_vertices]

        # Synchroniser les déplacements cibles (steel_in voit le déplacement de steel_o et inversement)
        for i, v_idx in enumerate(self.interface_vertices):
            start, end = v_idx * tdim, (v_idx + 1) * tdim
            self.u_target["steel_in"].x.array[start:end] = u_o[i]
            self.u_target["steel_o"].x.array[start:end] = u_in[i]

    def add_contact_residual(
        self,
        F_m: Dict[str, ufl.Form],
        u_global: fem.Function,
    ) -> Dict[str, ufl.Form]:
        """
        Ajoute un terme de pénalisation symbolique lissé.
        """
        penalty = self.contact_config["penalty"]
        initial_gap = self.contact_config["gap"]
        eps_smooth = 1e-7 # Paramètre de lissage pour aider Newton
        
        # Normale adaptée à la dimension du problème
        tdim = self.mesh.topology.dim
        normal = ufl.as_vector([1.0] + [0.0]*(tdim-1)) 
        
        labels = self.geometry_config.get("labels", {})
        interface_map = {"steel_in": "steel_in_right", "steel_o": "steel_o_left"}

        for name in self.subdomains:
            v_m = ufl.TestFunction(self.V_m[name])
            interface_tag = labels.get(interface_map[name])
            ds_contact = ufl.ds(domain=self.mesh, subdomain_data=self.facet_tags, subdomain_id=interface_tag)
            
            u_target = self.u_target[name]
            # Calcul du gap physique incluant le jeu initial
            # g = initial_gap + (u_target - u_global) . n
            gap_phys = initial_gap + ufl.dot(u_target - u_global, normal)
            
            # Pénétration lissée
            penetration_raw = 0.5 * (gap_phys - ufl.sqrt(gap_phys**2 + eps_smooth**2))
            
            # Sécurité additionnelle : limiter la pénétration prise en compte pour stabiliser Newton
            # Si la pénétration est trop forte (> 1mm), on sature la force pour éviter l'explosion
            penetration = ufl.conditional(ufl.lt(penetration_raw, -0.001), -0.001, penetration_raw)
            
            force_sign = 1.0 if name == "steel_in" else -1.0
            force = force_sign * penalty * penetration * normal
            
            F_m[name] += ufl.dot(force, v_m) * ds_contact

        return F_m

    def _get_facet_tags_for_subdomain(self, subdomain_name: str):
        """
        Récupère les tags des facettes pour un sous-domaine donné.
        """
        # À adapter selon comment les tags sont stockés dans Z3ST
        return self.facet_tags