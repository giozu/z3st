# import required model type
from gmshModel.Model import RandomInclusionRVE as RVE

# initialize new RVE
myRVE=RVE(size=[5,5,5], inclusionType="Sphere", inclusionSets=[1, 2])

# create Gmsh model
myRVE.createGmshModel()

# generate mesh
myRVE.createMesh()

# save resulting mesh to vtk
myRVE.saveMesh("myRVE.msh")

# visualize result
myRVE.visualizeMesh()

# finalize Gmsh-Python-API
myRVE.close()