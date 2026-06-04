Info    : Running '/home/baptiste/miniconda3/envs/z3st/bin/gmsh mesh.geo -3' [Gmsh 4.15.2, 1 node, max. 1 thread]
Info    : Started on Thu Jun  4 15:11:49 2026
Info    : Reading 'mesh.geo'...
Info    : Done reading 'mesh.geo'
Info    : Meshing 1D...
Info    : [ 10%] Meshing curve 14 (TrimmedCurve)
Info    : [ 30%] Meshing curve 16 (Line)
Info    : [ 30%] Meshing curve 17 (Line)
Info    : [ 40%] Meshing curve 18 (Line)
Info    : [ 50%] Meshing curve 19 (Line)
Info    : [ 50%] Meshing curve 20 (Line)
Info    : [ 60%] Meshing curve 21 (Line)
Info    : [ 70%] Meshing curve 22 (Line)
Info    : [ 70%] Meshing curve 23 (Line)
Info    : [ 80%] Meshing curve 24 (Line)
Info    : [ 90%] Meshing curve 25 (Line)
Info    : [ 90%] Meshing curve 26 (Line)
Info    : [100%] Meshing curve 27 (Line)
Info    : Done meshing 1D (Wall 0.0211306s, CPU 0.021131s)
Info    : Meshing 2D...
Info    : [  0%] Meshing surface 7 (BSpline surface, Frontal-Delaunay)
Info    : [  0%] Meshing surface 7 (BSpline surface, MeshAdapt)
Info    : [ 20%] Meshing surface 8 (Plane, Frontal-Delaunay)
Info    : [ 30%] Meshing surface 9 (Plane, Frontal-Delaunay)
Info    : [ 50%] Meshing surface 10 (Plane, Frontal-Delaunay)
Info    : [ 60%] Meshing surface 11 (Plane, Frontal-Delaunay)
Info    : [ 80%] Meshing surface 12 (Plane, Frontal-Delaunay)
Info    : [ 90%] Meshing surface 13 (Plane, Frontal-Delaunay)
Info    : Done meshing 2D (Wall 0.574791s, CPU 0.57314s)
Info    : Meshing 3D...
Info    : 3D Meshing 1 volume with 1 connected component
Info    : Tetrahedrizing 1796 nodes...
Info    : Done tetrahedrizing 1804 nodes (Wall 0.0273803s, CPU 0.027378s)
Info    : Reconstructing mesh...
Info    :  - Creating surface mesh
Info    :  - Identifying boundary edges
Info    :  - Recovering boundary
Info    : Done reconstructing mesh (Wall 0.0560154s, CPU 0.052137s)
Info    : Found void region
Info    : Found volume 1
Info    : It. 0 - 0 nodes created - worst tet radius 10.479 (nodes removed 0 0)
Info    : It. 500 - 500 nodes created - worst tet radius 1.80757 (nodes removed 0 0)
Info    : It. 1000 - 1000 nodes created - worst tet radius 1.40884 (nodes removed 0 0)
Info    : It. 1500 - 1500 nodes created - worst tet radius 1.24763 (nodes removed 0 0)
Info    : It. 2000 - 2000 nodes created - worst tet radius 1.12155 (nodes removed 0 0)
Info    : It. 2500 - 2500 nodes created - worst tet radius 1.04277 (nodes removed 0 0)
Info    : 3D refinement terminated (4645 nodes total):
Info    :  - 0 Delaunay cavities modified for star shapeness
Info    :  - 0 nodes could not be inserted
Info    :  - 23605 tetrahedra created in 0.190615 sec. (123836 tets/s)
Info    : 0 node relocations
Info    : Done meshing 3D (Wall 0.316044s, CPU 0.308803s)
Info    : Optimizing mesh...
Info    : Optimizing volume 1
Info    : Optimization starts (volume = 1.72597e-12) with worst = 0.00991723 / average = 0.755016:
Info    : 0.00 < quality < 0.10 :        54 elements
Info    : 0.10 < quality < 0.20 :       172 elements
Info    : 0.20 < quality < 0.30 :       297 elements
Info    : 0.30 < quality < 0.40 :       479 elements
Info    : 0.40 < quality < 0.50 :       805 elements
Info    : 0.50 < quality < 0.60 :      1539 elements
Info    : 0.60 < quality < 0.70 :      3339 elements
Info    : 0.70 < quality < 0.80 :      5772 elements
Info    : 0.80 < quality < 0.90 :      7506 elements
Info    : 0.90 < quality < 1.00 :      3639 elements
Info    : 509 edge swaps, 7 node relocations (volume = 1.72597e-12): worst = 0.272112 / average = 0.767719 (Wall 0.0130017s, CPU 0.013003s)
Info    : 511 edge swaps, 7 node relocations (volume = 1.72597e-12): worst = 0.272112 / average = 0.767745 (Wall 0.0151142s, CPU 0.015116s)
Info    : No ill-shaped tets in the mesh :-)
Info    : 0.00 < quality < 0.10 :         0 elements
Info    : 0.10 < quality < 0.20 :         0 elements
Info    : 0.20 < quality < 0.30 :         1 elements
Info    : 0.30 < quality < 0.40 :       477 elements
Info    : 0.40 < quality < 0.50 :       815 elements
Info    : 0.50 < quality < 0.60 :      1510 elements
Info    : 0.60 < quality < 0.70 :      3311 elements
Info    : 0.70 < quality < 0.80 :      5871 elements
Info    : 0.80 < quality < 0.90 :      7547 elements
Info    : 0.90 < quality < 1.00 :      3622 elements
Info    : Done optimizing mesh (Wall 0.0448061s, CPU 0.044802s)
Info    : 4646 nodes 26918 elements
Info    : Writing 'mesh.msh'...
Info    : Done writing 'mesh.msh'
Info    : Stopped on Thu Jun  4 15:11:51 2026 (From start: Wall 1.48271s, CPU 1.96578s)
