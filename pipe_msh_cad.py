import gmsh
import sys
import os
import math
import meshio
gmsh.initialize()
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
gmsh.model.add("pipe_CAD")
# merge STL, create surface patches that are reparametrizable (so we can remesh
# them) and compute the parametrizations
path = os.path.dirname(os.path.abspath(__file__))
print(path)
print(os.pardir)
gmsh.merge(os.path.join(path, 'pipe_CAD.step'))
#v = gmsh.model.occ.importShapes(os.path.join(path, os.pardir, 'pipe_CAD.step'))
#v=gmsh.model.occ.importShapes(os.path.join(path, 'pipe_CAD.step'))
#gmsh.model.mesh.classifySurfaces(math.pi, True, True)
#gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

#gmsh.model.mesh.createGeometry()




#gmsh.model.mesh.optimize('Relocate2D',True)

in_marker=1
out_marker=2
wall_marker=3
domain_marker=4
gmsh.model.geo.addPhysicalGroup(2,[4],in_marker)
gmsh.model.setPhysicalName(2,in_marker,'in')
gmsh.model.geo.addPhysicalGroup(2,[3],out_marker)
gmsh.model.setPhysicalName(2,out_marker,'out')
gmsh.model.geo.addPhysicalGroup(2,[1,2],wall_marker)
gmsh.model.setPhysicalName(2,wall_marker,'wall')
#gmsh.model.geo.addPhysicalGroup(2,[177],6)
#gmsh.model.setPhysicalName(2,6,'boundary')
gmsh.model.geo.addPhysicalGroup(3,[1],domain_marker)
gmsh.model.setPhysicalName(3,domain_marker,'domain')
    # create the inner volume
#gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(top_surf)])

gmsh.model.geo.synchronize()
#gmsh.model.mesh.setCompound(2,[2,3])
# use MeshAdapt for the resulting not-so-smooth parametrizations
gmsh.option.setNumber('Mesh.Algorithm', 1)
#gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.1)
gmsh.option.setNumber('Mesh.MeshSizeMin', 0.01)
gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
gmsh.model.mesh.generate(3)

gmsh.model.geo.synchronize()
gmsh.write('tipe_cad.msh')
#Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
msh=meshio.read("tipe_cad.msh")
meshio.write("pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
#boundary face
meshio.write("mf_pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
#boundary cell
meshio.write("cf_pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]},cell_data={"tetra": {"name_to_read":msh.cell_data["tetra"]["gmsh:physical"]}}))
