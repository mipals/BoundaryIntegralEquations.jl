using IntegralEquations
using MeshViz
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(600, 600))

path = "/Users/mpasc/Documents/testfiles"
# tri_mesh_file = path * "/test.obj"
# tri_mesh_file = path * "/waveguideDTU.stl"
# tri_mesh_file = path * "/binary.stl"
# tri_mesh_file = path * "/test2.off"
# tri_mesh_file = path * "/polygonal_face.obj"
tri_mesh_file = path * "/test_binary.ply"
# mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=:linear)
@time tri_mesh = IntegralEquations.load3dTriangularMesh(tri_mesh_file;geometry_order=:linear);
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)
