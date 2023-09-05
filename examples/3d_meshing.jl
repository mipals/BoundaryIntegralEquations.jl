using BoundaryIntegralEquations
using Meshes
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))

path = "/Users/mpasc/Documents/testfiles"
tri_mesh_file = path * "/test.obj"
# tri_mesh_file = path * "/waveguideDTU.stl"
tri_mesh_file = path * "/binary.stl"
# tri_mesh_file = path * "/test2.off"
# tri_mesh_file = path * "/polygonal_face.obj"
# tri_mesh_file = path * "/test_binary.ply"
@time tri_mesh  = load3dTriangularMesh(tri_mesh_file);
@time simple_tri_mesh = create_simple_mesh(tri_mesh);
viz(simple_tri_mesh;showfacets=true)

# quad_mesh_file = path * "/quad_sphere_binary.ply" # Can read but can not plot?
# quad_mesh_file = path * "/quad_sphere_text.ply"   # Can read but can not plot?
# quad_mesh_file = path * "/quad_sphere_binary.stl" # Can not read...
quad_mesh_file = path * "/quad_sphere_text.stl"
@time tri_mesh = load3dTriangularMesh(quad_mesh_file);
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)
