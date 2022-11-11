using IntegralEquations
using MeshViz
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(600, 600))

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


coords, topology, ents = IntegralEquations.load_mesh_file(quad_mesh_file)

import GeometryBasics
using FileIO, MeshIO
quad_mesh_file = path * "/quad_sphere_binary.ply"
mesh = load(quad_mesh_file)
mesh_coords = GeometryBasics.coordinates(mesh)
elements    = GeometryBasics.faces(mesh)
n_elements  = length(elements)
n_interpolations = length(elements[1])

topology = zeros(Int64,n_interpolations,n_elements)
coords   = zeros(3,length(mesh_coords))
for (index,element) in enumerate(elements)
    for i = 1:n_interpolations
        topology[i,index] = Int64(element[i].i) + 1
    end
end
for (index,coord) in enumerate(mesh_coords)
    coords[:,index] = Float64.(coord)
end

ents = zeros(Int64, n_elements)
