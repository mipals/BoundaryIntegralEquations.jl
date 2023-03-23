#
using Meshes, BoundaryIntegralEquations
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))

# # Visualizing the `connected_sources`
# # Loading meshes
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes")
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_extra_coarse");
tri_mesh_file = joinpath(mesh_path,"sphere_1m");
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[4])

# Compute connections
depth = 2
cone,sene = BoundaryIntegralEquations.connected_sources(mesh,depth)
j = size(mesh.sources,2) # Probably an edge node
j = 1 # probably a corner node
T = reduce(hcat,[mesh.coordinates[:,top] for top in eachcol(mesh.topology[1:3,cone[j]])])
# Creating point data
points = Point.(T[1,:],T[2,:],T[3,:])
# Extracting the linear part of the discontinuous topology
new_topology = new_topology = reshape(1:length(points),3,length(cone[j]))
# Creating vector of connectivities
connectivities = [connect(Tuple(Float64.(face))) for face in eachcol(new_topology)]
# Combining connectivities and points to create a simple mesh
simple_mesh = SimpleMesh(points, connectivities)
elm_sources = mesh.sources[:,sene[j]]
full_mesh = create_simple_mesh(mesh)
viz(full_mesh;showfacets=true)
viz!(simple_mesh;showfacets=true)
wgl.scatter!(T[1,:],T[2,:],T[3,:],color=:red)
wgl.scatter!(elm_sources[1,:],elm_sources[2,:],elm_sources[3,:],color=:blue)
wgl.scatter!(mesh.sources[1,j],mesh.sources[2,j],mesh.sources[3,j],color=:green)
