"""
    create_simple_mesh(bem_mesh::Mesh3d)

Returns a `SimpleMesh` from the mesh. Can be used to plot with MeshViz.
"""
function create_simple_mesh()
end


"""
    create_bc_simple_mesh(mesh_file,mesh,boundary_condition_entities,boundary_conditions=true)
"""
function create_bc_simple_mesh()
end

"""
    create_vizualization_data(mesh,data)

Returns a linear `SimpleMesh` and interpolated `data_viz` vector.

Note that the interpolation of the data is linear even for quadratic elements.
"""
function create_vizualization_data()
end



# """
#     create_simple_mesh(bem_mesh::Mesh3d)

# Returns a `SimpleMesh` from the mesh. Can be used to plot with MeshViz.
# """
# function create_simple_mesh(bem_mesh::Mesh3d)
#     return create_simple_mesh(bem_mesh,bem_mesh.shape_function)
# end

# function create_simple_mesh(tri_mesh,shape_function::Triangular)
#     # MeshViz only support linear facets. We mimick this by remove 2nd order.
#     initial_topology = tri_mesh.topology[1:3,:]
#     # Removing unused nodes i.e. the 2nd order nodes
#     used_nodes   = sort(unique(initial_topology))
#     new_topology = remove_unused_nodes(initial_topology)
#     # Reordering coordinates
#     coordinates = tri_mesh.coordinates[:,used_nodes]
#     # Creating vector of points
#     points = Point.(coordinates[1,:], coordinates[2,:], coordinates[3,:])
#     # Creating vector of connectivities
#     connectivities = [connect(Tuple(Float64.(face))) for face in eachcol(new_topology)]
#     # Returning a SimpleMesh of the points and connectivities
#     return SimpleMesh(points, connectivities)
# end

# function create_simple_mesh(quad_mesh,shape_function::Quadrilateral)
#     # MeshViz only support linear facets. We mimick this by remove 2nd order.
#     # Furthermore we have to reorder the COMSOL layout to fit that of Meshes.jl
#     initial_topology = quad_mesh.topology[[1;2;4;3],:]
#     # Removing unused nodes i.e. the 2nd order nodes
#     used_nodes   = sort(unique(initial_topology))
#     new_topology = remove_unused_nodes(initial_topology)
#     # Reordering coordinates
#     coordinates = quad_mesh.coordinates[:,used_nodes]
#     # Creating vector of points
#     points = Point.(coordinates[1,:], coordinates[2,:], coordinates[3,:])
#     # Creating vector of connectivities
#     connectivities = [connect(Tuple(Float64.(face))) for face in eachcol(new_topology)]
#     # Returning a SimpleMesh of the points and connectivities
#     return SimpleMesh(points, connectivities)
# end


# """
#     create_bc_simple_mesh(mesh_file,mesh,boundary_condition_entities,boundary_conditions=true)
# """
# function create_bc_simple_mesh(mesh,boundary_condition_entities,boundary_conditions=true)
#     # For now re-read topology and entities
#     entities = mesh.entities
#     topology = mesh.topology
#     # Finding indicies of boundary conditions
#     boundary_condition_id = Bool.(sum(boundary_condition_entities' .∈ entities,dims=2))[:]
#     if boundary_conditions == false
#         boundary_condition_id = .!boundary_condition_id
#     end
#     #
#     boundary_condition_topology = topology[:,boundary_condition_id]
#     # We have to reorder the COMSOL layout to fit that of Meshes.jl
#     # initial_topology = boundary_condition_topology[[1;2;4;3],:]
#     if typeof(mesh.shape_function) <: Quadrilateral
#         initial_topology = boundary_condition_topology[[1;2;4;3],:]
#     elseif typeof(mesh.shape_function) <: Triangular
#         initial_topology = boundary_condition_topology[1:3,:]
#     end
#     # Remove unused nodes from topology
#     new_topology = remove_unused_nodes(initial_topology)
#     # Remove unused nodes from coordinate list
#     coordinates = mesh.coordinates[:,sort(unique(initial_topology))]
#     # Creating vector of points
#     points = Point.(coordinates[1,:], coordinates[2,:], coordinates[3,:])
#     # Creating vector of connectivities
#     connectivities = [connect(Tuple(Float64.(face))) for face in eachcol(new_topology)]
#     # Returning a SimpleMesh of the points and connectivities
#     return SimpleMesh(points, connectivities)
# end

# """
#     create_vizualization_data(mesh,data)

# Returns a linear `SimpleMesh` and interpolated `data_viz` vector.

# Note that the interpolation of the data is linear even for quadratic elements.
# """
# function create_vizualization_data(mesh,data)
#     return create_vizualization_data(mesh,data,mesh.physics_function)
# end
# function create_vizualization_data(mesh,data,physics_function::Triangular)
#     simple_mesh = create_simple_mesh(mesh)                        # Creating simple mesh
#     data_viz = data[sort!(unique(mesh.physics_topology[1:3,:]))]  # Removing quadratic parts
#     return simple_mesh, data_viz
# end
# function create_vizualization_data(mesh,data,physics_function::DiscontinuousTriangular)
#     # Repeating edge-points
#     T = reduce(hcat,[mesh.coordinates[:,top] for top in eachcol(mesh.topology[1:3,:])])
#     # Creating point data
#     points = Point.(T[1,:],T[2,:],T[3,:])
#     # Extracting the linear part of the discontinuous topology
#     new_topology = reshape(1:length(points),3,size(mesh.topology,2))
#     # Creating vector of connectivities
#     connectivities = [connect(Tuple(Float64.(face))) for face in eachcol(new_topology)]
#     # Combining connectivities and points to create a simple mesh
#     simple_mesh = SimpleMesh(points, connectivities)
#     # Pre-allocation of data_viz
#     data_viz = zeros(eltype(data), length(points))
#     # Making a physics element that interpolates on the corner nodes
#     physics_copy = deepcopy(physics_function)
#     set_interpolation_nodes!(physics_copy,[0.0;1.0;0.0],[0.0;0.0;1.0])
#     # Extrapolating data to corner points
#     for (data_top, phys_top) in zip(eachcol(new_topology),eachcol(mesh.physics_topology))
#         for (i,interp) in enumerate(eachcol(physics_copy.interpolation))
#             data_viz[data_top[i]] = dot(data[phys_top],interp)
#         end
#     end
#     return simple_mesh, data_viz
# end
