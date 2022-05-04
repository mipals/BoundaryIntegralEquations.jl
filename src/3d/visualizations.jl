function create_simple_mesh(bem_mesh::Mesh3d)
    return create_simple_mesh(bem_mesh,bem_mesh.shape_function)
end
function create_simple_mesh(tri_mesh,shape_function::Triangular)
    # Removing higher order terms
    initial_topology = tri_mesh.topology[1:3,:]
    used_nodes = sort(unique(initial_topology))
    topo       = remove_unused_nodes(initial_topology)
    x = tri_mesh.coordinates[1,used_nodes]
    y = tri_mesh.coordinates[2,used_nodes]
    z = tri_mesh.coordinates[3,used_nodes]
    points = Point.(x, y, z)
    connec = [connect(Tuple(Float64(c) for c in face)) for face in eachcol(topo[1:3,:])]
    return SimpleMesh(points, connec)
end

function create_simple_mesh(quad_mesh,shape_function::Quadrilateral)
    # We have to reorder the COMSOL layout to fit that of Meshes.jl
    initial_topology = quad_mesh.topology[[1;2;4;3],:]
    used_nodes = sort(unique(initial_topology))
    topo       = remove_unused_nodes(initial_topology)
    x = quad_mesh.coordinates[1,used_nodes]
    y = quad_mesh.coordinates[2,used_nodes]
    z = quad_mesh.coordinates[3,used_nodes]
    points = Point.(x, y, z)
    connec = [connect(Tuple(Float64(c) for c in face)) for face in eachcol(topo)]
    return SimpleMesh(points, connec)
end

function create_bc_simple_mesh(quad_mesh_file,quad_mesh,shape_function::Quadrilateral,bc_ents,bc=true)
    _,topo,ents = read_comsol_mesh(quad_mesh_file,quad_mesh.shape_function)
    if bc == true
        bc_id = Bool.(sum(bc_ents .∈ ents,dims=1))[:]
    else
        bc_id = .!Bool.(sum(bc_ents .∈ ents,dims=1))[:]
    end
    top_ents = topo[:,bc_id]
    # We have to reorder the COMSOL layout to fit that of Meshes.jl
    initial_topology = top_ents[[1;2;4;3],:]
    used_nodes   = sort(unique(initial_topology))
    topo        = IntegralEquations.remove_unused_nodes(initial_topology)
    x = quad_mesh.coordinates[1,used_nodes]
    y = quad_mesh.coordinates[2,used_nodes]
    z = quad_mesh.coordinates[3,used_nodes]
    points = Point.(x, y, z)
    connec = [connect(Tuple(Float64(c) for c in face)) for face in eachcol(topo)]
    return SimpleMesh(points, connec)
end

# function plot_mesh(bem_mesh::Mesh3d;showfacets=true)
#     plot_mesh(bem_mesh,bem_mesh.shape_function;showfacets=showfacets)
# end

# function plot_mesh(tri_mesh,geometry_type::Triangular;showfacets=true)
#     initial_topology = tri_mesh.topology[1:3,:]
#     used_nodes   = sort(unique(initial_topology))
#     topo        = IntegralEquations.remove_unused_nodes(initial_topology)
#     x = tri_mesh.coordinates[1,used_nodes]
#     y = tri_mesh.coordinates[2,used_nodes]
#     z = tri_mesh.coordinates[3,used_nodes]
#     points = Point.(x, y, z)
#     connec = [connect(Tuple(Float64(c) for c in face)) for face in eachcol(topo[1:3,:])]
#     m_lin = SimpleMesh(points, connec)
#     viz(m_lin, showfacets=showfacets)
# end
# function plot_mesh(quad_mesh,geometry_type::Quadrilateral;showfacets=true)
#     initial_topology = quad_mesh.topology[[1;2;4;3],:]
#     used_nodes   = sort(unique(initial_topology))
#     topo        = IntegralEquations.remove_unused_nodes(initial_topology)
#     x = quad_mesh.coordinates[1,used_nodes]
#     y = quad_mesh.coordinates[2,used_nodes]
#     z = quad_mesh.coordinates[3,used_nodes]
#     points = Point.(x, y, z)
#     connec = [connect(Tuple(Float64(c) for c in face)) for face in eachcol(topo)]
#     m_quad = SimpleMesh(points, connec)
#     viz(m_quad, showfacets=showfacets)
# end

