get_offset(physics_element::DiscontinuousTriangular)    = 0.10
get_offset(physics_element::ContinuousTriangular)       = 0.0
get_offset(physics_element::DiscontinuousQuadrilateral) = 0.0
get_offset(physics_element::ContinuousQuadrilateral)    = 0.1


get_beta(physicsElement::Triangular) = 0.0
get_beta(physicsElement::DiscontinuousTriangularConstant)  = 0.5
get_beta(physicsElement::DiscontinuousTriangularLinear)    = physicsElement.beta
get_beta(physicsElement::DiscontinuousTriangularQuadratic) = physicsElement.beta
#==========================================================================================
                                TriangularQuadratic Surface
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3
                                          6  5
                                          1  4  2
==========================================================================================#
number_of_corners(::Triangular)    = 3
number_of_corners(::Quadrilateral) = 4
shape_connections(::ShapeFunction) = []
shape_connections(::TriangularQuadratic) = [[1 2],[2 3],[3 1]]
shape_connections(::QuadrilateralQuadraticLagrange) = [[1 2],[1 3],[1 2 3 4],[2 4],[3 4]]
shape_connections(::QuadrilateralQuadratic) = [[1 2],[2 3],[3 4],[4 1]]



function connected_sources(mesh,depth=1)
    @assert depth ∈ [1, 2]
    n_sources = size(mesh.sources,2)
    source_connections  = [zeros(Int64,0) for i = 1:n_sources]
    element_connections = [zeros(Int64,0) for i = 1:n_sources]
    # Maybe we should use actual topology for this to work?
    n_corners = number_of_corners(mesh.physics_function)
    n_shape_functions = size(mesh.physics_topology,1)
    for element = 1:size(mesh.topology,2)
        for i = 1:n_shape_functions
            append!(element_connections[mesh.physics_topology[i,element]],element)
            append!(source_connections[mesh.physics_topology[i,element]],mesh.physics_topology[:,element])
        end
    end
    if depth == 1
        conns = shape_connections(mesh.physics_function)
        for element = 1:size(mesh.physics_topology,2)
            for idx = (n_corners + 1):n_shape_functions
                for j ∈ conns[idx-n_corners]
                    append!(element_connections[mesh.physics_topology[idx,element]],element_connections[mesh.physics_topology[j,element]])
                    append!(source_connections[mesh.physics_topology[idx,element]],  source_connections[mesh.physics_topology[j,element]])
                end
            end
        end
        return sort!.(unique.(element_connections)), sort!.(unique.(source_connections))
    elseif depth == 2
        element_connections2 = [zeros(Int64,0) for i = 1:n_sources]
        source_connections2  = [zeros(Int64,0) for i = 1:n_sources]
        for source=1:n_sources
            for idx ∈ source_connections[source]
                append!(element_connections2[source],element_connections[idx])
                append!(source_connections2[source],source_connections[idx])
            end
        end
        return sort!.(unique.(element_connections2)), sort!.(unique.(source_connections2))
    end

end

function create_row_indices(lengths,total_counts)
    row_indices =zeros(Int64,total_counts)
    row = 1
    idx = 0
    for i ∈ lengths
        row_indices[idx+1:idx+i] .= row
        row = row + 1
        idx = idx + i
    end
    return row_indices
end


function sparse_computing_integrals!(physics_interpolation,shape_function,
                                    normals, tangents,sangents,
                                    interpolation,jacobian_mul_weights,r,integrand,
                                    element_coordinates,
                                    fOn,gOn,submatrixF,submatrixG,k,source)

    my_mul!(interpolation,element_coordinates,shape_function.interpolation)
    compute_distances!(r,interpolation,source)
    jacobian!(shape_function,element_coordinates,
                                    normals,tangents,sangents,jacobian_mul_weights)
    integrand_mul!(jacobian_mul_weights,shape_function.weights)
    if gOn
        ### Evaluating the G-kernel (single-layer kernel) at the global nodes
        greens3d!(integrand,r,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating the integral and the adding the values to the BEM matrix
        mygemm_vec!(submatrixG,integrand,physics_interpolation)
    end
    if fOn
        ### Evaluating the F-kernel (double-layer kernel) at the global nodes
        freens3d!(integrand,r,interpolation,source,normals,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating integral and adding the value to the BEM matrix
        mygemm_vec!(submatrixF,integrand,physics_interpolation)
    end
end

function find_physics_nodes!(physics_nodes,idxs,di,physics_top)
    @inbounds for i = 1:length(physics_nodes)
        physics_nodes[i] = idxs + di[physics_top[i]]
    end
    # [(idx[source_node] + di[t]) for t in physics_topology[:,element]]
end

function sparse_assemble_parallel!(mesh::Mesh3d,k,sources,shape_function::TriangularLinear;
                                fOn=true,gOn=true,progress=true)
    topology    = get_topology(mesh)
    n_sources   = size(sources,2)
    n_nodes     = size(mesh.sources,2)
    coordinates = mesh.coordinates
    physics_topology = mesh.physics_topology
    physics_function = mesh.physics_function
    #======================================================================================
        Introducing three elements: (The hope here is to compute singular integrals)
        The numbers corresponds to the corner for which the GP-points are clustered
        As such choose the clustering that are closest to the source point
        ————————————————————————————————————  Grid  ————————————————————————————————————————
                                            3
                                            | \
                                            1 - 2
    ======================================================================================#
    # shape_function1 = create_rotated_element(shape_function,n,m,1)
    # shape_function2 = create_rotated_element(shape_function,n,m,2)
    # shape_function3 = create_rotated_element(shape_function,n,m,3)
    beta = get_beta(physics_function)
    tmp  = DiscontinuousTriangularLinear(physics_function,beta)
    offr = get_offset(physics_function)
    # Dealing with corners
    nodesX1,nodesY1,weights1 = singular_triangle_integration(tmp,3,[1.00;offr;offr],Diagonal(ones(3)),1e-6)
    nodesX2,nodesY2,weights2 = singular_triangle_integration(tmp,3,[offr;1.00;offr],Diagonal(ones(3)),1e-6)
    nodesX3,nodesY3,weights3 = singular_triangle_integration(tmp,3,[offr;offr;1.00],Diagonal(ones(3)),1e-6)
    shape_function1 = deepcopy(shape_function)
    shape_function2 = deepcopy(shape_function)
    shape_function3 = deepcopy(shape_function)
    set_interpolation_nodes!(shape_function1,nodesX1,nodesY1,weights1)
    set_interpolation_nodes!(shape_function2,nodesX2,nodesY2,weights2)
    set_interpolation_nodes!(shape_function3,nodesX3,nodesY3,weights3)
    physics_function1 = deepcopy(physics_function)
    physics_function2 = deepcopy(physics_function)
    physics_function3 = deepcopy(physics_function)
    set_interpolation_nodes!(physics_function1,nodesX1,nodesY1,weights1)
    set_interpolation_nodes!(physics_function2,nodesX2,nodesY2,weights2)
    set_interpolation_nodes!(physics_function3,nodesX3,nodesY3,weights3)

    # Copying interpolation of physics functions1
    physics_interpolation1 = copy(physics_function1.interpolation')
    physics_interpolation2 = copy(physics_function2.interpolation')
    physics_interpolation3 = copy(physics_function3.interpolation')

    # Connections
    element_connections, source_connections = connected_sources(mesh,2)
    dict = [Dict(zip(source_connections[i],1:length(source_connections[i]))) for i = 1:length(source_connections)]

    lengths = length.(source_connections)
    idx = [0; cumsum(lengths)]
    # Preallocation of return values
    F = zeros(ComplexF64, idx[end])
    G = zeros(ComplexF64, idx[end])

    # mn = m*n
    mn = length(nodesX1)
    # mn_threads = mn*Base.Threads.nthreads()
    # Normals       = zeros(3, mn_threads)
    # Tangents      = zeros(3, mn_threads)
    # Sangents      = zeros(3, mn_threads)
    # Interpolation = zeros(3, mn_threads)
    # Jacobian      = zeros(mn_threads)
    # R             = zeros(mn_threads)
    # Integrand_    = zeros(ComplexF64,mn_threads,3)

    n_physics_functions = number_of_shape_functions(physics_function)
    n_shape_functions = number_of_shape_functions(shape_function)

    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source     = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        normals       = zeros(3, mn)
        tangents      = zeros(3, mn)
        sangents      = zeros(3, mn)
        interpolation = zeros(3, mn)
        jacobian      = zeros(mn)
        r             = zeros(mn)
        integrand     = zeros(ComplexF64, mn)
        di = dict[source_node]
        physics_nodes = zeros(Int64, n_physics_functions)
        # physics_nodes = [1;2;3]
        element_coordinates = zeros(3,n_shape_functions)
        @inbounds for element ∈ element_connections[source_node]
            # Access element topology and coordinates
            element_coordinates .= coordinates[:,topology[:,element]]
            # physics_nodes        .= [(idx[source_node] + di[t]) for t in physics_topology[:,element]]
            find_physics_nodes!(physics_nodes,idx[source_node],di,physics_topology[:,element])
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[physics_nodes]
            submatrixG = @view G[physics_nodes]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,element_coordinates)
            if close_corner == 1
                sparse_computing_integrals!(physics_interpolation1,shape_function1,
                                            normals, tangents,sangents,
                                            interpolation,jacobian,r,integrand,
                                            element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            elseif close_corner == 2
                sparse_computing_integrals!(physics_interpolation2,shape_function2,
                                            normals, tangents,sangents,
                                            interpolation,jacobian,r,integrand,
                                            element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            else
                sparse_computing_integrals!(physics_interpolation3,shape_function3,
                                            normals, tangents,sangents,
                                            interpolation,jacobian,r,integrand,
                                            element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            end
        end
        if progress; next!(prog); end # For the progress meter
    end

    I = create_row_indices(lengths,idx[end])
    J = vcat(source_connections...)

    return sparse(I,J,F), sparse(I,J,G), F
end
