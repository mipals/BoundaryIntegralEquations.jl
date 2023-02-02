get_offset(physics_element::DiscontinuousTriangular)    = 0.1
get_offset(physics_element::ContinuousTriangular)       = 0.0
get_offset(physics_element::DiscontinuousQuadrilateral) = 0.0
get_offset(physics_element::ContinuousQuadrilateral)    = 0.1

get_beta(physicsElement::Triangular) = 0.0
get_beta(physicsElement::DiscontinuousTriangularConstant)  = 0.5
get_beta(physicsElement::DiscontinuousTriangularLinear)    = 0.05
get_beta(physicsElement::DiscontinuousTriangularQuadratic) = 0.05
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


function connected_topology(mesh)
    n_coordinates = size(mesh.coordinates,2)
    source_connections  = [zeros(Int64,0) for i = 1:n_coordinates]
    element_connections = [zeros(Int64,0) for i = 1:n_coordinates]
    # Maybe we should use actual topology for this to work?
    topology = mesh.topology
    n_physics_functions,n_elements = size(mesh.topology)
    for element = 1:n_elements
        for i = 1:3
            append!(element_connections[topology[i,element]],element)
        end
    end
    return sort!.(unique.(element_connections))
end

function connected_sources(mesh,depth,physics_function::T) where
        {T <: Union{DiscontinuousTriangular,DiscontinuousQuadrilateral}}
    cone = connected_topology(mesh)
    topology = mesh.topology
    physics_topology = mesh.physics_topology
    n_physics,n_elements = size(mesh.physics_topology)
    element_connections = [zeros(Int64,0) for i = 1:n_elements]
    n_corners = number_of_corners(physics_function)
    for element = 1:n_elements
        for i = 1:n_corners
            append!(element_connections[element], cone[topology[i,element]])
        end
    end

    n_sources = size(mesh.sources,2)
    source_connections  = [zeros(Int64,0) for i = 1:n_sources]
    for element = 1:n_elements
        for elm ∈ element_connections[element]
            for i = 1:n_physics
                append!(source_connections[physics_topology[i,element]],physics_topology[:,elm])
            end
        end
    end

    cone = inverse_rle(element_connections,n_physics*ones(Int64,length(element_connections)))

    return cone, sort!.(unique.(source_connections))
end


function connected_sources(mesh,depth=0)
    return connected_sources(mesh,depth,mesh.physics_function)
end

function connected_sources(mesh,depth,physics_function::T) where
    {T <: Union{ContinuousTriangular,ContinuousQuadrilateral}}
    @assert depth ∈ [0, 1, 2]
    n_sources = size(mesh.sources,2)
    source_connections  = [zeros(Int64,0) for i = 1:n_sources]
    element_connections = [zeros(Int64,0) for i = 1:n_sources]
    # Maybe we should use actual topology for this to work?
    topology = mesh.topology
    physics_topology = mesh.physics_topology
    n_corners = number_of_corners(mesh.physics_function)
    n_physics_functions,n_elements = size(mesh.physics_topology)
    for element = 1:n_elements
        for i = 1:n_physics_functions
            append!(element_connections[physics_topology[i,element]],element)
            append!(source_connections[physics_topology[i,element]],physics_topology[:,element])
        end
    end
    if depth == 1
        # conns = shape_connections(mesh.physics_function)
        # for element = 1:n_elements
        #     for idx = (n_corners + 1):n_physics_functions
        #         for j ∈ conns[idx-n_corners]
        #             append!(element_connections[physics_topology[idx,element]],element_connections[j,element])
        #             append!(source_connections[physics_topology[idx,element]],  source_connections[j,element])
        #         end
        #     end
        # end
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
    return sort!.(unique.(element_connections)), sort!.(unique.(source_connections))
end

function create_row_indices(lengths,total_counts)
    row_indices = zeros(Int64,total_counts)
    row = 1
    idx = 0
    for i ∈ lengths
        row_indices[idx+1:idx+i] .= row
        row = row + 1
        idx = idx + i
    end
    return row_indices
end

function mygemm_vec3!(C, B, A)
    @inbounds @fastmath for n ∈ eachindex(C)
        Cmn = zero(eltype(C))
        for k ∈ eachindex(A)
            Cmn += A[k] * B[n,k]
        end
        C[n] += Cmn
    end
end

function sparse_computing_integrals!(physics_interpolation,shape_function,
                                    normals, tangents,sangents,
                                    interpolation,jacobian_mul_weights,r,integrand,
                                    element_coordinates,
                                    fOn,gOn,submatrixF,submatrixG,k,source)

    mul!(interpolation,element_coordinates,shape_function.interpolation)
    jacobian!(shape_function,element_coordinates,
                normals,tangents,sangents,jacobian_mul_weights)
    integrand_mul!(jacobian_mul_weights,shape_function.weights)
    compute_distances!(r,interpolation,source)
    # jacobian_openbem!(shape_function,element_coordinates,normals,tangents,sangents,jacobian_mul_weights)
    # jacobian_mul_weights .*= shape_function.weights
    if fOn
        ### Evaluating the F-kernel (double-layer kernel) at the global nodes
        freens3d!(integrand,r,interpolation,source,normals,k)
        # F3d!(interpolation,source,k,normals,r,integrand,-1)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # integrand .= integrand .* jacobian_mul_weights
        # Approximating integral and adding the value to the BEM matrix
        mygemm_vec3!(submatrixF,physics_interpolation,integrand)
        # submatrixF .+= physics_interpolation*integrand
        # mul!(submatrixF,physics_interpolation,integrand,true,true)
    end
    if gOn
        ### Evaluating the G-kernel (single-layer kernel) at the global nodes
        greens3d!(integrand,r,k)
        integrand_mul!(integrand,jacobian_mul_weights)
        mygemm_vec3!(submatrixG,physics_interpolation,integrand)
        # G3d!(k,r,integrand,-1)
        # Multiplying integrand with jacobian and weights
        # integrand .= integrand .* jacobian_mul_weights
        # Approximating the integral and the adding the values to the BEM matrix
        # submatrixG .+= physics_interpolation*integrand
        # mul!(submatrixG,physics_interpolation,integrand,true,true)
    end
end

function find_physics_nodes!(physics_nodes,idxs,di,physics_top)
    @inbounds for i = 1:length(physics_nodes)
        physics_nodes[i] = idxs + di[physics_top[i]]
    end
    # [(idx[source_node] + di[t]) for t in physics_topology[:,element]]
end

function sparse_assemble_parallel!(mesh::Mesh3d,k,sources,shape_function::T;
        fOn=true,gOn=true,depth=1,offset=nothing,
        progress=true) where {T <: Union{TriangularLinear,DiscontinuousTriangularLinear}}
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
    beta = get_beta(physics_function)
    tmp  = DiscontinuousTriangularLinear(physics_function,beta)
    if typeof(offset) <: Real
        offr = 0.0
    else
        offr = get_offset(physics_function)
    end
    # offr = 0.0
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
    physics_interpolation1 = physics_function1.interpolation
    physics_interpolation2 = physics_function2.interpolation
    physics_interpolation3 = physics_function3.interpolation

    # Connections
    element_connections, source_connections = connected_sources(mesh,depth)
    lengths = length.(source_connections)
    dict = [Dict(zip(source_connections[i],1:lengths[i])) for i = 1:length(lengths)]

    idx = [0; cumsum(lengths)]
    # Preallocation of return values
    F = zeros(ComplexF64, idx[end])
    G = zeros(ComplexF64, idx[end])

    mn = length(nodesX1)
    n_physics_functions = number_of_shape_functions(physics_function)
    n_shape_functions   = number_of_shape_functions(shape_function)

    n_threads     = nthreads()
    Normals       = zeros(3, mn*n_threads)
    Tangents      = zeros(3, mn*n_threads)
    Sangents      = zeros(3, mn*n_threads)
    Interpolation = zeros(3, mn*n_threads)
    Jacobian      = zeros(mn,n_threads)
    R             = zeros(mn,n_threads)
    Integrand     = zeros(ComplexF64, mn, n_threads)
    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source     = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        normals       = @view Normals[:,((threadid()-1)*mn+1):threadid()*mn]
        tangents      = @view Tangents[:,((threadid()-1)*mn+1):threadid()*mn]
        sangents      = @view Sangents[:,((threadid()-1)*mn+1):threadid()*mn]
        interpolation = @view Interpolation[:,((threadid()-1)*mn+1):threadid()*mn]
        jacobian      = @view Jacobian[:,threadid()]
        r             = @view R[:,threadid()]
        integrand     = @view Integrand[:,threadid()]
        di = dict[source_node]
        physics_nodes = zeros(Int64, n_physics_functions)
        element_coordinates = zeros(3,n_shape_functions)
        @inbounds for element ∈ element_connections[source_node]
            # Access element topology and coordinates
            element_coordinates .= coordinates[:,topology[:,element]]
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

    return sparse(I,J,F), sparse(I,J,G)
end


function sparse_assemble_parallel!(mesh::Mesh3d,k,sources,shape_function::T;
    fOn=true,gOn=true,depth=1,offset=nothing,
    progress=true) where {T <: Union{TriangularQuadratic,DiscontinuousTriangularQuadratic}}
    topology    = get_topology(mesh)
    n_sources   = size(sources,2)
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
    beta = get_beta(physics_function)
    tmp  = DiscontinuousTriangularLinear(physics_function,beta)
    if typeof(offset) <: Real
        offr = offset
    else
        offr = get_offset(physics_function)
    end
    # Dealing with corners
    DO = Diagonal(ones(3))
    nodesX1,nodesY1,weights1 = singular_triangle_integration(tmp,3,[1.00;offr;offr],DO,1e-6)
    nodesX2,nodesY2,weights2 = singular_triangle_integration(tmp,3,[offr;1.00;offr],DO,1e-6)
    nodesX3,nodesY3,weights3 = singular_triangle_integration(tmp,3,[offr;offr;1.00],DO,1e-6)
    nodesX4,nodesY4,weights4 = singular_triangle_integration(tmp,3,[0.50;0.50;offr],DO,1e-6)
    nodesX5,nodesY5,weights5 = singular_triangle_integration(tmp,3,[offr;0.50;0.50],DO,1e-6)
    nodesX6,nodesY6,weights6 = singular_triangle_integration(tmp,3,[0.50;offr;0.50],DO,1e-6)
    shape_function1 = deepcopy(shape_function)
    shape_function2 = deepcopy(shape_function)
    shape_function3 = deepcopy(shape_function)
    shape_function4 = deepcopy(shape_function)
    shape_function5 = deepcopy(shape_function)
    shape_function6 = deepcopy(shape_function)
    set_interpolation_nodes!(shape_function1,nodesX1,nodesY1,weights1)
    set_interpolation_nodes!(shape_function2,nodesX2,nodesY2,weights2)
    set_interpolation_nodes!(shape_function3,nodesX3,nodesY3,weights3)
    set_interpolation_nodes!(shape_function4,nodesX4,nodesY4,weights4)
    set_interpolation_nodes!(shape_function5,nodesX5,nodesY5,weights5)
    set_interpolation_nodes!(shape_function6,nodesX6,nodesY6,weights6)
    physics_function1 = deepcopy(physics_function)
    physics_function2 = deepcopy(physics_function)
    physics_function3 = deepcopy(physics_function)
    physics_function4 = deepcopy(physics_function)
    physics_function5 = deepcopy(physics_function)
    physics_function6 = deepcopy(physics_function)
    set_interpolation_nodes!(physics_function1,nodesX1,nodesY1,weights1)
    set_interpolation_nodes!(physics_function2,nodesX2,nodesY2,weights2)
    set_interpolation_nodes!(physics_function3,nodesX3,nodesY3,weights3)
    set_interpolation_nodes!(physics_function4,nodesX4,nodesY4,weights4)
    set_interpolation_nodes!(physics_function5,nodesX5,nodesY5,weights5)
    set_interpolation_nodes!(physics_function6,nodesX6,nodesY6,weights6)

    # Copying interpolation of physics functions1
    physics_interpolation1 = physics_function1.interpolation
    physics_interpolation2 = physics_function2.interpolation
    physics_interpolation3 = physics_function3.interpolation
    physics_interpolation4 = physics_function4.interpolation
    physics_interpolation5 = physics_function5.interpolation
    physics_interpolation6 = physics_function6.interpolation

    # Connections
    element_connections, source_connections = connected_sources(mesh,depth)
    lengths = length.(source_connections)
    dict = [Dict(zip(source_connections[i],1:lengths[i])) for i = 1:length(lengths)]

    idx = [0; cumsum(lengths)]
    # Preallocation of return values
    F = zeros(ComplexF64, idx[end])
    G = zeros(ComplexF64, idx[end])

    # Preallocation of intermediate arrays
    n_threads = nthreads()
    mn = length(nodesX1)
    Corner_normals       = zeros(3, mn*n_threads)
    Corner_tangents      = zeros(3, mn*n_threads)
    Corner_sangents      = zeros(3, mn*n_threads)
    Corner_interpolation = zeros(3, mn*n_threads)
    Corner_jacobian      = zeros(mn,n_threads)
    Corner_r             = zeros(mn,n_threads)
    Corner_integrand     = zeros(ComplexF64, mn, n_threads)

    mp = length(nodesX4)
    Middle_normals       = zeros(3, mp*n_threads)
    Middle_tangents      = zeros(3, mp*n_threads)
    Middle_sangents      = zeros(3, mp*n_threads)
    Middle_interpolation = zeros(3, mp*n_threads)
    Middle_jacobian      = zeros(mp,n_threads)
    Middle_r             = zeros(mp,n_threads)
    Middle_integrand     = zeros(ComplexF64, mp, n_threads)

    n_physics_functions = number_of_shape_functions(physics_function)
    n_shape_functions   = number_of_shape_functions(shape_function)

    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source     = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        corner_normals       = @view Corner_normals[:,((threadid()-1)*mn+1):threadid()*mn]
        corner_tangents      = @view Corner_tangents[:,((threadid()-1)*mn+1):threadid()*mn]
        corner_sangents      = @view Corner_sangents[:,((threadid()-1)*mn+1):threadid()*mn]
        corner_interpolation = @view Corner_interpolation[:,((threadid()-1)*mn+1):threadid()*mn]
        corner_jacobian      = @view Corner_jacobian[:,threadid()]
        corner_r             = @view Corner_r[:,threadid()]
        corner_integrand     = @view Corner_integrand[:,threadid()]

        middle_normals       = @view Middle_normals[:,((threadid()-1)*mp+1):threadid()*mp]
        middle_tangents      = @view Middle_tangents[:,((threadid()-1)*mp+1):threadid()*mp]
        middle_sangents      = @view Middle_sangents[:,((threadid()-1)*mp+1):threadid()*mp]
        middle_interpolation = @view Middle_interpolation[:,((threadid()-1)*mp+1):threadid()*mp]
        middle_jacobian      = @view Middle_jacobian[:,threadid()]
        middle_r             = @view Middle_r[:,threadid()]
        middle_integrand     = @view Middle_integrand[:,threadid()]

        di = dict[source_node]
        physics_nodes = zeros(Int64, n_physics_functions)
        # physics_nodes = [1;2;3]
        element_coordinates = zeros(3,n_shape_functions)
        @inbounds for element ∈ element_connections[source_node]
            # Access element topology and coordinates
            element_coordinates .= coordinates[:,topology[:,element]]
            find_physics_nodes!(physics_nodes,idx[source_node],di,physics_topology[:,element])
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[physics_nodes]
            submatrixG = @view G[physics_nodes]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,element_coordinates)
            if close_corner == 1
                sparse_computing_integrals!(physics_interpolation1,shape_function1,
                                            corner_normals,corner_tangents,corner_sangents,
                                            corner_interpolation,corner_jacobian,corner_r,
                                            corner_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            elseif close_corner == 2
                sparse_computing_integrals!(physics_interpolation2,shape_function2,
                                            corner_normals,corner_tangents,corner_sangents,
                                            corner_interpolation,corner_jacobian,corner_r,
                                            corner_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            elseif close_corner == 3
                sparse_computing_integrals!(physics_interpolation3,shape_function3,
                                            corner_normals,corner_tangents,corner_sangents,
                                            corner_interpolation,corner_jacobian,corner_r,
                                            corner_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            elseif close_corner == 4
                sparse_computing_integrals!(physics_interpolation4,shape_function4,
                                            middle_normals,middle_tangents,middle_sangents,
                                            middle_interpolation,middle_jacobian,middle_r,
                                            middle_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            elseif close_corner == 5
                sparse_computing_integrals!(physics_interpolation5,shape_function5,
                                            middle_normals,middle_tangents,middle_sangents,
                                            middle_interpolation,middle_jacobian,middle_r,
                                            middle_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            else
                sparse_computing_integrals!(physics_interpolation6,shape_function6,
                                            middle_normals,middle_tangents,middle_sangents,
                                            middle_interpolation,middle_jacobian,middle_r,
                                            middle_integrand,element_coordinates,
                                            fOn,gOn,submatrixF,submatrixG,k,source)
            end
        end
        if progress; next!(prog); end # For the progress meter
    end

    I = create_row_indices(lengths,idx[end])
    J = vcat(source_connections...)

    return sparse(I,J,F), sparse(I,J,G)
end