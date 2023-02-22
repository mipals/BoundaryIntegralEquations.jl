#==========================================================================================
                                Utility functions
==========================================================================================#
"""
    count_elements_pr_node(mesh)

Counting the number of elements that each node is part of
(Cornes can have various of numbers. Midpoints only 2.)
"""
function count_elements_pr_node(mesh)
    topology = mesh.physics_topology
    T = zeros(eltype(topology), size(mesh.sources,2))
    for j=axes(topology,2), i=axes(topology,1)
        T[topology[i,j]] += 1
    end
    return T
end

#==========================================================================================
            Computing shape function derivatives at each node as the average
==========================================================================================#
"""
    interpolation_function_derivatives(mesh)

Returns interpolation function derivative matrics ``(\\mathbf{D}_x,\\mathbf{D}_y,\\mathbf{D}_z)`` of the 'physics_function'.
This means that given the BEM interpolation of p the derivatives are given by
```math
    \\frac{\\partial\\mathbf{p}^\\parallel}{\\partial x} = \\mathbf{D}_x\\mathbf{p}, \\newline
    \\frac{\\partial\\mathbf{p}^\\parallel}{\\partial y} = \\mathbf{D}_y\\mathbf{p}, \\newline
    \\frac{\\partial\\mathbf{p}^\\parallel}{\\partial z} = \\mathbf{D}_z\\mathbf{p}.
```
NOTE! This is only the tangential part of the gradient.
The reason being that the BEM interpolation on deal with surface values, and does
therefore not contain the derivative information orthogonal to it.
"""
function interpolation_function_derivatives(mesh;global_derivatives=true)
    topology    = mesh.topology
    coordinates = mesh.coordinates
    n_elements  = number_of_elements(mesh)

    # Making a copy of the element type
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    physics_topology = mesh.physics_topology
    # Setting interpolations to be at the nodal positions
    # (since this is where we want to compute the derivatives)
    set_nodal_interpolation!(physics_function)
    copy_interpolation_nodes!(shape_function,physics_function)
    # Finding the number of shape functions. Used to distribute derivatives/gradients
    n_physics = number_of_shape_functions(physics_function)
    n_shape   = number_of_shape_functions(shape_function)
    n_elements_pr_node = count_elements_pr_node(mesh)
    # Pre-Allocation
    global_gradients    = zeros(3,n_physics) # ∇ₓ
    local_gradients     = zeros(3,n_physics) # ∇ξ
    change_of_variables = zeros(3,3)         # COV-matrix
    dX = zeros(3,n_physics) # X-Derivatives
    dY = zeros(3,n_physics) # Y-Derivatives
    dZ = zeros(3,n_physics) # Z-Derivatives

    element_coordinates = zeros(3,n_shape)
    derivatives_u = shape_function.derivatives_u
    derivatives_v = shape_function.derivatives_v

    _, source_connections = connected_sources(mesh,1)
    lengths = length.(source_connections)
    dict    = [Dict(zip(source_connections[i],1:lengths[i])) for i = 1:length(lengths)]

    # Preallocation of return values
    idx = [0; cumsum(lengths)]
    Dx = zeros(idx[end])
    Dy = zeros(idx[end])
    Dz = zeros(idx[end])

    for element = 1:n_elements
        # Moving values for memory locality in multiplication
        element_coordinates  .= coordinates[:,topology[:,element]]
        # Interpolation the derivatives
        my_mul!(dX,element_coordinates,derivatives_u)
        my_mul!(dY,element_coordinates,derivatives_v)
        cross_product!(dZ,dX,dY) # Computing last row of change-of-variables

        # Add nodal gradients to the shape function derivatives
        physics_nodes = zeros(Int64, n_physics)
        for (node,source_node) ∈ enumerate(physics_topology[:,element])
            di = dict[source_node]
            find_physics_nodes!(physics_nodes,idx[source_node],di,physics_topology[:,element])
            # Creating change-of-variables matrix
            change_of_variables[1,:] .= dX[:,node]
            change_of_variables[2,:] .= dY[:,node]
            change_of_variables[3,:] .= dZ[:,node]
            # Inserting local derivatives into vector
            local_gradients[1,:]     .= physics_function.derivatives_u[:,node]
            local_gradients[2,:]     .= physics_function.derivatives_v[:,node]
            # Compute global TANGENTIAL gradient
            global_gradients         .= change_of_variables\local_gradients
            # Add result to output
            Dx[physics_nodes] += global_gradients[1,:] ./ n_elements_pr_node[source_node]
            Dy[physics_nodes] += global_gradients[2,:] ./ n_elements_pr_node[source_node]
            Dz[physics_nodes] += global_gradients[3,:] ./ n_elements_pr_node[source_node]
        end
    end

    # Compute indices
    I = inverse_rle(1:length(lengths),lengths)
    J = vcat(source_connections...)

    # Global derivatives is really what should be used. The local are only here for
    # backwards compatability
    if global_derivatives
        return sparse(I,J,Dx), sparse(I,J,Dy), sparse(I,J,Dz)
    else
        # Repeating coordinate `i` lengths[i] times using StatsBase function `inverse_rle`
        T1 = inverse_rle(mesh.tangents[1,:],lengths)
        T2 = inverse_rle(mesh.tangents[2,:],lengths)
        T3 = inverse_rle(mesh.tangents[3,:],lengths)
        S1 = inverse_rle(mesh.sangents[1,:],lengths)
        S2 = inverse_rle(mesh.sangents[2,:],lengths)
        S3 = inverse_rle(mesh.sangents[3,:],lengths)
        # From global to local coordinates (this is essentially directional derivatives)
        Dt = Dx .* T1 + Dy .* T2 + Dz .* T3
        Ds = Dx .* S1 + Dy .* S2 + Dz .* S3
        return sparse(I,J,Dt), sparse(I,J,Ds)
    end
end
