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
    N = size(mesh.sources,2)
    T = zeros(eltype(topology), N)
    m,n = size(topology)
    for j=1:n
        for i=1:m
            T[topology[i,j]] += 1
        end
    end
    return T
end

#==========================================================================================
            Computing shape function derivatives at each node as the average
==========================================================================================#
"""
    shape_function_derivatives(mesh)

Returns derivative matrics of the 'physics_function', menaing that such that
∂p∂x = Dx*p, ∂p∂y = Dy*p, ∂p∂z = Dz*p & ∂pₙ∂x = Dx*pₙ, ∂pₙ∂y = Dy*pₙ, ∂pₙ∂z = Dz*pₙ,
"""
function shape_function_derivatives(mesh;global_derivatives=false)
    topology    = mesh.topology
    coordinates = mesh.coordinates
    n_elements  = number_of_elements(mesh)
    n_sources   = size(mesh.sources,2)

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
    global_gradients    = zeros(3,n_physics)   # ∇ₓ
    local_gradients     = zeros(3,n_physics)   # ∇ξ
    change_of_variables = zeros(3,3)         # COV-matrix
    dX  = zeros(3,n_physics) # X-Derivatives
    dY  = zeros(3,n_physics) # Y-Derivatives
    dZ  = zeros(3,n_physics) # Z-Derivatives

    element_coordinates = zeros(3,n_shape)
    derivatives_u = shape_function.derivatives_u
    derivatives_v = shape_function.derivatives_v

    _, source_connections = connected_sources(mesh,0)
    lengths = length.(source_connections)
    dict    = [Dict(zip(source_connections[i],1:lengths[i])) for i = 1:length(lengths)]

    # Preallocation of return values
    idx = [0; cumsum(lengths)]
    Dx = zeros(idx[end])
    Dy = zeros(idx[end])
    Dz = zeros(idx[end])

    for element = 1:n_elements
        # Extract element and compute
        element_coordinates  .= coordinates[:,topology[:,element]]
        my_mul!(dX,element_coordinates,derivatives_u)
        my_mul!(dY,element_coordinates,derivatives_v)
        cross_product!(dZ,dX,dY)

        source_nodes = @view physics_topology[:,element]
        # Add nodal gradients to the shape function derivatives
        physics_nodes = zeros(Int64, n_physics)
        for node = 1:n_physics
            source_node = source_nodes[node]
            di = dict[source_node]
            find_physics_nodes!(physics_nodes,idx[source_node],di,physics_topology[:,element])
            change_of_variables[1,:] = dX[:,node]
            change_of_variables[2,:] = dY[:,node]
            change_of_variables[3,:] = dZ[:,node]
            local_gradients[1,:]    .= physics_function.derivatives_u[:,node]
            local_gradients[2,:]    .= physics_function.derivatives_v[:,node]
            global_gradients        .= change_of_variables\local_gradients

            Dx[physics_nodes] += global_gradients[1,:] ./ n_elements_pr_node[source_node]
            Dy[physics_nodes] += global_gradients[2,:] ./ n_elements_pr_node[source_node]
            Dz[physics_nodes] += global_gradients[3,:] ./ n_elements_pr_node[source_node]
            # Dx[physics_nodes] += global_gradients[1,:]
            # Dy[physics_nodes] += global_gradients[2,:]
            # Dz[physics_nodes] += global_gradients[3,:]
        end
    end

    I = create_row_indices(lengths,idx[end])
    J = vcat(source_connections...)

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
        #
        # averanging = inverse_rle(n_elements_pr_node,lengths)
        # return sparse(I,J,Dt./averanging), sparse(I,J,Ds./averanging)
        return sparse(I,J,Dt), sparse(I,J,Ds)
    end
end

function mesh_gradients(mesh,p)
    Dx,Dy,Dz = global_coordinate_shape_function_derivative(mesh)
    return [p'*Dx'; p'*Dy'; p'*Dz']
end

function mesh_normal_derivative(mesh,p)
    mesh_gradient = mesh_gradients(mesh,p)
    n_grad = zeros(length(p))
    for i = 1:length(p)
        n_grad[i] = dot(mesh_gradient[:,i],mesh.normals[:,i])
    end
    return n_grad
end

"""
    getTangentialDerivativeMatrix!(physics_function,dZ,dX,dY)

Returns a matrix whose columns are the tangential derivatives of the interpolating nodes.
NOTE: It also normalized "dX" which it uses as the tangential direction of choice
"""
function get_tangential_derivative_matrix!(physics_function,dZ,dX,dY)
    # Finding the number of shape functions. Used to distribute derivatives/gradients
    n_physics = number_of_shape_functions(physics_function)
    n_gauss_points = length(physics_function.weights)

    # Pre-Allocation
    global_gradients      = MMatrix{3,n_physics}(zeros(3,n_physics)) # ∇ₓ
    local_gradients       = MMatrix{3,n_physics}(zeros(3,n_physics)) # ∇ξ
    change_of_variables   = MMatrix{3,3}(zeros(3,3))           # COV-matrix
    gauss_point_gradients = zeros(n_gauss_points,n_physics)         # ∂x,∂y,∂z

    for gauss_node = 1:n_gauss_points
        dXg = @view dX[:,gauss_node]
        dYg = @view dY[:,gauss_node]
        dZg = @view dZ[:,gauss_node]
        # Add nodal gradients to the shape function derivatives
        change_of_variables[1,:] = dXg
        change_of_variables[2,:] = dYg
        change_of_variables[3,:] = dZg

        # Computing global gradient
        local_gradients[1,:] = physics_function.derivatives_u[:,gauss_node]
        local_gradients[2,:] = physics_function.derivatives_v[:,gauss_node]
        global_gradients     = change_of_variables\local_gradients

        # Getting gradient in tangent direction
        tangents!(dZg,dXg,dYg)
        gauss_point_gradients[gauss_node:gauss_node,:] = dZg'*global_gradients
    end
    # Return tangential derivative matrix at each gauss point
    return gauss_point_gradients
end
