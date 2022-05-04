#==========================================================================================
                                Utility functions
==========================================================================================#
"""
    count_elements_pr_node(mesh)

Counting the no. elements that each node is part of
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
    shapeFunctionDerivative(mesh)

Computes the tangential derivative using the shape function derivative.
"""
function shapeFunctionDerivative(mesh)
    Dx,Dy,Dz = global_coordinate_shape_function_derivative(mesh)

    Ts = mesh.tangentsX
    Tt = mesh.tangentsY

    return  Dx .* Ts[1,:] + Dy .* Ts[2,:] + Dz .* Ts[3,:],
            Dx .* Tt[1,:] + Dy .* Tt[2,:] + Dz .* Tt[3,:]
end

"""
    global_coordinate_shape_function_derivative(mesh)

Returns derivative matrics of the 'physics_function', menaing that such that
∂p∂x = Dx*p, ∂p∂y = Dy*p, ∂p∂z = Dz*p & ∂pₙ∂x = Dx*pₙ, ∂pₙ∂y = Dy*pₙ, ∂pₙ∂z = Dz*pₙ, 
"""

function global_coordinate_shape_function_derivative(mesh)
    topology    = mesh.topology
    coordinates = mesh.coordinates
    nElements   = number_of_elements(mesh)
    n_sources   = size(mesh.sources,2)

    # Making a copy of the element type
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)

    # Setting interpolations to be at the nodal positions 
    # (since this is where we want to compute the derivatives)
    set_nodal_interpolation!(physics_function)
    copy_interpolation_nodes!(shape_function,physics_function)
    # Finding the number of shape functions. Used to distribute derivatives/gradients
    n_shape = number_of_shape_functions(physics_function)
    
    # Pre-Allocation 
    global_gradients    = zeros(3,n_shape)   # ∇ₓ 
    local_gradients     = zeros(3,n_shape)   # ∇ξ
    change_of_variables = zeros(3,3)         # COV-matrix
    dX  = zeros(3,n_shape) # X-Derivatives
    dY  = zeros(3,n_shape) # Y-Derivatives
    dZ  = zeros(3,n_shape) # Z-Derivatives

    Dx  = spzeros(n_sources,n_sources)  # ∂x 
    Dy  = spzeros(n_sources,n_sources)  # ∂y 
    Dz  = spzeros(n_sources,n_sources)  # ∂z 
    n_elements_pr_node = count_elements_pr_node(mesh) 
    for element = 1:nElements
        # Extract element and compute 
        element_nodes        = @view topology[:,element]
        element_coordinates  = @view coordinates[:,element_nodes]
        mul!(dX,element_coordinates,shape_function.derivatives_u)
        mul!(dY,element_coordinates,shape_function.derivatives_v)
        cross_product!(dZ,dX,dY)

        physics_nodes = @view mesh.physics_topology[:,element]
        # Add nodal gradients to the shape function derivatives
        for node = 1:n_shape
            change_of_variables[1,:] = dX[:,node]
            change_of_variables[2,:] = dY[:,node]
            change_of_variables[3,:] = dZ[:,node]
            local_gradients[1,:] .= physics_function.derivatives_u[:,node]
            local_gradients[2,:] .= physics_function.derivatives_v[:,node]
            global_gradients     .= change_of_variables\local_gradients
            
            Dx[physics_nodes[node],physics_nodes] += global_gradients[1,:] ./ 
                                                     n_elements_pr_node[physics_nodes[node]]
            Dy[physics_nodes[node],physics_nodes] += global_gradients[2,:] ./ 
                                                     n_elements_pr_node[physics_nodes[node]]
            Dz[physics_nodes[node],physics_nodes] += global_gradients[3,:] ./ 
                                                     n_elements_pr_node[physics_nodes[node]]
        end
    end

    return Dx, Dy, Dz

end

"""
    getTangentialDerivativeMatrix(physics_function,dZ,dX,dY)

Returns a matrix whose columns are the tangential derivatives of the interpolating nodes.
NOTE: It also normalized "dX" which it uses as the tangential direction of choice
"""
function get_tangential_derivative_matrix!(physics_function,dZ,dX,dY)
    # Finding the number of shape functions. Used to distribute derivatives/gradients
    n_shape = number_of_shape_functions(physics_function) 
    n_gauss_points = length(physics_function.weights)
    
    # Pre-Allocation 
    global_gradients      = MMatrix{3,n_shape}(zeros(3,n_shape)) # ∇ₓ 
    local_gradients       = MMatrix{3,n_shape}(zeros(3,n_shape)) # ∇ξ
    change_of_variables   = MMatrix{3,3}(zeros(3,3))           # COV-matrix
    gauss_point_gradients = zeros(n_gauss_points,n_shape)         # ∂x,∂y,∂z 

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
