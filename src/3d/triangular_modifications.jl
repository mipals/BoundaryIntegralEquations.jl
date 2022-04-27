
"""
    create_rotated_element(basisElement::Triangular,m,n,clusterCorner)

Creates a rotated triangular element of the same as the input basisElement.
"""
function create_rotated_element(shape_function,n::Real,m::Real,clusterCorner)
    if clusterCorner == 1
        nodes_u,nodes_v,weights = rotated_triangular_quadpoints(n,m)
    elseif clusterCorner == 2
        nodes_u,nodes_v,weights = triangularQuadpoints(n,m)
    elseif clusterCorner == 3
        nodes_v,nodes_u,weights = triangularQuadpoints(n,m)
    end
    rotated_element = deepcopy(shape_function)
    rotated_element.gauss_u = nodes_u
    rotated_element.gauss_v = nodes_v
    rotated_element.weights = weights
    interpolate_on_nodes!(rotated_element)
    return rotated_element
end
