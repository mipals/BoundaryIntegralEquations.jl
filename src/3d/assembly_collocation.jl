struct SurfaceElement{T<:AbstractFloat}
    jacobians::AbstractArray{T,1}
    interpolation::AbstractArray{T,2}
    normals::AbstractArray{T,2}
end

function interpolate_elements(mesh::Mesh3d)
    nElements = size(mesh.topology,2)
    shape_function = mesh.shape_function

    element_interpolations = Array{SurfaceElement{eltype(shape_function)}}(undef,nElements)

    normals  = zeros(3,length(shape_function.weights))
    tangents = similar(normals)
    sangents = similar(normals)
    interpolations = similar(normals)
    jacobians = similar(shape_function.weights)

    for element = 1:nElements
        # Extracting element properties (connectivity and coordinates)
        elementCoordinates  = @view mesh.coordinates[:,mesh.topology[:,element]]
        # Computing interpolation
        mul!(interpolations,elementCoordinates,shape_function.interpolation)
        # Computing tangential directions as well a a normal at each node
        jacobian!(shape_function,elementCoordinates,normals,tangents,sangents,jacobians)
        # Save element interpolations
        element_interpolations[element] = SurfaceElement(deepcopy(jacobians),
                                                        deepcopy(interpolations),
                                                        deepcopy(normals))
    end
    return element_interpolations
end
