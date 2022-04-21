"""
    cross!(normal::AbstractVector,a::AbstractVector,b::AbstractVector)

Inplace computation of the cross product `a × b`. Result saved in normal.
"""
function cross!(normal::AbstractVector{T},a::AbstractVector{T},b::AbstractVector{T}) where {T}
    normal[1] = a[2]*b[3]-a[3]*b[2]
    normal[2] = a[3]*b[1]-a[1]*b[3]
    normal[3] = a[1]*b[2]-a[2]*b[1]
end
"""
    cross_product!(normals,tangent,sangent)

Inplace computation of the cross product between every column of `tangent` and `sangent`.
Results are saved in the corresponding column of `normals`.
"""
function cross_product!(normals,tangent,sangent)
    @inbounds for i = 1:size(tangent,2)
        normal  = @view normals[:,i]
        dX      = @view tangent[:,i]
        dY      = @view sangent[:,i]
        cross!(normal,dX,dY)
    end
end
"""
    column_norms!(jacobian,normals)

Computing norms of columns in `normals`. Saved in `jacobian`. Everything assumed to be 3d.
"""
function column_norms!(jacobian,normals)
    @inbounds for i = 1:length(jacobian)
        jacobian[i] = hypot(normals[1,i],normals[2,i],normals[3,i])
    end
end
"""
    normalize!(normals,jacobian)

Dividing each column of `normals` with `jacobian`.
"""
function normalize!(normals,jacobian)
    @inbounds for i = 1:length(jacobian)
        normals[1,i] = -normals[1,i]/jacobian[i]
        normals[2,i] = -normals[2,i]/jacobian[i]
        normals[3,i] = -normals[3,i]/jacobian[i]
    end
end
"""
    jacobian!(basisElement::SurfaceFunction,coordinates,normals,tangent,sangent,jacobian)

Inplace computations of the jacobian of the basisElement::SurfaceFunction at the coordinates.
During the computations of the tangent- and normal-vectors is used and also saved. 
The results are saved in `tangent`, `sangent`, `normals` and `jacobian`.
"""
function jacobian!(basisElement::SurfaceFunction,
                   coordinates,normals,tangent,sangent,jacobian)
    mul!(tangent,coordinates,basisElement.derivatives_u)    # Computing tangent vector in X
    mul!(sangent,coordinates,basisElement.derivatives_v)    # Computing tangent vector in Y
    cross_product!(normals,tangent,sangent)                 # Computing normal vector
    column_norms!(jacobian,normals)                         # Jacobian = length of the normal
    normalize!(normals,jacobian)
end
