#==========================================================================================
                    Qudrature points for the quadrilateral elements                     
==========================================================================================#
"""
    quadrilateralQuadpoints(n,m)

Returns `m*n` integration points and `weights` for a quadrilateral element.
"""
function quadrilateralQuadpoints(n,m)
    nodex, wx = gausslegendre(n) # Nodes in the interval [-1, 1]
    nodey, wy = gausslegendre(m) # Nodes in the interval [-1, 1]

    ξ = kron(ones(m),nodex)
    η = kron(nodey,ones(n))
    weights = kron(wy,wx)
    return ξ,η,weights
end
function getQuadpoints(elementType::Quadrilateral,n=4,m=4)
    return quadrilateralQuadpoints(n,m)
end


#==========================================================================================
                     Qudrature points for triangular elements                        
==========================================================================================#
"""
    duffyBasis(u,v)

Computes the duffy transformation of `u` and `v`.
"""
function duffyBasis(u,v)
    return u, v.*(1.0 .- u)
end

"""
    duffyBasisDerivative(u,v)

Computes the derivative of the duffy transformation of `u` and `v`.
"""
function duffyBasisDerivative(u,v)
    return [ones(1,length(v)); -v],[zeros(1,length(u)); 1.0 .- u]
end

"""
    duffyJacobian(u,v)

Computes the Jacobian the duffy transformation of `u` and `v`.

This determiant can be calculated analytically
```math
                | 1       0 |
                  | -v  1 - u | = 1 - u  
```
"""
function duffyJacobian(u,v)
    return abs.(1.0 .- u)
end

"""
    triangularQuadpoints(n=4,m=4)

Returns `m*n` integration points and `weights` for a triangule using a duffy transformation.

The integration points are clustered around node "3" in the triangle:

                                        3          
                                        | \\
                                        1 - 2 
"""
function triangularQuadpoints(n,m)
    nodex, wx = curveLinearQuadpoints(n) # Nodes in the interval [0, 1]
    nodey, wy = curveLinearQuadpoints(m) # Nodes in the interval [0, 1]

    # Nodes in quadrilateral element (Not we can not )
    u = kron(ones(m),nodex)
    v = kron(nodey,ones(n))
    w = kron(wy,wx)

    # Transforming the quadrilateral element to an triangular element
    ξ,η = duffyBasis(u,v)
    weights = w .* duffyJacobian(u,v)

    return ξ,η,weights
end
function getQuadpoints(elementType::Triangular,n=4,m=4)
    return triangularQuadpoints(n,m)
end

#==========================================================================================
                                 Rotation of integration points 
==========================================================================================#
"""
    rotatedTriangularQuadpoints(n=4,m=4)

Returns `m*n` integration points and `weights` for a triangule using a duffy transformation.

The integration points are clustered around node "1" in the triangle:

                                    3          
                                    | \
                                    1 - 2 
"""
function rotated_triangular_quadpoints(n=4,m=4)
    nodex, wx = curveLinearQuadpoints(n) # Nodes in the interval [0, 1]
    nodey, wy = curveLinearQuadpoints(m) # Nodes in the interval [0, 1]

    u = kron(ones(m),nodex)
    v = kron(nodey,ones(n))
    w = kron(wy,wx)

    # Transforming the quadrilateral element to an triangular element
    ξ = u .* v
    η = (1.0 .- u) .* v
    weights = w .* v

    return ξ,η,weights
end
