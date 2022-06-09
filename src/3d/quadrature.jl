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
                                    | \\
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

#==========================================================================================
                                Spherical Coordinates
==========================================================================================#
"""
    linear_quadrature_transformation(gp,w,a,b)

Transforms Guassian quadrature points `gp` and weights `w` from the interval [-1,1] to [a,b]
"""
function linear_quadrature_transformation(gp,w,a,b)
    gp = (a .+ b)/2.0 .+ (b .- a)/2.0.*gp
    w  = (b .- a)/2.0.*w
    return gp,w
end
max_rad1(x) = 1.0 ./ (cos.(x) + sin.(x))             #
max_rad2(x) = 0.5 * sec.(x)                          # Per. definition of secant: 1/cos(x)
max_rad3(x) = sin.(π/4.0) ./ sin.(x .- π/4.0)*0.5    # Sinus relation
"""
    polar_quadrature

Returns the `n_rad` radii and `n_ang` angles of polar coordinates of a triangle defined by
the values `minAngle`, `maxAngle`, `minRadius` and the `maxRadius` function.
"""
function polar_quadrature(n_rad,n_ang=0,minAngle=0.0,maxAngle=π/2.0,minRadius=0.0,maxRad=max_rad1)
    # Getting standard Gaussian Points
    gp_ang, w_ang = gausslegendre(n_ang)
    gp_rad, w_rad = gausslegendre(n_rad)

    # Angles from 0 to π/2
    theta,thetaW = linear_quadrature_transformation(gp_ang,w_ang,minAngle,maxAngle)
    thetaW = kron(thetaW,ones(n_rad))

    # The corresponding radius depends on the angles.
    maxRadius = maxRad(theta)
    radius  = zeros(n_ang*n_rad)
    radiusW = zeros(n_ang*n_rad)

    for i = 1:n_ang
        rad,radW = linear_quadrature_transformation(gp_rad,w_rad,minRadius,maxRadius[i])
        radius[(i-1)*n_rad + 1:i*n_rad]  = rad
        radiusW[(i-1)*n_rad + 1:i*n_rad] = radW .* rad # rad is the jacobian
    end
    angles = kron(theta,ones(n_rad))
    w      = thetaW .* radiusW
    return angles,radius,w
end
"""
    polar_gaussian(n)

Returns `n^2` x,y-coordinates and weights for quadrature on a triangle with a singularity
on vertex 1.
"""
function polar_gaussian(n_rad)
    n_rad = Int(ceil(n_rad/2.0))*2
    n_ang = n_rad

    # NOTE: Things go anti-clockwise here
    angles,radius,w = polar_quadrature(n_rad,n_ang,0.0,π/2,0,max_rad1)
    x = radius .* cos.(angles)
    y = radius .* sin.(angles)

    return x,y,w

end
"""
    polar_gaussianMidpoint(n)

Returns `n^2` x,y-coordinates and weights for quadrature on a triangle with a singularity
on vertex 4 (i.e. vertex between vertex 1 and 2).
"""
function polar_gaussianMidpoint(n_rad)
    n_rad = Int(ceil(n_rad/2.0))*2
    n_ang = Int(n_rad/2)

    # NOTE: The angles go clockwise here
    angles1,radius1,w1 = polar_quadrature(n_rad,n_ang,0.0,atan(2),0.0,max_rad2)
    x1 = radius1 .* cos.(π .- angles1)
    y1 = radius1 .* sin.(π .- angles1)

    angles2,radius2,w2 = polar_quadrature(n_rad,n_ang,atan(2),π,0.0,max_rad3)
    x2 = radius2 .* cos.(π .- angles2)
    y2 = radius2 .* sin.(π .- angles2)

    return [x1;x2] .+ 0.5, [y1;y2], [w1; w2]
end

"""
    getpolar_gaussian(n,vertexNumber)

Returns `n^2` x,y-coordinates and weights for quadrature on a triangle with a singularity
on vertex number `vertexNumber`. Currently only vertex numbers between 1 and 6 is supported.
"""
function getpolar_gaussian(n,vertexNumber)
    if vertexNumber == 1
        return polar_gaussian(n)
    elseif vertexNumber == 2
        X,Y,W = polar_gaussian(n)
        return 1.0 .- X - Y, X, W
    elseif vertexNumber == 3
        X,Y,W = polar_gaussian(n)
        return Y, 1.0 .- X .- Y,W
    elseif vertexNumber == 4
        return polar_gaussianMidpoint(n)
    elseif vertexNumber == 5
        X,Y,W = polar_gaussianMidpoint(n)
        return 1.0 .- X - Y, X, W
    elseif vertexNumber == 6
        X,Y,W = polar_gaussianMidpoint(n)
        return Y, 1.0 .- X .- Y,W
    else
        error("Node number does not lie in range 1-6")
    end

end
