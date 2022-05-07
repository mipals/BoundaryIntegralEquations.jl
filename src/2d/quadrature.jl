#==========================================================================================
                                Computing Quadrature Points 
==========================================================================================#
# Scaling of legendre points to the interval [0,1]
function curveLinearQuadpoints(n=4)
    node, w = legendre(n)
    w = 0.5*w;
    node = 0.5*(node .+ 1.0)

    return node, w
end
function getQuadpoints(elementType::CurveLinear,n=4)
    return curveLinearQuadpoints(n)
end
# Standard legendre quadrature (interval [-1,1])
function curveQuadraticQuadpoints(n=4)
    node, w = legendre(n)
    return node, w
end
function getQuadpoints(elementType::CurveQuadratic,n=4)
    return curveQuadraticQuadpoints(n)
end
