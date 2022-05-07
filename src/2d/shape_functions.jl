#==========================================================================================
                                Jacobian Computations
==========================================================================================#
# For Linear Curve Shape Functions
function jacobian(curveFunction::CurveLinear,coordinates,ξ)
    len = length(ξ)
    if len == 1
        return norm(coordinates[:,1] - coordinates[:,2]);
    else
        return (ones(1,length(ξ))) * norm(coordinates[:,1] - coordinates[:,2]);
    end
end
function jacobian(curveFunction::CurveLinear,coordinates)
    ξ = curveFunction.nodes;
    return (ones(1,length(ξ))) * norm(coordinates[:,1] - coordinates[:,2]);
end 
# For Quadratic Curve Shape Functions
function jacobian(curveFunction::CurveQuadratic,coordinates,ξ)
    len = length(ξ)
    dX = coordinates * curveQuadraticDerivative(ξ);
    if len == 1
        return sqrt(dot(dX,dX));
    else
        return sqrt.(sum(dX.^2, dims=1))
    end
end
function jacobian(curveFunction::CurveQuadratic,coordinates)
    return sqrt.(sum((coordinates * curveFunction.derivatives).^2, dims=1))
end
#==========================================================================================
                                Continuous Quadratic
                                   1 ---------- 2
==========================================================================================#
curveLinear(ξ) = [(1.0 .- ξ); ξ];
function basisFunction(curveFunction::CurveLinear,ξ)
    return curveLinear(ξ)
end
curveLinearDerivative(ξ) = [-1.0*ones(eltype(ξ),length(ξ))';
                             1.0*ones(eltype(ξ),length(ξ))']
function basisFunctionDerivative(curveFunction::CurveLinear,ξ)
    return curveLinearDerivative(ξ);
end
#==========================================================================================
                                Continuous Quadratic
                                   1 --- 2 --- 3
==========================================================================================#
curveQuadratic(ξ) = [0.5 * ξ .* (ξ .- 1.0); 1.0 .- ξ .^2; 0.5 * ξ .* (ξ .+ 1.0)];
function basisFunction(curveFunction::CurveQuadratic,ξ)
    return curveQuadratic(ξ)
end
curveQuadraticDerivative(ξ) = [ξ .- 0.5; -2.0*ξ; ξ .+ 0.5]
function basisFunctionDerivative(curveFunction::CurveQuadratic,ξ)
    return curveQuadraticDerivative(ξ)
end
#==========================================================================================
                        Discontinuous Continuous Constant
                                -------- 1 --------
==========================================================================================#
curveConstant(ξ) = ones(eltype(ξ),1,length(ξ))
function basisFunction(curveFunction::DiscontinuousCurveLinear,ξ)
    return curveConstant(ξ)
end
function basisFunctionDerivative(curveFunction::DiscontinuousCurveLinear,ξ)
    return zeros(eltype(ξ),1,length(ξ))
end
#==========================================================================================
                        Discontinuous Continuous Linear
                                -- 1 ---------- 2 --
==========================================================================================#
function basisFunction(curveFunction::DiscontinuousCurveLinear,ξ)
    return curveLinear(ξ/curveFunction.alpha)
end
function basisFunctionDerivative(curveFunction::DiscontinuousCurveLinear,ξ)
    return curveLinearDerivative(ξ/curveFunction.alpha)/curveFunction.alpha
end
#==========================================================================================
                        Discontinuous Continuous Quadratic
                                -- 1 --- 2 --- 3 --
==========================================================================================#
function basisFunction(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return curveQuadratic(ξ/curveFunction.alpha);
end
function basisFunctionDerivative(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return curveQuadraticDerivative(ξ/curveFunction.alpha)/curveFunction.alpha
end
#==========================================================================================
                            Defining relevant constructors 
==========================================================================================#
### For Linear Curve Shape Functions
function CurveLinear(ξ::AbstractArray)
    weights = similar(ξ)
    interpolation = curveLinear(ξ')
    derivatives   = curveLinearDerivative(ξ')
    return CurveLinear(weights,ξ,derivatives,interpolation)
end
function CurveLinear(nodes::AbstractArray,weights::AbstractArray)
    interpolation = curveLinear(nodes')
    derivatives   = curveLinearDerivative(nodes')
    return CurveLinear(weights,nodes,derivatives,interpolation)
end
function CurveLinear(n::Real)
    nodes, weights = curveLinearQuadpoints(n)
    interpolation  = curveLinear(nodes')
    derivatives    = curveLinearDerivative(nodes')
    return CurveLinear(weights,nodes,derivatives,interpolation)
end
### FOr Quadrat
function CurveQuadratic(ξ::AbstractArray)
    weights = similar(ξ)
    interpolation = curveQuadratic(ξ)
    derivatives   = curveQuadraticDerivative(ξ)
    return CurveQuadratic(weights,ξ,derivatives,interpolation)
end
function CurveQuadratic(nodes::AbstractArray,weights::AbstractArray)
    interpolation = curveQuadratic(nodes')
    derivatives   = curveQuadraticDerivative(nodes')
    return CurveQuadratic(weights,nodes,derivatives,interpolation)
end
function CurveQuadratic(n::Real)
    nodes, weights = curveQuadraticQuadpoints(n);
    interpolation  = curveQuadratic(nodes')
    derivatives    = curveQuadraticDerivative(nodes')
    return CurveQuadratic(weights,nodes,derivatives,interpolation);
end
