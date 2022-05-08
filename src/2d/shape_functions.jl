#==========================================================================================
                                Jacobian Computations
==========================================================================================#
# For Linear Curve Shape Functions
function jacobian(curveFunction::ContinuousCurveLinear,coordinates,ξ)
    len = length(ξ)
    if len == 1
        return norm(coordinates[:,1] - coordinates[:,2]);
    else
        return (ones(1,length(ξ))) * norm(coordinates[:,1] - coordinates[:,2]);
    end
end
function jacobian(curveFunction::ContinuousCurveLinear,coordinates)
    ξ = curveFunction.nodes;
    return (ones(1,length(ξ))) * norm(coordinates[:,1] - coordinates[:,2]);
end 
# For Quadratic Curve Shape Functions
function jacobian(curveFunction::ContinuousCurveQuadratic,coordinates,ξ)
    len = length(ξ)
    dX = coordinates * curveQuadraticDerivative(ξ);
    if len == 1
        return sqrt(dot(dX,dX));
    else
        return sqrt.(sum(dX.^2, dims=1))
    end
end
function jacobian(curveFunction::ContinuousCurveQuadratic,coordinates)
    return sqrt.(sum((coordinates * curveFunction.derivatives).^2, dims=1))
end
#==========================================================================================
                                Continuous Quadratic
                                   1 ---------- 2
==========================================================================================#
curveLinear(ξ) = [(1.0 .- ξ); ξ];
function basisFunction(curveFunction::ContinuousCurveLinear,ξ)
    return curveLinear(ξ)
end
curveLinearDerivative(ξ) = [-1.0*ones(eltype(ξ),length(ξ))';
                             1.0*ones(eltype(ξ),length(ξ))']
function basisFunctionDerivative(curveFunction::ContinuousCurveLinear,ξ)
    return curveLinearDerivative(ξ);
end
#==========================================================================================
                                Continuous Quadratic
                                   1 --- 2 --- 3
==========================================================================================#
curveQuadratic(ξ) = [0.5 * ξ .* (ξ .- 1.0); 1.0 .- ξ .^2; 0.5 * ξ .* (ξ .+ 1.0)];
function basisFunction(curveFunction::ContinuousCurveQuadratic,ξ)
    return curveQuadratic(ξ)
end
curveQuadraticDerivative(ξ) = [ξ .- 0.5; -2.0*ξ; ξ .+ 0.5]
function basisFunctionDerivative(curveFunction::ContinuousCurveQuadratic,ξ)
    return curveQuadraticDerivative(ξ)
end
#==========================================================================================
                        Discontinuous Continuous Constant
                                -------- 1 --------
==========================================================================================#
discontinuousCurveConstant(ξ) = ones(eltype(ξ),1,length(ξ))
function basisFunction(curveFunction::DiscontinuousCurveConstant,ξ)
    return discontinuousCurveConstant(ξ)
end
discontinuousCurveConstantDerivative(ξ) = zeros(eltype(ξ),1,length(ξ))
function basisFunctionDerivative(curveFunction::DiscontinuousCurveConstant,ξ)
    return discontinuousCurveConstantDerivative(ξ)
end
#==========================================================================================
                        Discontinuous Continuous Linear
                                -- 1 ---------- 2 --
==========================================================================================#
discontinuousCurveLinear(ξ) = curveLinear(ξ/curveFunction.alpha)
function basisFunction(curveFunction::DiscontinuousCurveLinear,ξ)
    return discontinuousCurveLinear(ξ)
end
# Just a simple chain rule
discontinuousCurveLinearDerivative(ξ) = curveLinearDerivative(ξ/curveFunction.alpha)/
                                                                curveFunction.alpha
function basisFunctionDerivative(curveFunction::DiscontinuousCurveLinear,ξ)
    return discontinuousCurveLinearDerivative(ξ)
end
#==========================================================================================
                        Discontinuous Continuous Quadratic
                                -- 1 --- 2 --- 3 --
==========================================================================================#
discontinuousCurveQuadratic(ξ) = curveQuadratic(ξ/curveFunction.alpha);
function basisFunction(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return discontinuousCurveQuadratic(ξ)
end
# Just a simple chain rule
discontinuousCurveQuadraticDerivative(ξ) = curveQuadraticDerivative(ξ/curveFunction.alpha)/
                                                                      curveFunction.alpha
function basisFunctionDerivative(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return discontinuousCurveQuadraticDerivative(ξ)
end
#==========================================================================================
                            Defining relevant constructors 
==========================================================================================#
### For Linear Curve Shape Functions
function ContinuousCurveLinear(ξ::AbstractArray)
    weights = similar(ξ)
    interpolation = curveLinear(ξ')
    derivatives   = curveLinearDerivative(ξ')
    return ContinuousCurveLinear(weights,ξ,derivatives,interpolation)
end
function ContinuousCurveLinear(nodes::AbstractArray,weights::AbstractArray)
    interpolation = curveLinear(nodes')
    derivatives   = curveLinearDerivative(nodes')
    return ContinuousCurveLinear(weights,nodes,derivatives,interpolation)
end
function ContinuousCurveLinear(n::Real)
    nodes, weights = curveLinearQuadpoints(n)
    interpolation  = curveLinear(nodes')
    derivatives    = curveLinearDerivative(nodes')
    return ContinuousCurveLinear(weights,nodes,derivatives,interpolation)
end
### For Quadratic
function ContinuousCurveQuadratic(ξ::AbstractArray)
    weights = similar(ξ)
    interpolation = curveQuadratic(ξ)
    derivatives   = curveQuadraticDerivative(ξ)
    return ContinuousCurveQuadratic(weights,ξ,derivatives,interpolation)
end
function ContinuousCurveQuadratic(nodes::AbstractArray,weights::AbstractArray)
    interpolation = curveQuadratic(nodes')
    derivatives   = curveQuadraticDerivative(nodes')
    return ContinuousCurveQuadratic(weights,nodes,derivatives,interpolation)
end
function ContinuousCurveQuadratic(n::Real)
    nodes, weights = curveQuadraticQuadpoints(n);
    interpolation  = curveQuadratic(nodes')
    derivatives    = curveQuadraticDerivative(nodes')
    return ContinuousCurveQuadratic(weights,nodes,derivatives,interpolation)
end
### For Discontinuous Constant
function DiscontinuousCurveConstant(ξ::AbstractArray)
    weights       = similar(ξ)
    interpolation = discontinuousCurveConstant(ξ)
    derivatives   = discontinuousCurveConstantDerivative(ξ)
    return DiscontinuousCurveConstant(weights,ξ,derivatives,interpolation)
end
function DiscontinuousCurveConstant(nodes::AbstractArray,weights::AbstractArray)
    interpolation = discontinuousCurveConstant(nodes')
    derivatives   = discontinuousCurveConstantDerivative(nodes')
    return DiscontinuousCurveConstant(weights,nodes,derivatives,interpolation)
end
function DiscontinuousCurveConstant(n::Real)
    nodes, weights = curveQuadraticQuadpoints(n)
    interpolation  = discontinuousCurveConstant(nodes')
    derivatives    = discontinuousCurveConstantDerivative(nodes')
    return DiscontinuousCurveConstant(weights,nodes,derivatives,interpolation)
end
