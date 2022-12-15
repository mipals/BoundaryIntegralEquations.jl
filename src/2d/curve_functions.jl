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
curveLinear(ξ) = [-0.5*ξ .+ 0.5; 0.5*ξ .+ 0.5];
function basisFunction(curveFunction::ContinuousCurveLinear,ξ)
    return curveLinear(ξ)
end
curveLinearDerivative(ξ) = [-0.5*ones(eltype(ξ),length(ξ))';
                             0.5*ones(eltype(ξ),length(ξ))']
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
discontinuousCurveLinear(ξ,beta) = curveLinear(ξ/beta)
function basisFunction(curveFunction::DiscontinuousCurveLinear,ξ)
    return discontinuousCurveLinear(ξ,curveFunction.beta)
end
# Just a simple chain rule
discontinuousCurveLinearDerivative(ξ,beta) = curveLinearDerivative(ξ/beta)/beta
function basisFunctionDerivative(curveFunction::DiscontinuousCurveLinear,ξ)
    return discontinuousCurveLinearDerivative(ξ,curveFunction.beta)
end
#==========================================================================================
                        Discontinuous Continuous Quadratic
                                -- 1 --- 2 --- 3 --
==========================================================================================#
discontinuousCurveQuadratic(ξ,beta) = curveQuadratic(ξ/beta);
function basisFunction(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return discontinuousCurveQuadratic(ξ,curveFunction.beta)
end
# Just a simple chain rule
discontinuousCurveQuadraticDerivative(ξ,beta) = curveQuadraticDerivative(ξ/beta)/beta
function basisFunctionDerivative(curveFunction::DiscontinuousCurveQuadratic,ξ)
    return discontinuousCurveQuadraticDerivative(ξ,curveFunction.beta)
end
#==========================================================================================
                            Defining relevant constructors
==========================================================================================#
### For Linear Curve Shape Functions
function ContinuousCurveLinear(n::Int)
    nodes, weights = curveQuadraticQuadpoints(n)
    interpolation  = curveLinear(nodes')
    derivatives    = curveLinearDerivative(nodes')
    return ContinuousCurveLinear(weights,nodes,derivatives,interpolation)
end
### For Quadratic
function ContinuousCurveQuadratic(n::Int)
    nodes, weights = curveQuadraticQuadpoints(n);
    interpolation  = curveQuadratic(nodes')
    derivatives    = curveQuadraticDerivative(nodes')
    return ContinuousCurveQuadratic(weights,nodes,derivatives,interpolation)
end
### For Discontinuous Constant
function DiscontinuousCurveConstant(n::Int)
    nodes, weights = curveQuadraticQuadpoints(n)
    interpolation  = discontinuousCurveConstant(nodes')
    derivatives    = discontinuousCurveConstantDerivative(nodes')
    return DiscontinuousCurveConstant(weights,nodes,derivatives,interpolation)
end
### For Discontinuous Linear
function DiscontinuousCurveLinear(n::Int,beta=0.5773502691896258)
    nodes, weights = curveQuadraticQuadpoints(n)
    interpolation  = discontinuousCurveLinear(nodes',beta)
    derivatives    = discontinuousCurveLinearDerivative(nodes',beta)
    return DiscontinuousCurveLinear(weights,nodes,derivatives,interpolation,beta)
end
### For Discontinuous Quadratic
function DiscontinuousCurveQuadratic(n::Int,beta=0.7745966692414834)
    nodes, weights = curveQuadraticQuadpoints(n)
    interpolation  = discontinuousCurveQuadratic(nodes',beta)
    derivatives    = discontinuousCurveQuadraticDerivative(nodes',beta)
    return DiscontinuousCurveQuadratic(weights,nodes,derivatives,interpolation,beta)
end
