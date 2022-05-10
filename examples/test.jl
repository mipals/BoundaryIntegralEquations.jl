using SparseArrays

function fill_matrix!(M)
    for i = 1:size(M,1), j = 1:size(M,2)
        M[i,j] = 1.0
    end
end
function fill_matrix_partially!(M,p=0.01)
    for i = 1:size(M,1), j = 1:Int(ceil(size(M,2)*p))
        M[i,j] = 1.0
    end
end

A = zeros(3000,3000)
@time fill_matrix!(A)
@time fill_matrix_partially!(A)

B = spzeros(size(A))
@time fill_matrix_partially!(B,0.05)
# We definately can not just start with a sparse matrix!
# We need to compute Coordinates I,J,V and build sparse(I,J,V)...


using StaticArrays
mutable struct MutQuadElement{T<:AbstractFloat}
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
struct QuadElement{T<:AbstractFloat}
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end

n = 100
weights = rand(n)
gauss_u = rand(n)
gauss_v = rand(n)
derivatives_u = rand(4,n)
derivatives_v = rand(4,n)
interpolation = rand(4,n)
MQ = MutQuadElement(weights,gauss_u,gauss_v,derivatives_u,derivatives_v,interpolation)
Q = QuadElement(weights,gauss_u,gauss_v,derivatives_u,derivatives_v,interpolation)

X = rand(3,4)
Y = X*interpolation

using LinearAlgebra
@time mul!(Y,X,MQ.interpolation);
@time mul!(Y,X,Q.interpolation);





using IntegralEquations
import IntegralEquations: basisFunctionSecondOrderDerivative
TL = TriangularLinear(3,3)
TQ = TriangularQuadratic(3,3)
DTC = DiscontinuousTriangularConstant(TL)
DTL = DiscontinuousTriangularLinear(TL,0.2)
DTQ = DiscontinuousTriangularQuadratic(TL,0.2)

QL = QuadrilateralLinear4(3,3)
QQ = QuadrilateralQuadraticLagrange(3,3)
DQC = DiscontinuousQuadrilateralConstant(3,3)
DQL = DiscontinuousQuadrilateralLinear4(3,3)
DQQ = DiscontinuousQuadrilateralQuadraticLagrange(3,3)

x = rand(5)
y = rand(5)

d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(TL,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(TQ,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DTC,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DTL,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DTQ,x,y)

d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(QL,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(QQ,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DQC,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DQL,x,y)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DQQ,x,y)


CL = ContinuousCurveLinear(3)
CQ = ContinuousCurveQuadratic(3)
DCC = DiscontinuousCurveConstant(3)
DCL = DiscontinuousCurveLinear(3)
DCQ = DiscontinuousCurveQuadratic(3)

d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(CL,x)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(CQ,x)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DCC,x)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DCL,x)
d2x,d2y,d2xy=basisFunctionSecondOrderDerivative(DCQ,x)
