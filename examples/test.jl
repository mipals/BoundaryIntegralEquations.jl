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



