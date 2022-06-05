using LinearAlgebra, BenchmarkTools, StaticArrays

Ar = zeros(10,10)
Br = zeros(10,10)
Ac = zeros(ComplexF64,10,10)
Bc = zeros(ComplexF64,10,10)
Cc = zeros(ComplexF64,10,10)

@benchmark mul!(Cc,Bc,Ar,true,true) # Ok.
@benchmark mul!(Cc,Br,Ac,true,true) # Slow.
@benchmark Cc .+= Bc*Ar
@benchmark Cc .+= Br*Ac
@benchmark Cc .+= Bc*Ac
@benchmark mul!(Cc,Bc,Ac,true,true)


C = zeros(ComplexF64,100,100)
Cc = @view C[1:10,1:10]

@benchmark matmul!(Cc,Bc,Ar,Octavian.StaticInt(1),Octavian.StaticInt(1))
@benchmark matmul!(Cc,Br,Ac,Octavian.StaticInt(1),Octavian.StaticInt(1))
@benchmark matmul!(Cc,Bc,Ac,Octavian.StaticInt(1),Octavian.StaticInt(1))

using LoopVectorization
using BenchmarkTools

function mygemm2!(C, A, B)
    @inbounds @fastmath for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] += Cmn
    end
end
function mygemm3!(C, A, B)
    @inbounds @fastmath for n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(B,1)
            Cmn += A[k] * B[k,n]
        end
        C[n] += Cmn
    end
end
using BenchmarkTools
physics_nodes = [1; 3; 5]
m = 3
n = 50
C = zeros(ComplexF64,100,100);
Cc = @view C[1:1,physics_nodes]
A = zeros(ComplexF64,1,n)
B = zeros(n,m)
Cc - A*B

C3 = @view C[1,physics_nodes]
A3 = zeros(ComplexF64, n)
B3 = zeros(n, m)

@benchmark mygemm2!($Cc,$A,$B)
@benchmark mygemm3!($C3,$A3,$B3)


###
using IntegralEquations

C = zeros(ComplexF64,100,100);
physics_nodes = [1; 3; 5]
m = 3
n = 50
C3 = @view C[1,physics_nodes]
A3 = zeros(ComplexF64, n)
B3 = zeros(n, m)

normals = ones(3,n)
integrand = A3
r = ones(n)
interpolation = ones(3,n)
sources = rand(3)
k = 1.0
@time IntegralEquations.freens3d!(integrand,r,interpolation,sources,normals,k)

C = zeros(ComplexF64,10)
c = @view C[1:1]
weights = ones(size(integrand))
@time IntegralEquations.add_to_c!(c,integrand,weights)


function add_to(c,j,w)
    @inbounds for i = 1:length(j)
        c += w[i] * j[i]
    end
end
c = 0.0*im
j = ones(ComplexF64,10)
w = ones(10)
add_to(c,j,w)
c
