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


function mygemm!(C,A,B,v)
    Cre = reinterpret(reshape, Float64, C)
    Are = reinterpret(reshape, Float64, A)
    @turbo for n = eachindex(v)
      Cmn_re = zero(ComplexF64)
      Cmn_im = zero(ComplexF64)
    #   for k = indices((Are,Bre),(3,2))
      for k = indices(Are,3)
        Cmn_re += Are[1, 1, k] * B[k, n]
        Cmn_im += Are[2, 1, k] * B[k, n]
      end
      Cre[1,1,v[n]] += Cmn_re
      Cre[2,1,v[n]] += Cmn_im
    end
end

function mygemm2!(C, A, B)
    @inbounds @fastmath for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] += Cmn
    end
end
using BenchmarkTools
physics_nodes = [1; 3; 5]
m = 3
n = 10
C = zeros(ComplexF64,100,100);
Cc = @view C[1:1,physics_nodes]
A = zeros(ComplexF64,1,10)
B = zeros(ComplexF64,n,m)
Cc - A*B
# @benchmark mygemm!($parent(Cc),$A,$B,$Cc.indices[2])
@benchmark mygemm2!($Cc,$A,$B)


Cv = view(zeros(ComplexF64, 10, 10), 1:1, collect(1:6));
A = zeros(ComplexF64, 1, 25);
B = zeros(25, 6);

@benchmark mygemm!($parent(Cv),$A,$B,$Cv.indices[2])
@benchmark mygemm2!($Cv,$A,$B)


function mygemm3!(C, A, B)
    @inbounds @fastmath for n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A,2)
            Cmn += A[k] * B[k,n]
        end
        C[n] += Cmn
    end
end

C3 = @view C[1,physics_nodes]
A3 = zeros(ComplexF64, n);
B3 = zeros(ComplexF64, n, m);

@benchmark mygemm2!($Cc,$A,$B)
@benchmark mygemm3!($C3,$A3,$B3)
