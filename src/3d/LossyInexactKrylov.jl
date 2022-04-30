#===========================================================================================
                                Helper Functions
===========================================================================================#
l_mul_r0!(y,lx,tmp) = @inbounds for i=1:length(lx) y[i]  = lx[i]*tmp[i] end # y = lx*tmp
l_mul_r!(y,lx,tmp)  = @inbounds for i=1:length(lx) y[i] += lx[i]*tmp[i] end # y = ly*tmp + y
"""
    l_GH_r!(G,luG,H,lx,ly,lz,rx,ry,rz,tmp1,tmp2,y,x,setzero=true,uselu=false)

Computes the product
           ((lx*rx' + ly*ry' + lz*rz')∘(G^{-1}H))*x
Without performing the hadamard product by applying the following three times
           ((l*r')∘(G^{-1}H))*x = (diag(l)*G^{-1}*H*diag(r))*x
"""
function l_GH_r!(G,luG,H,lx,ly,lz,rx,ry,rz,tmp1,tmp2,y,x,setzero=true,uselu=false)
    # Compute ((rx*lx')∘(G^{-1}H))*x
    tmp1 .= rx .* x   # Vectorized version of: Diagonal(rx)*x
    mul!(tmp2,H,tmp1) # Inplace Mat-Vec product: tmp2 = H*tmp1
    if uselu
        ldiv!(luG,tmp2)      # Solve "G*x = tmp2" using lu-factorization of G
    else
        tmp1 .= tmp2
        gmres!(tmp2,G,tmp1)  # Solve "G*x = tmp2" gmres
    end
    setzero ? l_mul_r0!(y,lx,tmp2) : l_mul_r!(y,lx,tmp2)  # Add to sum
    # Compute ((ry*ly')∘(G^{-1}H))*x
    tmp1 .= ry.*x     # Vectorized version of: Diagonal(ry)*x
    mul!(tmp2,H,tmp1) # Inplace Mat-Vec product: tmp2 = H*tmp1
    if uselu
        ldiv!(luG,tmp2)     # Solve "G*x = tmp2" using lu-factorization of G
    else
        tmp1 .= tmp2
        gmres!(tmp2,G,tmp1) # Solve "G*x = tmp2" gmres
    end
    l_mul_r!(y,ly,tmp2) # Add to sum
    # Compute ((rz*lz')∘(G^{-1}H))*x
    tmp1 .= rz .* x   # Vectorized version of: Diagonal(rz)*x
    mul!(tmp2,H,tmp1) # Inplace Mat-Vec product: tmp2 = H*tmp1
    if uselu
        ldiv!(luG,tmp2)     # Solve "G*x = tmp2" using lu-factorization of G
    else
        tmp1 .= tmp2
        gmres!(tmp2,G,tmp1) # Solve "G*x = tmp2" gmres 
    end
    l_mul_r!(y,lz,tmp2) # Add to sum
end
"""
    l_D_r!(D,lx,ly,lz,rx,ry,rz,tmp1,tmp2,y,x,setzero=true) 

Computing the product
           ((lx*rx' + ly*ry' + lz*rz')∘D)*x
Without performing the hadamard product by applying the following three times
          ((l*r')∘D)*x = (diag(l)*D*diag(r))*x
"""
function l_D_r!(D,lx,ly,lz,rx,ry,rz,tmp1,tmp2,y,x,setzero=true) 
    # Compute ((lx*rx')∘D)*x
    tmp1 .= rx.*x
    mul!(tmp2,D,tmp1)
    setzero ? l_mul_r0!(y,lx,tmp2) : l_mul_r!(y,lx,tmp2) # Overwrite or add to tmp2
    # Compute ((ly*ry')∘D)*x
    tmp1 .= ry.*x
    mul!(tmp2,D,tmp1)
    l_mul_r!(y,ly,tmp2)
    # Compute ((lz*rz')∘D)*x
    tmp1 .= rz .* x
    mul!(tmp2,D,tmp1)
    l_mul_r!(y,lz,tmp2)
end

#===========================================================================================
                            Defining Inner Struct
===========================================================================================#
"""
    LossyOneVariableInner{T}

Memory efficient representation of the matrix
(n₁n₁ᵀ + n₂n₂ᵀ + n₃n₃ᵀ)∘(Gᵥ⁻Hᵥ) + (t₁n₁ᵀ + t₂n₂ᵀ + t₃n₃ᵀ)∘Dt₁ + (s₁n₁ᵀ + s₂n₂ᵀ + s₃n₃ᵀ)∘Dt₂,
which is part of the LossyOneVariableOuter{T}. 

Multiplication with the struct is implemented and can therefore be used memory efficiently 
solve linear systems using an iterative solver. 
"""
struct LossyOneVariableInner{T} <: LinearMaps.LinearMap{T}
    N::Int64                # The matrix is of size N × N
    # Viscosity Matrices
    Hv::AbstractArray{T}    # Viscous BEM H
    Gv::AbstractArray{T}    # Viscous BEM G
    luGv::Factorization{T}  # Factorization of Viscous BEM H-matrix
    # Derivative matrices
    Dt1::AbstractArray{T}   # Tangential derivative in tangent direction 1
    Dt2::AbstractArray{T}   # Tangential derivative in tangent direction 2
    # Normals components
    nx::AbstractArray       # x-coordinates of normal direction
    ny::AbstractArray       # y-coordinates of normal direction      
    nz::AbstractArray       # z-coordinates of normal direction      
    # Tangent components
    tx::AbstractArray       # x-coordinates of 1st tangent direction
    ty::AbstractArray       # y-coordinates of 1st tangent direction
    tz::AbstractArray       # z-coordinates of 1st tangent direction
    # Tagent components
    sx::AbstractArray       # x-coordinates of 2nd tangent direction
    sy::AbstractArray       # y-coordinates of 2nd tangent direction
    sz::AbstractArray       # z-coordinates of 2nd tangent direction
    # Temporary vectors - To avoid allocation when multiplying
    tmp1::AbstractArray{T} 
    tmp2::AbstractArray{T} 
    # How to solve 
    lu_on::Bool
end

# Define size operator - Used when calling LinearMaps.check_dim_mul
Base.size(A::LossyOneVariableInner) = (A.N, A.N)

# Standard Multiplication
function LinearAlgebra.mul!(y::AbstractVecOrMat{T}, 
                              A::LossyOneVariableInner{T}, 
                              x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    l_GH_r!(A.Gv,A.luGv, A.Hv, A.nx, A.ny, A.nz, A.nx, A.ny, A.nz,   A.tmp1, A.tmp2, y, x, true, A.lu_on)
    l_D_r!(A.Dt1,        A.tx, A.ty, A.tz, A.nx, A.ny, A.nz, A.tmp1, A.tmp2,         y, x, false)
    l_D_r!(A.Dt2,        A.sx, A.sy, A.sz, A.nx, A.ny, A.nz, A.tmp1, A.tmp2,         y, x, false)
end
#===========================================================================================
                            Defining Outer Struct
===========================================================================================#
"""
    LossyOneVariableOuter{T}

Struct representing the system matrix (A) in the 1-variable system: A*pₐ = b. 
"""
struct LossyOneVariableOuter{T} <: LinearMaps.LinearMap{T}
    N::Int64                            # Total system size (10n
    # Acoustic Matrices
    Ha::AbstractArray{T}                # Acoustic BEM H
    Ga::AbstractArray{T}                # Acoustic BEM G
    # Thermal Matrices
    Hh::AbstractArray{T}                # Thermal BEM H
    Gh::AbstractArray{T}                # Thermal BEM G
    luGh::Factorization{T}              # Factorization of Thermal BEM G
    # Viscosity Matrices
    Hv::AbstractArray{T}                # Viscous BEM H
    Gv::AbstractArray{T}                # Viscous BEM G
    luGv::Factorization{T}              # Factorization of Viscous BEM G
    # The inner matrix 
    inner::LossyOneVariableInner{T}
    # Derivative matrices
    Dt1::AbstractArray{T}               # Tangential derivative in direction 1
    Dt2::AbstractArray{T}               # Tangential derivative in direction 2
    # Normals components
    nx::AbstractArray                   # x-coordinates of normal direction      
    ny::AbstractArray                   # y-coordinates of normal direction
    nz::AbstractArray                   # z-coordinates of normal direction
    # Tangent components
    tx::AbstractArray                   # x-coordinates of 1st tangent direction
    ty::AbstractArray                   # y-coordinates of 1st tangent direction
    tz::AbstractArray                   # z-coordinates of 1st tangent direction
    # Tangent components
    sx::AbstractArray                   # x-coordinates of 2nd tangent direction
    sy::AbstractArray                   # y-coordinates of 2nd tangent direction
    sz::AbstractArray                   # z-coordinates of 2nd tangent direction
    ### Constants
    ϕₐ::T
    ϕₕ::T
    τₐ::T
    τₕ::T
    # Temporary vectors - To avoid allocation when multiplying
    tmp1::AbstractArray{T}
    tmp2::AbstractArray{T}
    res::AbstractArray{T}
    x1::AbstractArray{T}
    x2::AbstractArray{T}
    # Solving 
    lu_on::Bool
    # Saving frequency
    freq
end

# Define size operator - Used when calling LinearMaps.check_dim_mul
Base.size(A::LossyOneVariableOuter) = (A.N, A.N)

# Standard Multiplication
function LinearAlgebra.mul!(y::AbstractVecOrMat{T}, 
                              A::LossyOneVariableOuter{T}, 
                              x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    # First mulitplication
    mul!(A.x1,A.Dt1,x)
    l_GH_r!(A.Gv,A.luGv, A.Hv, A.nx, A.ny, A.nz, A.tx, A.ty, A.tz,   A.tmp1, A.tmp2, A.res, A.x1, true, A.lu_on)  
    l_D_r!(A.Dt1,        A.tx, A.ty, A.tz, A.tx, A.ty, A.tz, A.tmp1, A.tmp2,         A.res, A.x1, false)
    l_D_r!(A.Dt2,        A.sx, A.sy, A.sz, A.tx, A.ty, A.tz, A.tmp1, A.tmp2,         A.res, A.x1, false)
    # Second multiplication
    mul!(A.x2,A.Dt2,x)
    l_GH_r!(A.Gv, A.luGv, A.Hv, A.nx, A.ny, A.nz, A.sx, A.sy, A.sz,   A.tmp1, A.tmp2, A.res, A.x2, false, A.lu_on)  
    l_D_r!(A.Dt1,         A.tx, A.ty, A.tz, A.sx, A.sy, A.sz, A.tmp1, A.tmp2,         A.res, A.x2, false)
    l_D_r!(A.Dt2,         A.sx, A.sy, A.sz, A.sx, A.sy, A.sz, A.tmp1, A.tmp2,         A.res, A.x2, false)
    # Solving the inner system
    gmres!(A.x1,A.inner,A.res)
    # Reuse temporary vectors to avoid allocations
    mul!(A.tmp2,A.Hh,x)
    if A.lu_on
        ldiv!(A.luGh,A.tmp2)
    else
        gmres!(A.tmp1,A.Gh,A.tmp2)
    end
    A.x2 .= -A.ϕₕ*A.τₐ/A.τₕ*A.tmp1 + (A.ϕₐ - A.ϕₕ*A.τₐ/A.τₕ)*A.x1
    # Adding things to y
    mul!(y,A.Ha,x,A.ϕₐ,false)
    mul!(y,A.Ga,A.x2,true,true)
end
#===========================================================================================
                            Computing Right Hand Side
===========================================================================================#
l_GH_r(G,H,lx,ly,lz,rx,ry,rz,x) = lx.*(G\(H*(rx.*x))) + ly.*(G\(H*(ry.*x))) + lz.*(G\(H*(rz.*x)))
l_D_r(D,lx,ly,lz,rx,ry,rz,x)    = lx.*(D*   (rx.*x))  + ly.*(   D*(ry.*x))  + lz.*(   D*(rz.*x))
"""
    compute_lossy_rhs(A::LossyOneVariableOuter,vbn,vt1,vt2)

Computes the right-hand-side of the final 1-variable lossy system.
Requires normal and tangential velocities. 
"""
function compute_lossy_rhs(A::LossyOneVariableOuter,vbn,vt1,vt2)
    # First multiplication
    l_GH_r!(A.Gv,A.luGv, A.Hv, A.nx, A.ny, A.nz, A.tx, A.ty, A.tz,   A.tmp1, A.tmp2, A.x1, vt1, true, A.lu_on)  
    l_D_r!(A.Dt1,        A.tx, A.ty, A.tz, A.tx, A.ty, A.tz, A.tmp1, A.tmp2,         A.x1, vt1, false)
    l_D_r!(A.Dt2,        A.sx, A.sy, A.sz, A.tx, A.ty, A.tz, A.tmp1, A.tmp2,         A.x1, vt1, false)
    # Second multiplication
    l_GH_r!(A.Gv, A.luGv, A.Hv, A.nx, A.ny, A.nz, A.sx, A.sy, A.sz,   A.tmp1, A.tmp2, A.x2, vt2, true, A.lu_on)  
    l_D_r!(A.Dt1,         A.tx, A.ty, A.tz, A.sx, A.sy, A.sz, A.tmp1, A.tmp2,         A.x2, vt2, false)
    l_D_r!(A.Dt2,         A.sx, A.sy, A.sz, A.sx, A.sy, A.sz, A.tmp1, A.tmp2,         A.x2, vt2, false)
    return A.Ga*(vbn + gmres(A.inner,A.x1 + A.x2;verbose=true))
end
