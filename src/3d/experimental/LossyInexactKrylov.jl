#===========================================================================================
                                Helper Functions
===========================================================================================#
l_mul_r0!(y,lx,tmp) = @inbounds for i=1:length(lx) y[i]  = lx[i]*tmp[i] end # y = lx*tmp
l_mul_r!(y,lx,tmp)  = @inbounds for i=1:length(lx) y[i] += lx[i]*tmp[i] end # y = ly*tmp + y
"""
    l_GH_r!(G,luG,H,lx,ly,lz,rx,ry,rz,tmp1,tmp2,y,x,setzero=true,uselu=false)

Computes the product
```math
    ((\\mathbf{l}_x\\mathbf{r}_x^\\top + \\mathbf{l}_y\\mathbf{r}_y^\\top + \\mathbf{l}_z\\mathbf{r}_z^\\top)\\circ (\\mathbf{G}^{-1}\\mathbf{H}))\\mathbf{x}
```
Without performing the Hadamard product by applying the following three times
```math
    ((\\mathbf{l}\\mathbf{r}^\\top)\\circ (\\mathbf{G}^{-1}\\mathbf{H}))\\mathbf{x} =
    (\\text{diag}(\\mathbf{l})\\mathbf{G}^{-1}\\mathbf{H}\\text{diag}(\\mathbf{r}))\\mathbf{x}
```
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
```math
    ((\\mathbf{l}_x\\mathbf{r}_x^\\top + \\mathbf{l}_y\\mathbf{r}_y^\\top + \\mathbf{l}_z\\mathbf{r}_z^\\top)\\circ \\mathbf{D})\\mathbf{x}
```
Without performing the Hadamard product by applying the following three times
```math
    ((\\mathbf{l}\\mathbf{r}^\\top)\\circ \\mathbf{D})\\mathbf{x} = (\\text{diag}(\\mathbf{l})\\mathbf{D}\\text{diag}(\\mathbf{r}))\\mathbf{x}
```
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
```math
(\\mathbf{n}_1\\mathbf{n}_1^\\top +
 \\mathbf{n}_2\\mathbf{n}_2^\\top +
 \\mathbf{n}_3\\mathbf{n}_3^\\top)\\circ(\\mathbf{G}_v^{-1}\\mathbf{H}_v) +
(\\mathbf{t}_1\\mathbf{n}_1^\\top +
 \\mathbf{t}_2\\mathbf{n}_2^\\top +
 \\mathbf{t}_3\\mathbf{n}_3^\\top)\\circ\\mathbf{D}_{t_1} +
(\\mathbf{s}_1\\mathbf{n}_1^\\top +
 \\mathbf{s}_2\\mathbf{n}_2^\\top +
 \\mathbf{s}_3\\mathbf{n}_3^\\top )\\circ\\mathbf{D}_{t_2},
```
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

Struct representing the system matrix in the 1-variable system:
```math
    \\mathbf{A}\\mathbf{p}_a = \\mathbf{b}.
```
"""
struct LossyOneVariableOuter{T} <: LinearMaps.LinearMap{T}
    N::Int64                            # Total system size (10n
    # Acoustic Matrices
    # Ha::AbstractArray{T}                # Acoustic BEM H
    # Ga::AbstractArray{T}                # Acoustic BEM G
    Ha                                  # Acoustic BEM H
    Ga                                  # Acoustic BEM G
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

#===========================================================================================
                                Constructors
===========================================================================================#
function LossyOneVariableOuter(mesh::Mesh3d,freq;S=1,depth=1,exterior=true,
                            nearfield=true,n=3,thres=1e-6,offset=0.2)
    # Extracting local (n,t,s) coordinate systems
    nx = mesh.normals[1,:]
    ny = mesh.normals[2,:]
    nz = mesh.normals[3,:]
    tx = mesh.tangents[1,:]
    ty = mesh.tangents[2,:]
    tz = mesh.tangents[3,:]
    sx = mesh.sangents[1,:]
    sy = mesh.sangents[2,:]
    sz = mesh.sangents[3,:]
    N = length(nx)

    # Allocating memory used for saving intermediate values when multiplying with the struct
    inner_tmp1 = zeros(ComplexF64,N)
    inner_tmp2 = similar(inner_tmp1)
    outer_tmp1 = similar(inner_tmp1)
    outer_tmp2 = similar(inner_tmp1)
    outer_res  = similar(inner_tmp1)
    outer_x1   = similar(inner_tmp1)
    outer_x2   = similar(inner_tmp1)
    # Computing relevant BEM-matrices
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=S)

    ### Extracting our sources
    sources = mesh.sources

    ### Extracting the number of nodes (and hence also the total matrix size)

    ### Assembling the 3 BEM systems
    # Thermal matrices
    @info "Thermal Matrices:"
    Fₕ,Bₕ = assemble_parallel!(mesh,kₕ,sources;sparse=true,depth=depth);
    Aₕ = (exterior ?  Fₕ + I/2 : -Fₕ + I/2)
    # Viscous matrices
    @info "Viscous matrices:"
    Fᵥ,Bᵥ  = assemble_parallel!(mesh,kᵥ,sources;sparse=true,depth=depth);
    Aᵥ = (exterior ?  Fᵥ + I/2 : -Fᵥ + I/2)

    ### Computing tangential derivatives
    Dt₁, Dt₂ = interpolation_function_derivatives(mesh;global_derivatives=false)
    Dt₁ = convert.(eltype(Bᵥ),Dt₁)
    Dt₂ = convert.(eltype(Bᵥ),Dt₂)

    luGv = lu(rand(ComplexF64,2,2)) # Junk factorization
    luGh = lu(rand(ComplexF64,2,2)) # Junk factorization

    inner = LossyOneVariableInner(N,Aᵥ,Bᵥ,luGv,Dt₁,Dt₂,
                                nx,ny,nz,tx,ty,tz,sx,sy,sz,
                                inner_tmp1, inner_tmp2,false)

    # Computing relevant constants
    @info "Acoustic Matrices:"
    Ga = FMMGOperator(mesh,kₐ;n_gauss=n,tol=thres,offset=offset,nearfield=nearfield)
    Fa = FMMFOperator(mesh,kₐ;n_gauss=n,tol=thres,offset=offset,nearfield=nearfield)
    Ha = (exterior ? Fa + I/2 : -Fa + I/2)
    return LossyOneVariableOuter(N,Ha,Ga,Aₕ,Bₕ,luGh,Aᵥ,Bᵥ,
                            luGv,inner,Dt₁,Dt₂,
                            nx,ny,nz,tx,ty,tz,sx,sy,sz,ϕₐ,ϕₕ,τₐ,τₕ,
                            outer_tmp1,outer_tmp2,outer_res,outer_x1,outer_x2,false,freq)
end

function LossyOneVariableOuter(mesh::Mesh3d,BB::LossyBlockMatrix,freq;depth=1,
                        lu_on=false,fmm_on=false,nearfield=true,n=3,thres=1e-6,offset=0.2)
    τₐ = BB.τₐ
    τₕ = BB.τₕ
    ϕₐ = BB.ϕₐ
    ϕₕ = BB.ϕₕ
    # Extracting local (n,t,s) coordinate systems
    nx = mesh.normals[1,:]
    ny = mesh.normals[2,:]
    nz = mesh.normals[3,:]
    tx = mesh.tangents[1,:]
    ty = mesh.tangents[2,:]
    tz = mesh.tangents[3,:]
    sx = mesh.sangents[1,:]
    sy = mesh.sangents[2,:]
    sz = mesh.sangents[3,:]
    N = length(nx)

    # Allocating memory used for saving intermediate values when multiplying with the struct
    inner_tmp1 = zeros(ComplexF64,N)
    inner_tmp2 = similar(inner_tmp1)
    outer_tmp1 = similar(inner_tmp1)
    outer_tmp2 = similar(inner_tmp1)
    outer_res  = similar(inner_tmp1)
    outer_x1   = similar(inner_tmp1)
    outer_x2   = similar(inner_tmp1)

    # Given that Bᵥ and Bₕ in some limited tets has been found to be VERY well conditioned
    # the lu-factorization is disabled by default. Instead gmres is used.
    if lu_on
        luGv = lu(BB.Bᵥ)
        luGh = lu(BB.Bₕ)
    else
        luGv = lu(rand(ComplexF64,2,2)) # Junk factorization
        luGh = lu(rand(ComplexF64,2,2)) # Junk factorization
    end
    inner = LossyOneVariableInner(N,BB.Aᵥ,BB.Bᵥ,luGv,BB.Dt₁,BB.Dt₂,
                                    nx,ny,nz,tx,ty,tz,sx,sy,sz,
                                    inner_tmp1, inner_tmp2,lu_on)
    if fmm_on
        _,_,_,ka,_,_,_,_,_,_,_,_ = visco_thermal_constants(;freq=freq,S=1)
        Ga = FMMGOperator(mesh,ka;n_gauss=n,tol=thres,offset=offset,nearfield=nearfield,depth=depth)
        Ha = FMMFOperator(mesh,ka;n_gauss=n,tol=thres,offset=offset,nearfield=nearfield,depth=depth) + I/2
        outer = LossyOneVariableOuter(N,Ha,Ga,BB.Aₕ,BB.Bₕ,luGh,BB.Aᵥ,BB.Bᵥ,
                                luGv,inner,BB.Dt₁,BB.Dt₂,
                                nx,ny,nz,tx,ty,tz,sx,sy,sz,ϕₐ,ϕₕ,τₐ,τₕ,
                                outer_tmp1,outer_tmp2,outer_res,outer_x1,outer_x2,lu_on,freq)
    else
        outer = LossyOneVariableOuter(N,BB.Aₐ,BB.Bₐ,BB.Aₕ,BB.Bₕ,luGh,BB.Aᵥ,BB.Bᵥ,
                                luGv,inner,BB.Dt₁,BB.Dt₂,
                                nx,ny,nz,tx,ty,tz,sx,sy,sz,ϕₐ,ϕₕ,τₐ,τₕ,
                                outer_tmp1,outer_tmp2,outer_res,outer_x1,outer_x2,lu_on,freq)
    end
    return outer
end
