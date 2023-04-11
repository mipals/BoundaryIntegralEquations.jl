#==========================================================================================
                            Defining LossyGlobalInner
==========================================================================================#
"""
    LossyGlobalInner

A `LinearMap` corresponding to the inner system of the reduced lossy system.
"""
struct LossyGlobalInner{T} <: LinearMaps.LinearMap{T}
    n::Int64                # Number of nodes
    # Acoustic Matrices
    Hv::AbstractArray{T}    # Viscous BEM A
    Gv                      # Viscous BEM B
    # Normal transformation
    Nd::AbstractArray       # Collection of normal tr
    # Gradients
    Dr::AbstractArray       # [Dx  Dy  Dz]
end

Base.size(A::LossyGlobalInner) = (A.n, A.n)

# Defining multiplication with `LossyGlobalInner`
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::LossyGlobalInner{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions of input
    LinearMaps.check_dim_mul(y, A, x)
    # Extracting relevant values
    if typeof(A.Gv) <: Factorization
        y .= A.Dr*(A.Nd*x) - A.Nd'*(A.Gv\(A.Hv*(A.Nd*x)))
    else
        y .= A.Dr*(A.Nd*x) - A.Nd'*(gmres(A.Gv,A.Hv*(A.Nd*x)))
    end
    return y
end

function Base.Matrix(A::LossyGlobalInner)
    return A.Dr*A.Nd - A.Nd'*(Matrix(A.Gv)\(A.Hv*A.Nd))
end
#==========================================================================================
                            Defining LossyGlobalOuter
==========================================================================================#
struct LossyGlobalOuter{T} <: LinearMaps.LinearMap{T}
    n::Int64                    # Number of nodes
    Ha                          # Acoustical BEM H
    Ga                          # Acoustical BEM G
    Hh::AbstractArray{T}        # Thermal BEM H
    Gh                          # Thermal BEM G
    Hv::AbstractArray{T}        # Viscous BEM H
    Gv                          # Viscous BEM GH
    Nd::AbstractArray           # [diag(nx);diag(ny);diag(nz)]
    Dc::AbstractArray           # [Dx; Dy; Dc]
    Dr::AbstractArray           # [Dx  Dy  Dz]
    inner::LossyGlobalInner{T}  # Inner matrix
    ### Constants
    # For the constraint: vᵦ = ϕₐ∇pₐ + ϕₕ∇pₕ + vᵥ on Γ
    phi_a::T
    phi_h::T
    # For the constraint: T = τₐpₐ + τₕpₕ = 0 on Γ
    tau_a::T
    tau_h::T
    # Combinations
    mu_a::T
    mu_h::T
end
#==========================================================================================
                    Constructor (Assembling) of a LossyBlockMatrix
==========================================================================================#
"""
    LossyGlobalOuter(mesh::Mesh, freq;
                m=3,n=3,l=90,p=90,S=-1,sparsity=20.0,
                exterior=true,adaptive=false)

A `LinearMap` corresponding to the reduced lossy system.
"""
function LossyGlobalOuter(mesh::Mesh,freq;
                            progress=true,integral_free_term=[],
                            hmatrix_on=false,sparse_offset=nothing,
                            depth=1,sparse_assembly=true,exterior=true,sparse_lu=false,
                            m=3,n=3,S=1,fmm_on=false,nearfield=true,thres=1e-6,fmm_offset=0.2)
    if fmm_on == false && size(mesh.normals,2) > 20000
        @warn "Using a dense formulation with a mesh of this size can be problematic"
    end
    if (typeof(mesh.physics_function) <: DiscontinuousTriangularConstant)
        ArgumentError("Constant elements will have a tangential derivative equal to zero.")
    end
    # Computing physical constants
    # rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=S)

    ### Extracting our sources
    sources = mesh.sources

    ### Extracting the number of nodes (and hence also the total matrix size)
    nSource = size(sources,2)

    # Defining Diagonal Entries
    if isempty(integral_free_term)
        C0 = Diagonal(ones(eltype(kₕ),nSource)/2)
    elseif length(integral_free_term) == nSource
        C0 = Diagonal(integral_free_term)
    elseif !(length(integral_free_term) == nSource)
        throw(DimensionMismatch("Length of user-specified integral free terms,"*
        "$(length(integral_free_term)), is not equal to the number of sources, $(nSource)"))
    end
    ### Assembling the 3 BEM systems
    if progress; @info("Acoustic Matrices:"); end
    if fmm_on && hmatrix_on
        throw(ArgumentError("You can not both use the FMM and H-matrices"))
    elseif fmm_on && !hmatrix_on
        Ga = FMMGOperator(mesh,kₐ;
                        n_gauss=n,tol=thres,offset=fmm_offset,nearfield=nearfield,depth=depth)
        Ha = FMMFOperator(mesh,kₐ;
                        n_gauss=n,tol=thres,offset=fmm_offset,nearfield=nearfield,depth=depth) + I/2
    elseif !fmm_on && hmatrix_on
        Ga = HGOperator(mesh,kₐ;n_gauss=n,tol=thres,offset=fmm_offset,nearfield=nearfield,depth=depth)
        Ha = HFOperator(mesh,kₐ;n_gauss=n,tol=thres,offset=fmm_offset,nearfield=nearfield,depth=depth)
    else
        Ha,Ga,C = assemble_parallel!(mesh,kₐ,sources;m=m,n=n,progress=progress)
        C0 = (exterior ? Diagonal(C) : Diagonal(1.0 - C))
        Ha = (exterior ? Ha + C0 : -Ha + C0)
    end
    # Thermal matrices
    if progress; @info("Thermal Matrices:"); end
    Hh,Gh = assemble_parallel!(mesh,kₕ,sources;offset=sparse_offset,
                        sparse=sparse_assembly,depth=depth,progress=progress);
    Hh = (exterior ?  Hh + C0 : -Hh + C0)
    # Viscous matrices
    if progress; @info("Viscous matrices:"); end
    Fᵥ,Bᵥ  = assemble_parallel!(mesh,kᵥ,sources;offset=sparse_offset,
                        sparse=sparse_assembly,depth=depth,progress=progress);
    Aᵥ = (exterior ?  Fᵥ + C0 : -Fᵥ + C0)
    ### Computing tangential derivatives
    Dx,Dy,Dz = interpolation_function_derivatives(mesh)

    nx = mesh.normals[1,:]
    ny = mesh.normals[2,:]
    nz = mesh.normals[3,:]

    ## Creating NullDivergence
    Dc = [Dx;Dy;Dz]
    Dr = [Dx Dy Dz]
    # Plus or minus?
    Nd = -[sparse(Diagonal(nx)); sparse(Diagonal(ny)); sparse(Diagonal(nz))]
    # Nd = [sparse(Diagonal(nx)); sparse(Diagonal(ny)); sparse(Diagonal(nz))]
    Gv = blockdiag(Bᵥ, Bᵥ, Bᵥ)
    Hv = blockdiag(Aᵥ, Aᵥ, Aᵥ)

    mu_a = ϕₐ - τₐ*ϕₕ/τₕ
    mu_h =      τₐ*ϕₕ/τₕ # ϕₐ - mu_a

    if sparse_lu
        vlu = lu(Gv)
        inner = LossyGlobalInner(nSource,Hv,vlu,Nd,Dr)
        return LossyGlobalOuter(nSource,
                                Ha,Ga,
                                Hh,lu(Gh),
                                Hv,vlu,
                                Nd,Dc,Dr,
                                inner,
                                ϕₐ,ϕₕ,τₐ,τₕ,mu_a,mu_h)
    else
        inner = LossyGlobalInner(nSource,Hv,Gv,Nd,Dr)
        return LossyGlobalOuter(nSource,
                            Ha,Ga,
                            Hh,Gh,
                            Hv,Gv,
                            Nd,Dc,Dr,
                            inner,
                            ϕₐ,ϕₕ,τₐ,τₕ,mu_a,mu_h)
    end
end
#==========================================================================================
    Defining relevant routines for LinearMaps.jl to work on the LossyBlockMatrix format
==========================================================================================#
# Size. Required for the LinearMaps.jl (and IterativeSolvers.jl package)
Base.size(A::LossyGlobalOuter) = (A.n, A.n)
# Defining multiplication with `LossyGlobalOuter`
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::LossyGlobalOuter{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions of input
    LinearMaps.check_dim_mul(y, A, x)
    # Adding the contribution from Ha
    y .= -A.phi_a*(A.Ha*x)
    # Create temporary
    if typeof(A.Gh) <: Factorization
        y .+= A.Ga*(A.mu_h*(A.Gh\(A.Hh*x)) +
                        A.mu_a*gmres(A.inner, A.Dr*(A.Dc*x) -
                        A.Nd'*(A.Gv\(A.Hv*(A.Dc*x)))))
    else
        y .+= A.Ga*(A.mu_h*gmres(A.Gh,A.Hh*x) +
                    A.mu_a*gmres(A.inner, A.Dr*(A.Dc*x) -
                    A.Nd'*gmres(A.Gv,A.Hv*(A.Dc*x))))
    end
    return y
end

function _full10(A::LossyGlobalOuter)

    n = size(A,1)
    F = zeros(eltype(A), 10n,10n)

    # Acoustical Mode
    F[0n+1:1n,0n+1:1n] = A.Ga
    F[0n+1:1n,1n+1:2n] = A.Ha
    # Thermal Mode
    F[1n+1:2n,2n+1:3n] = A.Gh
    F[1n+1:2n,3n+1:4n] = A.Hh
    # Viscous Mode
    F[2n+1:5n,4n+1:7n]  = A.Gv
    F[2n+1:5n,7n+1:10n] = A.Hv
    # Divergence
    F[5n+1:6n,4n+1:7n] = A.Dr
    F[5n+1:6n,7n+1:10n] = A.Nd'
    # Isothermal BC
    F[6n+1:7n,0n+1:1n] = A.tau_a*Diagonal(ones(n))
    F[6n+1:7n,2n+1:3n] = A.tau_h*Diagonal(ones(n))
    # No-slip
    F[7n+1:10n,0n+1:1n] = A.phi_a*A.Dc
    F[7n+1:10n,1n+1:2n] = A.phi_a*A.Nd
    F[7n+1:10n,2n+1:3n] = A.phi_h*A.Dc
    F[7n+1:10n,3n+1:4n] = A.phi_h*A.Nd
    F[7n+1:10n,4n+1:7n] = Diagonal(ones(3n))

    return F
end

function _full4(A::LossyGlobalOuter)

    n = size(A,1)
    F = zeros(eltype(A), 4n, 4n)

    R = A.Dr - A.Nd'*(A.Gv\Matrix(A.Hv))

    F[0n+1:1n,1n+1:4n] = R

    F[1n+1:4n,0n+1:1n] = A.mu_a*A.Dc + A.Nd*(A.mu_h*(A.Gh\Matrix(A.Hh)) - A.phi_a*(A.Ga\A.Ha))

    F[1n+1:4n,1n+1:4n] = Diagonal(ones(3n))


    return F
end

function _full1(A::LossyGlobalOuter)

    R = A.Dr - A.Nd'*(A.Gv\Matrix(A.Hv))

    F = A.Ga*(A.mu_a*((R*A.Nd)\Matrix(R*A.Dc)) + A.mu_h*((A.Gh)\Matrix(A.Hh))) - A.phi_a*A.Ha

    return F
end
