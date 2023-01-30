#==========================================================================================
                            Defining LossyGlobalInner
==========================================================================================#
struct LossyGlobalInner{T} <: LinearMaps.LinearMap{T}
    n::Int64                # Number of nodes
    # Acoustic Matrices
    Hv::AbstractArray{T}    # Viscous BEM A
    Gv::AbstractArray{T}    # Viscous BEM B
    # Normal transformation
    Nd::AbstractArray       # Collection of normal tr
    # Gradients
    Dr::AbstractArray       # [Dx  Dy  Dz]
    # Temporary allocations. For testing purposes
    tmp1::AbstractArray{T}
    tmp2::AbstractArray{T}
end
Base.size(A::LossyGlobalInner) = (A.n, A.n)
# Defining multiplication with `LossyGlobalInner`
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::LossyGlobalInner{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions of input
    LinearMaps.check_dim_mul(y, A, x)
    # Extracting relevant values
    y .= A.Dr*(A.Nd*x) - A.Nd'*(gmres(A.Gv,A.Hv*(A.Nd*x)))
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
    Gh::AbstractArray{T}        # Thermal BEM G
    Hv::AbstractArray{T}        # Viscous BEM H
    Gv::AbstractArray{T}        # Viscous BEM GH
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
                exterior=true,adaptive=false,blockoutput=false)

Computes the Block matrix corresponding to the reduced lossy system.
If `blockoutput=false` returns sparse matrix.
If `blockoutput=true` returns a `LossyBlockMatrix` struct used for iterative solvers
"""
function LossyGlobalOuter(mesh::Mesh,freq;
                            progress=true,integral_free_term=[],
                            depth=1,sparse_assembly=true,exterior=true,
                            m=3,n=3,S=1,fmm_on=false,nearfield=true,thres=1e-6,offset=0.2)
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

    ### Assembling the 3 BEM systems
    if progress; @info("Acoustic Matrices:"); end
    if fmm_on
        # Defining Diagonal Entries
        if isempty(integral_free_term)
            C0 = Diagonal(ones(eltype(kₕ),nSource)/2)
        elseif length(integral_free_term) == nSource
            C0 = Diagonal(integral_free_term)
        elseif !(length(integral_free_term) == nSource)
            throw(DimensionMismatch("Length of user-specified integral free terms,"*
            "$(length(integral_free_term)), is not equal to the number of sources, $(nSource)"))
        end
        Ga = FMMGOperator(mesh,kₐ;
                        n=n,eps=thres,offset=offset,nearfield=nearfield,depth=depth)
        Ha = FMMHOperator(mesh,kₐ;
                        integral_free_term=integral_free_term,
                        n=n,eps=thres,offset=offset,nearfield=nearfield,depth=depth)
    else
        Ha,Ga,C = assemble_parallel!(mesh,kₐ,sources;m=m,n=n,progress=progress)
        C0 = (exterior ? Diagonal(C) : Diagonal(1.0 - C))
        Ha = (exterior ? Ha + C0 : Ha - C0)
    end
    # Thermal matrices
    if progress; @info("Thermal Matrices:"); end
    Hh,Gh = assemble_parallel!(mesh,kₕ,sources;
                        sparse=sparse_assembly,depth=depth,progress=progress);
    Hh = (exterior ?  Hh + C0 : -Hh + C0)
    # Viscous matrices
    if progress; @info("Viscous matrices:"); end
    Fᵥ,Bᵥ  = assemble_parallel!(mesh,kᵥ,sources;sparse=sparse_assembly,depth=depth,progress=progress);
    Aᵥ = (exterior ?  Fᵥ + C0 : -Fᵥ + C0)
    ### Computing tangential derivatives
    Dx,Dy,Dz = shape_function_derivatives(mesh;global_derivatives=true)

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

    tmp1 = zeros(eltype(Hv),3nSource) # Are not used right now
    tmp2 = zeros(eltype(Hv),3nSource) # Are not used right now
    inner = LossyGlobalInner(nSource,Hv,Gv,Nd,Dr,tmp1,tmp2)

    return LossyGlobalOuter(nSource,
                        Ha,Ga,
                        Hh,Gh,
                        Hv,Gv,
                        Nd,Dc,Dr,
                        inner,
                        ϕₐ,ϕₕ,τₐ,τₕ,mu_a,mu_h)
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
    # We only want to call the multiplication with Ga once pr. iteration
    y += A.Ga*(A.mu_h*gmres(A.Gh,A.Hh*x) +
               A.mu_a*gmres(A.inner, A.Dr*(A.Dc*x) -
               A.Nd'*gmres(A.Gv,A.Hv*(A.Dc*x))))
    return y
end
