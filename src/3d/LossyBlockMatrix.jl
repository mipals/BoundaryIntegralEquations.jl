#==========================================================================================
                            Defining LossyBlockMatrix
            Block matrix corresponding to 5 BEM systems and 5 constraints
==========================================================================================#
struct LossyBlockMatrix{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    N::Int64                            # Total system size (10n)
    # Acoustic Matrices
    Aₐ::AbstractArray{T}                # Acoustic BEM A
    Bₐ::AbstractArray{T}                # Acoustic BEM B
    # Thermal Matrices
    Aₕ::AbstractArray{T}                # Thermal BEM A
    Bₕ::AbstractArray{T}                # Thermal BEM B
    # Viscosity Matrices
    Aᵥ::AbstractArray{T}                # Viscous BEM A
    Bᵥ::AbstractArray{T}                # Viscous BEM B
    # Normals
    normals::AbstractArray{T}
    # Tangents
    tangents₁::AbstractArray{T}         # Tangents in direction 1
    tangents₂::AbstractArray{T}         # Tangents in direction 2
    # Tangential Derivatives
    Dt₁::AbstractArray{T}               # Tangential derivative in direction 1
    Dt₂::AbstractArray{T}               # Tangential derivative in direction 2
    # # Null divergence
    NullDivergence::AbstractArray{T}
    ### Constants
    # For the constraint: vᵦ = ϕₐ∇pₐ + ϕₕ∇pₕ + vᵥ on Γ      (<- This is a vector equation)
    ϕₐ::T                        # Acoustic gradient coupling parameter
    ϕₕ::T                        # Thermal gradient coupling parameter
    # For the constraint: T = τₐpₐ + τₕpₕ = 0 on Γ
    τₐ::T                        # Acoustic coupling parameter
    τₕ::T                        # Thermal coupling parameter
end

#==========================================================================================
                    Constructor (Assembling) of a LossyBlockMatrix
==========================================================================================#
"""
    LossyBlockMatrix(mesh::Mesh, freq;
                m=3,n=3,l=90,p=90,S=-1,sparsity=20.0,
                exterior=true,adaptive=false,blockoutput=false)

Computes the Block matrix corresponding to the reduced lossy system.
If `blockoutput=false` returns sparse matrix.
If `blockoutput=true` returns a `LossyBlockMatrix` struct used for iterative solvers
"""
function LossyBlockMatrix(mesh::Mesh,freq;depth=2,
                        m=3,n=3,S=-1,exterior=true,blockoutput=false)
    if (typeof(mesh.physics_function) <: DiscontinuousTriangularConstant)
        ArgumentError("Constant elements will have a tangential derivative equal to zero.")
    end
    # Computing physical constants
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=S)
    # Computing boundary layer thickness (approximate)
    dᵥ = 2.1e-3/sqrt(freq) # Viscous
    dₕ = 2.5e-3/sqrt(freq) # Thermal

    ### Extracting our sources
    sources = mesh.sources

    ### Extracting the number of nodes (and hence also the total matrix size)
    nSource = size(sources,2)
    N = 10*nSource

    ### Assembling the 3 BEM systems
    # Acoustic matrices
    println("Acoustic Matrices:")
    Fₐ,Bₐ,C₀ = assemble_parallel!(mesh,kₐ,sources;m=m,n=n)
    Aₐ = (exterior ? Fₐ + Diagonal(1.0 .- C₀) : Fₐ - Diagonal(C₀))
    Bₐ = (exterior ? Bₐ : Bₐ)
    # Thermal matrices
    println("Thermal Matrices:")
    Fₕ,Bₕ, = assemble_parallel!(mesh,kᵥ,sources;sparse=true,depth=depth);
    Aₕ = (exterior ? Fₕ + Diagonal(1.0 .- C₀) : Fₕ - Diagonal(C₀))
    Bₕ = (exterior ? Bₕ : Bₕ)
    # Viscous matrices
    println("Viscous matrices:")
    Fᵥ,Bᵥ  = assemble_parallel!(mesh,kₕ,sources;sparse=true,depth=depth);
    Aᵥ = (exterior ? Fᵥ + Diagonal(1.0 .- C₀) : Fᵥ - Diagonal(C₀))
    Bᵥ = (exterior ? Bᵥ : Bᵥ)

    #### Extracting the normal and tangent direction
    normals   = convert(typeof(Bᵥ),mesh.normals)
    tangents₁ = convert(typeof(Bᵥ),mesh.tangents)
    tangents₂ = convert(typeof(Bᵥ),mesh.sangents)

    ### Computing tangential derivatives
    Dt₁, Dt₂ = shape_function_derivatives(mesh)
    Dt₁ = convert.(eltype(Bᵥ),Dt₁)
    Dt₂ = convert.(eltype(Bᵥ),Dt₂)

    ## Creating NullDivergence
    ND = [(tangents₁[1,:].*Dt₁) (tangents₁[2,:].*Dt₁) (tangents₁[3,:].*Dt₁)] +
         [(tangents₂[1,:].*Dt₂) (tangents₂[2,:].*Dt₂) (tangents₂[3,:].*Dt₂)]

    LossyBlock = LossyBlockMatrix(nSource,N,
                        Aₐ,Bₐ,sparse(Aₕ),sparse(Bₕ),sparse(Aᵥ),sparse(Bᵥ),
                        normals,tangents₁,tangents₂,Dt₁,Dt₂,ND,
                        ϕₐ,ϕₕ,τₐ,τₕ)
    if blockoutput
        return LossyBlock
    else
        return full(LossyBlock)
    end
end

#==========================================================================================
    Defining relevant routines for LinearMaps.jl to work on the LossyBlockMatrix format
==========================================================================================#
# Size. Required for the LinearMaps.jl (and IterativeSolvers.jl package)
Base.size(A::LossyBlockMatrix) = (A.N, A.N)
# Standard Multiplication
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::LossyBlockMatrix{T},
                            x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    # Extracting relevant values
    n   = A.n
    ϕₐ  = A.ϕₐ
    ϕₕ  = A.ϕₕ
    τₐ  = A.τₐ
    τₕ  = A.τₕ
    Dt₁ = A.Dt₁
    Dt₂ = A.Dt₂
    # Normal vectors
    Nx = @view A.normals[1,:]
    Ny = @view A.normals[2,:]
    Nz = @view A.normals[3,:]
    # Tangents vectors 1
    Tx1 = @view A.tangents₁[1,:]
    Ty1 = @view A.tangents₁[2,:]
    Tz1 = @view A.tangents₁[3,:]
    # Tangent vectors 2
    Tx2 = @view A.tangents₂[1,:]
    Ty2 = @view A.tangents₂[2,:]
    Tz2 = @view A.tangents₂[3,:]
    # BEM-matrices
    Aₐ = A.Aₐ
    Bₐ = A.Bₐ
    Aₕ = A.Aₕ
    Bₕ = A.Bₕ
    Aᵥ = A.Aᵥ
    Bᵥ = A.Bᵥ

    # Acoustic variables
    pₐ  = @view x[0*n+1:1*n]
    vₐ  = @view x[1*n+1:2*n]
    # Thermal Variables
    pₕ  = @view x[2*n+1:3*n]
    vₕ  = @view x[3*n+1:4*n]
    # Viscous-Velocity x-components
    v₁  = @view x[4*n+1:5*n]
    ∂v₁ = @view x[5*n+1:6*n]
    # Viscous-Velocity y-components
    v₂  = @view x[6*n+1:7*n]
    ∂v₂ = @view x[7*n+1:8*n]
    # Viscous-Velocity z-components
    v₃  = @view x[8*n+1:9*n]
    ∂v₃ = @view x[9*n+1:end]

    # First 4 columns:             pₐ          vₐ                pₕ           vₕ
    y[0*n+1:1*n] .=           Aₐ * pₐ  -  Bₐ * vₐ
    y[1*n+1:2*n] .=                                         Aₕ * pₕ   -   Bₕ * vₕ
  # y[2*n+1:3*n] .= zeros(eltype(A),n)
  # y[3*n+1:4*n] .= zeros(eltype(A),n)
  # y[4*n+1:5*n] .= zeros(eltype(A),n)
    y[5*n+1:6*n] .=         τₐ * pₐ +                       τₕ * pₕ
    y[6*n+1:7*n] .=                       ϕₐ * vₐ +                      ϕₕ * vₕ
    y[7*n+1:8*n] .= ϕₐ * (Dt₁  * pₐ)  +           ϕₕ * (Dt₁  * pₕ)
    y[8*n+1:9*n] .= ϕₐ * (Dt₂  * pₐ)  +           ϕₕ * (Dt₂  * pₕ)
    # y[9*n+1:end] .= zeros(eltype(A),n)

    #  Last 6 columns:     v₁          ∂v₁         v₂          ∂v₂           v₃         ∂v₃
   #y[0*n+1:1*n] +=
   #y[1*n+1:2*n] +=
    y[2*n+1:3*n] .=   Aᵥ * v₁  -  Bᵥ * ∂v₁
    y[3*n+1:4*n] .=                           Aᵥ * v₂  -  Bᵥ * ∂v₂
    y[4*n+1:5*n] .=                                                     Aᵥ * v₃  -  Bᵥ * ∂v₃
   #y[5*n+1:6*n] +=
    y[6*n+1:7*n] += Nx  .* v₁         +      Ny .* v₂         +        Nz  .* v₃
    y[7*n+1:8*n] += Tx1 .* v₁         +      Ty1.* v₂         +        Tz1 .* v₃
    y[8*n+1:9*n] += Tx2 .* v₁         +      Ty2.* v₂         +        Tz2 .* v₃
    y[9*n+1:end] .= A.NullDivergence * [v₁; v₂; v₃] + Nx .* ∂v₁ +  Ny .* ∂v₂ + Nz .* ∂v₃
    return y
end
# Multiplication with adjoint
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            B::LinearMaps.AdjointMap{T,S},
                            x::AbstractVector) where {T <: ComplexF64, S <: LossyBlockMatrix}
    # Checking Dimenisons
    LinearMaps.check_dim_mul(y, B, x)
    # Extracting relevant values
    A   = B.lmap
    n   = A.n
    ϕₐ  = A.ϕₐ
    ϕₕ  = A.ϕₕ
    τₐ  = A.τₐ
    τₕ  = A.τₕ
    Dt₁ = A.Dt₁
    Dt₂ = A.Dt₂

    #
    x1  = @view x[0*n+1:1*n] # Acoustic variable
    x2  = @view x[1*n+1:2*n] # Acoustic variable
    x3  = @view x[2*n+1:3*n] # Thermal Variable
    x4  = @view x[3*n+1:4*n] # Thermal Variable
    x5  = @view x[4*n+1:5*n] # Viscous-Velocity x-components
    x6  = @view x[5*n+1:6*n] # Viscous-Velocity x-components
    x7  = @view x[6*n+1:7*n] # Viscous-Velocity y-components
    x8  = @view x[7*n+1:8*n] # Viscous-Velocity y-components
    x9  = @view x[8*n+1:9*n] # Viscous-Velocity z-components
    x10 = @view x[9*n+1:end] # Viscous-Velocity z-components

    # Extracing normal vectors
    Nx  = @view A.normals[1,:]
    Ny  = @view A.normals[2,:]
    Nz  = @view A.normals[3,:]
    # Extracing first tangential direction
    Tx1 = @view A.tangents₁[1,:]
    Ty1 = @view A.tangents₁[2,:]
    Tz1 = @view A.tangents₁[3,:]
    # Extracing second tangential direction
    Tx2 = @view A.tangents₂[1,:]
    Ty2 = @view A.tangents₂[2,:]
    Tz2 = @view A.tangents₂[3,:]
    # NullDivergence
    ND1 = @view (A.NullDivergence[:,0*n+1:1*n])
    ND2 = @view (A.NullDivergence[:,1*n+1:2*n])
    ND3 = @view (A.NullDivergence[:,2*n+1:3*n])

    # BEM-matrices
    Aₐ = A.Aₐ
    Bₐ = A.Bₐ
    Aₕ = A.Aₕ
    Bₕ = A.Bₕ
    Aᵥ = A.Aᵥ
    Bᵥ = A.Bᵥ

    # First 5 columns:        x1             x2            x3             x4           x5
    y[0*n+1:1*n] .=   (Aₐ' * x1)
    y[1*n+1:2*n] .=  -(Bₐ' * x1)
    y[2*n+1:3*n] .=                 (Aₕ' * x2)
    y[3*n+1:4*n] .=                -(Bₕ' * x2)
    y[4*n+1:5*n] .=                                (Aᵥ' * x3)
    y[5*n+1:6*n] .=                               -(Bᵥ' * x3)
    y[6*n+1:7*n] .=                                               (Aᵥ' * x4)
    y[7*n+1:8*n] .=                                              -(Bᵥ' * x4)
    y[8*n+1:9*n] .=                                                             (Aᵥ' * x5)
    y[9*n+1:end] .=                                                            -(Bᵥ' * x5)

    # Next 4 columns:         x6        x7                x8                 x9
    y[0*n+1:1*n] +=     τₐ' * x6 +          Dt₁' * (ϕₐ' * x8) + Dt₂' * (ϕₐ' * x9)
    y[1*n+1:2*n] +=               ϕₐ' * x7
    y[2*n+1:3*n] +=     τₕ' * x6 +          Dt₁' * (ϕₕ' * x8)  + Dt₂' * (ϕₕ' * x9)
    y[3*n+1:4*n] +=               ϕₕ' * x7
    y[4*n+1:5*n] +=               Nx .* x7    +    Tx1 .* x8    +     Tx2 .* x9
   #y[5*n+1:6*n] +=
    y[6*n+1:7*n] +=               Ny .* x7    +    Ty1 .* x8    +     Ty2 .* x9
   #y[7*n+1:8*n] +=
    y[8*n+1:9*n] +=               Nz .* x7    +    Tz1 .* x8    +     Tz2 .* x9
   #y[9*n+1:end] .=

    # Last column                                      x10
    y[4*n+1:5*n] += ND1' * x10
    y[5*n+1:6*n] += Nx .* x10
    y[6*n+1:7*n] += ND2' * x10
    y[7*n+1:8*n] += Ny .* x10
    y[8*n+1:9*n] += ND3' * x10
    y[9*n+1:end] += Nz .* x10
    return y
end

#==========================================================================================
                                Utility functions
==========================================================================================#
function full(B::LinearMaps.AdjointMap{T,S}) where {T <: ComplexF64, S <: LossyBlockMatrix}
    M = B.lmap
    n = M.n
    n = M.n
    # Prealloaction
    A = spzeros(eltype(M),M.N,M.N)

    # Column 1:
    A[0*n+1:1*n,0*n+1:1*n] =  M.Aₐ'
    A[1*n+1:2*n,0*n+1:1*n] = -M.Bₐ'
    # Column 2
    A[2*n+1:3*n,1*n+1:2*n] =  M.Aₕ'
    A[3*n+1:4*n,1*n+1:2*n] = -M.Bₕ'
    # Column 3
    A[4*n+1:5*n,2*n+1:3*n] =  M.Aᵥ'
    A[5*n+1:6*n,2*n+1:3*n] = -M.Bᵥ'
    # Column 4
    A[6*n+1:7*n,3*n+1:4*n] =  M.Aᵥ'
    A[7*n+1:8*n,3*n+1:4*n] = -M.Bᵥ'
    # Column 5
    A[8*n+1:9*n,4*n+1:5*n] =  M.Aᵥ'
    A[9*n+1:end,4*n+1:5*n] = -M.Bᵥ'
    # Column 6
    A[0*n+1:1*n,5*n+1:6*n] = Diagonal(ones(eltype(M),n)*M.τₐ)'
    A[2*n+1:3*n,5*n+1:6*n] = Diagonal(ones(eltype(M),n)*M.τₕ)'
    # Column 7
    A[1*n+1:2*n,6*n+1:7*n] = Diagonal(ones(eltype(M),n)*M.ϕₐ)'
    A[3*n+1:4*n,6*n+1:7*n] = Diagonal(ones(eltype(M),n)*M.ϕₕ)'
    A[4*n+1:5*n,6*n+1:7*n] = Diagonal(M.normals[1,:])'
    A[6*n+1:7*n,6*n+1:7*n] = Diagonal(M.normals[2,:])'
    A[8*n+1:9*n,6*n+1:7*n] = Diagonal(M.normals[3,:])'
    # Column 8
    A[0*n+1:1*n,7*n+1:8*n] = (M.Dt₁*M.ϕₐ)'
    A[2*n+1:3*n,7*n+1:8*n] = (M.Dt₁*M.ϕₕ)'
    A[4*n+1:5*n,7*n+1:8*n] = Diagonal(M.tangents₁[1,:])'
    A[6*n+1:7*n,7*n+1:8*n] = Diagonal(M.tangents₁[2,:])'
    A[8*n+1:9*n,7*n+1:8*n] = Diagonal(M.tangents₁[3,:])'
    # Column 9
    A[0*n+1:1*n,8*n+1:9*n] = (M.Dt₂*M.ϕₐ)'
    A[2*n+1:3*n,8*n+1:9*n] = (M.Dt₂*M.ϕₕ)'
    A[4*n+1:5*n,8*n+1:9*n] = Diagonal(M.tangents₂[1,:])'
    A[6*n+1:7*n,8*n+1:9*n] = Diagonal(M.tangents₂[2,:])'
    A[8*n+1:9*n,8*n+1:9*n] = Diagonal(M.tangents₂[3,:])'
    # Column 10
    A[4*n+1:5*n,9*n+1:10*n] = M.NullDivergence[:,0*n+1:1*n]'
    A[6*n+1:7*n,9*n+1:10*n] = M.NullDivergence[:,1*n+1:2*n]'
    A[8*n+1:9*n,9*n+1:10*n] = M.NullDivergence[:,2*n+1:3*n]'
    A[5*n+1:6*n,9*n+1:10*n] = Diagonal(M.normals[1,:])'
    A[7*n+1:8*n,9*n+1:10*n] = Diagonal(M.normals[2,:])'
    A[9*n+1:end,9*n+1:10*n] = Diagonal(M.normals[3,:])'

    return A

end

"""
    full(M::LossyBlockMatrix)

Returns a sparse matrix representation of the `LossyBlockMatrix`.
"""
function full(M::LossyBlockMatrix)
    n = M.n
    A = spzeros(eltype(M),M.N,M.N)
    # Row 1: Acoustics BEM matrices
    A[0*n+1:1*n,0*n+1:1*n] =  M.Aₐ
    A[0*n+1:1*n,1*n+1:2*n] = -M.Bₐ
    # Row 2: Thermal BEM matrices
    A[1*n+1:2*n,2*n+1:3*n] =  sparse(M.Aₕ)
    A[1*n+1:2*n,3*n+1:4*n] = -sparse(M.Bₕ)
    # Row 3-5: Viscous BEM matrices
    A[2*n+1:3*n,4*n+1:5*n] =  sparse(M.Aᵥ)
    A[2*n+1:3*n,5*n+1:6*n] = -sparse(M.Bᵥ)
    A[3*n+1:4*n,6*n+1:7*n] =  sparse(M.Aᵥ)
    A[3*n+1:4*n,7*n+1:8*n] = -sparse(M.Bᵥ)
    A[4*n+1:5*n,8*n+1:9*n] =  sparse(M.Aᵥ)
    A[4*n+1:5*n,9*n+1:end] = -sparse(M.Bᵥ)
    # Row 6: Pressure constraints
    A[5*n+1:6*n,0*n+1:1*n] = Diagonal(ones(eltype(M),n)*M.τₐ)
    A[5*n+1:6*n,2*n+1:3*n] = Diagonal(ones(eltype(M),n)*M.τₕ)
    # Row 7: Normal velocity constraint
    A[6*n+1:7*n,1*n+1:2*n] = Diagonal(ones(eltype(M),n)*M.ϕₐ)
    A[6*n+1:7*n,3*n+1:4*n] = Diagonal(ones(eltype(M),n)*M.ϕₕ)
    A[6*n+1:7*n,4*n+1:5*n] = Diagonal(M.normals[1,:])
    A[6*n+1:7*n,6*n+1:7*n] = Diagonal(M.normals[2,:])
    A[6*n+1:7*n,8*n+1:9*n] = Diagonal(M.normals[3,:])
    # Row 8: Tangent-1 velocity constraints
    A[7*n+1:8*n,0*n+1:1*n] = M.Dt₁*M.ϕₐ
    A[7*n+1:8*n,2*n+1:3*n] = M.Dt₁*M.ϕₕ
    A[7*n+1:8*n,4*n+1:5*n] = Diagonal(M.tangents₁[1,:])
    A[7*n+1:8*n,6*n+1:7*n] = Diagonal(M.tangents₁[2,:])
    A[7*n+1:8*n,8*n+1:9*n] = Diagonal(M.tangents₁[3,:])
    # Row 9: Tangent-2 velocity constraints
    A[8*n+1:9*n,0*n+1:1*n] = M.Dt₂*M.ϕₐ
    A[8*n+1:9*n,2*n+1:3*n] = M.Dt₂*M.ϕₕ
    A[8*n+1:9*n,4*n+1:5*n] = Diagonal(M.tangents₂[1,:])
    A[8*n+1:9*n,6*n+1:7*n] = Diagonal(M.tangents₂[2,:])
    A[8*n+1:9*n,8*n+1:9*n] = Diagonal(M.tangents₂[3,:])
    # Row 10: Null-divergence
    A[9*n+1:10*n,5*n+1:6*n] = Diagonal(M.normals[1,:])
    A[9*n+1:10*n,7*n+1:8*n] = Diagonal(M.normals[2,:])
    A[9*n+1:10*n,9*n+1:end] = Diagonal(M.normals[3,:])
    A[9*n+1:10*n,4*n+1:5*n] = M.NullDivergence[:,0*n+1:1*n]
    A[9*n+1:10*n,6*n+1:7*n] = M.NullDivergence[:,1*n+1:2*n]
    A[9*n+1:10*n,8*n+1:9*n] = M.NullDivergence[:,2*n+1:3*n]

    return A
end

"""
    fullGlobal(M::LossyBlockMatrix,mesh)

Returns a sparse matrix representation of the `LossyBlockMatrix`.
The caveat being that the shape function derivatives are computed from the global gradient.
"""
function fullGlobal(M::LossyBlockMatrix,mesh)
    n = M.n
    A = spzeros(eltype(M),M.N,M.N)
    Dx,Dy,Dz = globalCoordinateShapeFunctionDerivative(mesh)
    nvect0 = mesh.normals
    tvect1 = mesh.tangentsX
    tvect2 = mesh.tangentsY
    Dt1 = Dx .* tvect1[1,:] + Dy .* tvect1[2,:] + Dz .* tvect1[3,:]
    Dt2 = Dx .* tvect2[1,:] + Dy .* tvect2[2,:] + Dz .* tvect2[3,:]
    one = Diagonal(ones(eltype(M),n))
    # Row 1: Acoustics BEM matrices
    A[0*n+1:1*n,0*n+1:1*n] =  M.Aₐ
    A[0*n+1:1*n,1*n+1:2*n] = -M.Bₐ
    # Row 2: Thermal BEM matrices
    A[1*n+1:2*n,2*n+1:3*n] =  sparse(M.Aₕ)
    A[1*n+1:2*n,3*n+1:4*n] = -sparse(M.Bₕ)
    # Row 3-5: Viscous BEM matrices
    A[2*n+1:3*n,4*n+1:5*n] =  sparse(M.Aᵥ)
    A[2*n+1:3*n,5*n+1:6*n] = -sparse(M.Bᵥ)
    A[3*n+1:4*n,6*n+1:7*n] =  sparse(M.Aᵥ)
    A[3*n+1:4*n,7*n+1:8*n] = -sparse(M.Bᵥ)
    A[4*n+1:5*n,8*n+1:9*n] =  sparse(M.Aᵥ)
    A[4*n+1:5*n,9*n+1:end] = -sparse(M.Bᵥ)
    # Row 6: Isothermal boundary coupling
    A[5*n+1:6*n,0*n+1:1*n] = one*M.τₐ
    A[5*n+1:6*n,2*n+1:3*n] = one*M.τₕ
    # Row 7: x-velocity coupling
    A[6*n+1:7*n,0*n+1:1*n] = Dx*M.ϕₐ
    A[6*n+1:7*n,2*n+1:3*n] = Dx*M.ϕₕ
    A[6*n+1:7*n,4*n+1:5*n] = one
    # A[6*n+1:7*n,6*n+1:7*n] = one
    # A[6*n+1:7*n,8*n+1:9*n] = one
    # Row 8: y-velocity coupling
    A[7*n+1:8*n,0*n+1:1*n] = Dy*M.ϕₐ
    A[7*n+1:8*n,2*n+1:3*n] = Dy*M.ϕₕ
    # A[7*n+1:8*n,4*n+1:5*n] = one
    A[7*n+1:8*n,6*n+1:7*n] = one
    # A[7*n+1:8*n,8*n+1:9*n] = one
    # Row 9: z-velocity coupling
    A[8*n+1:9*n,0*n+1:1*n] = Dz*M.ϕₐ
    A[8*n+1:9*n,2*n+1:3*n] = Dz*M.ϕₕ
    # A[8*n+1:9*n,4*n+1:5*n] = one
    # A[8*n+1:9*n,6*n+1:7*n] = one
    A[8*n+1:9*n,8*n+1:9*n] = one
    # Row 10: Null-divergence constraint
    A[9*n+1:10*n,4*n+1:5*n] = Dt1 + Dt2
    A[9*n+1:10*n,5*n+1:6*n] = one
    A[9*n+1:10*n,6*n+1:7*n] = Dt1 + Dt2
    A[9*n+1:10*n,7*n+1:8*n] = one
    A[9*n+1:10*n,8*n+1:9*n] = Dt1 + Dt2
    A[9*n+1:10*n,9*n+1:end] = one

    return A
end

"""
    fullGlobal2(M::LossyBlockMatrix,mesh)

Returns a sparse matrix representation of the `LossyBlockMatrix`.
NB! This is experimental. Trying to see if the global gradient can be used directly.
Does not seem to be the case...
"""
function fullGlobal2(M::LossyBlockMatrix,mesh)
    n = M.n
    A = spzeros(eltype(M),M.N,M.N)
    Dx,Dy,Dz = globalCoordinateShapeFunctionDerivative(mesh)
    nvect0 = mesh.normals
    tvect1 = mesh.tangentsX
    tvect2 = mesh.tangentsY
    Dt1 = Dx .* tvect1[1,:] + Dy .* tvect1[2,:] + Dz .* tvect1[3,:]
    Dt2 = Dx .* tvect2[1,:] + Dy .* tvect2[2,:] + Dz .* tvect2[3,:]
    # Row 1: Acoustics BEM matrices
    A[0*n+1:1*n,0*n+1:1*n] =  M.Aₐ
    A[0*n+1:1*n,1*n+1:2*n] = -M.Bₐ
    # Row 2: Thermal BEM matrices
    A[1*n+1:2*n,2*n+1:3*n] =  sparse(M.Aₕ)
    A[1*n+1:2*n,3*n+1:4*n] = -sparse(M.Bₕ)
    # Row 3-5: Viscous BEM matrices
    A[2*n+1:3*n,4*n+1:5*n] =  sparse(M.Aᵥ)
    A[2*n+1:3*n,5*n+1:6*n] = -sparse(M.Bᵥ)
    A[3*n+1:4*n,6*n+1:7*n] =  sparse(M.Aᵥ)
    A[3*n+1:4*n,7*n+1:8*n] = -sparse(M.Bᵥ)
    A[4*n+1:5*n,8*n+1:9*n] =  sparse(M.Aᵥ)
    A[4*n+1:5*n,9*n+1:end] = -sparse(M.Bᵥ)
    # Row 6: Pressure constraints
    A[5*n+1:6*n,0*n+1:1*n] = Diagonal(ones(eltype(M),n)*M.τₐ)
    A[5*n+1:6*n,2*n+1:3*n] = Diagonal(ones(eltype(M),n)*M.τₕ)
    # Row 7: Normal velocity constraint
    A[6*n+1:7*n,1*n+1:2*n] = Diagonal(ones(eltype(M),n)*M.ϕₐ)
    A[6*n+1:7*n,3*n+1:4*n] = Diagonal(ones(eltype(M),n)*M.ϕₕ)
    A[6*n+1:7*n,4*n+1:5*n] = Diagonal(nvect0[1,:])
    A[6*n+1:7*n,6*n+1:7*n] = Diagonal(nvect0[2,:])
    A[6*n+1:7*n,8*n+1:9*n] = Diagonal(nvect0[3,:])
    # Row 8: Tangent-1 velocity constraints
    A[7*n+1:8*n,0*n+1:1*n] = Dt1*M.ϕₐ
    A[7*n+1:8*n,2*n+1:3*n] = Dt1*M.ϕₕ
    A[7*n+1:8*n,4*n+1:5*n] = Diagonal(tvect1[1,:])
    A[7*n+1:8*n,6*n+1:7*n] = Diagonal(tvect1[2,:])
    A[7*n+1:8*n,8*n+1:9*n] = Diagonal(tvect1[3,:])
    # Row 9: Tangent-2 velocity constraints
    A[8*n+1:9*n,0*n+1:1*n] = Dt2*M.ϕₐ
    A[8*n+1:9*n,2*n+1:3*n] = Dt2*M.ϕₕ
    A[8*n+1:9*n,4*n+1:5*n] = Diagonal(tvect2[1,:])
    A[8*n+1:9*n,6*n+1:7*n] = Diagonal(tvect2[2,:])
    A[8*n+1:9*n,8*n+1:9*n] = Diagonal(tvect2[3,:])
    # Row 10: Null-divergence
    A[9*n+1:10*n,4*n+1:5*n] = Dt1.*tvect1[1,:] + Dt2.*tvect2[1,:]
    A[9*n+1:10*n,5*n+1:6*n] = Diagonal(M.normals[1,:])
    A[9*n+1:10*n,6*n+1:7*n] = Dt1.*tvect1[2,:] + Dt2.*tvect2[2,:]
    A[9*n+1:10*n,7*n+1:8*n] = Diagonal(M.normals[2,:])
    A[9*n+1:10*n,8*n+1:9*n] = Dt1.*tvect1[3,:] + Dt2.*tvect2[3,:]
    A[9*n+1:10*n,9*n+1:end] = Diagonal(M.normals[3,:])

    return A
end
