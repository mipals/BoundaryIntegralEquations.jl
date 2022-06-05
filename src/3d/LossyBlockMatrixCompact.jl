#==========================================================================================
                        Defining LossyBlockMatrixCompact
    Block matrix corresponding to the reduced lossy system (Schur Complemented system)
==========================================================================================#
struct LossyBlockMatrixCompact{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    N::Int64                            # Total system size
    # Row1
    ABah::Matrix{T}                     # BEM matrix
    diagNX::Diagonal{T}                 # x-component of normal vector
    diagNY::Diagonal{T}                 # y-component of normal vector
    diagNZ::Diagonal{T}                 # z-component of normal vector
    # Row2
    NullDivergence::AbstractArray{T}    # Null-divergence constraint
    # Row3
    NoSlip1::AbstractArray{T}           # Contains the first NoSlip condition
    diagTX1::Diagonal{T}                # x-component of 1st tangent vector direction
    diagTY1::Diagonal{T}                # y-component of 1st tangent vector direction
    diagTZ1::Diagonal{T}                # z-component of 1st tangent vector direction
    # Row4
    NoSlip2::AbstractArray{T}           # Contains the second NoSlip condition
    diagTX2::Diagonal{T}                # x-component of 2nd tangent vector direction
    diagTY2::Diagonal{T}                # y-component of 2nd tangent vector direction
    diagTZ2::Diagonal{T}                # z-component of 2nd tangent vector direction
end
#==========================================================================================
                    Constructor/assembly of the LossyBlockMatrixDense
==========================================================================================#
"""
    LossyBlockMatrixCompact(mesh,freq;
            m=3,n=3,l=90,p=90,S=-1,sparsity=20.0,
            exterior=true,adaptive=false,blockoutput=false)

Computes the Block matrix corresponding to the reduced lossy system.
If `blockoutput=false` returns matrix.
If `blockoutput=true` returns a `LossyBlockMatrixCompact` struct used for iterative solvers
"""
function LossyBlockMatrixCompact(mesh,freq;
                m=3,n=3,l=90,p=90,S=-1,exterior=true,blockoutput=false)
    if (typeof(mesh.physics_function) <: DiscontinuousTriangularConstant)
        ArgumentError("Constant elements will have a tangential derivative equal to zero.")
    end
    # Computing physical constants
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=S)
    # Computing boundary layer thickness
    dᵥ = 2.1e-3/sqrt(freq) # Viscous boundary layer thickness
    dₕ = 2.5e-3/sqrt(freq) # Thermal boundary layer thickness
    # Extracting our sources
    sources = mesh.sources

    # Extracting the number of nodes (and hence also the total matrix size)
    nSource = size(sources,2)
    N = 4*nSource

    ### Assembling the 3 BEM systems
    # Acoustic matrices
    println("Acoustic Matrices:")
    Fₐ,Bₐ,C₀ = assemble_parallel!(mesh,kₐ,sources;m=m,n=n)
    Aₐ = (exterior ? Fₐ - Diagonal(1.0 .+ C₀) : Fₐ + Diagonal(C₀))
    # Thermal matrices
    println("Thermal Matrices:")
    println("Thermal Matrices:")
    Fₕ,Bₕ, = assemble_parallel!(mesh,kᵥ,mesh.sources;sparse=true);
    Aₕ =  exterior ? Fₕ - Diagonal(1.0 .+ C₀) : Fₕ - Diagonal(C₀)
    Bₕ = (exterior ? Bₕ : Bₕ)
    # Viscous matrices
    println("Viscous matrices:")
    Fᵥ,Bᵥ  = assemble_parallel!(mesh,kₕ,mesh.sources;sparse=true);
    Aᵥ = exterior ? Fᵥ - Diagonal(1.0 .+ C₀) : Fᵥ - Diagonal(C₀)
    Bᵥ = (exterior ? Bᵥ : Bᵥ)

    ### Extracting the normal and tangent direction
    normals   = convert(typeof(Aᵥ),(exterior ? mesh.normals : -mesh.normals))
    tangents₁ = convert(typeof(Aᵥ),mesh.tangents)
    tangents₂ = convert(typeof(Aᵥ),mesh.sangents)

    ### Computing tangential derivatives
    Dt₁, Dt₂ = shapeFunctionDerivative(mesh)
    Dt₁ = convert.(eltype(Aᵥ),Dt₁)
    Dt₂ = convert.(eltype(Aᵥ),Dt₂)

    # Row1
    ABah   = (ϕₐ*(Bₐ\Aₐ) - (ϕₕ*τₐ/τₕ)*(Matrix(Bₕ)\Matrix(Aₕ)))
    diagNX = Diagonal(normals[1,:])
    diagNY = Diagonal(normals[2,:])
    diagNZ = Diagonal(normals[3,:])
    # Row2: NullDivergence
    ABᵥ = Matrix(Bᵥ)\Matrix(Aᵥ)
    ND = [(  normals[1,:].*ABᵥ) (  normals[2,:].*ABᵥ)   (  normals[3,:].*ABᵥ)] +
         [(tangents₁[1,:].*Dt₁) (tangents₁[2,:].*Dt₁)   (tangents₁[3,:].*Dt₁)] +
         [(tangents₂[1,:].*Dt₂) (tangents₂[2,:].*Dt₂)   (tangents₂[3,:].*Dt₂)]
    # Row3
    NoSlip1 = (ϕₐ - (ϕₕ*τₐ)/τₕ)*Dt₁
    diagTX1 = Diagonal(tangents₁[1,:])
    diagTY1 = Diagonal(tangents₁[2,:])
    diagTZ1 = Diagonal(tangents₁[3,:])
    # Row4
    NoSlip2 = (ϕₐ - (ϕₕ*τₐ)/τₕ)*Dt₂
    diagTX2 = Diagonal(tangents₂[1,:])
    diagTY2 = Diagonal(tangents₂[2,:])
    diagTZ2 = Diagonal(tangents₂[3,:])

    LossyBlockCompact = LossyBlockMatrixCompact(nSource,N,
                                ABah,   diagNX, diagNY, diagNZ,
                                                            ND,
                                sparse(NoSlip1),diagTX1,diagTY1,diagTZ1,
                                sparse(NoSlip2),diagTX2,diagTY2,diagTZ2)
    if blockoutput
        return LossyBlockCompact
    else
        return full(LossyBlockCompact)
    end
end
#==========================================================================================
Defining relevant routines for LinearMaps.jl to work on the LossyBlockMatrixCompact format:
==========================================================================================#
# Size. Required for the LinearMaps.jl (and IterativeSolvers.jl package)
Base.size(A::LossyBlockMatrixCompact) = (A.N, A.N)
# Standard multiplication
function LinearAlgebra.mul!(y::AbstractVecOrMat,
                            A::LossyBlockMatrixCompact,
                            x::AbstractVector)
    # Checking Dimension
    LinearMaps.check_dim_mul(y, A, x)
    # Extract relevant informatinon from A
    n = A.n
    # For convenience
    x1 = @view x[1:n]
    x2 = @view x[n+1:2*n]
    x3 = @view x[2*n+1:3*n]
    x4 = @view x[3*n+1:4*n]
    # Extracting matrices
    ABah = A.ABah
    diagNX = A.diagNX
    diagNY = A.diagNY
    diagNZ = A.diagNZ
    diagTX1 = A.diagTX1
    diagTY1 = A.diagTY1
    diagTZ1 = A.diagTZ1
    diagTX2 = A.diagTX2
    diagTY2 = A.diagTY2
    diagTZ2 = A.diagTZ2
    NoSlip1 = A.NoSlip1
    NoSlip2 = A.NoSlip2
    NullDivergence = A.NullDivergence
    # Multiplication with A
    y[1:n]       .=    ABah*x1 +  diagNX*x2 +  diagNY*x3 +  diagNZ*x4
    y[n+1:2*n]   .=                       NullDivergence * x[n+1:end]
    y[2*n+1:3*n] .= NoSlip1*x1 + diagTX1*x2 + diagTY1*x3 + diagTZ1*x4
    y[3*n+1:4*n] .= NoSlip2*x1 + diagTX2*x2 + diagTY2*x3 + diagTZ2*x4

    return y
end
# Multiplication with the adjoint
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            B::LinearMaps.AdjointMap{T,S},
                            x::AbstractVector) where {T <: ComplexF64, S <: LossyBlockMatrixCompact}

    # Checking Dimension
    LinearMaps.check_dim_mul(y, B, x)
    # Extract relevant informatinon from B
    A = B.lmap
    n = A.n
    # For convenience
    x1 = @view x[1:n]
    x2 = @view x[n+1:2*n]
    x3 = @view x[2*n+1:3*n]
    x4 = @view x[3*n+1:4*n]
    #
    ABah = A.ABah
    diagNX = A.diagNX
    diagNY = A.diagNY
    diagNZ = A.diagNZ
    diagTX1 = A.diagTX1
    diagTY1 = A.diagTY1
    diagTZ1 = A.diagTZ1
    diagTX2 = A.diagTX2
    diagTY2 = A.diagTY2
    diagTZ2 = A.diagTZ2
    NoSlip1 = A.NoSlip1
    NoSlip2 = A.NoSlip2
    NullDivergence = A.NullDivergence

    # Multiplication with A' (conjugate transpose/adjoint)
    y[1:n]       .=   ABah'*x1            +           NoSlip1'*x3 + NoSlip2'*x4
    y[n+1:2*n]   .= diagNX'*x1            +           diagTX1'*x3 + diagTX2'*x4
    y[2*n+1:3*n] .= diagNY'*x1            +           diagTY1'*x3 + diagTY2'*x4
    y[3*n+1:4*n] .= diagNZ'*x1            +           diagTZ1'*x3 + diagTZ2'*x4
    y[n+1:end]   +=               NullDivergence'*x2     # Prettier to take x2 on its own

    return y
end
#==========================================================================================
                                    Utility function
==========================================================================================#
"""
    full(M::LossyBlockMatrixCompact)

Returns a matrix representation of the `LossyBlockMatrixCompact`.
"""
function full(M::LossyBlockMatrixCompact)
    N = M.N
    n = M.n

    A = zeros(eltype(M),N,N)
    A[1:n,       1:N] = [M.ABah                 M.diagNX    M.diagNY    M.diagNZ ]
    A[  n+1:2*n, 1:N] = [zeros(eltype(M),n,n)                    M.NullDivergence]
    A[2*n+1:3*n, 1:N] = [M.NoSlip1              M.diagTX1  M.diagTY1    M.diagTZ1]
    A[3*n+1:4*n, 1:N] = [M.NoSlip2              M.diagTX2  M.diagTY2    M.diagTZ2]

    return A
end
# Populating a the full matrix
function full(B::LinearMaps.AdjointMap{T,S}) where {T <: ComplexF64, S <: LossyBlockMatrixCompact}
    N = M.N
    n = M.n
    M = B.lmap

    A = zeros(eltype(M),N,N)
    A[1:N, 0*n+1:1*n] = [M.ABah';                 M.diagNX';    M.diagNY';    M.diagNZ']
    A[1:N, 1*n+1:2*n] = [zeros(eltype(M),n,n);                         M.NullDivergence']
    A[1:N, 2*n+1:3*n] = [M.NoSlip1';              M.diagTX1';  M.diagTY1';    M.diagTZ1']
    A[1:N, 3*n+1:4*n] = [M.NoSlip2';              M.diagTX2';  M.diagTY2';    M.diagTZ2']

    return A
end
# Creating the RHS from the "Shape-Function-Derivative"-Paper
function createRHS(mesh,u0)
    n = size(mesh.sources,2)
    vn0  = u0*mesh.normals[3,:]
    vt10 = u0*mesh.tangents[3,:]
    vt20 = u0*mesh.sangents[3,:]
    return [vn0; zeros(n); vt10; vt20]
end

function LossyBlockMatrixCompact(T::LossyBlockMatrix;blockoutput=false)
    normals   = T.normals
    tangents₁ = T.tangents₁
    tangents₂ = T.tangents₂
    Aₐ = T.Aₐ
    Bₐ = T.Bₐ
    ϕₐ = T.ϕₐ
    τₐ = T.τₐ
    Aᵥ = T.Aᵥ
    Bᵥ = T.Bᵥ
    Aₕ = T.Aₕ
    Bₕ = T.Bₕ
    ϕₕ = T.ϕₕ
    τₕ = T.τₕ
    Dt₁ = T.Dt₁
    Dt₂ = T.Dt₂
    nSource = T.n
    N = 4*nSource

    println("Reducing system")
    ABah   = (ϕₐ*(Bₐ\Aₐ) - (ϕₕ*τₐ/τₕ)*(Matrix(Bₕ)\Matrix(Aₕ)))
    diagNX = Diagonal(normals[1,:])
    diagNY = Diagonal(normals[2,:])
    diagNZ = Diagonal(normals[3,:])
    # Row2: NullDivergence
    println("Computing null divergence")
    ABᵥ = Matrix(Bᵥ)\Matrix(Aᵥ)
    ND = [(Diagonal(  normals[1,:])*ABᵥ) (Diagonal(  normals[2,:])*ABᵥ) (Diagonal(  normals[3,:])*ABᵥ)] +
         [(Diagonal(tangents₁[1,:])*Dt₁) (Diagonal(tangents₁[2,:])*Dt₁) (Diagonal(tangents₁[3,:])*Dt₁)] +
         [(Diagonal(tangents₂[1,:])*Dt₂) (Diagonal(tangents₂[2,:])*Dt₂) (Diagonal(tangents₂[3,:])*Dt₂)]
    # Row3
    println("Computing no-slip conditions")
    NoSlip1 = (ϕₐ - (ϕₕ*τₐ)/τₕ)*Dt₁
    diagTX1 = Diagonal(tangents₁[1,:])
    diagTY1 = Diagonal(tangents₁[2,:])
    diagTZ1 = Diagonal(tangents₁[3,:])
    # Row4
    NoSlip2 = (ϕₐ - (ϕₕ*τₐ)/τₕ)*Dt₂
    diagTX2 = Diagonal(tangents₂[1,:])
    diagTY2 = Diagonal(tangents₂[2,:])
    diagTZ2 = Diagonal(tangents₂[3,:])
    println("Creating LossyBlockMatrixCompact")
    LossyBlockCompact = LossyBlockMatrixCompact(nSource,N,
                                    ABah,   diagNX, diagNY, diagNZ,
                                                                ND,
                                    sparse(NoSlip1),diagTX1,diagTY1,diagTZ1,
                                    sparse(NoSlip2),diagTX2,diagTY2,diagTZ2)
    println("Returning output")
    if blockoutput
        return LossyBlockCompact
    else
        return full(LossyBlockCompact)
    end
end
