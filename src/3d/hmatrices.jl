#==========================================================================================
                                    Helper functions
==========================================================================================#
"""
    create_single_layer_matrix(X, Y, k)

Creates a `KernelMatrix` representing the single layer potential.
"""
function create_single_layer_matrix(X, Y, k)
    f = (x, y) -> begin
        # EPS = 1e-8 # fudge factor to avoid division by zero
        d = norm(x - y)
        exp(im * k * d) / (4π * d)
    end
    return KernelMatrix(f, X, Y)
end
#==========================================================================================
                        Defining G-operator (single-layer potential)
==========================================================================================#
"""
    HGOperator(k,G,C,nearfield_correction,coefficients)
    HGOperator(mesh,k;tol=1e-4,n_gauss=3,nearfield=true,offset=0.2,depth=1)

A `LinearMap` that represents the BEM ``\\mathbf{G}`` matrix through the H-matrix approach.
This matrix has ``k``th row given by ``\\mathbf{z}=\\mathbf{z}_k`` in the following
```math
\\begin{aligned}
    \\left(\\int_{\\Gamma} G(\\mathbf{x},\\mathbf{z})\\mathbf{T}(\\mathbf{x}) \\mathrm{d}S_\\mathbf{x}\\right)\\mathbf{y}
    &\\approx \\left(\\sum_{e=1}^{N}\\left(\\sum_{i=1}^{Q}G(\\mathbf{x}^e(\\mathbf{u}_i),\\mathbf{z})\\text{jacobian}(\\mathbf{u}_i)w_i\\mathbf{T}^e(\\mathbf{u}_i)\\right)\\mathbf{L}^e\\right)\\mathbf{y}         \\newline
    &= \\left(\\sum_{j=1}^{NQ}G(\\mathbf{x}_j,\\mathbf{z})\\underbrace{\\text{jacobian}(\\mathbf{u}_j)w_j\\mathbf{T}^{e(j)}(\\mathbf{u}_j)\\mathbf{L}^{e(j)}}_{j\\text{th row of } \\mathbf{C}}\\right)\\mathbf{y}   \\newline
    &=
    \\begin{bmatrix}
        G(\\mathbf{x}_1,\\mathbf{z}) & G(\\mathbf{x}_2,\\mathbf{z}) & \\dots & G(\\mathbf{x}_{NQ},\\mathbf{z})
    \\end{bmatrix}
    \\mathbf{C}\\mathbf{y},
\\end{aligned}
```
where the subscript ``j`` refers to an ordering of the collection of Gaussian points from all elements and ``e(j)`` is a function that returns the element number that Gaussian point ``j`` is located on. Note that ``\\mathbf{C}`` is the same same for all ``k``'s.

The remaining multiplication with the Green's functions is done utilizing the `HMatrices.jl` package.
"""
struct HGOperator{T} <: LinearMaps.LinearMap{T}
    k
    G
    C::AbstractMatrix{Float64}
    nearfield_correction::AbstractMatrix{T}
    coefficients::AbstractVecOrMat{T}
end
Base.size(A::HGOperator) = (size(A.G,1),size(A.G,1))
function LinearMaps._unsafe_mul!(y, A::HGOperator, x::AbstractVector)
    mul!(A.coefficients,A.C,x)
    mul!(y,A.G,A.coefficients);
    mul!(y,A.nearfield_correction,x,true,true);
    return y
end
function HGOperator(mesh,k;tol=1e-4,n_gauss=3,nearfield=true,offset=0.2,depth=1)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Setup operator
    sources,_,C_map,nearfield_correction = setup_fast_operator(mesh,zk,n_gauss,nearfield,
                                                                    offset,depth)
    # Extracting targets
    targets = mesh.sources
    # Creating temporary array
    coefficients = zeros(eltype(zk),size(sources,2))
    # Creating HMatrix
    X = [Point3D(x) for x in eachcol(targets)]
    Y = [Point3D(x) for x in eachcol(sources)]
    # Assembly of H-matrix
    H = assemble_hmat(create_single_layer_matrix(X, Y, k);rtol=tol)
    return HGOperator(k,H,C_map,nearfield_correction,coefficients)
end

#==========================================================================================
                            Defining H-operator (double-layer)
==========================================================================================#
"""
    HelmholtzDoubleLayer

```math
    \\frac{e^{ikr_j}}{4\\pi r_j^3}(ikr_j - 1) (x_j - y)\\cdot n
```
"""
struct HelmholtzDoubleLayer <: AbstractMatrix{ComplexF64}
    X::Vector{Point3D}
    Y::Vector{Point3D}
    NY::Vector{Point3D} # normals at Y coordinate
    k::Float64
end
function Base.getindex(K::HelmholtzDoubleLayer,i::Int,j::Int)
    # r = K.X[i] - K.Y[j]
    r = K.Y[i] - K.X[j]
    d = norm(r)
    return exp(im*K.k*d)/(4π*d^3) * (im*K.k*d - 1) * dot(r, K.NY[j])
end
Base.size(K::HelmholtzDoubleLayer) = length(K.X), length(K.Y)
"""
    HHOperator(k,H,C,nearfield_correction,coefficients)
    HHOperator(mesh,k;tol=1e-4,n_gauss=3,nearfield=true,offset=0.2,depth=1)

A `LinearMap` that represents the BEM ``\\mathbf{H}`` matrix through the H-matrix approach.
This matrix has ``k``th row given by ``\\mathbf{z}=\\mathbf{z}_k`` in the following
```math
\\begin{aligned}
    \\left(\\int_{\\Gamma}\\frac{\\partial G(\\mathbf{x}, \\mathbf{z})}{\\partial\\mathbf{n}(\\mathbf{x})}\\mathbf{T}(\\mathbf{x}) \\mathrm{d}S_{\\mathbf{x}}\\right)\\mathbf{y}
    &\\approx
    \\left(\\sum_{e=1}^{N}\\left(\\sum_{i=1}^{Q}\\frac{\\partial G(\\mathbf{x}^e(\\mathbf{u}_i), \\mathbf{z})}{\\partial\\mathbf{n}(\\mathbf{x})}\\text{jacobian}(\\mathbf{u}_i)w_i\\mathbf{T}^e(\\mathbf{u}_i)\\right)\\mathbf{L}^e\\right)\\mathbf{y}\\newline
    &= \\left(\\sum_{j=1}^{NQ}\\nabla G(\\mathbf{x}_j, \\mathbf{z})\\cdot \\mathbf{n}(\\mathbf{x}_j)\\underbrace{\\text{jacobian}(\\mathbf{u}_j)w_j\\mathbf{T}^{e(j)}(\\mathbf{u}_j)\\mathbf{L}^{e(j)}}_{j\\text{th row of }\\mathbf{C}}\\right)\\mathbf{y}\\newline
    &=
    \\begin{bmatrix}
        \\nabla G(\\mathbf{x}_1, \\mathbf{z}) \\cdot \\mathbf{n}(\\mathbf{x}_1) &
        \\dots &
        \\nabla G(\\mathbf{x}_{NQ}, \\mathbf{z})\\cdot \\mathbf{n}(\\mathbf{x}_{NQ})
    \\end{bmatrix}
    \\mathbf{C}\\mathbf{y},
\\end{aligned}
```
where ``\\mathbf{C}`` is coefficient mapping (the same for all ``k``).
The remaining multiplication with the Green's functions is performed utilizing the `HMatrices.jl` package.
"""
struct HHOperator{T} <: LinearMaps.LinearMap{T}
    k
    H
    C::AbstractMatrix{Float64}
    nearfield_correction::AbstractMatrix{T}
    coefficients::AbstractVecOrMat{T}
end
Base.size(H::HHOperator) = (size(H.H,1),size(H.H,1))
function LinearMaps._unsafe_mul!(y, A::HHOperator, x::AbstractVector)
    mul!(A.coefficients,A.C,x)
    mul!(y,A.H,A.coefficients);
    mul!(y,A.nearfield_correction,x,true,true);
    return y
end
function HHOperator(mesh,k;tol=1e-4,n_gauss=3,nearfield=true,offset=0.2,depth=1)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Setup operator
    sources,normals,C_map,nearfield_correction = setup_fast_operator(mesh,zk,n_gauss,
                                                nearfield,offset,depth;single_layer=false)
    # Extracting targets
    targets = mesh.sources
    # Creating temporary array
    coefficients = zeros(eltype(zk),size(sources,2))
    # Creating HMatrix
    X  = [Point3D(x) for x in eachcol(targets)]
    Y  = [Point3D(y) for y in eachcol(sources)]
    NY = [Point3D(n) for n in eachcol(normals)]
    # Creating cluster trees
    Xclt = ClusterTree(X)
    Yclt = ClusterTree(Y)
    # Assembling H-matrix representing the double-layer potential
    H = assemble_hmat(HelmholtzDoubleLayer(X,Y,NY,k),Xclt,Yclt;comp=PartialACA(;rtol=tol))
    return HHOperator(k,H,C_map,nearfield_correction + 0.5I,coefficients)
end
