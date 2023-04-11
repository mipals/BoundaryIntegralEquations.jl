using IFGF
using StaticArrays
using LinearMaps
import BoundaryIntegralEquations.setup_fast_operator
const Point3D = SVector{3,Float64}
#==========================================================================================
                                    Helper functions
==========================================================================================#
struct HelmholtzSingleLayer <: AbstractMatrix{ComplexF64}
    X::Vector{Point3D}
    Y::Vector{Point3D}
    k::Float64
end
IFGF.wavenumber(B::HelmholtzSingleLayer) = B.k
# Abstract matrix interface
Base.size(B::HelmholtzSingleLayer) = length(B.X), length(B.Y)
Base.getindex(B::HelmholtzSingleLayer,i::Int,j::Int) = A(B.X[i],B.Y[j])
# functor interface
function (B::HelmholtzSingleLayer)(x,y)
    k = IFGF.wavenumber(B)
    d = norm(x-y)
    return exp(im*k*d)/(4π*d)
end
#==========================================================================================
                        Defining G-operator (single-layer potential)
==========================================================================================#
"""
    IFGFGOperator(k,G,C,nearfield_correction,coefficients)
    IFGFGOperator(mesh,k;tol=1e-4,n_gauss=3,nearfield=true,offset=0.2,depth=1)

A `LinearMap` that represents the BEM ``\\mathbf{G}`` matrix through the IFGF-matrix approach.
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

The remaining multiplication with the Green's functions is done utilizing the `IFGF.jl` package.
"""
struct IFGFGOperator{T} <: LinearMaps.LinearMap{T}
    k
    G
    C::AbstractMatrix{Float64}
    nearfield_correction::AbstractMatrix{T}
    coefficients::AbstractVecOrMat{T}
end

Base.size(A::IFGFGOperator) = (size(A.G,1),size(A.G,1))

function LinearMaps._unsafe_mul!(y, A::IFGFGOperator, x::AbstractVector)
    mul!(A.coefficients,A.C,x)
    mul!(y,A.G,A.coefficients);
    mul!(y,A.nearfield_correction,x,true,true);
    return y
end

function IFGFGOperator(mesh,k;tol=1e-3,n_gauss=3,nearfield=true,offset=0.2,depth=1)
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
    H = assemble_ifgf(HelmholtzSingleLayer(X,Y,k),X,Y;tol=tol);
    return IFGFGOperator(k,H,C_map,nearfield_correction,coefficients)
end
