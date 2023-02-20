# Idea: Make a more general interface where both FMM and H-matrices can be specified to
# handle the multiplication with the structured matrix with $G[i,j] = K(x_i,x_j)$
# Reason: Currently the FMM- and H-operators is essentially a copy of the same code.
# would be nice to combine the code.
"""
    FastGOperator(k,tol,targets,sources,C,coefficients,nearfield_correction,G)
    FastGOperator(mesh,k;tol=1e-6,n_gauss=3,nearfield=true,offset=0.2,depth=1)

A `LinearMap` that represents the BEM ``\\mathbf{G}`` matrix through the FMM or H-matrix.
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

The remaining multiplication with the Green's functions is done utilizing either FMM or an H-matrix.
"""
struct FastGOperator{T} <: LinearMaps.LinearMap{T}
    # Physical Quantities
    k::T                                # Wavenumber
    tol::Float64                        # Precision
    # FMM setup
    targets::AbstractMatrix{Float64}    # FMM-targets (size = 3,n)
    sources::AbstractMatrix{Float64}    # FMM-sources (size = 3,m)
    # Mapping from global to coefficients
    C::AbstractMatrix{Float64}
    # For storing coefficients
    coefficients::AbstractVecOrMat{T}
    # Correting near field computations
    nearfield_correction::AbstractMatrix{T}
    # G-operator
    G
end



"""
    FastHOperator(k,tol,targets,sources,normals,C,coefficients,dipvecs,nearfield_correction)
    FastHOperator(mesh,k;n_gauss=3,tol=1e-6,nearfield=true,offset=0.2,depth=1,)

A `LinearMap` that represents the BEM ``\\mathbf{H}`` matrix through the FMM or a H-matarix.
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

The remaining multiplication with the Green's functions is performed utilizing either the FMM or a H-matrix.
"""
struct FastHOperator{T} <: LinearMaps.LinearMap{T}
    # Physical Quantities
    k::T                                # Wavenumber
    tol::Float64                        # Precision
    # Acoustic Matrices
    targets::AbstractMatrix{Float64}    # FMM-targets (size = 3,n)
    sources::AbstractMatrix{Float64}    # FMM-sources (size = 3,m)
    normals::AbstractMatrix{Float64}    # FMM-directions
    #
    C::AbstractMatrix{Float64}          # Coefficient mapping
    coefficients::AbstractVecOrMat{T}   # FMM-Coefficients
    dipvecs::AbstractVecOrMat{T}        # FMM-dipole vectors
    # Correting near field computations
    nearfield_correction::AbstractMatrix{T}
    # H-operator
    H
end
