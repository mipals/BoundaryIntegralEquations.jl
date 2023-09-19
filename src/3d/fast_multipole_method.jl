#==========================================================================================
                            Partial Assembly
==========================================================================================#
"""
    partial_assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                    fOn=true,gOn=true,n=3,progress=false,depth=1)

Assembles only local contributions in the BEM computations.
This is used for singularity extraction when using the Fast Mulitpole Method (FMM) for BEM.
"""
function partial_assemble_parallel!(mesh::Mesh3d, k, in_sources, shape_function::Triangular;
                                    fOn=true, gOn=true, n=3, progress=false, depth=1)
    # Extracting mesh information
    topology = get_topology(mesh)
    sources = convert.(eltype(shape_function), in_sources)
    n_sources = size(in_sources, 2)
    coordinates = convert.(eltype(shape_function), get_coordinates(mesh))
    physics_topology = mesh.physics_topology
    physics_function = mesh.physics_function
    n_physics_functions = number_of_shape_functions(physics_function)
    n_shape_functions = number_of_shape_functions(shape_function)
    # Creating copies of mesh shape and physics/interpolation functions
    shape_function1 = deepcopy(shape_function)
    physics_function1 = deepcopy(physics_function)
    # Setting interpolation to be at `n` gaussian points defined by `gauss_points_triangle`.
    set_interpolation_nodes!(shape_function1, gauss_points_triangle(n)...)
    set_interpolation_nodes!(physics_function1, gauss_points_triangle(n)...)
    # Copying interpolation of physics functions1
    physics_interpolation1 = physics_function1.interpolation
    # Creating dictionary of elements connected to the source/collocation point
    element_connections, source_connections = connected_sources(mesh, depth)
    lengths = length.(source_connections)
    dict = [Dict(zip(source_connections[i], 1:lengths[i])) for i in 1:length(lengths)]
    idx = [0; cumsum(lengths)]
    # Preallocation of return values
    F = zeros(ComplexF64, idx[end])
    G = zeros(ComplexF64, idx[end])
    # Assembly loop
    if progress
        prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50)
    end
    @inbounds @threads for source_node in 1:n_sources
        # Pre-allocate matrices
        integrand = zeros(ComplexF64, n)
        r = zeros(Float64, n)
        jacobian = similar(r)
        normals = zeros(Float64, 3, n)
        tangents = similar(normals)
        sangents = similar(normals)
        interpolation = similar(normals)
        element_coordinates = zeros(3, n_shape_functions)
        physics_nodes = zeros(Int64, n_physics_functions)
        # Accessing correct "source"
        source = sources[:, source_node]
        di = dict[source_node]
        # Only loop over elements close to sourge point
        @inbounds for element in element_connections[source_node]
            # Access element topology and coordinates
            element_coordinates .= coordinates[:, topology[:, element]]
            find_physics_nodes!(physics_nodes, idx[source_node], di,
                                physics_topology[:, element])
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[physics_nodes]
            submatrixG = @view G[physics_nodes]
            # Integrate using "physics_interpolation1"
            sparse_computing_integrals!(physics_interpolation1, shape_function1,
                                        normals, tangents, sangents,
                                        interpolation, jacobian, r, integrand,
                                        element_coordinates,
                                        fOn, gOn, submatrixF, submatrixG, k, source)
        end
        if progress
            next!(prog)
        end
    end
    # Computing row and column indices for the assembled contributions
    I = inverse_rle(1:n_sources, lengths)
    J = vcat(source_connections...)
    # Returning sparse matrices corresponding to the "wrong" integration done by the FMM.
    return sparse(I, J, F), sparse(I, J, G)
end

#==========================================================================================
                            Helper Functions
==========================================================================================#
"""
    unroll_interpolations(interpolations)

Merges the element interpolations from `interpolate_elements` into vectors.
"""
function unroll_interpolations(interpolations)
    # Number of elements
    n_elements = length(interpolations)
    # Number of gauss points pr. element
    Nweights = length(interpolations[1].jacobian_mul_weights)
    # Allocating output
    weights = zeros(n_elements * Nweights)
    interps = zeros(3, n_elements * Nweights)
    normals = zeros(3, n_elements * Nweights)
    # Extracting total number of weights, interpolations and normals
    for i in 1:n_elements
        weights[((i - 1) * Nweights + 1):(i * Nweights)] = interpolations[i].jacobian_mul_weights
        interps[:, ((i - 1) * Nweights + 1):(i * Nweights)] = interpolations[i].interpolation
        normals[:, ((i - 1) * Nweights + 1):(i * Nweights)] = interpolations[i].normals
    end
    return interps, weights, normals
end

"""
    create_coefficient_map(weights,physics_topology,physics_function,n_gauss)

Returns the mapping from nodal values to coefficients for the FMM/H-matrix.

This mapping is represented by a (sparse) matrix ``C`` with rows given as
```math
    \\underbrace{\\text{jacobian}(\\mathbf{u}_j)w_j\\mathbf{T}^{e(j)}(\\mathbf{u}_j)\\mathbf{L}^{e(j)}}_{j\\text{th row of } \\mathbf{C}},
```
where each of the ``j`` rows corresponds to the Gaussian points from all elements and ``e(j)`` is a function that returns the element number that Gaussian point ``j`` is located on.
"""
function create_coefficient_map(weights, physics_topology, physics_function, n_gauss)
    # Extracting number of shape functions and number of elements
    n_shape = number_of_shape_functions(physics_function)
    n_weights = length(weights)
    # Pre-Allocation of values
    V_mat = zeros(n_shape, n_weights)
    # Looping over each row. Note that weights[i] is equal to jacobian(u_j)*w_j.
    for i in 1:n_weights
        V_mat[:, i] = weights[i] *
                      physics_function.interpolation[:, mod(i - 1, n_gauss) + 1]
    end
    # Each row contains `n_shape` of non-zeros pr. row
    I = inverse_rle(collect(1:n_weights), n_shape * ones(Int, n_weights))
    # Each element topology is repeated `n_gauss` times and vectorized.
    J = repeat(physics_topology, n_gauss, 1)[:]
    # Return the sparse matrix with the rows vectorized
    return sparse(I, J, V_mat[:])
end

"""
    scale_columns!(dipvecs,normals,scalings)

Saves the ``i``th column of `normals` scaled by the ``i``th value of `scalings` in `dipvecs`.

Used when creating the "dipole-vectors" required for using the FMM.
"""
function scale_columns!(dipvecs, normals, scalings)
    @inbounds for i in eachindex(scalings)
        dipvecs[1, i] = normals[1, i] * scalings[i]
        dipvecs[2, i] = normals[2, i] * scalings[i]
        dipvecs[3, i] = normals[3, i] * scalings[i]
    end
    return dipvecs
end

"""
    setup_fast_operator(mesh,zk,n_gauss,nearfield,offset,depth;single_layer=true)

Computes `sources`, `normals`, `C_map` and the `nearfield_correction` required when setting up a fast operator.
"""
function setup_fast_operator(mesh, zk, n_gauss, nearfield, offset, depth; single_layer=true)
    # Setting up elements with the correct number of Gaussian points
    shape_function = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function, gauss_points_triangle(n_gauss)...)
    copy_interpolation_nodes!(physics_function, shape_function)
    # Interpolating on the mesh using the previous computed shape/physics functions
    interpolations = interpolate_elements(mesh, shape_function)
    sources, weights, normals = unroll_interpolations(interpolations)
    C_map = create_coefficient_map(weights, mesh.physics_topology, physics_function,
                                   n_gauss)
    # Extracting mesh information
    targets = mesh.sources
    N = size(targets, 2)
    # Computing near-field correction
    if nearfield && single_layer
        _, C = partial_assemble_parallel!(mesh, zk, targets, shape_function; fOn=false,
                                          depth=depth)
        _, S = assemble_parallel!(mesh, zk, targets; sparse=true, depth=depth, fOn=false,
                                  progress=false, offset=offset)
        nearfield_correction = -C + S
    elseif nearfield && !single_layer
        C, _ = partial_assemble_parallel!(mesh, zk, targets, shape_function; gOn=false,
                                          depth=depth)
        S, _ = assemble_parallel!(mesh, zk, targets; sparse=true, depth=depth, gOn=false,
                                  offset=offset, progress=false)
        nearfield_correction = -C + S
    else
        nearfield_correction = spzeros(ComplexF64, N, N)
    end
    return sources, normals, C_map, nearfield_correction
end
#==========================================================================================
                        Defining G-operator (single-layer potential)
==========================================================================================#
"""
    FMMGOperator(k,tol,targets,sources,C,coefficients,nearfield_correction)
    FMMGOperator(mesh,k;tol=1e-6,n_gauss=3,nearfield=true,offset=0.2,depth=1)

A `LinearMap` that represents the BEM ``\\mathbf{G}`` matrix through the FMM.
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

The remaining multiplication with the Green's functions is done utilizing the Flatiron Institute Fast Multipole libraries.
"""
struct FMMGOperator{T} <: LinearMaps.LinearMap{T}
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
end

# Overloading size. Used to check dimension of inputs
Base.size(A::FMMGOperator) = (size(A.targets, 2), size(A.targets, 2))

function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMGOperator{T},
                            x::AbstractVector) where {T<:ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Compute coefficients
    mul!(A.coefficients, A.C, x)
    # Computing the FMM sum
    vals = hfmm3d(A.tol, A.k, A.sources; charges=A.coefficients, targets=A.targets, pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # The FMM3D library does not divide by 4π so we do it manually
    return y .= vals.pottarg / 4π + A.nearfield_correction * x
end

function FMMGOperator(mesh, k; tol=1e-6, n_gauss=3, nearfield=true, offset=0.2, depth=1)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Setup operator
    sources, _, C_map, nearfield_correction = setup_fast_operator(mesh, zk, n_gauss,
                                                                  nearfield, offset, depth)
    # Extracting targets
    targets = mesh.sources
    # Allocating array for intermediate computations
    coefficients = zeros(eltype(zk), size(sources, 2))
    return FMMGOperator(zk, tol, targets, sources, C_map, coefficients,
                        nearfield_correction)
end

"""
    evaluate_targets(A::FMMGOperator,x,targets)

Evaluates the G-operator with coefficients `x` at the `targets` positions.
The wavenumber used is defined by the `FMMGOperaetor`.
"""
function evaluate_targets(A::FMMGOperator, x, targets)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients, A.C, x)
    # Computing the FMM sum
    return hfmm3d(A.tol, A.k, A.sources; charges=A.coefficients, targets=targets, pgt=1).pottarg / 4π
end

"""
    evaluate_targets(A::FMMGOperator,k,x,targets)

Evaluates the G-operator with coefficients `x` at the `targets` positions.
The wavenumber used, `k`, is supplied by the user.
"""
function evaluate_targets(A::FMMGOperator, k, x, targets)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients, A.C, x)
    # Computing the FMM sum
    return hfmm3d(A.tol, k, A.sources; charges=A.coefficients, targets=targets, pgt=1).pottarg / 4π
end

#==========================================================================================
                            Defining H-operator (double-layer)
==========================================================================================#
"""
    FMMFOperator(k,tol,targets,sources,normals,C,coefficients,dipvecs,nearfield_correction)
    FMMFOperator(mesh,k;n_gauss=3,tol=1e-6,nearfield=true,offset=0.2,depth=1,)

A `LinearMap` that represents the BEM ``\\mathbf{H}`` matrix through the FMM.
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

The remaining multiplication with the Green's functions is performed utilizing the Flatiron Institute Fast Multipole libraries.
"""
struct FMMFOperator{T} <: LinearMaps.LinearMap{T}
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
end

Base.size(A::FMMFOperator) = (size(A.targets, 2), size(A.targets, 2))

function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMFOperator{T},
                            x::AbstractVector) where {T<:ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients, A.C, x)
    # Scale "dipoles"
    scale_columns!(A.dipvecs, A.normals, A.coefficients)
    # Computing the FMM sum
    vals = hfmm3d(A.tol, A.k, A.sources; targets=A.targets, dipvecs=A.dipvecs, pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # The FMM3D library does not divide by 4π so we do it manually
    return y .= vals.pottarg / 4π + A.nearfield_correction * x
end

function FMMFOperator(mesh, k; n_gauss=3, tol=1e-6, nearfield=true, offset=0.2, depth=1)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Setup operator
    sources, normals, C_map, nearfield_correction = setup_fast_operator(mesh, zk, n_gauss,
                                                                        nearfield, offset,
                                                                        depth;
                                                                        single_layer=false)
    # Creating targets
    targets = mesh.sources
    # Allocating arrays for intermediate computations
    coefficients = zeros(eltype(zk), size(sources, 2))
    dipvecs = zeros(eltype(zk), 3, size(sources, 2))
    return FMMFOperator(zk, tol, targets, sources, normals, C_map, coefficients, dipvecs,
                        nearfield_correction)
end

### For evaluating the FMM at new target positions
"""
    evaluate_targets(A::FMMFOperator,x,targets)

Evaluates the G-operator with coefficients `x` at the `targets` positions.
The wavenumber used is defined by the `FMMGOperaetor`.
"""
function evaluate_targets(A::FMMFOperator, x, targets)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients, A.C, x)
    # Scale "dipoles"
    scale_columns!(A.dipvecs, A.normals, A.coefficients)
    # Computing the FMM sum
    return hfmm3d(A.tol, A.k, A.sources; targets=targets, dipvecs=A.dipvecs, pgt=1).pottarg / 4π
end

"""
    evaluate_targets(A::FMMFOperator,k,x,targets)

Evaluates the F-operator with coefficients `x` at the `targets` positions.
The wavenumber used, `k`, is supplied by the user.
"""
function evaluate_targets(A::FMMFOperator, k, x, targets)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients, A.C, x)
    # Scale "dipoles"
    scale_columns!(A.dipvecs, A.normals, A.coefficients)
    # Computing the FMM sum
    return hfmm3d(A.tol, k, A.sources; targets=targets, dipvecs=A.dipvecs, pgt=1).pottarg / 4π
end
