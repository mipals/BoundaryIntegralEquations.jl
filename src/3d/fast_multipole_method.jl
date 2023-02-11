#==========================================================================================
                            Partial Assembly
==========================================================================================#
"""
    partial_assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                    fOn=true,gOn=true,n=3,progress=false,depth=1)

Assembles only local contributions in the BEM computations.
This is used for singularity extraction when using the Fast Mulitpole Method (FMM) for BEM.
"""
function partial_assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                    fOn=true,gOn=true,n=3,progress=false,depth=1)
    topology    = get_topology(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    coordinates = convert.(eltype(shape_function),get_coordinates(mesh))
    physics_topology = mesh.physics_topology
    physics_function = mesh.physics_function
    #======================================================================================
        Introducing three elements: (The hope here is to compute singular integrals)
        The numbers corresponds to the corner for which the GP-points are clustered
        As such choose the clustering that are closest to the source point
     ————————————————————————————————————  Grid  ————————————————————————————————————————
                                           3
                                           | \
                                           1 - 2
    ======================================================================================#
    shape_function1   = deepcopy(shape_function)
    physics_function1 = deepcopy(physics_function)
    set_interpolation_nodes!(shape_function1,gauss_points_triangle(n)...)
    copy_interpolation_nodes!(physics_function1,shape_function1)

    # Copying interpolation of physics functions1
    physics_interpolation1 = physics_function1.interpolation

    # Connections
    element_connections, source_connections = connected_sources(mesh,depth)
    lengths = length.(source_connections)
    dict = [Dict(zip(source_connections[i],1:lengths[i])) for i = 1:length(lengths)]
    idx = [0; cumsum(lengths)]
    # Preallocation of return values
    F = zeros(ComplexF64, idx[end])
    G = zeros(ComplexF64, idx[end])

    n_physics_functions = number_of_shape_functions(physics_function)
    n_shape_functions   = number_of_shape_functions(shape_function)

    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        integrand  = zeros(ComplexF64,n)
        r          = zeros(Float64,n)
        jacobian   = similar(r)
        normals    = zeros(Float64,3,n)
        tangents   = similar(normals)
        sangents   = similar(normals)
        interpolation = similar(normals)

        di = dict[source_node]
        physics_nodes = zeros(Int64, n_physics_functions)
        element_coordinates = zeros(3,n_shape_functions)
        # element_coordinates = zeros(3,n_physics_functions)
        @inbounds for element ∈ element_connections[source_node]
            # Access element topology and coordinates
            element_coordinates .= coordinates[:,topology[:,element]]
            find_physics_nodes!(physics_nodes,idx[source_node],di,physics_topology[:,element])
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[physics_nodes]
            submatrixG = @view G[physics_nodes]
            # Use quadrature point clustered around the closest vertex
            sparse_computing_integrals!(physics_interpolation1,shape_function1,
                                        normals, tangents,sangents,
                                        interpolation,jacobian,r,integrand,
                                        element_coordinates,
                                        fOn,gOn,submatrixF,submatrixG,k,source)
        end
        if progress; next!(prog); end # For the progress meter
    end

    # Computing row and column indices for the assembled contributions
    I = create_row_indices(lengths,idx[end])
    J = vcat(source_connections...)

    return sparse(I,J,F), sparse(I,J,G)

end

#==========================================================================================
                            Helper Functions
==========================================================================================#
"""
    unroll_interpolations(interpolations)

Merges the element interpolations `interpolate_elements()` into vectors.
"""
function unroll_interpolations(interpolations)
    # Number of elements
    n_elements = length(interpolations)
    # Number of gauss points pr. element
    Nweights = length(interpolations[1].jacobian_mul_weights)
    # Allocating output
    weights = zeros(n_elements*Nweights)
    interps = zeros(3,n_elements*Nweights)
    normals = zeros(3,n_elements*Nweights)
    # Extracting total number of weights, interpolations and normals
    for i = 1:n_elements
        weights[  (i-1)*Nweights + 1:i*Nweights] = interpolations[i].jacobian_mul_weights
        interps[:,(i-1)*Nweights + 1:i*Nweights] = interpolations[i].interpolation
        normals[:,(i-1)*Nweights + 1:i*Nweights] = interpolations[i].normals
    end
    return interps, weights, normals
end

"""
    create_coefficient_map(weights,physics_topology,physics_function,n_gauss)

Returns the mappgin from nodal values to coefficients for the FMM.
This mapping is represented by a (sparse) matrix ``C`` with rows given as
```math
    \\underbrace{\\text{jacobian}(\\mathbf{u}_j)w_j\\mathbf{T}^{e(j)}(\\mathbf{u}_j)\\mathbf{L}^{e(j)}}_{j\\text{th row of } \\mathbf{C}},
```
where each of the ``j`` rows corresponds to the Gaussian points from all elements and ``e(j)`` is a function that returns the element number that Gaussian point ``j`` is located on.
"""
function create_coefficient_map(weights,physics_topology,physics_function,n_gauss)
    # Computing row and column indicies
    I = inverse_rle(collect(1:length(weights)),number_of_shape_functions(physics_function)*ones(Int,length(weights)))
    J = repeat(physics_topology,n_gauss,1)[:]
    # Extracting number of shape functions and number of elements
    n_shape = number_of_shape_functions(physics_function)
    n_elements = size(physics_topology,2)
    # Pre-Allocation of values
    V_mat = zeros(n_shape,n_elements*n_gauss)
    # Looping over each row of the coefficient mapping
    for i = 1:n_elements*n_gauss
        V_mat[:,i] = weights[i]*physics_function.interpolation[:,mod(i-1,n_gauss)+1]
    end
    return sparse(I,J,V_mat[:])
end

"""
    scale_columns!(weights,normals,scalings)

Saves the ``i``th column of `normals` scaled by the ``i``th value of `scalings` in `weights`.
"""
function scale_columns!(weights,normals,scalings)
    for i = eachindex(scalings)
        weights[1,i] = normals[1,i]*scalings[i]
        weights[2,i] = normals[2,i]*scalings[i]
        weights[3,i] = normals[3,i]*scalings[i]
    end
    return weights
end
#==========================================================================================
                        Defining G-operator (single-layer potential)
==========================================================================================#
"""
    FMMGOperator

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
where the subscript ``j`` refers to an ordering of the collection of Gaussian points from all elements and ``e(j)`` is a function that returns the element number that Gaussian point ``j`` is located on. Note that ``\\mathbf{C}`` is the same same for all ``k``'s. The remaining multiplication with the Green's functions is done utilizing the Flatiron Institute Fast Multipole libraries.
"""
struct FMMGOperator{T} <: LinearMaps.LinearMap{T}
    # Dimensions of operator
    n::Int64                            # Number of sources
    m::Int64                            # Number of targets
    # Physical Quantities
    k::T                                # Wavenumber
    eps::Float64                        # Precision
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
Base.size(A::FMMGOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMGOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Compute coefficients
    mul!(A.coefficients,A.C,x)
    # Computing the FMM sum
    vals = hfmm3d(A.eps,A.k,A.sources,charges=A.coefficients,targets=A.targets,pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # Note that the uses a Greens function that does not divide by 4π.
    y .= vals.pottarg/4π + A.nearfield_correction*x
end
function FMMGOperator(mesh,k;eps=1e-6,n_gauss=3,nearfield=true,offset=0.2,depth=1)
    # Setting up elements with the correct number of Gaussian points
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n_gauss)...)
    copy_interpolation_nodes!(physics_function,shape_function)
    # Interpolating on the mesh using the previous computed shape/physics functions
    interpolations    = interpolate_elements(mesh,shape_function)
    sources,weights,_ = unroll_interpolations(interpolations)
    # Creating map from nodes to FMM coefficients
    C_map = create_coefficient_map(weights,mesh.physics_topology,physics_function,n_gauss)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Extracting mesh information
    targets = mesh.sources
    N = size(targets,2)
    M = size(sources,2)
    # Computing near-field correction
    if nearfield
        _,C = partial_assemble_parallel!(mesh,zk,targets,shape_function;fOn=false,depth=depth)
        _,S = assemble_parallel!(mesh,zk,targets;
                                sparse=true,depth=depth,fOn=false,progress=false,offset=offset);
        nearfield_correction = - C + S
    else
        nearfield_correction = spzeros(ComplexF64,N,N)
    end
    # Creating temporary array
    coefficients = zeros(eltype(zk),M)
    return FMMGOperator(N,M,zk,eps,targets,sources,C_map,coefficients,nearfield_correction)
end
#==========================================================================================
                            Defining H-operator (double-layer)
==========================================================================================#
"""
    FMMHOperator

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
where ``\\mathbf{C}`` is coefficient mapping (the same for all ``k``). The remaining multiplication with the Green's functions is performed utilizing the Flatiron Institute Fast Multipole libraries.
"""
struct FMMHOperator{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    m::Int64                            # Number of gauss nodes
    # Physical Quantities
    k::T                                # Wavenumber
    eps::Float64                        # Precision
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
Base.size(A::FMMHOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMHOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Map the mesh nodes to the gauss nodes
    mul!(A.coefficients,A.C,x)
    # Scale "dipoles"
    scale_columns!(A.dipvecs,A.normals,A.coefficients)
    # Computing the FMM sum
    vals = hfmm3d(A.eps,A.k,A.sources,targets=A.targets,dipvecs=A.dipvecs,pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # Note that FMM3D uses a Greens function that does not divide by 4π.
    y .= vals.pottarg/4π + A.nearfield_correction*x
end

function FMMHOperator(mesh,k;n_gauss=3,eps=1e-6,nearfield=true,offset=0.2,depth=1,
                                integral_free_term = [],progress=false)
    # Creating physics and geometry for the FMM operator
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n_gauss)...)
    copy_interpolation_nodes!(physics_function,shape_function)
    # Making sure that the wavenumber is complex (required for the FMM3D library)
    zk = Complex(k)
    # Interpolating on each element
    interpolations = interpolate_elements(mesh,shape_function)
    # Unrolling element interpolations
    sources,weights,normals = unroll_interpolations(interpolations)
    # Create coefficient map
    C_map = create_coefficient_map(weights,mesh.physics_topology,physics_function,n_gauss)
    # Set targets equal to sources
    targets = mesh.sources
    # Getting size of FMM matrix
    N = size(targets,2)
    M = size(sources,2)
    # Computing dipole directional strenght. Making sure they're complex (required by FMM3D)
    # Allocating arrays for intermediate computations
    coefficients = zeros(eltype(zk),length(weights))
    dipvecs = zeros(eltype(zk),3,length(weights))
    # If not set by the user set the integral free term equal to a half
    if isempty(integral_free_term)
        nearfield_correction = Diagonal(ones(eltype(zk),N)/2)
    elseif length(integral_free_term) == M
        nearfield_correction = Diagonal(integral_free_term)
    end
    # Computing near-field correction
    if nearfield
        C,_ = partial_assemble_parallel!(mesh,zk,targets,shape_function;gOn=false,depth=depth)
        S,_ = assemble_parallel!(mesh,zk,targets;sparse=true,depth=depth,gOn=false,offset=offset,progress=progress);
        nearfield_correction += - C + S
    end
    return FMMHOperator(N,M,zk,eps,targets,sources,normals,C_map,coefficients,dipvecs,nearfield_correction)
end
