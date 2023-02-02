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
    nodes_to_gauss!(gauss_points,elmement_interpolation,physics_topology,x)

Computes the global coordinates of all the Gauss points.
"""
function nodes_to_gauss!(gauss_points,elmement_interpolation,physics_topology,x)
    # Getting number of interpolations pr. element
    n_interps  = size(elmement_interpolation,2)
    # Copying the element interpolation
    elm_interp = copy(elmement_interpolation')
    # Compute the Gaussian points on each element
    @inbounds for i = axes(physics_topology,2)
        gauss_points[(i-1)*n_interps+1:i*n_interps] = elm_interp*x[physics_topology[:,i]]
    end
    return gauss_points
end
#==========================================================================================
                        Defining G-operator (single-layer potential)
==========================================================================================#
"""
    FMMGOperator

A `LinearMap` that represents the BEM G matrix through the FMM.
"""
struct FMMGOperator{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    m::Int64                            # Number of gauss nodes
    # Physical Quantities
    k::T                                # Wavenumber
    eps::Float64                        # Precision
    # Acoustic Matrices
    targets::AbstractMatrix{Float64}    # FMM-targets (size = 3,n)
    sources::AbstractMatrix{Float64}    # FMM-sources (size = 3,m)
    weights::AbstractVecOrMat{T}        # Nodal weights (size = m)
    # For mapping nodal-values to gauss-nodes
    element_interpolation::AbstractMatrix{Float64} # To go from p to
    physics_topology::AbstractMatrix{Int64}
    # Correting near field computations
    nearfield_correction::AbstractMatrix{T}
    # To avoid repeating allocation when multiplying
    tmp_weights::AbstractVecOrMat{T}
end
# Overloading size. Used to check dimension of inputs
Base.size(A::FMMGOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMGOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Map the mesh nodes to the gauss nodes
    nodes_to_gauss!(A.tmp_weights,A.element_interpolation,A.physics_topology,x)
    # Multiplying by the pre-computed weights
    integrand_mul!(A.tmp_weights,A.weights)
    # Computing the FMM sum
    vals = hfmm3d(A.eps,A.k,A.sources,charges=A.tmp_weights,targets=A.targets,pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # Note that the uses a Greens function that does not divide by 4π.
    y .= vals.pottarg/4π + A.nearfield_correction*x
end
function FMMGOperator(eps,k,targets,sources,weights,elm_interp,physics_topology)
    # Making sure the wavenumber is complex (required for the FMM3D library)
    zk = Complex(k)
    # Getting size of the FMM operator.
    n = size(targets,2)
    m = size(sources,2)
    # Allocating array for intermediate computations
    tmp = zeros(eltype(zk),length(weights))
    # The near-field correction is here set to the zero matrix. Warn the user.
    @warn "No-near field correction computed"
    nearfield_correction = 0*I
    return FMMGOperator(n,m,zk,eps,targets,sources,weights,elm_interp,
                            physics_topology,nearfield_correction,tmp)
end
function FMMGOperator(mesh,k;eps=1e-6,n=3,nearfield=true,offset=0.2,depth=1)
    # Creating physics and geometry for the FMM operator
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n)...)
    copy_interpolation_nodes!(physics_function,shape_function)
    # Making sure the wave number is complex
    zk = Complex(k)
    # Interpolating on the mesh using the previous computed shape/physics functions
    interpolations = interpolate_elements(mesh,shape_function)
    sources,weights,_ = unroll_interpolations(interpolations)
    # Extracting mesh information
    targets = mesh.sources
    physics_topology = mesh.physics_topology
    element_interpolation = physics_function.interpolation
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
    tmp = zeros(eltype(zk),length(weights))
    # FMM3D also requires the weights to be complex
    return FMMGOperator(N,M,zk,eps,targets,sources,Complex.(weights),element_interpolation,
                            physics_topology,nearfield_correction,tmp)
end
#==========================================================================================
                            Defining H-operator (double-layer)
==========================================================================================#
"""
    FMMHOperator

A `LinearMap` that represents the BEM H matrix through the FMM.
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
    weights::AbstractMatrix{T}          # Nodal normals.*weights (size = m)
    #
    element_interpolation::AbstractMatrix{Float64} # To go from p to
    physics_topology::AbstractMatrix{Int64}
    # Correting near field computations
    nearfield_correction::AbstractMatrix{T}
    #
    tmp::AbstractVecOrMat{T}
    tmp_weights::AbstractVecOrMat{T}
    #
    # C                                   # C-constants
end
Base.size(A::FMMHOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMHOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    # Checking dimensions
    LinearMaps.check_dim_mul(y, A, x)
    # Map the mesh nodes to the gauss nodes
    nodes_to_gauss!(A.tmp,A.element_interpolation,A.physics_topology,x)
    # Scale "dipoles"
    scale_columns!(A.tmp_weights,A.weights,A.tmp)
    # Computing the FMM sum
    vals = hfmm3d(A.eps,A.k,A.sources,targets=A.targets,dipvecs=A.tmp_weights,pgt=1)
    # Ouput equal to the FMM contributions + near field corrections
    # Note that FMM3D uses a Greens function that does not divide by 4π.
    y .= vals.pottarg/4π + A.nearfield_correction*x
end
function scale_columns!(weights,normals,tmp)
    for i = eachindex(tmp)
        weights[1,i] = normals[1,i]*tmp[i]
        weights[2,i] = normals[2,i]*tmp[i]
        weights[3,i] = normals[3,i]*tmp[i]
    end
    return weights
end
function FMMHOperator(eps,k,targets,sources,normals,weights,elm_interp,physics_topology,
                    integral_free_term = [])
    # Making sure the wavenumber is complex (required for the FMM3D library)
    zk = Complex(k)
    # Getting size of the FMM operator.
    N = size(targets,2)
    M = size(sources,2)
    # Compute dipole direction times dipole strengths
    dipvecs = normals .* weights'
    # Create tempory arrays
    tmp = zeros(eltype(zk),length(weights))
    tmp_weights = zeros(eltype(zk),3,length(weights))
    # The near-field correction is here set to the zero matrix. Warn the user.
    @warn "No-near field correction computed"
    if isempty(integral_free_term)
        nearfield_correction = Diagonal(ones(eltype(zk),N)/2)
    elseif length(integral_free_term) == N
        nearfield_correction = Diagonal(integral_free_term)
    end
    # FMM3D uses a Greens function that does not divide by 4π.
    return FMMHOperator(N,M,zk,eps,targets,sources,dipvecs,elm_interp,physics_topology,
                            nearfield_correction,tmp,tmp_weights)
end
function FMMHOperator(mesh,k;n=3,eps=1e-6,nearfield=true,offset=0.2,depth=1,
                                integral_free_term = [],progress=false)
    # Creating physics and geometry for the FMM operator
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n)...)
    copy_interpolation_nodes!(physics_function,shape_function)
    # Making sure that the wavenumber is complex (required for the FMM3D library)
    zk = Complex(k)
    # Interpolating on each element
    interpolations = interpolate_elements(mesh,shape_function)
    # Unrolling element interpolations
    sources,weights,normals = unroll_interpolations(interpolations)
    # Set targets equal to sources
    targets = mesh.sources
    # Getting size of FMM matrix
    N = size(targets,2)
    M = size(sources,2)
    # Computing dipole directional strenght. Making sure they're complex (required by FMm3D)
    dipvecs = Complex.(normals .* weights')
    # Allocating arrays for intermediate computations
    tmp = zeros(eltype(zk),length(weights))
    tmp_weights = zeros(eltype(zk),3,length(weights))
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
    return FMMHOperator(N,M,zk,eps,targets,sources,dipvecs,physics_function.interpolation,
                            mesh.physics_topology,nearfield_correction,tmp,tmp_weights)
end
