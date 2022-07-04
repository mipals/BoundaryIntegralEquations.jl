#==========================================================================================
                            Partial Assembly
==========================================================================================#
function partial_assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                    fOn=true,gOn=true,n=3,progress=true,depth=1)
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
    physics_interpolation1 = copy(physics_function1.interpolation')

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

    I = create_row_indices(lengths,idx[end])
    J = vcat(source_connections...)

    return sparse(I,J,F), sparse(I,J,G)

end

#==========================================================================================
                            Helper Functions
==========================================================================================#
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

# function mygemm_vec2!(C, A, B)
#     @inbounds @fastmath for n ∈ eachindex(C)
#         Cmn = zero(eltype(C))
#         for k ∈ eachindex(A)
#             Cmn += A[k] * B[k,n]
#         end
#         C[n] = Cmn
#     end
# end
function nodes_to_gauss!(tmp,elmement_interpolation,physics_topology,x)
    # Getting sizess
    n_elements = size(physics_topology,2)
    n_interps  = size(elmement_interpolation,2)
    elm_interp = copy(elmement_interpolation')
    # Interpolating on each element
    @inbounds for i = 1:n_elements
        tmp[(i-1)*n_interps+1:i*n_interps] = elm_interp*x[physics_topology[:,i]]
        # mul!(tmp[(i-1)*n_interps+1:i*n_interps],elm_interp,x[physics_topology[:,i]])
    end
    return tmp
end
#==========================================================================================
                            Defining G-operator (single-layer)
==========================================================================================#
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
Base.size(A::FMMGOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMGOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    nodes_to_gauss!(A.tmp_weights,A.element_interpolation,A.physics_topology,x)
    integrand_mul!(A.tmp_weights,A.weights)
    vals = hfmm3d(A.eps,A.k,A.sources,charges=A.weights,targets=A.targets,pgt=1)
    y .= vals.pottarg + A.nearfield_correction*x
end
function FMMGOperator(eps,k,targets,sources,weights,elm_interp,physics_topology)
    zk = Complex(k)
    n = size(targets,2)
    m = size(sources,2)
    tmp = zeros(eltype(zk),length(weights))
    @warn "No-near field correction computed"
    nearfield_correction = 0*I
    # FMM3D uses a greens functiont that does not divide by 4π.
    return FMMGOperator(n,m,zk,eps,targets,sources,weights/(4π + 0im),elm_interp,
                            physics_topology,nearfield_correction,tmp)
end
function FMMGOperator(mesh,k;eps=1e-6,n=3)
    zk = Complex(k)
    interpolations = interpolate_elements(mesh,mesh.shape_function)
    sources,weights,_ = unroll_interpolations(interpolations)
    targets = mesh.sources
    physics_topology = mesh.physics_topology
    n = size(targets,2)
    m = size(sources,2)
    tmp = zeros(eltype(zk),length(weights))
    # Computing near-field correction
    _,C = partial_assemble_parallel!(mesh,zk,mesh.sources,mesh.shape_function;fOn=false)
    _,S = assemble_parallel!(mesh,zk,mesh.sources;sparse=true,depth=1,fOn=false,n=n);
    nearfield_correction = - C + S
    #
    element_interpolation = mesh.physics_function.interpolation
    # FMM3D uses a greens functiont that does not divide by 4π.
    return FMMGOperator(n,m,zk,eps,targets,sources,weights/(4π + 0im),element_interpolation,
                            physics_topology,nearfield_correction,tmp)
end
#==========================================================================================
                            Defining F-operator (double-layer)
==========================================================================================#
struct FMMFOperator{T} <: LinearMaps.LinearMap{T}
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
end
Base.size(A::FMMFOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMFOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    nodes_to_gauss!(A.tmp,A.element_interpolation,A.physics_topology,x)
    # integrand_mul!(A.tmp,A.weights)
    scale_columns!(A.tmp_weights,A.weights,A.tmp)
    vals = hfmm3d(A.eps,A.k,A.sources,targets=A.targets,dipvecs=A.tmp_weights,pgt=1)
    y .= vals.pottarg + A.nearfield_correction*x
end
function scale_columns!(weights,normals,tmp)
    for i = 1:length(tmp)
        weights[1,i] = normals[1,i]*tmp[i]
        weights[2,i] = normals[2,i]*tmp[i]
        weights[3,i] = normals[3,i]*tmp[i]
    end
    return weights
end
function FMMFOperator(eps,k,targets,sources,normals,weights,elm_interp,physics_topology)
    zk = Complex(k)
    n = size(targets,2)
    m = size(sources,2)
    dipvecs = normals .* weights'
    tmp = zeros(eltype(zk),length(weights))
    tmp_weights = zeros(eltype(zk),3,length(weights))
    @warn "No-near field correction computed"
    nearfield_correction = 0.5*I
    # FMM3D uses a greens functiont that does not divide by 4π.
    return FMMFOperator(n,m,zk,eps,targets,sources,dipvecs/(4π + 0im),elm_interp,physics_topology,
                            nearfield_correction,tmp,tmp_weights)
end
function FMMFOperator(mesh,k;n=3,eps=1e-6)
    shape_function   = deepcopy(mesh.shape_function)
    physics_function = deepcopy(mesh.physics_function)
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n)...)
    copy_interpolation_nodes!(physics_function,shape_function)

    zk = Complex(k)
    interpolations = interpolate_elements(mesh,shape_function)
    sources,weights,normals = unroll_interpolations(interpolations)
    targets = mesh.sources
    physics_topology = mesh.physics_topology
    n = size(targets,2)
    m = size(sources,2)
    dipvecs = normals .* weights'
    tmp = zeros(eltype(zk),length(weights))
    tmp_weights = zeros(eltype(zk),3,length(weights))
    # Computing near-field correction
    # _,C = partial_assemble_parallel!(mesh,zk,mesh.sources,shape_function;fOn=false)
    # _,S = assemble_parallel!(mesh,zk,mesh.sources;n=n,sparse=true,depth=1,fOn=false);
    # nearfield_correction = - C + S + 0.5*I
    nearfield_correction = spzeros(ComplexF64,n,n) + 0.5*I
    # Element interpolation
    element_interpolation = physics_function.interpolation
    # FMM3D uses a greens functiont that does not divide by 4π.
    return FMMFOperator(n,m,zk,eps,targets,sources,dipvecs/(4π + 0im),element_interpolation,
                        physics_topology,nearfield_correction,tmp,tmp_weights)
end
