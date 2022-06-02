#==========================================================================================
                Saving the element interpolations to avoid recomputations
==========================================================================================#
struct SurfaceElement{T<:AbstractFloat}
    jacobian_mul_weights::AbstractArray{T,1}
    interpolation::AbstractArray{T,2}
    normals::AbstractArray{T,2}
end
function create_shape_function(shape_function::SurfaceFunction;n=4,m=4)
    nodes_u,nodes_v,weights = getQuadpoints(shape_function,m,n)
    new_shape_function      = deepcopy(shape_function)
    new_shape_function.gauss_u = nodes_u
    new_shape_function.gauss_v = nodes_v
    new_shape_function.weights = weights
    interpolate_on_nodes!(new_shape_function)
    return new_shape_function
end
function interpolate_elements(mesh::Mesh3d;n=4,m=4)
    shape_function = create_shape_function(mesh.shape_function;n=n,m=m)
    interpolate_elements(mesh,shape_function)
end
function interpolate_elements(mesh::Mesh3d,shape_function::SurfaceFunction)

    @assert typeof(mesh.shape_function) <: typeof(shape_function)

    n_elements = size(mesh.topology,2)

    element_interpolations = Array{SurfaceElement{eltype(shape_function)}}(undef,n_elements)

    jacobians = similar(shape_function.weights)
    normals   = zeros(3,length(jacobians))
    tangents  = similar(normals)
    sangents  = similar(normals)
    interps   = similar(normals)

    for element = 1:n_elements
        # Extracting element properties (connectivity and coordinates)
        element_coordinates  = @view mesh.coordinates[:,mesh.topology[:,element]]
        # Computing interpolation
        mul!(interps,element_coordinates,shape_function.interpolation)
        # Computing tangential directions as well a a normal at each node
        jacobian!(shape_function,element_coordinates,normals,tangents,sangents,jacobians)
        # Save element interpolations
        element_interpolations[element] = SurfaceElement(deepcopy(jacobians.*shape_function.weights),
                                                         deepcopy(interps),
                                                         deepcopy(normals))
    end
    return element_interpolations
end

"""
    find_closest_corner(source,element_coordinates)

Computes the element corner which the source is closest to.
"""
function find_closest_corner(source,element_coordinates)
    max_distance   = Inf
    closest_corner = 0
    @inbounds for i = 1:size(element_coordinates,2)
        tmp_dist = hypot(element_coordinates[1,i] - source[1],
                         element_coordinates[2,i] - source[2],
                         element_coordinates[3,i] - source[3])
        if tmp_dist < max_distance
            max_distance   = tmp_dist
            closest_corner = i
        end
    end
    return closest_corner
end

#==========================================================================================
                                Utility functions
==========================================================================================#
function compute_distances!(r,interpolation,source)
    @inbounds for i = 1:size(r,1), j = 1:size(r,2)
        r[i,j] = hypot(interpolation[1,i] - source[1,j],
                       interpolation[2,i] - source[2,j],
                       interpolation[3,i] - source[3,j])
    end
end
function integrand_mul!(integrand,jacobian)
    @inbounds for i = 1:length(integrand)
        integrand[i] = integrand[i]*jacobian[i]
    end
end
function dotC!(C,integrand,jacobian)
    @inbounds for i = 1:length(integrand)
        C[1] += integrand[i] * jacobian[i]
    end
end

#==========================================================================================
                                Assembly Functions
==========================================================================================#
"""
    computing_integrals!

Approximates the integral of `freens3d!`, `greens3d!` and `freens3dk0!` multiplied by the
shapefunction over single `shape_function`, defined by `coordinates`.
"""
function computing_integrals!(physics_interpolation,interpolation_element,
                                submatrixF,submatrixG,subvectorC,k,
                                source,integrand,r)
    interpolation        = interpolation_element.interpolation
    normals              = interpolation_element.normals
    jacobian_mul_weights = interpolation_element.jacobian_mul_weights
    # Note that everything here has been written to avoid re-allocating things.
    # As such we re-use the memory of "integrand" for all computations
    compute_distances!(r,interpolation,source)

    ### Evaluating the F-kernel (double-layer kernel) at the global nodes
    freens3d!(integrand,r,interpolation,source,normals,k)
    # Computing the integrand
    integrand_mul!(integrand,jacobian_mul_weights)
    # Approximating integral and adding the value to the BEM matrix
    # mygemm!(submatrixF,integrand,Transpose(physics_interpolation))
    mygemm_vec!(submatrixF,integrand,physics_interpolation)

    ### Evaluating the G-kernel (single-layer kernel) at the global nodes
    greens3d!(integrand,r,k)
    # onefunction!(integrand,r,k)
    # Computing the integrand
    integrand_mul!(integrand,jacobian_mul_weights)
    # Approximating the integral and the adding the values to the BEM matrix
    # mygemm!(submatrixG,integrand,Transpose(physics_interpolation))
    mygemm_vec!(submatrixG,integrand,physics_interpolation)


    ### Evaluating the G0-kernel (used for computing the C-constant)
    # Recomputing integrand
    freens3dk0!(integrand,r,interpolation,source,normals)
    # Approximating the integral and adding the value to the C-vector
    add_to_c!(subvectorC,integrand,jacobian_mul_weights)
end


"""
    assemble_parallel!(mesh::Mesh3d,k,sources;m=5,n=5,progress=true)

Assembles the BEM matrices for F, G and G0 kernels over the elements on the mesh.
"""
function assemble_parallel!(mesh::Mesh3d,k,in_sources;m=4,n=4,progress=true)
    return assemble_parallel!(mesh::Mesh3d,k,in_sources,mesh.shape_function;
                                    m=m,n=n,progress=progress)
end
function assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                            m=3,n=3,progress=true)
    topology    = get_topology(mesh)
    n_elements  = number_of_elements(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    n_nodes     = size(mesh.sources,2)
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
    shape_function1 = create_rotated_element(shape_function,n,m,1)
    shape_function2 = create_rotated_element(shape_function,n,m,2)
    shape_function3 = create_rotated_element(shape_function,n,m,3)
    physics_function1 = deepcopy(physics_function)
    physics_function2 = deepcopy(physics_function)
    physics_function3 = deepcopy(physics_function)
    copy_interpolation_nodes!(physics_function1,shape_function1)
    copy_interpolation_nodes!(physics_function2,shape_function2)
    copy_interpolation_nodes!(physics_function3,shape_function3)

    # Computing interpolation on each element
    interpolation_list1 = interpolate_elements(mesh,shape_function1)
    interpolation_list2 = interpolate_elements(mesh,shape_function2)
    interpolation_list3 = interpolate_elements(mesh,shape_function3)

    # Converting interpolation to complex numbers
    # physics_interpolation1 = convert.(Complex,physics_function1.interpolation)
    # physics_interpolation2 = convert.(Complex,physics_function2.interpolation)
    # physics_interpolation3 = convert.(Complex,physics_function3.interpolation)
    physics_interpolation1 = copy(physics_function1.interpolation')
    physics_interpolation2 = copy(physics_function2.interpolation')
    physics_interpolation3 = copy(physics_function3.interpolation')

    # Preallocation of return values
    F = zeros(ComplexF64, n_sources, n_nodes)
    G = zeros(ComplexF64, n_sources, n_nodes)
    C = zeros(ComplexF64, n_sources)

    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        integrand = zeros(ComplexF64,n*m)
        r         = zeros(Float64,n*m)
        @inbounds for element = 1:n_elements
            # Access element topology and coordinates
            element_coordinates = @view coordinates[:,topology[:,element]]
            physics_nodes       = @view physics_topology[:,element]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[source_node,physics_nodes]
            submatrixG = @view G[source_node,physics_nodes]
            subvectorC = @view C[source_node]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,element_coordinates)
            if close_corner == 1
                computing_integrals!(physics_interpolation1,interpolation_list1[element],
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            elseif close_corner == 2
                computing_integrals!(physics_interpolation2,interpolation_list2[element],
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            else
                computing_integrals!(physics_interpolation3,interpolation_list3[element],
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            end
        end
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end

function assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::SurfaceFunction;
                            m=4,n=4,progress=true)
    n_elements   = number_of_elements(mesh)
    sources      = convert.(eltype(shape_function),in_sources)
    n_sources    = size(in_sources,2)
    n_nodes      = size(mesh.sources,2)
    physics_topology = mesh.physics_topology
    physics_function = mesh.physics_function
    #======================================================================================

    ======================================================================================#
    shape_function1    = create_shape_function(shape_function;n=n,m=m)
    physics_function1  = deepcopy(physics_function)
    copy_interpolation_nodes!(physics_function1,shape_function1)
    # mesh.physics_function = physics_function1
    interpolation_list = interpolate_elements(mesh,shape_function1)
    # Avoiding to have gauss-node on a singularity. Should be handled differently,
    # but, this is good for now.
    if typeof(shape_function) <: QuadrilateralQuadraticLagrange
        for (x,y) in zip(shape_function1.gauss_u,shape_function1.gauss_v)
            if isapprox.(x, 0.0, atol=1e-15) && isapprox.(y, 0.0, atol=1e-15)
                error("Gauss Node On Singularity.")
            end
        end
    end
    #======================================================================================
                        Preallocation of return values & Intermediate values
    ======================================================================================#
    F = zeros(ComplexF64, n_sources, n_nodes)
    G = zeros(ComplexF64, n_sources, n_nodes)
    C = zeros(ComplexF64, n_sources)

    #======================================================================================
                                    Assembly
    ======================================================================================#
    physics_interpolation = copy(physics_function1.interpolation')
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source    = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        integrand = zeros(ComplexF64,n*m)
        r         = zeros(Float64,n*m)
        @inbounds for element = 1:n_elements
            physics_nodes = @view physics_topology[:,element]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[source_node,physics_nodes]
            submatrixG = @view G[source_node,physics_nodes]
            subvectorC = @view C[source_node]
            # Interpolating on the mesh.
            computing_integrals!(physics_interpolation,interpolation_list[element],
                                submatrixF,submatrixG,subvectorC,k,
                                source,integrand,r)
        end
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end
