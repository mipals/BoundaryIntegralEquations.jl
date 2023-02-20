#==========================================================================================
                Saving the element interpolations to avoid recomputations
==========================================================================================#
struct SurfaceElement{T<:AbstractFloat}
    jacobian_mul_weights::AbstractArray{T,1}
    interpolation::AbstractArray{T,2}
    normals::AbstractArray{T,2}
    center::AbstractArray{T,1}
    max_side_length::T
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

    center_element = define_center_element(shape_function)

    n_shape_functions = number_of_shape_functions(shape_function)
    element_coordinates = zeros(3,n_shape_functions)
    for element = 1:n_elements
        # Extracting element properties (connectivity and coordinates)
        element_coordinates .= mesh.coordinates[:,mesh.topology[:,element]]
        # Computing interpolation
        my_mul!(interps,element_coordinates,shape_function.interpolation)
        # Computing tangential directions as well a a normal at each node
        jacobian!(shape_function,element_coordinates,normals,tangents,sangents,jacobians)
        # Computing center of element
        center = element_coordinates * center_element.interpolation
        # Computing maximum sidelength
        max_side = maximum_sidelength(element_coordinates[:,1],
                                      element_coordinates[:,2],
                                      element_coordinates[:,3])
        # Save element interpolations
        element_interpolations[element] = SurfaceElement(deepcopy(jacobians.*shape_function.weights),
                                                         deepcopy(interps),
                                                         deepcopy(normals),
                                                         deepcopy(center[:]),
                                                         deepcopy(max_side))
    end
    return element_interpolations
end

"""
    define_center_element(surface_function)

Returns a surface function with interpolation point equal to the center of the element.
"""
function define_center_element(surface_function::Triangular)
    center_element = deepcopy(surface_function)
    set_interpolation_nodes!(center_element,[1.0/3.0],[1.0/3.0])
    return center_element
end

function define_center_element(surface_function::Quadrilateral)
    center_element = deepcopy(surface_function)
    set_interpolation_nodes!(center_element,[0.0],[0.0])
    return center_element
end

"""
    maximum_sidelength(x,y,z)

Computes the maximum distance between ``\\mathbf{x}, \\mathbf{y}, \\mathbf{z} \\in \\mathbb{R}^3``.
"""
function maximum_sidelength(x,y,z)
    d1 = hypot(x[1] - y[1],x[2] - y[2], x[3] - y[3])
    d2 = hypot(y[1] - z[1],y[2] - z[2], y[3] - z[3])
    d3 = hypot(z[1] - x[1],z[2] - x[2], z[3] - x[3])
    return max(d1,d2,d3)
end

function compute_distance(source,center)
    return hypot(source[1] - center[1], source[2] - center[2], source[3] - center[3])
end

"""
    find_closest_corner(source,element_coordinates)

Computes the element corner which the source is closest to.
"""
function find_closest_corner(source,element_coordinates)
    max_distance   = Inf
    closest_corner = 0
    @inbounds for i = axes(element_coordinates,2)
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

"""
    create_rotated_element(basisElement::Triangular,m,n,clusterCorner)

Creates a rotated triangular element of the same as the input basisElement.
"""
function create_rotated_element(shape_function,n::Real,m::Real,clusterCorner)
    if clusterCorner == 1
        nodes_u,nodes_v,weights = rotated_triangular_quadpoints(n,m)
    elseif clusterCorner == 2
        nodes_u,nodes_v,weights = triangularQuadpoints(n,m)
    elseif clusterCorner == 3
        nodes_v,nodes_u,weights = triangularQuadpoints(n,m)
    end
    rotated_element = deepcopy(shape_function)
    rotated_element.gauss_u = nodes_u
    rotated_element.gauss_v = nodes_v
    rotated_element.weights = weights
    interpolate_on_nodes!(rotated_element)
    return rotated_element
end
#==========================================================================================
                                Utility functions
==========================================================================================#
"""
    compute_distances!(r,interpolation,sources)

Computes the euclidean distance between the ``j``th column of `sources` and the ``i``th columns of `interpolation` and saves them in `r[i,j]`.
In many case the sources only have a single column.
"""
function compute_distances!(r,interpolation,sources)
    @inbounds for i = 1:size(r,1), j = 1:size(r,2)
        r[i,j] = hypot(interpolation[1,i] - sources[1,j],
                       interpolation[2,i] - sources[2,j],
                       interpolation[3,i] - sources[3,j])
    end
    return r
end

"""
    integrand_mul!(integrand,jacobian)

Inplace multiplication of integrand[i] and jacobian[i] saved in integrand[i].
Used to e.g. scale the jacobian at the Gaussian points with the Gaussian weights.
"""
function integrand_mul!(integrand,jacobian)
    @inbounds @fastmath for i = eachindex(integrand)
        integrand[i] = integrand[i]*jacobian[i]
    end
    return integrand
end

"""
    dotC!(C,integrand,jacobian)

Adds the innerproduct of `integrand` and `jacobian` to C[1].
"""
function dotC!(C,integrand,jacobian)
    @inbounds @fastmath for i = eachindex(integrand)
        C[1] = C[1] + integrand[i] * jacobian[i]
    end
    return C
end

"""
    dotC!(C,integrand,jacobian)

Adds the product of `integrand[i]` and `jacobian[i]` to c[i].
Used when e.g. computing the integral free term.
"""
function c_integrand_mul!(c,integrand,jacobian)
    @inbounds @fastmath for i = eachindex(c)
        c[i] += integrand[i]*jacobian[i]
    end
    return c
end

"""
    sum_to_c!(y,integrand)

Adds all elements of `integrand` to c.
"""
function sum_to_c!(c,integrand)
    @inbounds @fastmath for i = eachindex(integrand)
        c[1] += integrand[i]
    end
    return c
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
                                fOn,gOn,cOn,
                                submatrixF,submatrixG,subvectorC,k,
                                source,integrand,r)
    interpolation        = interpolation_element.interpolation
    normals              = interpolation_element.normals
    jacobian_mul_weights = interpolation_element.jacobian_mul_weights
    # Note that everything here has been written to avoid re-allocating things.
    # As such we re-use the memory of "integrand" for all computations
    compute_distances!(r,interpolation,source)
    if gOn
        ### Evaluating the G-kernel (single-layer kernel) at the global nodes
        greens3d!(integrand,r,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating the integral and the adding the values to the BEM matrix
        mygemm_vec!(submatrixG,integrand,physics_interpolation)
    end
    if cOn
        ### Evaluating the G0-kernel (used for computing the C-constant)
        freens3dk0!(integrand,r,interpolation,source,interpolation_element.normals)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating the integral and adding the value to the C-vector
        sum_to_c!(subvectorC,integrand)
    end
    if fOn && cOn && false
        # Converting the integrand of the integral-free-form to the double-layer integrand
        freens3dk0_to_freens3d!(integrand,r,k)
        # Computing the integrand
        mygemm_vec!(submatrixF,integrand,physics_interpolation)
    elseif fOn
        ### Evaluating the F-kernel (double-layer kernel) at the global nodes
        freens3d!(integrand,r,interpolation,source,normals,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating integral and adding the value to the BEM matrix
        mygemm_vec!(submatrixF,integrand,physics_interpolation)
    end
end

"""
    assemble_parallel!(mesh::Mesh3d,k,sources;m=4,n=4,progress=true)

Assembles the BEM matrices for F, G and G0 kernels over the elements on the mesh.
"""
function assemble_parallel!(mesh::Mesh3d,k,in_sources;fOn=true,gOn=true,cOn=true,
                            sparse=false,m=4,n=4,progress=true,depth=2,offset=nothing)
    if sparse
        return sparse_assemble_parallel!(mesh,k,in_sources,mesh.shape_function;
                                fOn=fOn,gOn=gOn,progress=progress,depth=depth,offset=offset)
    else
        return assemble_parallel!(mesh::Mesh3d,k,in_sources,mesh.shape_function;
                                fOn=fOn,gOn=gOn,cOn=cOn,m=m,n=n,progress=progress)
    end
end

function assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                fOn=true,gOn=true,cOn=true,m=3,n=3,progress=true)
    n_elements  = number_of_elements(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    n_nodes     = size(mesh.sources,2)
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

    # Copying interpolation of physics functions1
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
        integrand  = zeros(ComplexF64,n*m)
        r          = zeros(Float64,n*m)
        subvectorC = @view C[source_node]
        # element_coordinates = zeros(3,n_physics_functions)
        @inbounds for element = 1:n_elements
            # Access element topology and coordinates
            physics_nodes        = @view physics_topology[:,element]
            physics_coordinates  = @view sources[:,physics_nodes]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[source_node,physics_nodes]
            submatrixG = @view G[source_node,physics_nodes]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,physics_coordinates)
            if close_corner == 1
                computing_integrals!(physics_interpolation1,interpolation_list1[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            elseif close_corner == 2
                computing_integrals!(physics_interpolation2,interpolation_list2[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            else
                computing_integrals!(physics_interpolation3,interpolation_list3[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            end
        end
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end

function assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::TriangularQuadratic;
                                fOn=true,gOn=true,cOn=true,m=3,n=3,progress=true)
    n_elements  = number_of_elements(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    n_nodes     = size(mesh.sources,2)
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

    # nodesX4,nodesY4,weights4 = getpolar_gaussian(n,4)
    # nodesX5,nodesY5,weights5 = getpolar_gaussian(n,5)
    # nodesX6,nodesY6,weights6 = getpolar_gaussian(n,6)
    nodesX4,nodesY4,weights4 = rotated_midpoint_triangular_quadpoints(n,m,4)
    nodesX5,nodesY5,weights5 = rotated_midpoint_triangular_quadpoints(n,m,5)
    nodesX6,nodesY6,weights6 = rotated_midpoint_triangular_quadpoints(n,m,6)
    physics_function4 = deepcopy(physics_function)
    physics_function5 = deepcopy(physics_function)
    physics_function6 = deepcopy(physics_function)
    shape_function4 = deepcopy(shape_function)
    shape_function5 = deepcopy(shape_function)
    shape_function6 = deepcopy(shape_function)
    set_interpolation_nodes!(shape_function4,nodesX4,nodesY4,weights4)
    set_interpolation_nodes!(shape_function5,nodesX5,nodesY5,weights5)
    set_interpolation_nodes!(shape_function6,nodesX6,nodesY6,weights6)
    copy_interpolation_nodes!(physics_function4,shape_function4)
    copy_interpolation_nodes!(physics_function5,shape_function5)
    copy_interpolation_nodes!(physics_function6,shape_function6)

    # Computing interpolation on each element
    interpolation_list1 = interpolate_elements(mesh,shape_function1)
    interpolation_list2 = interpolate_elements(mesh,shape_function2)
    interpolation_list3 = interpolate_elements(mesh,shape_function3)
    interpolation_list4 = interpolate_elements(mesh,shape_function4)
    interpolation_list5 = interpolate_elements(mesh,shape_function5)
    interpolation_list6 = interpolate_elements(mesh,shape_function6)

    # Copying interpolation of physics functions1
    physics_interpolation1 = copy(physics_function1.interpolation')
    physics_interpolation2 = copy(physics_function2.interpolation')
    physics_interpolation3 = copy(physics_function3.interpolation')
    physics_interpolation4 = copy(physics_function4.interpolation')
    physics_interpolation5 = copy(physics_function5.interpolation')
    physics_interpolation6 = copy(physics_function6.interpolation')

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
        integrand  = zeros(ComplexF64,n*m)
        r          = zeros(Float64,n*m)
        mid_integrand = zeros(ComplexF64,length(nodesX4))
        mid_r         = zeros(Float64,length(nodesX4))
        subvectorC    = @view C[source_node]
        @inbounds for element = 1:n_elements
            # Access element topology and coordinates
            physics_nodes        = @view physics_topology[:,element]
            physics_coordinates  = @view sources[:,physics_nodes]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[source_node,physics_nodes]
            submatrixG = @view G[source_node,physics_nodes]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,physics_coordinates)
            if close_corner == 1
                computing_integrals!(physics_interpolation1,interpolation_list1[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            elseif close_corner == 2
                computing_integrals!(physics_interpolation2,interpolation_list2[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            elseif close_corner == 3
                computing_integrals!(physics_interpolation3,interpolation_list3[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r)
            elseif close_corner == 4
                computing_integrals!(physics_interpolation4,interpolation_list4[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,mid_integrand,mid_r)
            elseif close_corner == 5
                computing_integrals!(physics_interpolation5,interpolation_list5[element],
                                    fOn,gOn,cOn,
                                    submatrixF,submatrixG,subvectorC,k,
                                    source,mid_integrand,mid_r)
            else
                computing_integrals!(physics_interpolation6,interpolation_list6[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,mid_integrand,mid_r)
            end
        end
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end

function assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::SurfaceFunction;
                            fOn=true,gOn=true,cOn=true,m=4,n=4,progress=true)
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
    interpolation_list = interpolate_elements(mesh,shape_function1)
    # Avoiding to have gauss-node on a singularity. Should be handled differently.
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
        integrand  = zeros(ComplexF64,n*m)
        r          = zeros(Float64,n*m)
        subvectorC = @view C[source_node]
        @inbounds for element = 1:n_elements
            physics_nodes = @view physics_topology[:,element]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[source_node,physics_nodes]
            submatrixG = @view G[source_node,physics_nodes]
            # Interpolating on the mesh.
            computing_integrals!(physics_interpolation,interpolation_list[element],
                                    fOn,gOn,cOn,
                                    submatrixF,submatrixG,subvectorC,k,
                                    source,integrand,r)
        end
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end
