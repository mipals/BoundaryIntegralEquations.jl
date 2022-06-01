#==========================================================================================
                Saving the element interpolations to avoid recomputations
==========================================================================================#
struct CurveElement{T<:AbstractFloat}
    jacobian_mul_weights::AbstractArray{T,1}
    interpolation::AbstractArray{T,2}
    normals::AbstractArray{T,2}
end
function create_shape_function(shape_function::CurveFunction;n=4)
    nodes,weights               = gausslegendre(n)
    new_shape_function          = deepcopy(shape_function)
    new_shape_function.gauss    = nodes
    new_shape_function.weights  = weights
    interpolate_on_nodes!(new_shape_function)
    return new_shape_function
end
function interpolate_elements(mesh::Mesh2d;n=4)
    shape_function = create_shape_function(mesh.shape_function;n=n)
    interpolate_elements(mesh,shape_function)
end
function interpolate_elements(mesh::Mesh2d,shape_function::CurveFunction)

    @assert typeof(mesh.shape_function) <: typeof(shape_function)

    n_elements = size(mesh.topology,2)

    element_interpolations = Array{CurveElement{eltype(shape_function)}}(undef,n_elements)

    jacobians = similar(shape_function.weights)
    normals   = zeros(2,length(jacobians))
    tangents  = similar(normals)
    interps   = similar(normals)

    for element = 1:n_elements
        # Extracting element properties (connectivity and coordinates)
        element_coordinates  = @view mesh.coordinates[:,mesh.topology[:,element]]
        # Computing interpolation
        mul!(interps,element_coordinates,shape_function.interpolation)
        # Computing tangential directions as well a a normal at each node
        jacobian!(shape_function,element_coordinates,normals,tangents,jacobians)
        # Save element interpolations
        element_interpolations[element] = CurveElement(deepcopy(jacobians.*shape_function.weights),
                                                       deepcopy(interps),
                                                       deepcopy(normals))
    end
    return element_interpolations
end

"""
    compute_distances_2d!(r,interpolation,source)

Saves distances in `r`.
"""
function compute_distances_2d!(r,interpolation,source)
    @turbo for i = 1:size(r,1), j = 1:size(r,2)
        r[i,j] = hypot(interpolation[1,j] - source[1,i],
                       interpolation[2,j] - source[2,i])
    end
end
"""
    normalize_2d!(normals,jacobian)

Dividing each column of `normals` with `jacobian`.
"""
function normalize_2d!(normals,jacobian)
    @turbo for i = 1:length(jacobian)
        normals[1,i] = -normals[1,i]/jacobian[i]
        normals[2,i] = -normals[2,i]/jacobian[i]
    end
end
"""
    column_norms_2d!(jacobian,normals)

Computing norms of columns in `normals`. Saved in `jacobian`.
"""
function column_norms_2d!(jacobian,normals)
    @turbo for i = 1:length(jacobian)
        jacobian[i] = hypot(normals[1,i],normals[2,i])
    end
end
"""
    tangent_to_normal!(normals,tangent)

Computes the normal by rotating of the normal 90 degrees in the plane.
"""
function tangent_to_normal!(normals,tangent;clockwise=true)
    if clockwise
        @inbounds for i = 1:size(normals,2)
            normals[1,i] =  tangent[2,i]
            normals[2,i] = -tangent[1,i]
        end
    else
        @inbounds for i = 1:size(normals,2)
            normals[1,i] = -tangent[2,i]
            normals[2,i] =  tangent[1,i]
        end
    end
end
"""
    jacobian!(basisElement::CurveFunction,coordinates,normals,tangent,jacobian)

Inplace computations of the jacobian of the `basisElement::CurveFunction` at the coordinates.
The results are saved in `tangent`, `normals` and `jacobian`.
"""
function jacobian!(basisElement::CurveFunction,coordinates,normals,tangent,jacobian)
    mul!(tangent,coordinates,basisElement.derivatives)  # Computing tangent vector in X
    tangent_to_normal!(normals,tangent)                 # Computing normal vector
    column_norms_2d!(jacobian,normals)                  # Jacobian = length of the normal
    normalize_2d!(normals,jacobian)
end
"""
    integrand_mul_weights!(integrand,weights)

Inplace scaling of each value of `integrand` by the corresponding value in `weights`.
"""
function integrand_mul_weights!(integrand,weights)
    @fastmath @inbounds for i = 1:length(weights)
        integrand[i] = integrand[i]*weights[i]
    end
end
"""
    add_to_c!(C,integrand,weights)

Inplace adding of `dot(integrand,weights)` into `C[1]`.
"""
function add_to_c!(C,integrand,weights)
    @fastmath @inbounds for i = 1:length(weights)
        C[1] += integrand[i]*weights[i]
    end
end

function mygemm!(C, A, B)
    # @fastmath @inbounds for m ∈ axes(A,1), n ∈ axes(B,2)
        for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] += Cmn
    end
end

function compute_integrands!(Fslice,Gslice,Cslice,fOn,gOn,cOn,
                        physics_interpolation,interpolation_element,source,integrand,r,k)
    # Extracting element data
    normals         = interpolation_element.normals
    interpolation   = interpolation_element.interpolation
    jac_mul_weights = interpolation_element.jacobian_mul_weights
    # Computing distances
    compute_distances_2d!(r,interpolation,source)

    if gOn # Only compute G if you need it
        # Evaluate Greens function
        G!(integrand,k,r)
        # Compute integrand = integrand .* jacobian .* weights
        integrand_mul_weights!(integrand,jac_mul_weights)
        # Integration with physics interpolation
        # mul!(Gslice,integrand,physics_interpolation,true,true)
        # mul!(Gslice,integrand,physics_interpolation)
        # Gslice .+= integrand*physics_interpolation
        mygemm!(Gslice,integrand,physics_interpolation)
    end
    if fOn # Only compute F if you need it
        # Evaluate ∂ₙGreens function
        F!(integrand,interpolation,source,k,normals,r)
        # Compute integrand = integrand .* jacobian .* weights
        integrand_mul_weights!(integrand,jac_mul_weights)
        # Integration with basis function
        # mul!(Fslice,integrand,physics_interpolation,true,true)
        # mul!(Fslice,integrand,physics_interpolation)
        # Fslice .+= integrand*physics_interpolation
        mygemm!(Fslice,integrand,physics_interpolation)
    end
    if cOn # Only compute c if you need it
        C!(integrand,interpolation,source,normals,r)      # Evaluate G₀
        add_to_c!(Cslice,integrand,jac_mul_weights)
        # Cslice .+= mydotavx(integrand,jac_mul_weights)
        # Cslice[1] += dot(integrand,jac_mul_weights)  # Integration of G₀
        # @muladd Cslice = Cslice + integrand*jac_mul_weights
    end
end


function assemble_parallel!(mesh::Mesh2d,k,in_sources;
                    gOn=true,fOn=true,cOn=true,interior=true,n=4,progress=true)
    return assemble_parallel!(mesh::Mesh2d,k,in_sources,mesh.shape_function;
                        gOn=gOn,fOn=fOn,cOn=cOn,interior=interior,n=n,progress=progress)
end

function assemble_parallel!(mesh::Mesh2d,k,in_sources,shape_function::CurveFunction;
                        gOn=true,fOn=true,cOn=true,interior=true,n=4,progress=true)
    # Extracting physics
    physics_topology      = mesh.physics_topology
    physics_function      = create_shape_function(mesh.physics_function;n=n)
    physics_interpolation = convert.(ComplexF64,physics_function.interpolation')
    # Grabbing sizes
    n_physics_functions   = number_of_shape_functions(physics_function)
    n_nodes     = size(mesh.sources,2)
    n_sources   = size(in_sources,2)
    n_elements  = size(physics_topology,2)
    # Interpolating on full mesh
    interpolation_list = interpolate_elements(mesh,n=n)
    # Preallocation
    F,G = (zeros(ComplexF64, n_sources, n_nodes) for _ in 1:2)
    C   =  zeros(ComplexF64, n_sources)
    # Assembly
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    Threads.@threads for source = 1:n_sources
        Cslice      = @view C[source:source]
        source_node = in_sources[:,source]
        r = zeros(1,n)
        integrand     = zeros(ComplexF64,1,n)
        physics_nodes = zeros(Int64,n_physics_functions)
        # Fslice = zeros(ComplexF64,1,n_physics_functions)
        # Gslice = zeros(ComplexF64,1,n_physics_functions)
        for element = 1:n_elements
            physics_nodes .= physics_topology[:,element]
            Fslice         = @view F[source:source,physics_nodes]
            Gslice         = @view G[source:source,physics_nodes]
            compute_integrands!(Fslice,Gslice,Cslice,
                                fOn,gOn,cOn,
                                physics_interpolation,interpolation_list[element],
                                source_node,integrand,r,k)
            # F[source:source,physics_nodes] .+= Fslice
            # G[source:source,physics_nodes] .+= Gslice
        end
        if progress; next!(prog); end
    end
    # return F,G,(interior ? C : 1.0 .+ C)
    return F,G,C
end
