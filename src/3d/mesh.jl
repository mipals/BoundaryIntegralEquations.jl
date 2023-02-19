"""
    remove_unused_nodes(topology)

Reindexing topology so that the indexing has no gaps.
"""
function remove_unused_nodes(topology)
    # Find used nodes
    used_nodes = sort(unique(topology))
    # Find the maximum node value
    maxium_node_number = maximum(used_nodes)
    # If all nodes are used return the topology (i.e. remove nothing)
    if length(used_nodes) == maxium_node_number
        return topology
    end
    # If some nodes are not used make "master" list of all possible node numbers
    offset = collect(1:maxium_node_number)
    # Finding the nodes not used
    offset[used_nodes] .= 0
    nonzero = findall(x -> x != 0,offset)
    # Remake master list
    offset = collect(1:maxium_node_number)
    # For every node not used subtract values from the nodes with higher index
    # this essentially makes a list of the "new" nodal indicies
    for nz in nonzero[end:-1:1]
        offset[nz:end] .-= 1
    end
    # Return the topology with the new nodal indicies
    return offset[topology]
end

"""
    normalize_vector!(normals)

Normalizes the columns of `normals`.
"""
function normalize_vector!(normals)
    @inbounds for i = axes(normals,2)
        scale = hypot(normals[1,i],normals[2,i],normals[3,i])
        normals[1,i] = normals[1,i]/scale
        normals[2,i] = normals[2,i]/scale
        normals[3,i] = normals[3,i]/scale
    end
    return normals
end

"""
    tangents!(normals, tangent1, tangent2)

Inplace computation of two tagents of the columns of `normals`.
Results are saved in `tangent1` and `tangent2` respectively.
"""
function tangents!(normals, tangent1, tangent2)
    # Defining standard basis vectors
    e2 = [0.0;1.0;0.0]
    e1 = [1.0;0.0;0.0]
    for i = 1:size(tangent1,2)
        normal  = @view normals[:,i]
        dX      = @view tangent1[:,i]
        dY      = @view tangent2[:,i]
        cross!(dX,e2,normal)
        cross!(dY,normal,dX)
        if abs.(dot(e2,normal) - 1.0) <= 1e-10
            cross!(dY,normal,e1)
            cross!(dX,dY,normal)
        end
    end
    normalize_vector!(tangent1)
    normalize_vector!(tangent2)
    return tangent1,tangent2
end

"""
    get_element_normals(shape_function::Triangular,coordinates, topology)

Compute the normal as the average normal of all normals from the different elements
that a `coordinate` is connected to via `topology`.
"""
function get_element_normals(shape_function::SurfaceFunction,coordinates,topology)
    # Copying geometry element
    surface_function = deepcopy(shape_function)
    # Setting the interpolation to be on the nodal values
    set_nodal_interpolation!(surface_function)
    # Getting number shape functions as well as number of elements
    n_shape_functions = number_of_shape_functions(surface_function)
    n_elements        = size(topology,2)
    # Pre-Allocations
    normals  = zeros(size(coordinates))
    average  = zeros(1,size(coordinates,2))
    tangent1 = zeros(3,n_shape_functions)
    tangent2 = zeros(3,n_shape_functions)
    normal   = zeros(3,n_shape_functions)
    for element = 1:n_elements
        # Extracting element properties (connectivity and coordinates)
        element_nodes       = @view topology[:,element]
        element_coordinates = @view coordinates[:,element_nodes]
        # Computing nodal normals
        mul!(tangent1,element_coordinates,surface_function.derivatives_u)
        mul!(tangent2,element_coordinates,surface_function.derivatives_v)
        cross_product!(normal,tangent1,tangent2)
        # Normalizing the normal before adding it to the list
        normals[:,element_nodes] .+= normal ./ sqrt.(sum(normal.^2,dims=1))
        # Add one to average (Midpoints are always connected to 2 elements. Corners varies.)
        average[element_nodes]       .+= 1.0
    end
    # Setting normal equal to the average normal
    normals = normals ./ average
    return -normals ./ sqrt.(sum(abs2,normals,dims=1)) # Normalizing + fixing orientation
end

"""
    compute_sources(shape_function,physics_function,topology,coordinates)

Computes the nodal position for discontinuous elements.
"""
function compute_sources(shape_function,physics_function,topology,coordinates)
    # Getting the number of functions for the discontinuous element
    n_shape_functions = number_of_shape_functions(physics_function)
    # Get the number of elemenets
    n_elements = size(topology,2)
    # The number of sources are computed
    n_sources = n_shape_functions * n_elements
    # Making the physics topology
    physics_topology = reshape(collect(1:n_sources),n_shape_functions,n_elements)
    # Setting up the interpolation to be the same as the physics functions
    surface_function = deepcopy(shape_function)
    set_interpolation_nodes!(surface_function,physics_function)
    # Preallocation
    sources   = zeros(3,n_sources)
    normals   = similar(sources)
    tangent1  = zeros(3,n_shape_functions)
    tangent2  = zeros(3,n_shape_functions)
    normal    = zeros(3,n_shape_functions)
    @inbounds for element = 1:n_elements
        # Access element topology and coordinates
        element_nodes       = @view topology[:,element]
        element_coordinates = @view coordinates[:,element_nodes]
        element_sources     = @view physics_topology[:,element]
        # Interpolate on element to find source position
        sources[:,element_sources] = element_coordinates * surface_function
        # Compute source normals
        mul!(tangent1,element_coordinates,surface_function.derivatives_u)
        mul!(tangent2,element_coordinates,surface_function.derivatives_v)
        cross_product!(normal,tangent1,tangent2)
        normals[:,element_sources] .= -normal ./ sqrt.(sum(normal.^2,dims=1))
    end
    return sources,normals,physics_topology
end

"""
    set_physics_element(physics_order,shape_function,beta_type)

Returns element of type `physics_order` with the same interpolation as `shape_function`.
For discontinuous elements `beta_type` defines the position of the internal nodes.
"""
function set_physics_element(physics_order,shape_function,beta_type)
    if physics_order == :geometry
        return deepcopy(shape_function)
    elseif physics_order == :linear
        return TriangularLinear(shape_function)
    elseif physics_order == :disctriconstant
        return DiscontinuousTriangularConstant(shape_function)
    elseif physics_order == :disctrilinear
        beta = get_beta_tri_linear(beta_type)
        return DiscontinuousTriangularLinear(shape_function,beta)
    elseif physics_order == :disctriquadratic
        beta = get_beta_tri_quadratic(beta_type)
        return DiscontinuousTriangularQuadratic(shape_function,beta)
    elseif physics_order == :discquadconstant
        return DiscontinuousQuadrilateralConstant(shape_function)
    elseif physics_order == :discquadlinear
        beta = get_beta_tri_linear(beta_type)
        return DiscontinuousQuadrilateralLinear4(shape_function,beta)
    elseif physics_order == :discquadquadratic
        beta = get_beta_tri_quadratic(beta_type)
        return DiscontinuousQuadrilateralQuadraticLagrange(shape_function,beta)
    else
        return deepcopy(shape_function)
    end
end

function get_beta_tri_linear(beta_type)
    if typeof(beta_type) <: Number
        return beta_type
    elseif beta_type == :legendre
        return 0.1667
    elseif beta_type == :equidistant
        return 0.25
    end
end
function get_beta_tri_quadratic(beta_type)
    if typeof(beta_type) <: Number
        return beta_type
    elseif beta_type == :legendre
        return 0.0916
    elseif beta_type == :equidistant
        return 0.1667
    end
end
function get_beta_quad_linear(beta_type)
    if typeof(beta_type) <: Number
        return beta_type
    elseif beta_type == :legendre
        # Zeros of the Legendre polynomials
        nodes,_ = gausslegendre(2)
        return 1 + nodes[1]
    elseif beta_type == :equidistant
        return 0.5
    end
end
function get_beta_quad_quadratic(beta_type)
    if typeof(beta_type) <: Number
        return beta_type
    elseif beta_type == :legendre
        # Zeros of the Legendre polynomials
        nodes,_ = gausslegendre(3)
        return 1 + nodes[1]
    elseif beta_type == :equidistant
        return 1.0/3.0
    end
end


"""
    get_hmax(mesh::Mesh3d)

Computes the `h_max` of the mesh.
"""
function get_hmax(mesh::Mesh3d)
    # Interpolating on each elements
    interpolations = interpolate_elements(mesh)
    # Preallocation
    areas = zeros(length(interpolations))
    # Extracting the areas of each elements as the sum of jacobian*weights
    for (index,interpolation) in enumerate(interpolations)
        areas[index] = sum(interpolation.jacobian_mul_weights)
    end
    # Return correct hmax depending on element type
    if typeof(mesh.shape_function) <: Triangular
        return sqrt(2*maximum(areas))
    elseif typeof(mesh.shape_function) <: Quadrilateral
        return sqrt(maxium(areas))
    end
end
