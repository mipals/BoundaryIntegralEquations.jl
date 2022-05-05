"""
    get_element_normals(shape_function::Triangular,coordinates, topology)

Compute the normal as the average normal of all normals from the different elements
that a `coordinate` is connected to via `topology`.
"""
function get_element_normals(shape_function::SurfaceFunction,coordinates, topology)
    # Copying geometry element and setting interpolation nodes to be at corners & midpoints
    geometryElement = deepcopy(shape_function)
    set_nodal_interpolation!(geometryElement)

    # Getting number shape functions as well as number of elements
    nShape = number_of_shape_functions(geometryElement)
    nElements = size(topology,2)

    # Pre-Allocations
    normals  = zeros(size(coordinates))
    avg      = zeros(1,size(coordinates,2))
    tangentX = zeros(3,nShape)
    tangentY = zeros(3,nShape)
    normal   = zeros(3,nShape)

    for element = 1:nElements
        # Extracting element properties (connectivity and coordinates)
        elementNodes        = @view topology[:,element]
        elementCoordinates  = @view coordinates[:,elementNodes]
        # Computing tangential directions as well a a normal at each node
        tangentX .= elementCoordinates * geometryElement.derivatives_u
        tangentY .= elementCoordinates * geometryElement.derivatives_v
        cross_product!(normal,tangentX,tangentY)
        # Normalizing the normal before adding it to the list
        normals[:,elementNodes] .+= normal ./ sqrt.(sum(normal.^2,dims=1))
        # Add one to average (Middle points are connected with maximum 2 elements. Corners more.)
        avg[elementNodes]       .+= 1.0
    end

    # Averaging
    normals = normals ./ avg
 
    return -normals ./ sqrt.(sum(abs2,normals,dims=1)) # Normalizing + fixing orientation
end

function normalize_vector!(normals)
    @inbounds for i = 1:size(normals,2)
        scale = hypot(normals[1,i],normals[2,i],normals[3,i])
        normals[1,i] = normals[1,i]/scale
        normals[2,i] = normals[2,i]/scale
        normals[3,i] = normals[3,i]/scale
    end
end

"""
    tangents!(normals, tangentX, tangentY)

Inplace computation of two tagents of the columns of `normals`.
Results are saved in `tangentX` and `tangentY` respectively.
"""
function tangents!(normals, tangentX, tangentY)
    e2 = [0.0;1.0;0.0]
    e1 = [1.0;0.0;0.0]
    for i = 1:size(tangentX,2)
        normal  = @view normals[:,i]
        dX      = @view tangentX[:,i]
        dY      = @view tangentY[:,i]
        cross!(dX,e2,normal)
        cross!(dY,normal,dX)
        if abs.(dot(e2,normal) - 1.0) <= 1e-10
            cross!(dY,normal,e1)
            cross!(dX,dY,normal)
        end
    end
    normalize_vector!(tangentX)
    normalize_vector!(tangentY)
end

function compute_sources(shape_function,physicsElement,topology,coordinates)
    nShapeFunctions = number_of_shape_functions(physicsElement)
    nElements = size(topology,2)
    nSources  = nShapeFunctions * nElements
    sources   = zeros(3,nSources)
    normals   = similar(sources)
    geometryElement = deepcopy(shape_function)
    set_interpolation_nodes!(geometryElement,physicsElement)
    dX = geometryElement.derivatives_u
    dY = geometryElement.derivatives_v
    
    tangentX            = zeros(3,nShapeFunctions)
    tangentY            = zeros(3,nShapeFunctions)
    normal              = zeros(3,nShapeFunctions)
    physicsTopology     = reshape(collect(1:nSources),nShapeFunctions,nElements)
    @inbounds for element = 1:nElements
        # Access element topology and coordinates
        elementNodes        = @view topology[:,element]
        elementCoordinates  = @view coordinates[:,elementNodes]
        elementSources      = @view physicsTopology[:,element]
        sources[:,elementSources] = elementCoordinates * geometryElement
        tangentX .= elementCoordinates * dX
        tangentY .= elementCoordinates * dY
        cross_product!(normal,tangentX,tangentY)
        normals[:,elementSources] .= -normal ./ sqrt.(sum(normal.^2,dims=1))
    end

    return sources,normals,physicsTopology
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
        # nodes,w = gausslegendre(2), return 1 + nodes[1]
        return 0.4226
    elseif beta_type == :equidistant
        return 0.5
    end
end
function get_beta_quad_quadratic(beta_type)
    if typeof(beta_type) <: Number
        return beta_type
    elseif beta_type == :legendre
        # Zeros of the Legendre polynomials
        # nodes,w = gausslegendre(3), return 1 + nodes[1]
        return 0.2254
    elseif beta_type == :equidistant
        return 1.0/3.0
    end
end

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

function remove_unused_nodes(topology)
    used_nodes = sort(unique(topology))
    maxium_node_number = maximum(used_nodes)
    if length(used_nodes) == maxium_node_number
        return topology
    end
    offset = collect(1:maxium_node_number)
    offset[used_nodes] .= 0
    nonzero = findall(x -> x != 0,offset)
    offset = collect(1:maxium_node_number)
    for nz in nonzero[end:-1:1]
        offset[nz:end] .-= 1
    end
    new_topology = similar(topology)
    nElementNodes, nElements = size(topology)
    for i = 1:nElementNodes
        for j = 1:nElements
            new_topology[i,j] = offset[topology[i,j]]
        end
    end
    return new_topology
end
