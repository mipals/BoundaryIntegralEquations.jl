"""
    get_element_normals(basisElement::Triangular,coordinates, topology)

Compute the normal as the average normal of all normals from the different elements
that a `coordinate` is connected to via `topology`.
"""
function get_element_normals(basisElement::SurfaceFunction,coordinates, topology)
    # Copying geometry element and setting interpolation nodes to be at corners & midpoints
    geometryElement = deepcopy(basisElement)
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

function compute_sources(basisElement,physicsElement,topology,coordinates)
    nShapeFunctions = number_of_shape_functions(physicsElement)
    nElements = size(topology,2)
    nSources  = nShapeFunctions * nElements
    sources   = zeros(3,nSources)
    normals   = similar(sources)
    geometryElement = deepcopy(basisElement)
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

function get_beta_tri_linear(betaType)
    if lowercase(betaType) == "legendre"
        return 0.1667
    elseif lowercase(betaType) == "equidistant"
        return 0.25
    end
end
function get_beta_tri_quadratic(betaType)
    if typeof(betaType) <: Number
        return betaType
    elseif lowercase(betaType) == "legendre"
        return 0.0916
    elseif lowercase(betaType) == "equidistant"
        return 0.1667
    end
end
function get_beta_quad_linear(betaType)
    if lowercase(betaType) == "legendre"
        return 0.4226
    elseif lowercase(betaType) == "equidistant"
        return 0.5
    end
end
function get_beta_quad_quadratic(betaType)
    if typeof(betaType) <: Number
        return betaType
    elseif lowercase(betaType) == "legendre"
        return 0.2254
    elseif lowercase(betaType) == "equidistant"
        return 1.0/3.0
    end
end

function set_physics_element(physicsType,basisElement,betaType)
    if lowercase(physicsType) == "geometry"
        return deepcopy(basisElement)
    elseif lowercase(physicsType) == "linear"
        return TriangularLinear(basisElement)
    elseif lowercase(physicsType) == "disctriconstant"
        return DiscontinuousTriangularConstant(basisElement)
    elseif lowercase(physicsType) == "disctrilinear"
        beta = get_beta_tri_linear(betaType)
        return DiscontinuousTriangularLinear(basisElement,beta)
    elseif lowercase(physicsType) == "disctriquadratic"
        beta = get_beta_tri_quadratic(betaType)
        return DiscontinuousTriangularQuadratic(basisElement,beta)
    elseif lowercase(physicsType) == "discquadconstant"
        error("Discontinuous Constant Quadrilaterals are not implemented yet")
    elseif lowercase(physicsType) == "discquadlinear"
        beta = get_beta_tri_linear(betaType)
        error("Discontinuous Linear Quadrilaterals are not implemented yet")
        # return DiscontinuousQuadrilateralLinear4(basisElement,beta)
    elseif lowercase(physicsType) == "discquadquadratic"
        beta = get_beta_tri_quadratic(betaType)
        error("Discontinuous Quadratic Quadrilaterals are not implemented yet")
        # return DiscontinuousQuadrilateralLinear9(basisElement,beta)
    else 
        return deepcopy(basisElement)
    end
end

function remove_unused_nodes(topology)
    usedNodes = sort(unique(topology))
    maxNodeNumber = maximum(usedNodes)
    if length(usedNodes) == maxNodeNumber
        return topology
    end
    offset = collect(1:maxNodeNumber)
    offset[usedNodes] .= 0
    nonzero = findall(x -> x != 0,offset)
    offset = collect(1:maxNodeNumber)
    for nz in nonzero[end:-1:1]
        offset[nz:end] .-= 1
    end
    newTopology = similar(topology)
    nElementNodes, nElements = size(topology)
    for i = 1:nElementNodes
        for j = 1:nElements
            newTopology[i,j] = offset[topology[i,j]]
        end
    end
    return newTopology
end
