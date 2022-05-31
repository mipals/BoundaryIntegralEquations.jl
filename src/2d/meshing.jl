function compute_sources_2d(shape_function,physics_function,topology,coordinates)
    n_shape_functions = number_of_shape_functions(physics_function)
    n_elements = size(topology,2)
    n_sources  = n_shape_functions * n_elements
    sources    = zeros(2,n_sources)
    normals    = similar(sources)
    surface_function = deepcopy(shape_function)
    set_interpolation_nodes!(surface_function,physics_function)
    dX = surface_function.derivatives

    tangent  = zeros(2,n_shape_functions)
    normal   = zeros(2,n_shape_functions)
    physics_topology = reshape(collect(1:n_sources),n_shape_functions,n_elements)
    @inbounds for element = 1:n_elements
        # Access element topology and coordinates
        element_nodes       = @view topology[:,element]
        element_coordinates = @view coordinates[:,element_nodes]
        element_sources     = @view physics_topology[:,element]
        sources[:,element_sources] = element_coordinates * surface_function.interpolation
        tangent .= element_coordinates * dX
        tangent_to_normal!(normal,tangent)
        normals[:,element_sources] .= -normal ./ sqrt.(sum(normal.^2,dims=1))
    end

    return sources,normals,physics_topology
end

function create_pointmat(segments)
    nnodes = 2*Int(sum(segments[:,5]))
    pointmat = fill(1.0,nnodes, 3)
    counter = 1
    for segment in eachrow(segments)
        nel  = Int(segment[5])
        xtmp = range(segment[1],segment[3],length=2*nel+1)
        ytmp = range(segment[2],segment[4],length=2*nel+1)
        pointmat[counter:counter+2*nel-1,1] = xtmp[1:end-1]
        pointmat[counter:counter+2*nel-1,2] = ytmp[1:end-1]
        counter += 2*nel
    end
    return pointmat
end

function create_topology(segments)
    n_elements = Int(sum(segments[:,5]));
    topology = fill(1, n_elements, 4)
    idx = 1
    for i = 1:n_elements
        for j = 1:2
            topology[i,j] = idx
            idx += 1;
        end
        topology[i,3] = idx
    end
    topology[end,3] = 1
    return topology
end

function create_circle(shape_function::ContinuousCurveLinear,n_elements,radius=1.0, angle_output=0)
    θ = reverse(collect(range(-pi,pi,length=n_elements+1)))
    θ = θ[1:end-1]'

    if angle_output == 0
        return radius*[cos.(θ); sin.(θ)]
    elseif angle_output == 1
        return radius*[cos.(θ); sin.(θ)], θ'
    end
end

function create_circle(shape_function::ContinuousCurveQuadratic,n_elements,radius=1.0, angle_output=0)
    θ = reverse(collect(range(-pi,pi,length=2*n_elements+1)))
    θ = θ[1:end-1]'

    if angle_output == 0
        return radius*[cos.(θ); sin.(θ)]
    elseif angle_output == 1
        return radius*[cos.(θ); sin.(θ)], θ'
    end

end

function create_top(shape_function::ContinuousCurveLinear,n_elements)
    topology = ones(Int64, 2, n_elements)
    topology[1,:] = 1:1:n_elements
    topology[2,1:end-1] = 2:1:n_elements
    return topology
end

function create_top(shape_function::ContinuousCurveQuadratic,n_elements)
    topology = zeros(Int64, 3, n_elements)
    idx = 1
    for i = 1:n_elements
        for j = 1:2
            topology[j,i] = idx
            idx += 1
        end
        topology[3,i] = idx
    end
    topology[3,end] = 1
    return topology
end

function circle_nodegen(shape_function,n_elements,radius=1.0)
    coordinates = create_circle(shape_function,n_elements, radius)
    topology    = create_top(shape_function,n_elements)
    return coordinates, topology
end

function normals_to_tangents!(tangents,normals)
    @inbounds for i = 1:size(normals,2)
        tangents[1,i] = -normals[2,i]
        tangents[2,i] =  normals[1,i]
    end
end


function get_element_normals_2d(shape_function,coordinates,topology)
        # Copying geometry element and setting interpolation nodes to be at corners & midpoints
        surface_function = deepcopy(shape_function)
        set_nodal_interpolation!(surface_function)

        # Getting number shape functions as well as number of elements
        n_shape = number_of_shape_functions(surface_function)
        n_elements = size(topology,2)

        # Pre-Allocations
        normals = zeros(size(coordinates))
        avg     = zeros(1,size(coordinates,2))
        tangent = zeros(2,n_shape)
        normal  = zeros(2,n_shape)

        for element = 1:n_elements
            # Extracting element properties (connectivity and coordinates)
            element_nodes        = @view topology[:,element]
            element_coordinates  = @view coordinates[:,element_nodes]
            # Computing tangential directions as well a a normal at each node
            tangent .= element_coordinates * surface_function.derivatives
            tangent_to_normal!(normal,tangent)
            # Normalizing the normal before adding it to the list
            normals[:,element_nodes] .+= normal ./ sqrt.(sum(normal.^2,dims=1))
            # Add one to average (Middle points are connected with maximum 2 elements. Corners more.)
            avg[element_nodes]       .+= 1.0
        end

        # Averaging
        normals = normals ./ avg

        return -normals ./ sqrt.(sum(abs2,normals,dims=1))
end

function mesh_circle(shape_function,n_elements;radius=1.0)
    coordinates, topology = circle_nodegen(shape_function,n_elements,radius)
    normals  = get_element_normals_2d(shape_function,coordinates,topology)
    sources  = coordinates
    tangents = similar(normals)
    normals_to_tangents!(tangents,normals)
    physics_function = shape_function
    physics_topology = topology
    return Mesh2d(sources,coordinates,topology,normals,tangents,
                shape_function,physics_function,physics_topology)
end

function mesh_circle(shape_function,physics_function::DiscontinuousCurveFunction,n_elements;
                    radius = 1.0)
    coordinates, topology            = circle_nodegen(shape_function,n_elements,radius)
    sources,normals,physics_topology = compute_sources_2d(shape_function,physics_function,
                                                            topology,coordinates)
    tangents = similar(normals)
    normals_to_tangents!(tangents,normals)
    return Mesh2d(sources,coordinates,topology,normals,tangents,
                    shape_function,physics_function,physics_topology)
end
