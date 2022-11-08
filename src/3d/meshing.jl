import GeometryBasics

function load_mesh_file(file)
    mesh = load(file)
    mesh_coords = GeometryBasics.coordinates(mesh)
    elements    = GeometryBasics.faces(mesh)
    n_elements  = length(elements)
    n_interpolations = length(elements[1])

    topology = zeros(Int64,n_interpolations,n_elements)
    coords   = zeros(3,length(mesh_coords))
    for (index,element) in enumerate(elements)
        for i = 1:n_interpolations
            topology[i,index] = Int64(element[i].i) + 1
        end
    end
    for (index,coord) in enumerate(mesh_coords)
        coords[:,index] = Float64.(coord)
    end

    ents = zeros(Int64, n_elements)
    return coords, topology, ents
end

function load3dTriangularMesh(meshFile;m=3,n=3,
                geometry_order=:linear,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removedEntites=[-1])

    if geometry_order !== :linear
        error("Currently only linear element supported")
    end

    # Figuring out the element type
    initialCoordinates,initialTopology,ents = load_mesh_file(meshFile)
    if geometry_order == :quadratic
        shape_function = TriangularQuadratic(m,n)
    elseif geometry_order == :linear
        shape_function = TriangularLinear(m,n)
        initialTopology = initialTopology[1:3,:]
    else
        error("Only quadratic and linear geometries are currently supported")
    end

    mask = .!convert.(Bool,sum(ents' .∈ removedEntites,dims=1))[:]
    usedNodes   = sort(unique(initialTopology[:,mask]))
    topology    = remove_unused_nodes(initialTopology[:,mask])
    coordinates = initialCoordinates[:,usedNodes]
    sources = coordinates[:,sort(unique(topology))]

    if physics_order == :geometry
        physics_function = set_physics_element(physics_order,shape_function,beta_type)
        normals = get_element_normals(shape_function,coordinates,topology)
        physics_topology = topology
    elseif physics_order == :linear
        physics_function = TriangularLinear(3,3)
        normals = get_element_normals(shape_function,coordinates,topology)
        physics_topology = topology[1:3,:]
        used_nodes = sort(unique(physics_topology))
        normals = normals[:,used_nodes]
        sources = coordinates[:,used_nodes]
    else
        physics_function = set_physics_element(physics_order,shape_function,beta_type)
        sources,normals,physics_topology = compute_sources(shape_function,
                                                        physics_function,
                                                        topology,
                                                        coordinates)
    end

    # Allocate and compute tangent directions from normal
    # This is to backwards-compatible to the old implementation of losses.
    tangents = similar(normals)
    sangents = similar(normals)
    tangents!(normals,tangents,sangents)

    # Setting the local-nodes of elements to be that of a regular triangle
    set_interpolation_nodes!(shape_function,gauss_points_triangle(n)...)
    copy_interpolation_nodes!(physics_function,shape_function)
    # Finalizing mesh
    mesh = Mesh3d(sources,coordinates,topology,normals,tangents,sangents,shape_function,
                 physics_function,physics_topology,ents[mask])
    if entites
        return mesh,ents[mask]
    else
        return mesh
    end
end
