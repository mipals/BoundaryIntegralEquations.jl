import GeometryBasics

"""
    load_mesh_file(file)

Loads 3D mesh files.
"""
function load_mesh_file(file)
    # Read mesh
    mesh = load(file)
    # Extract coordinates and connectivities
    mesh_coords = GeometryBasics.coordinates(mesh)
    elements    = GeometryBasics.faces(mesh)
    # Get number of elements and number of nodes on the element (element order)
    n_elements  = length(elements)
    n_interpolations = length(elements[1])
    # Preallocate topology and coordinates
    topology = zeros(Int64,n_interpolations,n_elements)
    coords   = zeros(3,length(mesh_coords))
    # Insert connectivities into topology matrix
    for (index,element) in enumerate(elements)
        for i = 1:n_interpolations
            topology[i,index] = Int64(element[i].i) + 1
        end
    end
    # Extract coordinates
    for (index,coord) in enumerate(mesh_coords)
        coords[:,index] = Float64.(coord)
    end
    # Make dummy entities
    ents = zeros(Int64, n_elements)
    # Return coordiates, topology and ents
    return coords, topology, ents
end

"""
    load3dTriangularMesh(mesh_file;m=3,n=3,geometry_order=:linear,physics_order=:geometry,
                                    beta_type=:legendre,entites=false,removed_entites=[-1])

Returns a `Mesh3d` from the geoemtry described in `mesh_file`.

Supports .obj, .ply, .stl, .off, .2DM files using `MeshIO.jl` and `FileIO.jl`.
"""
function load3dTriangularMesh(mesh_file;m=3,n=3,
                geometry_order=:linear,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removed_entites=[-1])

    if geometry_order !== :linear
        error("Currently only linear element supported")
    end

    # Figuring out the element type
    initialCoordinates,initialTopology,ents = load_mesh_file(mesh_file)
    if geometry_order == :quadratic
        shape_function = TriangularQuadratic(m,n)
    elseif geometry_order == :linear
        shape_function = TriangularLinear(m,n)
        initialTopology = initialTopology[1:3,:]
    else
        error("Only quadratic and linear geometries are currently supported")
    end

    mask = .!convert.(Bool,sum(ents' .∈ removed_entites,dims=1))[:]
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


"""
    load3dTriangularGmshMesh(gmshfile;m=3,n=3, geometry_order=:quadratic,
        physics_order=:geometry,beta_type=:legendre,entites=false,removed_entites=[-1])

"""
function load3dTriangularGmshMesh(gmshfile;m=3,n=3,
                renumber=false,set_order = false,refine_multiple=0,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removed_entites=[-1])

    nDim = 3
    gmsh.initialize()
    gmsh.open(gmshfile)
    # Renumber nodes
    # Generate (if only a .geo file is passed)
    gmsh.model.mesh.generate(2)
    renumber && gmsh.model.mesh.renumber_nodes()
    renumber && gmsh.model.mesh.renumber_elements()
    for i = 1:refine_multiple
        gmsh.model.mesh.refine()
    end
    if set_order ∈ [1,2]
        gmsh.model.mesh.set_order(set_order)
    end
    gmsh.model.mesh.optimize("HighOrderElastic")
    gmsh.model.mesh.optimize("HighOrder")
    # Load elements
    element_types, elementTags, nodeTags = gmsh.model.mesh.get_elements(2)
    if length(element_types) == 0
        error("No elements found.")
    end
    element_type = element_types[1]
    if !(element_type in [2,9])
        error("Element type not Triangular.")
    end
    if element_type == 2
        geometry_order = :linear
        nodes_per_element = 3
    elseif element_type == 9
        geometry_order = :quadratic
        nodes_per_element = 6
    end
    # Load all nodes
    _, nodes, _ = gmsh.model.mesh.get_nodes();
    # Compute number of nodes
    number_of_nodes = div(length(nodes),nDim);
    # Reshaping into the desired format
    coordinates = reshape(nodes,nDim,number_of_nodes);
    intialTopology    = Int.(reshape(nodeTags[1],nodes_per_element,
                        div(length(nodeTags[1]),nodes_per_element)))
    initial_topology = remove_unused_nodes(intialTopology)

    # Finding unused nodes
    usedNodes = [used for used in 1:number_of_nodes if (used in nodeTags[1])]
    initial_coordinates = coordinates[:,usedNodes]
    # Make dummy ents
    ents = zeros(Int64, size(initial_topology,2))

    gmsh.finalize()

    if geometry_order == :quadratic
        shape_function = TriangularQuadratic(m,n)
    elseif geometry_order == :linear
        shape_function = TriangularLinear(m,n)
        initial_topology = initial_topology[1:3,:]
    else
        error("Only quadratic and linear geometries are currently supported")
    end

    mask = .!convert.(Bool,sum(ents' .∈ removed_entites,dims=1))[:]
    used_nodes  = sort(unique(initial_topology[:,mask]))
    topology    = remove_unused_nodes(initial_topology[:,mask])
    coordinates = initial_coordinates[:,used_nodes]
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
