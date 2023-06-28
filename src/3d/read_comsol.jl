#==========================================================================================
                            Triangular Meshes
==========================================================================================#
"""
    read_comsol_mesh(mesh_file)

Loads the coordinates, topology and entites from a ".mphtxt" file.
"""
function read_comsol_mesh(mesh_file)
    # Opening COMSOL mesh file
    fid = open(mesh_file * ".mphtxt")
    # Initializing output values
    coordinates = Vector{Vector{Float64}}(undef,0)
    topology    = Vector{Vector{Int64}}(undef,0)
    entities    = Vector{Int64}(undef,0)
    for line in eachline(fid)
        if length(line) < 4
            continue
        end
        if line == "# Mesh vertex coordinates"
            new_line = readline(fid)
            while !isempty(new_line)
                coordinates = push!(coordinates, parse.(Float64,split(new_line)))
                new_line = readline(fid)
            end
        end
        if line == "# Elements"
            new_line  = readline(fid)
            while !isempty(new_line)
                topology = push!(topology, parse.(Int64,split(new_line)))
                new_line = readline(fid)
            end
        end
        if line == "# Geometric entity indices"
            new_line = readline(fid)
            while !isempty(new_line)
                entities = push!(entities, parse.(Int64,new_line))
                new_line = readline(fid)
            end
        end
    end

    close(fid)
    # Concatination the coordinates and topology
    coordinates = reduce(hcat,coordinates)
    topology    = reduce(hcat,topology)
    # Transform COMSOLs 0-indexing to Julias 1-indexing
    topology = topology .+ 1
    if size(topology,1) == 6
        # Fixing COMSOLs weird triangle layout
        topology[5:6,:] = topology[6:-1:5,:]
    end
    return coordinates,topology,entities'
end

"""
    load3dTriangularComsolMesh(mesh_file;m=3,n=3, geometry_order=:quadratic,
        physics_order=:geometry,beta_type=:legendre,entites=false,removed_entites=[-1])

"""
function load3dTriangularComsolMesh(mesh_file;m=3,n=3,
                geometry_order=:quadratic,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removed_entites=[-1])

    # Figuring out the element type
    initial_coordinates,initial_topology,ents = read_comsol_mesh(mesh_file)
    if geometry_order == :quadratic
        shape_function = TriangularQuadratic(m,n)
    elseif geometry_order == :linear
        shape_function = TriangularLinear(m,n)
        initial_topology = initial_topology[1:3,:]
    else
        error("Only quadratic and linear geometries are currently supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removed_entites,dims=1))[:]
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

"""
    load3dQuadComsolMesh(mesh_file;m=4,n=4,geometry_order=:quadratic,
            physics_order=:geometry,beta_type=:legendre,entites=false,removed_entites=[-1])

"""
function load3dQuadComsolMesh(mesh_file;m=4,n=4,
                geometry_order=:quadratic,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removed_entites=[-1])
    # Figuring out the element type
    initial_coordinates,initial_topology,ents = read_comsol_mesh(mesh_file)
    # Checking input
    if geometry_order == :quadratic
        shape_function = QuadrilateralQuadraticLagrange(m,n)
    elseif geometry_order == :linear
        shape_function = QuadrilateralLinear4(m,n)
        initial_topology = initial_topology[1:4,:]
    else
        error("Only QuadQuadratic and QuadLinear geometry types supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removed_entites,dims=1))[:]
    used_nodes  = sort(unique(initial_topology[:,mask]))
    topology    = remove_unused_nodes(initial_topology[:,mask])
    coordinates = initial_coordinates[:,used_nodes]
    sources = coordinates[:,sort(unique(topology))]

    if physics_order == :geometry
        physics_function = set_physics_element(physics_order,shape_function,beta_type)
        normals = get_element_normals(shape_function,coordinates,topology)
        physics_topology = topology
    elseif physics_order == :linear
        physics_function = QuadrilateralLinear4(3,3)
        copy_interpolation_nodes!(physics_function,shape_function)
        normals = get_element_normals(shape_function,coordinates,topology)
        physics_topology = topology[1:4,:]
        used_nodes = sort(unique(topology[1:4,:]))
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

    # Finalizing mesh
    mesh = Mesh3d(sources,coordinates,topology,normals,tangents,sangents,shape_function,
                 physics_function,physics_topology,ents[mask])
    if entites
        return mesh,ents[mask]
    else
        return mesh
    end
end
