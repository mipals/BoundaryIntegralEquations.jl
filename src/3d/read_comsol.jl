#==========================================================================================
                            Triangular Meshes
==========================================================================================#
"""
    read_comsol_mesh(meshName)

Loads the coordinates, topology and entites from a ".mphtxt" file.
"""
function read_comsol_mesh(meshName,elementType)
    # Number of shape functions. Needed for the size of the topology
    nShape = number_of_shape_functions(elementType)
    # Opening COMSOL mesh file
    fid = open(meshName * ".mphtxt")
    # Initializing output values
    coordinates = zeros(3,0)
    topology    = zeros(Int64,nShape,0)
    entities    = zeros(Int64,1,0)
    for line in eachline(fid)
        if line == "# Mesh vertex coordinates"
            newLine = readline(fid)
            while !isempty(newLine)
                coordinates = [coordinates parse.(Float64,split(newLine))]
                newLine = readline(fid)
            end
        end
        if line == "# Elements"
            newLine  = readline(fid)
            while !isempty(newLine)
                topology = [topology parse.(Int64,split(newLine))]
                newLine = readline(fid)
            end
        end
        if line == "# Geometric entity indices"
            newLine = readline(fid)
            while !isempty(newLine)
                entities = [entities parse.(Int64,newLine)]
                newLine = readline(fid)
            end
        end
    end

    close(fid)
    # Transform COMSOLs 0-indexing to Julias 1-indexing
    topology = topology .+ 1
    if typeof(elementType) <: TriangularQuadratic
        # Fixing COMSOLs weird triangle layout
        topology[5:6,:] = topology[6:-1:5,:]
    end
    return coordinates,topology,entities
end

function load3dTriangularComsolMesh(meshFile;m=3,n=3,
                geometry_order=:quadratic,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removedEntites=[-1])

    # Figuring out the element type
    initialCoordinates,initialTopology,ents = read_comsol_mesh(meshFile,TriangularQuadratic(2,2))
    if geometry_order == :quadratic
        shape_function = TriangularQuadratic(m,n)
    elseif geometry_order == :linear
        shape_function = TriangularLinear(m,n)
        initialTopology = initialTopology[1:3,:]
    else
        error("Only quadratic and linear geometries are currently supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removedEntites,dims=1))[:]
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

function load3dQuadComsolMesh(meshFile;m=4,n=4,
                geometry_order=:quadratic,
                physics_order=:geometry,beta_type=:legendre,
                entites=false,removedEntites=[-1])
    # Figuring out the element type
    initialCoordinates,initialTopology,ents = read_comsol_mesh(meshFile,QuadrilateralQuadraticLagrange(2,2))
    # Checking input
    if geometry_order == :quadratic
        shape_function = QuadrilateralQuadraticLagrange(m,n)
    elseif geometry_order == :linear
        shape_function = QuadrilateralLinear4(m,n)
        initialTopology = initialTopology[1:4,:]
    else
        error("Only QuadQuadratic and QuadLinear geometry types supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removedEntites,dims=1))[:]
    usedNodes   = sort(unique(initialTopology[:,mask]))
    topology    = remove_unused_nodes(initialTopology[:,mask])
    coordinates = initialCoordinates[:,usedNodes]
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

#==========================================================================================
                                 Utility functions
==========================================================================================#
"""
    load_comsol_results(file)

Reads comsol results from file.
"""
function load_comsol_results(file)
    fid = open(file)
    coordinates = zeros(3,0)
    pLossy = zeros(ComplexF64,0,1)
    for line in eachline(fid)
            splitLine = split(line)
            coordinates = [coordinates parse.(Float64,splitLine[1:3])]
            pLossy      = [pLossy; parse.(ComplexF64,splitLine[end])]
    end
    close(fid)
    return coordinates,pLossy[:]
end
"""
    write_sources(file,sources)

Writes `sources` to file.
"""
function write_sources(file,sources)
    open(file, "w") do io
        writedlm(io, sources)
    end
end


function find_boundary_conditions(mesh,ents,boundaryConditionEntityRes)
    entsBC = ents .== boundaryConditionEntityRes
    return entsBC
end
