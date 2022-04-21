#==========================================================================================
                            Triangular Meshes
==========================================================================================#
"""
    read_comsol_mesh(meshName)

Loads the coordinates, topology and entites from a ".mphtxt" file.
"""
function read_comsol_mesh(meshName,elementType)
    nShape = number_of_shape_functions(elementType)
    fid = open(meshName * ".mphtxt")
    coordinates = zeros(3,0)
    topology = zeros(Int64,nShape,0)
    entities = zeros(Int64,1,0)
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
    topology = topology .+ 1             # Fixing 0-index
    if typeof(elementType) <: TriangularQuadratic
        topology[5:6,:] = topology[6:-1:5,:] # Fixing COMSOLs weird triangle layout
        topology[2:3,:] = topology[3:-1:2,:]
        topology[4:6,:] = topology[6:-1:4,:]
    end
    return coordinates,topology,entities
end

function load3dTriangularComsolMesh(meshFile;m=4,n=4,
                geometryType="TriQuadratic",
                physicsType="geometry",betaType="Legendre",
                entites=false,removedEntites=[-1])
    # Figuring out the element type
    initialCoordinates,initialTopology,ents = read_comsol_mesh(meshFile,TriangularQuadratic(2,2))
    if lowercase(geometryType) == "triquadratic"
        basisElement = TriangularQuadratic(m,n)
    elseif lowercase(geometryType) == "trilinear"
        basisElement = TriangularLinear(m,n)
        initialTopology = initialTopology[1:3,:]
    else
        error("Only TriQuadratic and TriLinear geometry types supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removedEntites,dims=1))[:]
    usedNodes   = sort(unique(initialTopology[:,mask]))
    topology    = remove_unused_nodes(initialTopology[:,mask])
    coordinates = initialCoordinates[:,usedNodes]
    sources = coordinates[:,sort(unique(topology))]

    println("")
    # Compute normals
    # physicsElement = set_physics_element(physicsType,basisElement,betaType)
    if lowercase(physicsType) == "geometry"
        physicsElement = set_physics_element(physicsType,basisElement,betaType)
        normals = get_element_normals(basisElement,coordinates,topology)
        physicsTopology = topology
    elseif lowercase(physicsType) == "linear"
        physicsElement = TriangularLinear(3,3)
        normals = get_element_normals(basisElement,coordinates,topology)
        physicsTopology = topology[1:3,:]
        normals = normals[:,unique(topology[1:3,:])]
        sources = coordinates[:,unique(topology[1:3,:])]
    else
        physicsElement = set_physics_element(physicsType,basisElement,betaType)
        sources,normals,physicsTopology = compute_sources(basisElement,
                                                        physicsElement,
                                                        topology,
                                                        coordinates)
    end

    # Allocate and compute tangent directions from normal
    tangentX = similar(normals)
    tangentY = similar(normals)
    tangents!(normals,tangentX,tangentY)
    
    # Create mesh
    mesh = Mesh3d(sources,coordinates,topology,normals,tangentX,tangentY,basisElement,
                 physicsElement,physicsTopology)
    if entites
        return mesh,ents[mask]
    else
        return mesh
    end
end

function load3dQuadComsolMesh(meshFile;m=4,n=4,
                geometryType="QuadQuadratic",
                physicsType="geometry",betaType="Legendre",
                entites=false,removedEntites=[-1])
    # Figuring out the element type
    initialCoordinates,initialTopology,ents = read_comsol_mesh(meshFile,QuadrilateralQuadratic9(2,2))
    if lowercase(geometryType) == "quadquadratic"
        basisElement = QuadrilateralQuadratic9(m,n)
    elseif lowercase(geometryType) == "quadlinear"
        basisElement = QuadrilateralLinear4(m,n)
        initialTopology = initialTopology[1:4,:]
    else
        error("Only QuadQuadratic and QuadLinear geometry types supported")
    end

    mask = .!convert.(Bool,sum(ents .∈ removedEntites,dims=1))[:]
    usedNodes   = sort(unique(initialTopology[:,mask]))
    topology    = remove_unused_nodes(initialTopology[:,mask])
    coordinates = initialCoordinates[:,usedNodes]
    sources = coordinates[:,sort(unique(topology))]

    println("")
    # Compute normals
    # physicsElement = set_physics_element(physicsType,basisElement,betaType)
    if lowercase(physicsType) == "geometry"
        physicsElement = set_physics_element(physicsType,basisElement,betaType)
        normals = get_element_normals(basisElement,coordinates,topology)
        physicsTopology = topology
    elseif lowercase(physicsType) == "linear"
        physicsElement = QuadrilateralLinear4(3,3)
        normals = get_element_normals(basisElement,coordinates,topology)
        physicsTopology = topology[1:4,:]
        normals = normals[:,unique(topology[1:4,:])]
        sources = coordinates[:,unique(topology[1:4,:])]
    else
        physicsElement = set_physics_element(physicsType,basisElement,betaType)
        sources,normals,physicsTopology = compute_sources(basisElement,
                                                        physicsElement,
                                                        topology,
                                                        coordinates)
    end

    # Allocate and compute tangent directions from normal
    tangentX = similar(normals)
    tangentY = similar(normals)
    tangents!(normals,tangentX,tangentY)
    
    # Create mesh
    mesh = Mesh3d(sources,coordinates,topology,normals,tangentX,tangentY,basisElement,
                 physicsElement,physicsTopology)
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
