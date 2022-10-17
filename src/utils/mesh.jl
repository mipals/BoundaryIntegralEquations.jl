#==========================================================================================
                            Defining mesh structs
==========================================================================================#
abstract type Mesh end

struct Mesh2d{T} <: Mesh where {T <: AbstractArray}
    sources::T                                  # Collocation nodes
    coordinates::T                              # Mesh coordinates
    topology::AbstractArray{Int64,2}            # Element connectivities
    normals::T                                  # Normal at collocation nodes
    tangents::T                                 # Tangent at collocation nodes
    shape_function::ShapeFunction               # Defines Element Interpolation
    physics_function::ShapeFunction             # Defines Physics Interpolation
    physics_topology::AbstractArray{Int64,2}    # Physics connectivities
    # entites                                    # should probably be added at some point
end
struct Mesh3d{T} <: Mesh where {T <: AbstractArray}
    sources::T                                  # Collocation nodes
    coordinates::T                              # Mesh coordinates
    topology::AbstractArray{Int64,2}            # Element connectivities
    normals::T                                  # Normal at collocation nodes
    tangents::T                                 # First tangent at collocation nodes
    sangents::T                                 # Second tangent at collocation nodes
    shape_function::ShapeFunction               # Defines Element Interpolation
    physics_function::ShapeFunction             # Defines Physics Interpolation
    physics_topology::AbstractArray{Int64,2}    # Physics connectivities
    entities                                    # COMSOL numbering of faces (for BCs)
end
#==========================================================================================
                            Helper Functions
==========================================================================================#
number_of_elements(mesh::Mesh)::Int64   = size(mesh.topology,2)
number_of_nodes(mesh::Mesh)::Int64      = size(mesh.coordinates,2)
get_element_type(mesh::Mesh)            = mesh.elementType
get_coordinates(mesh::Mesh)             = mesh.coordinates
get_topology(mesh::Mesh)                = mesh.topology
get_sources(mesh::Mesh3d)               = mesh.sources
get_normals(mesh::Mesh3d)               = mesh.normals
get_tangents(mesh::Mesh3d)              = mesh.tangents,mesh.sangents
get_entities(mesh::Mesh3d)              = mesh.entities
#==========================================================================================
                            Printing symmary of mesh
==========================================================================================#
function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh2d)
    println(io, "Number of elements: \t $(size(mesh.topology,2))")
    println(io, "Number of unkowns:  \t $(size(mesh.sources,2))")
    println(io, "Geometry defined by:\t $(typeof(mesh.shape_function))")
    println(io, "Physics defined by: \t $(typeof(mesh.physics_function))")
end
function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh3d)
    println(io, "Number of elements: \t $(size(mesh.topology,2))")
    println(io, "Number of unkowns:  \t $(size(mesh.sources,2))")
    println(io, "Geometry defined by:\t $(typeof(mesh.shape_function))")
    println(io, "Physics defined by: \t $(typeof(mesh.physics_function))")
end
