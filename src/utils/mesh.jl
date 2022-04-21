#==========================================================================================
                            Defining mesh structs
==========================================================================================#
abstract type Mesh end

mutable struct Mesh2d{T} <: Mesh where {T <: AbstractArray}
    coordinates::T
    topology::AbstractArray{Int64,2}
    shape_function::ShapeFunction
end
mutable struct Mesh3d{T} <: Mesh where {T <: AbstractArray}
    sources::T
    coordinates::T
    topology::AbstractArray{Int64,2}
    normals::T
    tangents::T
    sangents::T
    shape_function::ShapeFunction
    physics_function::ShapeFunction
    physics_topology::AbstractArray{Int64,2}
end
#==========================================================================================
                            Helper Functions
==========================================================================================#
number_of_elements(mesh::Mesh)::Int64   = size(mesh.topology,2)
number_if_nodes(mesh::Mesh)::Int64      = size(mesh.coordinates,2)
get_element_type(mesh::Mesh)            = mesh.elementType
get_coordinates(mesh::Mesh)             = mesh.coordinates
get_topology(mesh::Mesh)                = mesh.topology
get_sources(mesh::Mesh3d)               = mesh.sources
get_normals(mesh::Mesh3d)               = mesh.normals
get_tangents(mesh::Mesh3d)              = mesh.tangents,mesh.sangents
#==========================================================================================
                            Printing symmary of mesh
==========================================================================================#
function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh2d)
    println(io, "Number of elements:\t $(size(mesh.topology,2))")
    println(io, "Number of unkowns:\t $(size(mesh.coordinates,2))")
    println(io, "Geometry defined by:\t $(typeof(mesh.shape_function))")
end
function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh3d)
    println(io, "Number of elements:\t $(size(mesh.topology,2))")
    println(io, "Number of unkowns:\t $(size(mesh.sources,2))")
    println(io, "Geometry defined by:\t $(typeof(mesh.shape_function))")
    println(io, "Physics defined by: \t $(typeof(mesh.physics_function))")
end