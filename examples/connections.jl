#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, IntegralEquations, Plots
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/sphere_1m"
tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])

# Standard output: Return all elements that are connected to each sources
cone,cons = IntegralEquations.connected_sources(mesh)
source = 3
sum(mesh.topology[:,cone[source]] .== source) == length(cone[source])

#
cone,cons = IntegralEquations.connected_sources(mesh,1)
source = 3
sum(mesh.topology[:,cone[source]] .== source) == length(cone[source])

##
interpolations = IntegralEquations.interpolate_elements(mesh;n=2,m=2)


function unroll_interpolations(interpolations)
    N = length(interpolations)
    Nweights = length(interpolations[1].jacobian_mul_weights)

    weights = zeros(N*Nweights)
    interps = zeros(3,N*Nweights)

    for i = 1:N
        weights[(i-1)*Nweights + 1:i*Nweights]   = interpolations[i].jacobian_mul_weights
        interps[:,(i-1)*Nweights + 1:i*Nweights] = interpolations[i].interpolation
    end

    return interps, weights
end

interp,weights = unroll_interpolations(interpolations)
