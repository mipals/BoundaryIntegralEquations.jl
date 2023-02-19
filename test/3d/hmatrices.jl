#==========================================================================================
                                Using relevant packages
==========================================================================================#
using Test
using BoundaryIntegralEquations
using LinearAlgebra
using IterativeSolvers
#==========================================================================================
                            Testing mesh interpolation schemes
==========================================================================================#
@testset "Triangular Mesh" begin
    mesh_file = "../examples/meshes/sphere_1m_extremely_fine"
    geometry_orders = [:linear,:quadratic]
    physics_orders  = [:linear,:quadratic,:disctriconstant,:disctrilinear,:disctriquadratic]
    for go in geometry_orders, po in physics_orders
        if go == :linear && po == :quadratic
            continue
        end
        mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        # TODO
    end
end
