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
# Hmm, the linear physics fails. Here, but not if run seperately
@testset "Triangular Mesh" begin
    mesh_file = "../examples/meshes/sphere_1m"
    geometry_orders = [:linear,:quadratic]
    # physics_orders  = [:linear,:quadratic,:disctriconstant,:disctrilinear,:disctriquadratic]
    physics_orders  = [:linear,:quadratic]
    for go in geometry_orders, po in physics_orders
        if go == :linear && po == :quadratic
            continue
        end
        mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        k = 1.0
        HG = BoundaryIntegralEquations.HGOperator(mesh,k)
        x = ones(size(HG,2))
        y = HG*x
        xg = gmres(HG,y)
        println(go)
        println(po)
        @test norm(xg - x)/norm(x) ≈ 0 atol=1e-4
    end
end

# mesh_file = "../examples/meshes/sphere_1m"
# geometry_orders = [:linear,:quadratic]
# physics_orders  = [:linear,:quadratic,:disctriconstant,:disctrilinear,:disctriquadratic]
# go = geometry_orders[1]
# po = physics_orders[1]
# # for go in geometry_orders, po in physics_orders
#     # if go == :linear && po == :quadratic
#         # continue
#     # end
#     mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
#     # TODO
# # end
# k = 1.0
# HG = BoundaryIntegralEquations.HGOperator(mesh,k)
# x = ones(size(HG,2))
# y = HG*x
# xg = gmres(HG,y)

# @test norm(xg - x)/norm(x) ≈ 0 atol=1e-4
