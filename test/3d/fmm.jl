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
        interpolations = BoundaryIntegralEquations.interpolate_elements(mesh;n=2,m=2);
        interp,weights,normals = BoundaryIntegralEquations.unroll_interpolations(interpolations);
        if go == :linear
            @test 4π ≈ sum(weights) atol=1e-1
        else
            @test 4π ≈ sum(weights) atol=1e-5
        end
    end
end
#==========================================================================================
            Testing the double-layer potentential for scattering of a sphere
==========================================================================================#
# @testset "BEM vs FMM scattering" begin
#     mesh_file = "../examples/meshes/sphere_1m"
#     geometry_orders = [:linear,:quadratic]
#     physics_orders  = [:linear,:quadratic,:disctriconstant,:disctrilinear,:disctriquadratic]
#     for go in geometry_orders, po in physics_orders
#         if go == :linear && po == :quadratic
#             continue
#         end
#         mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
#         zk = 1.0
#         radius = 1.0                                    # Radius of sphere_1m       [m]
#         # Computing incident pressure
#         angles = [π/2 0.0]
#         pI = BoundaryIntegralEquations.incoming_wave(angles,1.0,mesh.sources,zk)
#         Af = FMMHOperator(mesh,zk)
#         @time Fp,_,Cp = assemble_parallel!(mesh,zk,mesh.sources,n=2,m=2,gOn=false,sparse=false);
#         Ap = Fp + Diagonal(1.0 .- Cp);
#         # Testing if the multiplication is the same
#         @test maximum(abs.(Af*x - Ap*x)) atol=1e-3
#         # Checking if the solution of the system is the same
#         p_bem = gmres(Ap,pI)
#         p_fmm = gmres(Af,pI)
#         @test p_bem ≈ p_fmm atol=1e-2
#     end
# end
