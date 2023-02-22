#==========================================================================================
                                Using relevant packages
==========================================================================================#
using BoundaryIntegralEquations
using LinearAlgebra
using Test
import BoundaryIntegralEquations: tangents!, number_of_shape_functions
import BoundaryIntegralEquations: interpolation_function_derivatives
#==========================================================================================
                                    Creating Tests
==========================================================================================#
@testset "Triangular Mesh" begin
    mesh_file = "../examples/meshes/sphere_1m" #! Relative path from the "runtests.jl" file
    geometry_orders = [:linear,:quadratic]
    physics_orders  = [:linear,:quadratic,:disctriconstant,:disctrilinear,:disctriquadratic]
    for go in geometry_orders, po in physics_orders
        if go == :linear && po == :quadratic
            continue
        end
        mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        #? The elements described by 'topology' should match the number of shape functions
        @test size(mesh.topology)   == (number_of_shape_functions(mesh.shape_function),840)
        #? The number of sources should match the number of physics nodes
        @test size(mesh.sources)    == (3,length(unique(mesh.physics_topology)))
        #? Checking if the number of geometrical elements is the same as physics elements
        @test size(mesh.topology,2) == size(mesh.physics_topology,2)
        #? Checking if the sizes of the normals and tangents are the same
        @test size(mesh.normals)    == size(mesh.tangents)
        @test size(mesh.tangents)   == size(mesh.sangents)
        #? Checking for orthogonality between the normals and the tangents
        @test isapprox(sum(ones(1,3)*(mesh.normals.*mesh.tangents)),  0.0, atol=1e-14)
        @test isapprox(sum(ones(1,3)*(mesh.normals.*mesh.tangents)),  0.0, atol=1e-14)
        @test isapprox(sum(ones(1,3)*(mesh.tangents.*mesh.sangents)), 0.0, atol=1e-14)
        #? Checking of constant pressure will result in zero derivatives
        Dx,Dy,Dz = interpolation_function_derivatives(mesh)
        n_sources = size(mesh.sources,2)
        @test isapprox(sum(abs.(Dx*ones(n_sources)))/n_sources, 0.0, atol=1e-14)
        @test isapprox(sum(abs.(Dy*ones(n_sources)))/n_sources, 0.0, atol=1e-14)
        @test isapprox(sum(abs.(Dz*ones(n_sources)))/n_sources, 0.0, atol=1e-14)
        if go == :linear && po == :linear
            @test mesh.coordinates ≈ mesh.sources
        end
    end
end

@testset "Quadrilateral Mesh" begin
    mesh_file = "../examples/meshes/quad_sphere" #! Relative path from the "runtests.jl" file
    # Mesh with continuous quadrilateral quadratic geometry and physics
    geometry_orders = [:linear,:quadratic]
    physics_orders  = [:linear,:quadratic,:discquadconstant,:discquadlinear,:discquadquadratic]
    for go in geometry_orders, po in physics_orders
        if go == :linear && po == :quadratic
            continue
        end
        mesh = load3dQuadComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        #? The elements described by 'topology' should match the number of shape functions
        @test size(mesh.topology)   == (number_of_shape_functions(mesh.shape_function),315)
        #? The number of sources should match the number of physics nodes
        @test size(mesh.sources)    == (3,length(unique(mesh.physics_topology)))
        #? Checking if the number of geometrical elements is the same as physics elements
        @test size(mesh.topology,2) == size(mesh.physics_topology,2)
        #? Checking if the sizes of the normals and tangents are the same
        @test size(mesh.normals)    == size(mesh.tangents)
        @test size(mesh.tangents)   == size(mesh.sangents)
        #? Checking for orthogonality between the normals and the tangents
        @test isapprox(sum(ones(1,3)*(mesh.normals.*mesh.tangents)),  0.0, atol=1e-14)
        @test isapprox(sum(ones(1,3)*(mesh.normals.*mesh.tangents)),  0.0, atol=1e-14)
        @test isapprox(sum(ones(1,3)*(mesh.tangents.*mesh.sangents)), 0.0, atol=1e-14)
        #? Checking of constant pressure will result in zero derivatives
        Dx,Dy,Dz = interpolation_function_derivatives(mesh)
        n_sources = size(mesh.sources,2)
        @test isapprox(sum(abs.(Dx*ones(n_sources))/n_sources), 0.0, atol=1e-14)
        @test isapprox(sum(abs.(Dy*ones(n_sources))/n_sources), 0.0, atol=1e-14)
        @test isapprox(sum(abs.(Dz*ones(n_sources))/n_sources), 0.0, atol=1e-14)
    end
end
