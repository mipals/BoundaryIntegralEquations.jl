import IntegralEquations: tangents!

@testset "Mesh Functions" begin
    normals  = [[1.0 0.0 0.0]; [0.0 1.0 0.0]; [0.0 0.0 1.0]]
    tangents = similar(normals)
    sangents = similar(normals)
    tangents!(normals,tangents,sangents)
    # Testing Orthogonality of the nodal (normal,tangent,sangent)-system(s)
    @test ones(1,3)*(normals  .* tangents) ≈ zeros(1,3)
    @test ones(1,3)*(normals  .* sangents) ≈ zeros(1,3)
    @test ones(1,3)*(sangents .* tangents) ≈ zeros(1,3)
end

@testset "Triangular Mesh" begin
    mesh_file = "meshes/sphere_1m"
    # Mesh with continuous triangular quadratic geometry and physics
    mesh = load3dTriangularComsolMesh(mesh_file)
    @test size(mesh.topology) == (6,840)
    @test size(mesh.sources)  == (3,1682)
    @test size(mesh.normals)  == size(mesh.tangents)
    @test size(mesh.tangents) == size(mesh.sangents)
    # Mesh with continuous triangular linear geometry and physics
    mesh_tri_lin = load3dTriangularComsolMesh(mesh_file;geometryType="TriLinear")
    @test size(mesh_tri_lin.topology) == (3,840)
    @test size(mesh_tri_lin.sources)  == (3,422)
    @test size(mesh_tri_lin.normals)  == size(mesh_tri_lin.tangents)
    @test size(mesh_tri_lin.tangents) == size(mesh_tri_lin.sangents)
    @test size(mesh_tri_lin.topology) == size(mesh_tri_lin.physics_topology) 
    # Mesh with continuous triangular quadratic geometry and triangular linear physics
    mesh_phys_lin = load3dTriangularComsolMesh(mesh_file;physicsType="linear")
    @test size(mesh_phys_lin.topology)   == (6,840)
    @test size(mesh_phys_lin.sources)    == (3,422)
    @test size(mesh_phys_lin.normals)    == size(mesh_phys_lin.tangents)
    @test size(mesh_phys_lin.tangents)   == size(mesh_phys_lin.sangents)
    @test size(mesh_phys_lin.topology,2) == size(mesh_phys_lin.physics_topology,2) 
end

@testset "Quadrilateral Mesh" begin
    mesh_file = "meshes/quad_cylinder"
    # Mesh with continuous quadrilateral quadratic geometry and physics
    mesh = load3dQuadComsolMesh(mesh_file)
    @test size(mesh.topology) == (9,2360)
    @test size(mesh.sources)  == (3,9442)
    @test size(mesh.normals)  == size(mesh.tangents)
    @test size(mesh.tangents) == size(mesh.sangents)
    # Mesh with continuous quadrilateral linear geometry and physics
    mesh_quad_lin = load3dQuadComsolMesh(mesh_file;geometryType="quadlinear")
    @test size(mesh_quad_lin.topology) == (4,2360)
    @test size(mesh_quad_lin.sources)  == (3,2362)
    @test size(mesh_quad_lin.normals)  == size(mesh_quad_lin.tangents)
    @test size(mesh_quad_lin.tangents) == size(mesh_quad_lin.sangents)
    @test size(mesh_quad_lin.topology) == size(mesh_quad_lin.physics_topology) 
    # Mesh with continuous quadrilateral quadratic geometry and triangular linear physics
    mesh_phys_lin = load3dQuadComsolMesh(mesh_file;physicsType="linear")
    @test size(mesh_phys_lin.topology)   == (9,2360)
    @test size(mesh_phys_lin.sources)    == (3,2362)
    @test size(mesh_phys_lin.normals)    == size(mesh_phys_lin.tangents)
    @test size(mesh_phys_lin.tangents)   == size(mesh_phys_lin.sangents)
    @test size(mesh_phys_lin.topology,2) == size(mesh_phys_lin.physics_topology,2) 
end

