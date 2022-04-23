@testset "Linear Triangular Shape Functions" begin
    linear_triangular = TriangularLinear(3,3)
    # Weights should add to area of triangle
    @test sum(linear_triangular.weights) ≈ 0.5 
    # Constructing a linear element with same gauss nodes as linear_triangular
    tri_lin  = TriangularLinear(linear_triangular)
    @test tri_lin.weights == linear_triangular.weights
    @test tri_lin.gauss_u == linear_triangular.gauss_u
    @test tri_lin.gauss_v == linear_triangular.gauss_v
    @test tri_lin.derivatives_u == linear_triangular.derivatives_u
    @test tri_lin.derivatives_v == linear_triangular.derivatives_v
    # Constructing a quadratic element with same gauss nodes as linear_triangular
    tri_quad = TriangularQuadratic(linear_triangular)
    @test tri_quad.weights == linear_triangular.weights
    @test tri_quad.gauss_u == linear_triangular.gauss_u
    @test tri_quad.gauss_v == linear_triangular.gauss_v
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(linear_triangular)
    @test linear_triangular.interpolation ≈ I 
end

@testset "Quadratic Triangular Shape Functions" begin
    quadratic_triangular = TriangularQuadratic(3,3)
    # Weights should add to area of triangle
    @test sum(quadratic_triangular.weights) ≈ 0.5 
    # Constructing a linear element with same gauss nodes as linear_triangular
    tri_lin  = TriangularLinear(quadratic_triangular)
    @test tri_lin.weights == quadratic_triangular.weights
    @test tri_lin.gauss_u == quadratic_triangular.gauss_u
    @test tri_lin.gauss_v == quadratic_triangular.gauss_v
    # Constructing a quadratic element with same gauss nodes as linear_triangular
    tri_quad = TriangularQuadratic(quadratic_triangular)
    @test tri_quad.weights == quadratic_triangular.weights
    @test tri_quad.gauss_u == quadratic_triangular.gauss_u
    @test tri_quad.gauss_v == quadratic_triangular.gauss_v
    @test tri_quad.derivatives_u == quadratic_triangular.derivatives_u
    @test tri_quad.derivatives_v == quadratic_triangular.derivatives_v
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(quadratic_triangular)
    @test quadratic_triangular.interpolation ≈ I
end

@testset "Linear Quadrilateral Shape Functions" begin
    linear_quad = QuadrilateralLinear4(3,3)
    # Weights should add to area of quad
    @test sum(linear_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(linear_quad)
    # Checking if nodal set_nodal_interpolation! works as intended
    @test linear_quad.interpolation ≈ I
end

@testset "Quadratic Quadrilateral Shape Functions" begin
    quadratic_quad = QuadrilateralQuadratic9(3,3)
    # Weights should add to area of quad
    @test sum(quadratic_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(quadratic_quad)
    # Checking if nodal set_nodal_interpolation! works as intended 
    @test quadratic_quad.interpolation ≈ I
end
