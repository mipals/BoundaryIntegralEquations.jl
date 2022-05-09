import IntegralEquations: get_beta_tri_linear, get_beta_tri_quadratic

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

@testset "Discontinuous Triangular Constant Shape Functions" begin
    quadratic_triangular = TriangularQuadratic(3,3)
    disc_constant_tri = DiscontinuousTriangularConstant(quadratic_triangular)
    # Weights should add to area of quad
    @test sum(disc_constant_tri.weights) ≈ 0.5 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_constant_tri)
    # Checking if nodal set_nodal_interpolation! works as intended
    @test disc_constant_tri.interpolation ≈ I
end

@testset "Discontinuous Triangular Linear Shape Functions" begin
    quadratic_triangular = TriangularQuadratic(3,3)
    beta = get_beta_tri_linear(:legendre)
    disc_linear_tri = DiscontinuousTriangularLinear(quadratic_triangular,beta)
    # Weights should add to area of quad
    @test sum(disc_linear_tri.weights) ≈ 0.5 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_linear_tri)
    # Checking if nodal set_nodal_interpolation! works as intended
    @test disc_linear_tri.interpolation ≈ I
end

@testset "Discontinuous Triangular Quadratic Shape Functions" begin
    quadratic_triangular = TriangularQuadratic(3,3)
    beta = get_beta_tri_quadratic(:legendre)
    disc_linear_tri = DiscontinuousTriangularQuadratic(quadratic_triangular,beta)
    # Weights should add to area of quad
    @test sum(disc_linear_tri.weights) ≈ 0.5 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_linear_tri)
    # Checking if nodal set_nodal_interpolation! works as intended
    @test disc_linear_tri.interpolation ≈ I
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
    quadratic_quad = QuadrilateralQuadraticLagrange(3,3)
    # Weights should add to area of quad
    @test sum(quadratic_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(quadratic_quad)
    # Checking if nodal set_nodal_interpolation! works as intended 
    @test quadratic_quad.interpolation ≈ I
end

@testset "Discontinuous Quadrilateral Constant Shape Functions" begin
    disc_constant_quad = DiscontinuousQuadrilateralConstant(4,4)
    # Weights should add to area of quad
    @test sum(disc_constant_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_constant_quad)
    # Checking if nodal set_nodal_interpolation! works as intended 
    @test disc_constant_quad.interpolation ≈ I
end
@testset "Discontinuous Quadrilateral Linear Shape Functions" begin
    disc_linear_quad = DiscontinuousQuadrilateralLinear4(4,4)
    # Weights should add to area of quad
    @test sum(disc_linear_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_linear_quad)
    # Checking if nodal set_nodal_interpolation! works as intended 
    @test disc_linear_quad.interpolation ≈ I
end
@testset "Discontinuous Quadrilateral Quadratic Shape Functions" begin
    disc_quadratic_quad = DiscontinuousQuadrilateralQuadraticLagrange(4,4)
    # Weights should add to area of quad
    @test sum(disc_quadratic_quad.weights) ≈ 4.0 
    # Set interpolation to be on the interpolating nodes
    set_nodal_interpolation!(disc_quadratic_quad)
    # Checking if nodal set_nodal_interpolation! works as intended 
    @test disc_quadratic_quad.interpolation ≈ I
end
