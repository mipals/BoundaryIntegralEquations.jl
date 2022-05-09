#==========================================================================================
                                Continuous Segments
==========================================================================================#
@testset "Continuous Linear Segments" begin
    con_curve_linear = ContinuousCurveLinear(3)
    # Weights should add to length of line segment of triangle
    @test sum(con_curve_linear.weights) ≈ 2.0
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(con_curve_linear)
    @test con_curve_linear.interpolation ≈ I 
end
@testset "Continuous Quadratic Segments" begin
    con_curve_quadratic = ContinuousCurveQuadratic(3)
    # Weights should add to length of line segment of triangle
    @test sum(con_curve_quadratic.weights) ≈ 2.0
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(con_curve_quadratic)
    @test con_curve_quadratic.interpolation ≈ I 
end
#==========================================================================================
                            Discontinuous Segments
==========================================================================================#
@testset "Discontinuous Constant Segments" begin
    discon_curve_constant = DiscontinuousCurveConstant(3)
    # Weights should add to length of line segment of triangle
    @test sum(discon_curve_constant.weights) ≈ 2.0
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(discon_curve_constant)
    @test discon_curve_constant.interpolation ≈ I 
end
@testset "Discontinuous Linear Segments" begin
    discon_curve_linear = DiscontinuousCurveLinear(3)
    # Weights should add to length of line segment of triangle
    @test sum(discon_curve_linear.weights) ≈ 2.0
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(discon_curve_linear)
    @test discon_curve_linear.interpolation ≈ I 
end
@testset "Discontinuous Quadratic Segments" begin
    discon_curve_quadratic = DiscontinuousCurveQuadratic(3)
    # Weights should add to length of line segment of triangle
    @test sum(discon_curve_quadratic.weights) ≈ 2.0
    # Checking if set_nodal_interpolation! works as intended
    set_nodal_interpolation!(discon_curve_quadratic)
    @test discon_curve_quadratic.interpolation ≈ I 
end
