using IntegralEquations
using LinearAlgebra
using Test
import IntegralEquations: set_nodal_interpolation!

@testset "Testing 2d" begin
    # Write your tests here.
end


@testset "Testing 3d" begin
    include("3d/shape_functions.jl")
    include("3d/meshing.jl")
end

@testset "Testing utility function" begin
    # Write your tests here.
end

