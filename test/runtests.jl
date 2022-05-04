using IntegralEquations
using LinearAlgebra
using IterativeSolvers
using Test

# For testing if interpolation nodes get set correct
import IntegralEquations: set_nodal_interpolation!
# For testing assembly against analytical scattering of sphere
import IntegralEquations: incoming_wave, plane_wave_scattering_sphere

# @testset "Testing 2d" begin
#     # Write your tests here.
# end

@testset "Testing 3d" begin
    include("3d/shape_functions.jl")
    include("3d/meshing.jl")
    # include("3d/assembly.jl") # A pretty computationally heavy test
end

# @testset "Testing utility function" begin
#     # Write your tests here.
# end
