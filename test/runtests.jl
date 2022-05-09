using IntegralEquations
using LinearAlgebra
using IterativeSolvers
using Test

# For testing if interpolation nodes get set correct
import IntegralEquations: set_nodal_interpolation!
# For testing assembly against analytical scattering of sphere
import IntegralEquations: incoming_wave, plane_wave_scattering_sphere

#==========================================================================================
                                        2D
==========================================================================================#
@testset "Curve Functions" begin
    include("2d/curve_functions.jl")
end
# @testset "Meshing - 2D" begin
    # include("2d/meshing.jl")
# end
# @testset "Assembly - 2D" begin
#     # include("2d/assembly.jl") # A pretty computationally heavy test
# end
#==========================================================================================
                                        3D
==========================================================================================#
@testset "Functions" begin
    include("3d/surface_functions.jl")
end
@testset "Meshing - 3D" begin
    include("3d/meshing.jl")
end
# @testset "Assembly - 3D" begin
#     # include("3d/assembly.jl") # A pretty computationally heavy test
# end
