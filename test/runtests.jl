using SafeTestsets, Test

@safetestset "2D: Curve Functions  " begin include("2d/curve_functions.jl")     end
@safetestset "2D: Assembly         " begin include("2d/assembly.jl")            end
@safetestset "3D: Surface Functions" begin include("3d/surface_functions.jl")   end
@safetestset "3D: Meshing          " begin include("3d/meshing.jl")             end
# @safetestset "3D-Assembly          " begin include("3d/assembly.jl")            end
