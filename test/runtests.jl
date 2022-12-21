using SafeTestsets, Test

# 2D Tests
@safetestset "2D: Curve Functions  " begin include("2d/curve_functions.jl")     end
@safetestset "2D: Assembly         " begin include("2d/assembly.jl")            end
# @safetestset "2D: FMM              " begin include("2d/fmm.jl")                 end
# 3D Tests
@safetestset "3D: Surface Functions" begin include("3d/surface_functions.jl")   end
@safetestset "3D: Meshing          " begin include("3d/meshing.jl")             end
@safetestset "3D: Assembly         " begin include("3d/assembly.jl")            end
@safetestset "3D: FMM              " begin include("3d/fmm.jl")                 end
@safetestset "3D: Losses           " begin include("3d/losses.jl")              end
