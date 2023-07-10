using SafeTestsets

# Formatting
@safetestset "Aqua testing         " begin include("aqua_test.jl") end
# 2D Tests
@safetestset "2D: Curve Functions  " begin include("2d/curve_functions_test.jl")     end
@safetestset "2D: Assembly         " begin include("2d/assembly_test.jl")            end
# @safetestset "2D: FMM              " begin include("2d/fmm_test.jl")                 end
# 3D Tests
@safetestset "3D: Surface Functions" begin include("3d/surface_functions_test.jl")   end
@safetestset "3D: Meshing          " begin include("3d/meshing_test.jl")             end
@safetestset "3D: Assembly         " begin include("3d/assembly_test.jl")            end
@safetestset "3D: FMM              " begin include("3d/fmm_test.jl")                 end
@safetestset "3D: Losses           " begin include("3d/losses_test.jl")              end
