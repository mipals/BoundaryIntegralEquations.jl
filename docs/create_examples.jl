using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/examples")

examples = [
    # "2d_element_usage.jl",
    # "3d_element_usage.jl",
    # "3d_rosebem.jl",
    # "3d_chebyshev.jl",
    # "3d_lossy_sphere.jl",
    # "3d_visualization.jl",
    "3d_rigid_sphere.jl",
    "3d_oscilating_sphere.jl",
    "3d_pulsating_sphere.jl",
    "3d_cube_wave.jl",
    "3d_cube_wave_anechoic.jl",
    # "3d_fmm.jl",
    # "3d_hmatrix.jl",
    # "2d_infinite_cylinder.jl",
]

function uncomment_objects(str)
    str = replace(str, "###```@raw" => "```\n\n```@raw")
    str = replace(str, "###<object" => "<object")
    str = replace(str, "###```\n```" => "```")
    str
end

for example in examples
    example_filepath = normpath(joinpath(EXAMPLES_DIR, example))
    Literate.markdown(example_filepath, OUTPUT_DIR; execute=true, postprocess = uncomment_objects)
end
