using Documenter
using DocumenterCitations
using BoundaryIntegralEquations
using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/examples")

examples = [
    # "2d_element_usage.jl",
    # "3d_element_usage.jl",
    # "3d_visualization.jl",
    # "3d_rigid_sphere.jl",
    # "3d_oscilating_sphere.jl",
    # "3d_pulsating_sphere.jl",
    # "3d_cube_wave.jl",
    # "3d_cube_wave_anechoic.jl",
    # "3d_lossy_sphere.jl",
    "3d_rosebem.jl",
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

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric  # default
)

makedocs(bib,
    doctest=false,
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    modules = [BoundaryIntegralEquations],
    sitename="BoundaryIntegralEquations.jl",
    authors = "Mikkel Paltorp",
    pages = Any[
        "Home" => "index.md",
        "Background Theory"          => "theory_background.md",
        "ROSEBEM Theory"             => "theory_rosebem.md",
        "Viscothermal Losses Theory" => "theory_lossy.md",
        "Examples" => [
            "3D" => [
                "examples/3d_rosebem.md",
                "examples/3d_rigid_sphere.md",
                "examples/3d_oscilating_sphere.md",
                "examples/3d_pulsating_sphere.md",
                "examples/3d_cube_wave.md",
                "examples/3d_cube_wave_anechoic.md",
                "examples/3d_lossy_sphere.md",
                "examples/3d_element_usage.md",
            ],
            "2D" => [
                "examples/2d_infinite_cylinder.md",
                "examples/2d_element_usage.md",
            ],
        ],
        "Internals" => [
            "3D" => [
                "Kernels"                       => "3d_kernels.md",
                "Element Types"                 => "3d_elements.md",
                "Meshes"                        => "3d_meshes.md",
                "Fast Operators"                => "3d_fast_operators.md",
                "Viscous and Thermal Losses"    => "3d_losses.md",
                "ROSEBEM"                       => "3d_rosebem.md",
                # "Experimental Features"         => "3d_experimental.md",
                # "Analytical Expressions"        => "3d_analytical.md"
            ],
            "2D" => [
                "Kernels"       => "2d_kernels.md",
                "Element Types" => "2d_elements.md",
                "Meshes"        => "2d_meshes.md",
            ],
        ]
    ]
    )

# servedocs(skip_dir=joinpath("docs","src","examples"))

deploydocs(
    repo = "github.com/mipals/BoundaryIntegralEquations.jl.git",
)
