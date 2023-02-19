using Documenter, BoundaryIntegralEquations, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/examples")

examples = [
    # "2d_element_usage.jl",
    # "3d_element_usage.jl",
    # "3d_visualization.jl",
    "3d_fmm.jl",
    "3d_hmatrix.jl",
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

makedocs(
    doctest=false,
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    modules = [BoundaryIntegralEquations],
    sitename="BoundaryIntegralEquations.jl",
    authors = "Mikkel Paltorp",
    pages = Any[
        "Home" => "index.md",
        "Background Theory" => "background_theory.md",
        # [
        #     "Computers and functions" => "computers_and_functions.md",
        #     "Elements" => "elements.md",
        #     "Quadrature" => "quadrature.md",
        #     "The Boundary Element Method" => "bem.md",
        #     "The Fast Multipole Method" => "fmm.md",
        # ],

        "2D" => [
            "Kernels"       => "2d_kernels.md",
            "Element Types" => "2d_elements.md",
            "Meshes"        => "2d_meshes.md",
        ],
        "3D" => [
            "Kernels"                       => "3d_kernels.md",
            "Element Types"                 => "3d_elements.md",
            "Meshes"                        => "3d_meshes.md",
            "Fast Operators"                => "3d_fast_operators.md",
            "Viscous and Thermal Losses"    => "3d_losses.md",
            "Experimental Features"         => "3d_experimental.md",
            ],
        "Examples" => [
            "examples/2d_element_usage.md",
            "examples/2d_infinite_cylinder.md",
            "examples/3d_element_usage.md",
            "examples/3d_visualization.md",
            "examples/3d_fmm.md",
            "examples/3d_hmatrix.md",
            ]
        ]
    )

    # servedocs(skip_dir=joinpath("docs","src","examples"))

deploydocs(
    repo = "github.com/mipals/BoundaryIntegralEquations.jl.git",
)
