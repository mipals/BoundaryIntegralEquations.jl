module IntegralEquations

#==========================================================================================
                        Importing Base functions to overload
==========================================================================================#
using Base: Number
import Base: adjoint, *
#==========================================================================================
                            Using relevant packagess
==========================================================================================#
using LinearAlgebra
using SpecialFunctions
using GaussQuadrature
using ProgressMeter
using DelimitedFiles
using SparseArrays
using ForwardDiff
using LegendrePolynomials
using Meshes
using Base.Threads
#==========================================================================================
                            General Properties of Shape Functions
==========================================================================================#
include("ShapeFunction.jl")
#==========================================================================================
                                    Utility functions
==========================================================================================#
include("utils/mesh.jl")
include("utils/ambient_to_properties.jl")
include("utils/lossy_constants.jl")
#==========================================================================================
                                    2D Routines
==========================================================================================#
include("2d/CurveFunction.jl")
include("2d/shape_functions.jl")
include("2d/quadrature.jl")
include("2d/kernels.jl")
include("2d/meshing.jl")
#==========================================================================================
                                    3D Routines
==========================================================================================#
include("3d/SurfaceFunction.jl")
include("3d/kernels.jl")
include("3d/quadrature.jl")
include("3d/shape_functions.jl")
include("3d/shape_function_derivatives.jl")
include("3d/jacobian.jl")
include("3d/mesh.jl")
include("3d/read_comsol.jl")
include("3d/adaptive_integration.jl")
include("3d/triangular_modifications.jl")
include("3d/assembly_collocation.jl")
include("3d/assembly_galerkin.jl")
include("3d/visualizations.jl")
#==========================================================================================
                                Analytical solutions
==========================================================================================#
include("analytical/scattering.jl")
#==========================================================================================
                                Exporting relevant function
==========================================================================================#
# Helper functions
export adjoint, *

# Utility functions
export visco_thermal_constants, ambient_to_properties
# Mesh-related functions
export Mesh, Mesh3d, Mesh2d
export load3dTriangularComsolMesh,load3dQuadComsolMesh,read_comsol_mesh

# 1D element types (for 2D)

# 2D element types (for 3D)
export TriangularLinear, TriangularQuadratic, Triangular,
       DiscontinuousTriangularConstant, DiscontinuousTriangularLinear,
       DiscontinuousTriangularQuadratic, DiscontinuousTriangular,
       Quadrilateral,QuadrilateralLinear,QuadrilateralLinear4,
       QuadrilateralQuadratic,QuadrilateralQuadraticLagrange,
       DiscontinuousQuadrilateralConstant,DiscontinuousQuadrilateralLinear4,
       DiscontinuousQuadrilateralQuadraticLagrange
# Assembly
export assemble_parallel!
# visualizations
export create_simple_mesh, create_bc_simple_mesh

end
