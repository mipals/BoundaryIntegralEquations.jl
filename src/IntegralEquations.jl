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
using ForwardDiff
using LegendrePolynomials
using Meshes
using MeshViz
using Base.Threads
import GLMakie as Mke
#==========================================================================================
                                    Utility functions
==========================================================================================#
include("utils/basis.jl")
include("utils/mesh.jl")
include("utils/ambient_to_properties.jl")
include("utils/lossy_constants.jl")
#==========================================================================================
                                    2D Routines
==========================================================================================#
#==========================================================================================
                                    3D Routines
==========================================================================================#
include("3d/kernels.jl")
include("3d/quadrature.jl")
include("3d/shape_functions.jl")
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
       QuadrilateralQuadratic,QuadrilateralQuadratic9
# Assembly
export assemble_parallel!
# visualizations
export create_simple_mesh, create_bc_simple_mesh, viz

end
