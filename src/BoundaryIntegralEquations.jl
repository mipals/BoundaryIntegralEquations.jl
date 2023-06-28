module BoundaryIntegralEquations

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
using FastGaussQuadrature
using ProgressMeter
using DelimitedFiles
using SparseArrays
using LegendrePolynomials
using Meshes
using LoopVectorization
using LinearMaps
using IterativeSolvers
using StatsBase
using FMM3D
using FileIO
using MeshIO
using Base.Threads
using StaticArrays
# The gives out an @Info from the HMatrices when loaded
# using Gmsh: Gmsh, gmsh
using ForwardDiff
# using HMatrices: HMatrices
import HMatrices: assemble_hmat, PartialACA, ClusterTree
#==========================================================================================
                                        Constants
==========================================================================================#
const Point3D = SVector{3,Float64}
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
include("2d/curve_functions.jl")
include("2d/quadrature.jl")
include("2d/kernels.jl")
include("2d/meshing.jl")
include("2d/assembly_collocation.jl")
#==========================================================================================
                                    3D Routines
==========================================================================================#
include("3d/SurfaceFunction.jl")
include("3d/kernels.jl")
include("3d/quadrature.jl")
include("3d/surface_functions.jl")
include("3d/interpolation_function_derivatives.jl")
include("3d/jacobian.jl")
include("3d/mesh.jl")
include("3d/meshing.jl")
include("3d/read_comsol.jl")
include("3d/adaptive_integration.jl")
include("3d/assembly_collocation.jl")
include("3d/assembly_collocation_losses.jl")
include("3d/LossyGlobalMatrix.jl")
include("3d/visualizations.jl")
include("3d/fast_multipole_method.jl")
include("3d/hmatrices.jl")
include("3d/rosebem.jl")
# Old/Experimental features
include("3d/experimental/LossyBlockMatrix.jl")
include("3d/experimental/LossyBlockMatrixCompact.jl")
include("3d/experimental/LossyInexactKrylov.jl")
include("3d/experimental/assembly_galerkin.jl")
#==========================================================================================
                                Analytical solutions
==========================================================================================#
include("analytical/2d_scattering.jl")
include("analytical/3d_scattering.jl")
include("analytical/sphere_first_order.jl")
#==========================================================================================
                                Exporting relevant function
==========================================================================================#
# Helper functions
export adjoint,*

# Abstract Shape Function Types
export CurveFunction,ShapeFunction,SurfaceFunction

# Utility functions
export visco_thermal_constants,ambient_to_properties

# Mesh-related functions
export Mesh,Mesh3d,Mesh2d
export load3dTriangularComsolMesh,load3dQuadComsolMesh,read_comsol_mesh,load3dTriangularMesh

# 1D element types (for 2D)
export ContinuousCurveLinear,ContinuousCurveQuadratic,
       DiscontinuousCurveConstant,DiscontinuousCurveLinear,DiscontinuousCurveQuadratic

# 2D element types (for 3D)
export TriangularLinear,TriangularQuadratic,Triangular,
       DiscontinuousTriangularConstant,DiscontinuousTriangularLinear,
       DiscontinuousTriangularQuadratic,DiscontinuousTriangular
export Quadrilateral,QuadrilateralLinear,QuadrilateralLinear4,
       QuadrilateralQuadratic,QuadrilateralQuadraticLagrange,
       DiscontinuousQuadrilateralConstant,DiscontinuousQuadrilateralLinear4,
       DiscontinuousQuadrilateralQuadraticLagrange

# Assembly
export assemble_parallel!

# Fast mulitpole Method
export FMMGOperator, FMMFOperator, evaluate_targets

# H-Matrices
export HGOperator, HFOperator

# Visualizations
export create_simple_mesh, create_bc_simple_mesh, viz, create_vizualization_data

# Losses
export LossyBlockMatrix,LossyBlockMatrixCompact,LossyOneVariableInner,LossyOneVariableOuter
export compute_lossy_rhs,LossyGlobalOuter

end
