#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, IntegralEquations, Plots, IterativeSolvers
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/meta_boundary"
tri_mesh_file = "examples/meshes/meta_new"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[2])
#==========================================================================================
                                3d Visualization
==========================================================================================#
using MeshViz
import WGLMakie as wgl
wgl.set_theme!(resolution=(1500, 1500))
tri_bc_ents = [358] .- 1
# Plotting entities in `bc_ents` red
tri_simple_mesh = create_bc_simple_mesh(tri_mesh,tri_bc_ents,false)
viz(tri_simple_mesh;showfacets=true)
tri_simple_bc   = create_bc_simple_mesh(tri_mesh,tri_bc_ents)
viz!(tri_simple_bc;showfacets=true,color=:red)
