#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using IntegralEquations
using Plots
using IterativeSolvers
#==========================================================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finest"
# tri_mesh_file = "examples/meshes/sphere_1m_35k"
# tri_mesh_file = "examples/meshes/sphere_1m_77k"
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[2])

# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/test_binary.ply");
#==========================================================================================
                                    3d Visualization
==========================================================================================#
# using MeshViz
# import WGLMakie as wgl
# simple_mesh = create_simple_mesh(mesh)
# wgl.set_theme!(resolution=(800, 800))
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
ω = 2π*freq
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
xyzb = mesh.sources
n = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
LGM_on  = IntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true, depth=1,n=3,exterior=true)
LGM_off = IntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,exterior=true)
#===========================================================================================
                                        Tests
===========================================================================================#
using Test
x = ones(ComplexF64,size(LGM_on,1))
LGM_on.Ha*x - LGM_off.Ha*x
LGM_on.Ga*x - LGM_off.Ga*x
(LGM_on.Ha*x)./ (LGM_off.Ha*x)

y = rand(ComplexF64,size(LGM_on,1))
sources = mesh.sources
n = 3
m = 3
Ha,Ga,C = assemble_parallel!(mesh,kₐ,sources;m=m,n=n,progress=false)
sHa,sGa = assemble_parallel!(mesh,kₐ,sources;m=m,n=n,progress=false,sparse=true,depth=2)
sHa[1:4,1:4] ./ Ha[1:4,1:4]
