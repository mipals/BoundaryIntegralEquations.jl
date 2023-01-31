#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, BoundaryIntegralEquations, Plots, IterativeSolvers
#==========================================================================================
                                    Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/meta_new"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])
#==========================================================================================
                                3d Visualization
==========================================================================================#
using MeshViz
import WGLMakie as wgl
wgl.set_theme!(resolution=(2000, 1500))
tri_bc_ents = [1,358] .- 1
# Plotting boundary conditions
tri_simple_mesh = create_bc_simple_mesh(tri_mesh,tri_bc_ents,false)
viz(tri_simple_mesh;showfacets=true,alpha=0.5)
# Termination end plotted in red
tri_simple_bc   = create_bc_simple_mesh(tri_mesh,tri_bc_ents[1])
viz!(tri_simple_bc;showfacets=true,color=:red)
# Source end plotted in blue
tri_simple_bc   = create_bc_simple_mesh(tri_mesh,tri_bc_ents[2])
viz!(tri_simple_bc;showfacets=true,color=:blue)

# Tests
normals  = tri_mesh.normals
tangent1 = tri_mesh.tangents
tangent2 = tri_mesh.sangents
# Creating RHS
u₀ = 1e-2
vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]
n    = length(vn0)
v0   = [zeros(2n); u₀*ones(n)]
#==========================================================================================
                                Defining frequency grid
==========================================================================================#
freqs = collect(1000:100.0:5000) # Grid from paper
n_freqs = length(freqs)
results = zeros(n_freqs)
# prog = Progress(n_freqs, 0.2, "Frequency sweep: \t", 50)
freq = freqs[1]
# for freq in freqs
    rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
    k   = 2*π*freq/c                              # Wavenumber                [1/m]
    #======================================================================================
                            Iterative Solution of the 1-variable system
    ======================================================================================#
    LGM = BoundaryIntegralEquations.LossyGlobalOuter(tri_mesh,freq;fmm_on=false,depth=1,n=3)
    rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0));verbose=true)
    @time pa = gmres(LGM,rhs;verbose=true);
