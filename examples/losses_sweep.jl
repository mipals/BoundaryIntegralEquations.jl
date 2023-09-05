#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, BoundaryIntegralEquations, Plots, IterativeSolvers
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/meta_boundary"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                physics_order=tri_physics_orders[2])
#==========================================================================================
                                3d Visualization
==========================================================================================#
using Meshes
import WGLMakie as wgl
wgl.set_theme!(resolution=(900, 900))
tri_bc_ents = [358] .- 1
# Plotting entities in `bc_ents` red
tri_simple_mesh = create_bc_simple_mesh(tri_mesh_file,tri_mesh,tri_bc_ents,false)
viz(tri_simple_mesh;showfacets=true)
tri_simple_bc   = create_bc_simple_mesh(tri_mesh_file,tri_mesh,tri_bc_ents)
viz!(tri_simple_bc;showfacets=true,color=:red)
#==========================================================================================
                                Setting up constants
==========================================================================================#
radius = 1.0                                     # Radius of sphere_1m       [m]
xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents
# Creating RHS
vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]
n    = length(vn0)
v0 = [zeros(2n); u₀*ones(n)]
#==========================================================================================
                                Defining frequency grid
==========================================================================================#
freqs = collect(1000:100.0:5000) # Grid from paper
n_freqs = length(freqs)
results = zeros(n_freqs)
prog = Progress(n_freqs, 0.2, "Frequency sweep: \t", 50)
freq = freqs[1]
# for freq in freqs
    rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
    k   = 2*π*freq/c                              # Wavenumber                [1/m]
    #======================================================================================
                            Iterative Solution of the 1-variable system
    ======================================================================================#
    LGM = BoundaryIntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true,depth=1,n=3)
    rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0));verbose=true)
    @time pa = gmres(LGM,rhs;verbose=true);
    #======================================================================================
                                    Checking results
    ======================================================================================#
    ph   = -LGM.tau_a/LGM.tau_h * pa
    dph  = -gmres(LGM.Gh,LGM.Hh*ph)
    tmp1 =  gmres(LGM.Gh,LGM.Hh*pa;verbose=true)
    dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true) # <- This is the bottleneck...

    v   = v0 - (LGM.mu_a*LGM.Dc*pa + LGM.mu_h*LGM.Nd*tmp1 + LGM.phi_a*LGM.Nd*dpa)
    dvn = -gmres(LGM.Gv,LGM.Hv*v) # Here you solve that it works out.

    vx = v[0n+1:1n]
    vy = v[1n+1:2n]
    vz = v[2n+1:3n]
    dvx = dvn[0n+1:1n]
    dvy = dvn[1n+1:2n]
    dvz = dvn[2n+1:3n]

    # using Statistics
    # Why are we off by this constant???
    scaling  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v))
    scalings = (-LGM.Nd'*dvn) ./ (LGM.Dr*v)

    # Local components of the viscous velocity on the boundary
    v_n0   = LGM.Nd'*v
    v_t1   = tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
    v_t2   = tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
    vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

    # next!(prog)
# end
