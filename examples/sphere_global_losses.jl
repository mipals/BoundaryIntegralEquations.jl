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
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finest"
# tri_mesh_file = "examples/meshes/sphere_1m_35k"
tri_mesh_file = "examples/meshes/sphere_1m_77k"
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])
#==========================================================================================
                                    3d Visualization
==========================================================================================#
using MeshViz
import WGLMakie as wgl
simple_mesh = create_simple_mesh(mesh)
wgl.set_theme!(resolution=(1200, 1200))
viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 1000.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
ω = 2π*freq
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)

# Creating RHS
vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]

n  = length(vn0)
v0 = [zeros(2n); u₀*ones(n)]
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
LGM = IntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true,depth=1,n=3)
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0));verbose=true)
@time pa = gmres(LGM,rhs;verbose=true);

# Generating analytical solution
coordinates = [radius*ones(n,1) acos.(xyzb[3,:]/radius)]
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                IntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)

# Plotting
K = 1
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM-global",marker=:cross,markersize=2,color=:black)
# scatter!(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:red)
ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
# savefig("/Users/mpasc/Dropbox/Apps/ShareLaTeX/JTCA_Iterative_losses/figures/notes/pa$(n).png")
#===========================================================================================
                                Checking results
===========================================================================================#
ph   = -LGM.tau_a/LGM.tau_h * pa
dph  = -gmres(LGM.Gh,LGM.Hh*ph)
tmp1 =  gmres(LGM.Gh,LGM.Hh*pa;verbose=true)
# 517 iterations for 69k DOF at freq=1000
dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true) # <- This is the bottleneck...
# The bottleneck could be the reason for the extended system not working?

v   = v0 - (LGM.mu_a*LGM.Dc*pa + LGM.mu_h*LGM.Nd*tmp1 + LGM.phi_a*LGM.Nd*dpa)
dvn = -gmres(LGM.Gv,LGM.Hv*v) # Here you solve that it works out.

vx = v[0n+1:1n]
vy = v[1n+1:2n]
vz = v[2n+1:3n]
dvx = dvn[0n+1:1n]
dvy = dvn[1n+1:2n]
dvz = dvn[2n+1:3n]

using Statistics
# Why are we off by this constant???
scaling  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v))
scalings = (-LGM.Nd'*dvn) ./ (LGM.Dr*v)

# Local components of the viscous velocity on the boundary
v_n0   = LGM.Nd'*v
v_t1   = tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   = tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

#===========================================================================================
                                    Tests!
                        Null-divergence is not satisfied?
===========================================================================================#
using Test
@test LGM.Hh*ph ≈ -LGM.Gh*dph       # Checking Acoustical BEM System
@test LGM.Ha*pa ≈ -LGM.Ga*dpa       # Checking Thermal BEM system
@test LGM.Hv*v  ≈ -LGM.Gv*dvn       # Checking viscous BEM system
@test LGM.tau_a*pa ≈ -LGM.tau_h*ph  # Checking Isothermal condition (kind of redundant)
@test LGM.Dr*v ≈ -LGM.Nd'*dvn       # Checking null-divergence (off by a factor??)
# Checking no-slip condition
@test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0

#===========================================================================================
                                Scaling of model
===========================================================================================#
nDOF = n
nDOF = 154_000
# A single BEM matrix
3(nDOF)^2*2*8/(2^30)
# Full dense system (100 times a single BEM Matrix)
(10nDOF)^2*2*8/(2^30)
# 6 BEM Systems, 2 inverse products (Gv^{-1}Hv, Gh^{-1}Hv) and 1 collection
(9)*(nDOF)^2*2*8/(2^30) # Computing inverse would also kill performance...

Base.summarysize(LGM)/(2^30)
#===========================================================================================
                                Memory tracking
===========================================================================================#
Ga = FMMGOperator(mesh,kₐ;n=1,eps=1e-6,offset=0.2,nearfield=true,depth=1)
@time Ga*pa;
Ha = FMMFOperator(mesh,kₐ;n=3,eps=1e-6,offset=0.2,nearfield=true,depth=1)
@time Ha*pa;

mem_matrix = Base.summarysize(LGM)/(2^20)
mem_multiply = @allocated begin
    LGM*vn0
end
@time LGM*vn0;
@time LGM.Ga*vn0;
@time LGM.Ha*vn0;

@time LGM.inner*vn0;

Base.summarysize(LGM.Gv)/(2^20)
Base.summarysize(LGM.Hv)/(2^20)
Base.summarysize(LGM.Gh)/(2^20)
Base.summarysize(LGM.Hh)/(2^20)
Base.summarysize(LGM.Ga)/(2^20)
Base.summarysize(LGM.Ha)/(2^20)
Base.summarysize(LGM.Nd)/(2^20)
# Note that Dc requires less memory as we use a CSC format (and Dc has less columns)
Base.summarysize(LGM.Dr)/(2^20)
Base.summarysize(LGM.Dc)/(2^20)
Base.summarysize(LGM.inner)/(2^20)

(Base.summarysize(LGM.Gv) + Base.summarysize(LGM.Hv) +
 Base.summarysize(LGM.Gh) + Base.summarysize(LGM.Hh) +
 Base.summarysize(LGM.Ga) + Base.summarysize(LGM.Ha) +
 Base.summarysize(LGM.Dr) + Base.summarysize(LGM.Dc) +
 Base.summarysize(LGM.Nd))/(2^20)
Base.summarysize(LGM)/(2^20) # Why is this less?

LGM.Dr[1:n,0n+1:1n] - LGM.Dc[0n+1:1n,1:n]
LGM.Dr[1:n,1n+1:2n] - LGM.Dc[1n+1:2n,1:n]
LGM.Dr[1:n,2n+1:3n] - LGM.Dc[2n+1:3n,1:n]
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Plotting
K = 100
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|p|"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
# plt2 = scatter(ang_axis[1:K:end],abs.(v_n0[1:K:end])/scaling,label="BEM",marker=:cross,markersize=2,color=:black)
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]/mean(scalings)),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|Vn|"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
# ylims!((-5e-6,5e-6))
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|Vt|")
plot(plt1,plt2,plt3,layout=(3,1),dpi=500)
# savefig("/Users/mpasc/Dropbox/Apps/ShareLaTeX/JTCA_Iterative_losses/figures/notes/test.png")

#===========================================================================================
                            Full Solution
===========================================================================================#
B*LGM.Nd

## FULL
B = Matrix(LGM.Dr) - Matrix(LGM.Nd')*(Matrix(LGM.Gv)\Matrix(LGM.Hv))
C = LGM.mu_a*LGM.Dc + LGM.Nd*(LGM.mu_h*Matrix(LGM.Gh)\Matrix(LGM.Hh) - LGM.phi_a*LGM.Ga\LGM.Ha)

N = zeros(ComplexF64,3n,3n)
N[0n+1:1n,0n+1:1n] = LGM.Ga
N[1n+1:2n,1n+1:2n] = LGM.Ga
N[2n+1:3n,2n+1:3n] = LGM.Ga
M = [zeros(ComplexF64,n,n) B;
     C         I]

rhs2 = [zeros(ComplexF64,n); v0]

sol = gmres(M,rhs2;verbose=true,maxiter=100)
pa  = gmres(LGM,rhs;verbose=true);
# pa = sol[0n+1:1n]
vx = sol[1n+1:2n]
vy = sol[2n+1:3n]
vz = sol[3n+1:4n]
v = [vx;vy;vz]
dvn  = -gmres(LGM.Gv,LGM.Hv*[vx;vy;vz])
ph   = -LGM.tau_a/LGM.tau_h * pa
dph  =  LGM.tau_a/LGM.tau_h * gmres(LGM.Gh,LGM.Hh*pa)
tmp1 =  gmres(LGM.Gh,LGM.Hh*pa;verbose=true)
dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true)

@test LGM.Ha*pa ≈ -LGM.Ga*dpa          # Checking Acoustical BEM System
@test LGM.Hh*ph ≈ -LGM.Gh*dph          # Checking Thermal BEM system
@test LGM.Hv*v  ≈ -LGM.Gv*dvn          # Checking viscous BEM system
@test LGM.tau_a*pa ≈ -LGM.tau_h*ph
@test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0 atol=1e-1
@test LGM.Dr*v ≈ -LGM.Nd'*dvn      # Checking null-divergence

v_n0   = LGM.Nd'*[vx;vy;vz]
v_t1   = tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   = tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Plotting
K = 1
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|p|"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|Vn|"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
# ylims!((-5e-6,5e-6))
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|Vt|")
plot(plt1,plt2,plt3,layout=(3,1))
# plt = plot(plt1,plt2,plt3,layout=(3,1),dpi=300)
