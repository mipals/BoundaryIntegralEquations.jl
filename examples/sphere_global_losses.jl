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
tri_mesh_file = "examples/meshes/sphere_1m_fine"
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
using MeshViz
import WGLMakie as wgl
simple_mesh = create_simple_mesh(mesh)
wgl.set_theme!(resolution=(1200, 1200))
viz(simple_mesh, showfacets = true)
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
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
LGM = IntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3)
# Creating RHS
n  = length(normals[3,:])
id = 4
if id == 1
    v0 = [u₀*ones(n); zeros(2n)]
    coordinates = [radius*ones(n,1) acos.(xyzb[id,:]/radius)]
elseif id == 2
    v0 = [zeros(n); u₀*ones(n); zeros(n)]
    coordinates = [radius*ones(n,1) acos.(xyzb[id,:]/radius)]
elseif id == 3
    v0 = [zeros(2n); u₀*ones(n)]
    coordinates = [radius*ones(n,1) acos.(xyzb[id,:]/radius)]
else
    # Rotate coordinates so that the
    using Rotations
    tmp_coordinates = RotX(π/4)*mesh.coordinates
    id = 3
    v0 = [zeros(n); u₀/sqrt(2)*ones(n); u₀/sqrt(2)*ones(n)]
    coordinates = [radius*ones(n,1) acos.(tmp_coordinates[3,:]/radius)]
end
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0));verbose=true);
@time pa = gmres(LGM,rhs;verbose=true);
# 323695 Dofs: 2438.138033 seconds (136.03 M allocations: 91.786 GiB, 0.29% gc time, 0.00% compilation time)
# 2*20 + 14 = 54 iterations

# Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                IntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = coordinates[:,2]*180.0/pi
perm = sortperm(ang_axis)

# Plotting
K = 1
gr(size=(600,300))
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM-global",marker=:cross,
        markersize=2,color=:black,dpi=600,xtickfontsize=15,ytickfontsize=15,
        legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
xlabel!("Angle")
# savefig("pa$(n).png")
# savefig("/Users/mpasc/Dropbox/Apps/ShareLaTeX/JTCA_Iterative_losses/figures/notes/pa$(n).png")
#===========================================================================================
                                Checking results
===========================================================================================#
ph   = -LGM.tau_a/LGM.tau_h * pa
dph  = -gmres(LGM.Gh,LGM.Hh*ph)
# 517 iterations for 69k DOF at freq=1000
dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true) # <- This is the bottleneck...
# The bottleneck could be the reason for the extended system not working?
v   = v0 - (LGM.phi_a*LGM.Dc*pa +
            LGM.phi_a*LGM.Nd*dpa +
            LGM.phi_h*LGM.Dc*ph +
            LGM.phi_h*LGM.Nd*dph)

dvn = -gmres(LGM.Gv,LGM.Hv*v;verbose=true)

# The scaling
using Statistics
# Why are we off by this constant???
scaling  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v))
scalings = (-LGM.Nd'*dvn) ./ (LGM.Dr*v)
plot(ang_axis[perm],abs.(scalings[perm]),yaxis=:log)
scatter!(ang_axis[perm],abs.(scalings[perm]),yaxis=:log,markersize=1)
#===========================================================================================
                                    Tests!
                        Null-divergence is not satisfied?
===========================================================================================#
using Test
# Checking constants
@test LGM.mu_a == LGM.phi_a - LGM.phi_h*LGM.tau_a/LGM.tau_h
@test LGM.mu_h == LGM.phi_h*LGM.tau_a/LGM.tau_h
# All of these are solved quantities - So they're satisfied by construction
@test LGM.Hh*ph ≈ -LGM.Gh*dph       # Checking Acoustical BEM System
@test LGM.Ha*pa ≈ -LGM.Ga*dpa       # Checking Thermal BEM system
@test LGM.Hv*v  ≈ -LGM.Gv*dvn       # Checking viscous BEM system
@test LGM.tau_a*pa ≈ -LGM.tau_h*ph  # Checking Isothermal condition
# These are the extra conditions
@test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0
@test LGM.phi_a*LGM.Dc*pa + LGM.phi_a*LGM.Nd*dpa + LGM.phi_h*LGM.Dc*ph + LGM.phi_h*LGM.Nd*dph + v ≈ v0
# Checking null-divergence (off by a factor??)
@test LGM.Dr*v ≈ -LGM.Nd'*dvn
# We can fix this by solving a rectangular system instead
rhs_v = -[LGM.Hv*v; LGM.Dr*v]
V   = [LGM.Gv; LGM.Nd']
dvn2 = lsqr(V,rhs_v;atol=1e-15,btol=1e-13,verbose=true);
@test LGM.Dr*v ≈ -LGM.Nd'*dvn2            # Now this is satisfied
@test LGM.Hv*v ≈ -LGM.Gv*dvn2 atol = 1e-4 # But the BEM system has a higher error

rhs_w = -[LGM.Gv*dvn; LGM.Nd'*dvn]
W  = [LGM.Hv; LGM.Dr]
v2  = lsqr(W,rhs_w;verbose=true);
@test LGM.Dr*v2 ≈ -LGM.Nd'*dvn            # Now this is satisfied
@test LGM.Hv*v2 ≈ -LGM.Gv*dvn atol = 1e-4 # But the BEM system has a higher error

# The rectangular system solves the scaling, but since "v" is the same we still have an error later on
scaling_v  =  median(abs.(-LGM.Nd'*dvn2) ./ abs.(LGM.Dr*v))
scalings_v = (-LGM.Nd'*dvn2) ./ (LGM.Dr*v)

scaling_w  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v2))
scalings_w = (-LGM.Nd'*dvn) ./ (LGM.Dr*v2)

## Can only be done on the
B1 = [zeros(ComplexF64,n,n)                  LGM.Dr - LGM.Nd'*(Matrix(LGM.Gv)\Matrix(LGM.Hv))]
B1 = [zeros(ComplexF64,n,n) LGM.Dr-LGM.Nd'*GvHv]
B2 = [LGM.mu_a*LGM.Dc + LGM.Nd*(LGM.mu_h*Matrix(LGM.Gh)\Matrix(LGM.Hh) - LGM.phi_a*LGM.Ga\LGM.Ha) I]
B = [zeros(ComplexF64,n,n)                  LGM.Dr - LGM.Nd'*(Matrix(LGM.Gv)\Matrix(LGM.Hv));
 LGM.mu_a*LGM.Dc + LGM.Nd*(LGM.mu_h*(Matrix(LGM.Gh)\Matrix(LGM.Hh)) - LGM.phi_a*(LGM.Ga\LGM.Ha)) I]
rhs_b = [zeros(n); v0]
solB = gmres(B,rhs_b;verbose=true, maxiter=100)
p_b = solB[0n+1:1n]
v_b = solB[1n+1:4n]
dvn_b =  -gmres(LGM.Gv,LGM.Hv*v_b;verbose=true)
scaling_b  =  median(abs.(-LGM.Nd'*dvn_b) ./ abs.(LGM.Dr*v_b))
scalings_b = (-LGM.Nd'*dvn_b) ./ (LGM.Dr*v_b)

plot(ang_axis[perm],abs.(pa./p_b)[perm],yaxis=:log)
scatter(ang_axis,abs.(pa./p_b),yaxis=:log,markersize=1)
# v = v_b
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
vx = v[0n+1:1n]
vy = v[1n+1:2n]
vz = v[2n+1:3n]
dvx = dvn[0n+1:1n]
dvy = dvn[1n+1:2n]
dvz = dvn[2n+1:3n]

# Local components of the viscous velocity on the boundary
v_n0   = LGM.Nd'*v
# v_n0   =  normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz
v_t1   = tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   = tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

# Plotting
K = 1
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]/mean(scalings)),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
# ylims!((-5e-6,5e-6))
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("Re(Vt)")
plot(plt1,plt2,plt3,layout=(3,1),dpi=500)
# savefig("/Users/mpasc/Dropbox/Apps/ShareLaTeX/JTCA_Iterative_losses/figures/notes/test.png")

#===========================================================================================
                                Scaling of model
===========================================================================================#
nDOF = n
nDOF = 154_000
# A single BEM matrix
2(nDOF)^2*2*8/(2^30)
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
Ha = FMMHOperator(mesh,kₐ;n=3,eps=1e-6,offset=0.2,nearfield=true,depth=1)
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
