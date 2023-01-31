#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
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
tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
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
wgl.set_theme!(resolution=(800, 800))
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
n = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
LGM = BoundaryIntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3)
#* Creating the right-hand side
id = 1
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
    #! Rotate coordinates so that the
    using Rotations
    tmp_coordinates = RotX(π/4)*mesh.coordinates
    id = 3
    v0 = [zeros(n); u₀/sqrt(2)*ones(n); u₀/sqrt(2)*ones(n)]
    coordinates = [radius*ones(n,1) acos.(tmp_coordinates[3,:]/radius)]
end
rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0));verbose=true);
@time pa = gmres(LGM,rhs;verbose=true,reltol=1e-17);
#? Try to solve using a direct solver / dense matrix

#* Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
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
# pa = p_b
ph   = -LGM.tau_a/LGM.tau_h * pa
dph  = -gmres(LGM.Gh,LGM.Hh*ph)
#* 517 iterations for 69k DOF at freq=1000
dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true) # <- This is the bottleneck...
#! The bottleneck could be the reason for the extended system not working?
v   = v0 - (LGM.phi_a*LGM.Dc*pa +
            LGM.phi_a*LGM.Nd*dpa +
            LGM.phi_h*LGM.Dc*ph +
            LGM.phi_h*LGM.Nd*dph)

dvn = -gmres(LGM.Gv,LGM.Hv*v;verbose=true)

#* The scaling
using Statistics
#? Why are we off by this constant?
scaling  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v))
#* LGM.Dr*v + LGM.Nd'*dvn = 0
scalings = (-LGM.Nd'*dvn) ./ (LGM.Dr*v)
plot(ang_axis[perm],abs.(scalings[perm]),yaxis=:log)
scatter!(ang_axis[perm],abs.(scalings[perm]),yaxis=:log,markersize=1)
#===========================================================================================
                                    Tests!
                        Null-divergence is not satisfied?
===========================================================================================#
using Test
#* Checking constants
@test LGM.mu_a == LGM.phi_a - LGM.phi_h*LGM.tau_a/LGM.tau_h
@test LGM.mu_h == LGM.phi_h*LGM.tau_a/LGM.tau_h
#* All of these are solved quantities - So they're satisfied by construction
@test LGM.Hh*ph ≈ -LGM.Gh*dph       #! Checking Thermal BEM System
@test LGM.Ha*pa ≈ -LGM.Ga*dpa       #! Checking Acoustics BEM system
@test LGM.Hv*v  ≈ -LGM.Gv*dvn       #! Checking viscous BEM system
@test LGM.tau_a*pa ≈ -LGM.tau_h*ph  #! Checking Isothermal condition
#* These are the extra conditions
@test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0
#! Checking null-divergence (off by a factor??)
@test LGM.Dr*v ≈ -LGM.Nd'*dvn
# #* We can fix this by solving a rectangular system instead
# rhs_v = -[LGM.Hv*v; LGM.Dr*v]
# V   = [LGM.Gv; LGM.Nd']
# dvn2 = lsqr(V,rhs_v;atol=1e-15,btol=1e-13,verbose=true);
# @test LGM.Dr*v ≈ -LGM.Nd'*dvn2            # Now this is satisfied
# @test LGM.Hv*v ≈ -LGM.Gv*dvn2 atol = 1e-4 # But the BEM system has a higher error

# rhs_w = -[LGM.Gv*dvn; LGM.Nd'*dvn]
# W  = [LGM.Hv; LGM.Dr]
# v2  = lsqr(W,rhs_w;verbose=true);
# @test LGM.Dr*v2 ≈ -LGM.Nd'*dvn            # Now this is satisfied
# @test LGM.Hv*v2 ≈ -LGM.Gv*dvn atol = 1e-4 # But the BEM system has a higher error

# #* The rectangular system solves the scaling, but since "v" is the same we still have problems
# scaling_v  =  median(abs.(-LGM.Nd'*dvn2) ./ abs.(LGM.Dr*v))
# scalings_v = (-LGM.Nd'*dvn2) ./ (LGM.Dr*v)

# scaling_w  =  median(abs.(-LGM.Nd'*dvn) ./ abs.(LGM.Dr*v2))
# scalings_w = (-LGM.Nd'*dvn) ./ (LGM.Dr*v2)

#! Pretty slow as all things are dense
B = [zeros(ComplexF64,n,n)                  LGM.Dr - LGM.Nd'*(Matrix(LGM.Gv)\Matrix(LGM.Hv));
 LGM.mu_a*LGM.Dc + LGM.Nd*(LGM.mu_h*(Matrix(LGM.Gh)\Matrix(LGM.Hh)) - LGM.phi_a*(LGM.Ga\LGM.Ha)) I]
rhs_b = [zeros(n); v0]
solB = gmres(B,rhs_b;verbose=true, maxiter=500) # Requires a lot iterations. Not promising
p_b = solB[0n+1:1n]
v_b = solB[1n+1:4n]
dvn_b =  -gmres(LGM.Gv,LGM.Hv*v_b;verbose=true)

scaling_b  =  median(abs.(-LGM.Nd'*dvn_b) ./ abs.(LGM.Dr*v_b))
scalings_b = (-LGM.Nd'*dvn_b) ./ (LGM.Dr*v_b)

@test LGM.Dr*v_b ≈ -LGM.Nd'*dvn_b atol=1e-3

plot(ang_axis[perm], abs.(pa./p_b)[perm])
scatter!(ang_axis[perm], abs.(pa./p_b)[perm],markersize=1)
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
vx = v[0n+1:1n]
vy = v[1n+1:2n]
vz = v[2n+1:3n]
dvx = dvn[0n+1:1n]
dvy = dvn[1n+1:2n]
dvz = dvn[2n+1:3n]
#* Normal compontent of the viscous flow
v_n0   = LGM.Nd'*v ./ scalings  #? why is the scaling needed
# v_n0   = LGM.Nd'*v
#* Computing the tangential velocity by substracting the normal information
v_t = v - LGM.Nd*v_n0
vt_sum = sqrt.(v_t[0n+1:1n].^2 + v_t[1n+1:2n].^2 + v_t[2n+1:3n].^2)
#! Skipping points. Useful when the mesh is large
K = 1
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",markersize=2)
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
ylabel!("Re(p)");  title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_rAN_V[perm]),label=false,linewidth=2)
ylabel!("Re(Vn)");
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label=false,linewidth=2)
xlabel!("Angle"); ylabel!("Re(Vt)")
plot(plt1,plt2,plt3,layout=(3,1),dpi=500)
# savefig("/Users/mpasc/Dropbox/Apps/ShareLaTeX/JTCA_Iterative_losses/figures/notes/test.png")

#===========================================================================================
                        Plotting solutions - Comparisons
===========================================================================================#
p_a = p_b
# p_a = pa
p_h   = -LGM.tau_a/LGM.tau_h * p_a
dp_h  = -gmres(LGM.Gh,LGM.Hh*p_h)
#* 517 iterations for 69k DOF at freq=1000
dp_a  = -gmres(LGM.Ga,LGM.Ha*p_a;verbose=true) # <- This is the bottleneck...
#! The bottleneck could be the reason for the extended system not working?
v   = v0 - (LGM.phi_a*LGM.Dc*p_a +
            LGM.phi_a*LGM.Nd*dp_a +
            LGM.phi_h*LGM.Dc*p_h +
            LGM.phi_h*LGM.Nd*dp_h)

vx = v[0n+1:1n]
vy = v[1n+1:2n]
vz = v[2n+1:3n]
dvx = dvn[0n+1:1n]
dvy = dvn[1n+1:2n]
dvz = dvn[2n+1:3n]

#* Normal compontent of the viscous flow
# v_n0   = LGM.Nd'*v ./ scalings  #? why is the scaling needed
v_n0   = LGM.Nd'*v
#* Computing the tangential velocity by substracting the normal information
v_t = v - LGM.Nd*v_n0
vt_sum = sqrt.(v_t[0n+1:1n].^2 + v_t[1n+1:2n].^2 + v_t[2n+1:3n].^2)
#
K = 1 #! Skipping points. Useful when the mesh is large
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",markersize=2)
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
ylabel!("Re(p)");  title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_rAN_V[perm]),label=false,linewidth=2)
ylabel!("Re(Vn)");
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label=false,linewidth=2)
xlabel!("Angle"); ylabel!("Re(Vt)")
plot(plt1,plt2,plt3,layout=(3,1),dpi=500)
