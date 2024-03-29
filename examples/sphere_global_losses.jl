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
tri_physics_orders  = [:linear,:geometry,:disctrilinear,:disctriquadratic]
# Triangular Meshes
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes")
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_coarser");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_fine");
tri_mesh_file = joinpath(mesh_path,"sphere_1m_finer");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_finest");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_35k");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_77k");
# tri_mesh_file = joinpath(mesh_path,"sphereTri6_ref3");
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[3])

# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/test_binary.ply");
#==========================================================================================
                                    3d Visualization
==========================================================================================#
# corners = sort!(unique(mesh.physics_topology[1:3,:]))
# edges   = sort!(unique(mesh.physics_topology[4:6,:]))
# using MeshViz
# import WGLMakie as wgl
# wgl.set_theme!(resolution=(800, 800))
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
# # wgl.scatter!(mesh.coordinates[1,:],mesh.coordinates[2,:],mesh.coordinates[3,:])
# wgl.scatter!(mesh.targets[1,corners],mesh.targets[2,corners],mesh.targets[3,corners],color=:red)
# wgl.scatter!(mesh.targets[1,edges],mesh.targets[2,edges],mesh.targets[3,edges],color=:blue)

# element = 1
# element_coordinates = mesh.coordinates[:,mesh.topology[:,element]]

# viz(simple_mesh)
# wgl.lines!(element_coordinates[1,:],element_coordinates[2,:],element_coordinates[3,:])
# wgl.scatter!(element_coordinates[1,:],element_coordinates[2,:],element_coordinates[3,:])

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
# LGM = BoundaryIntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true,depth=2,n=3)
# LGM = BoundaryIntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true,depth=1)
# LGM = BoundaryIntegralEquations.LossyGlobalOuter(mesh,freq;fmm_on=true,depth=2,sparse_lu=true)
LGM = LossyGlobalOuter(mesh,freq;fmm_on=true,depth=2,sparse_lu=true);
# LGM = LossyGlobalOuter(mesh,freq;fmm_on=true,depth=2);
# cond(Matrix(LGM.Gh))
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
if typeof(LGM.Gv) <: Factorization
    rhs = LGM.Ga*gmres(LGM.inner,LGM.Dr*v0 - LGM.Nd'*(LGM.Gv\(LGM.Hv*v0));verbose=true);
else
    rhs = LGM.Ga*gmres(LGM.inner,LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0);verbose=true);
end
# @time pa = gmres(LGM,rhs;verbose=true,maxiter=2);
@time pa = gmres(LGM,rhs;verbose=true);
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
ph   = -LGM.tau_a/LGM.tau_h * pa
dpa  = -gmres(LGM.Ga,LGM.Ha*pa;verbose=true) # <- This is the bottleneck...
if typeof(LGM.Gh) <: Factorization
    dph  = -(LGM.Gh\(LGM.Hh*ph))
else
    dph  = -gmres(LGM.Gh,LGM.Hh*ph)
end
v   = v0 - (LGM.phi_a*LGM.Dc*pa  +
            LGM.phi_a*LGM.Nd*dpa +
            LGM.phi_h*LGM.Dc*ph  +
            LGM.phi_h*LGM.Nd*dph)
if typeof(LGM.Gv) <: Factorization
    dvn = -(LGM.Gv\(LGM.Hv*v))
else
    dvn = -gmres(LGM.Gv,LGM.Hv*v;verbose=true)
end
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
v_n0   = LGM.Nd'*v
#* Computing the tangential velocity by substracting the normal information
v_t = v + LGM.Nd*v_n0
vt_sum = sqrt.(v_t[0n+1:1n].^2 + v_t[1n+1:2n].^2 + v_t[2n+1:3n].^2)
#! Skipping points. Useful when the mesh is large
K = 1
plt1 = scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM",markersize=2)
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
ylabel!("Re(p)");  title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],real.(v_n0[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_rAN_V[perm]),label=false,linewidth=2)
ylabel!("Re(Vn)")
plt3 = scatter(ang_axis[1:K:end],real.(vt_sum[1:K:end]),label=false,markersize=2)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label=false,linewidth=2)
xlabel!("Angle"); ylabel!("Re(Vt)")
plot(plt1,plt2,plt3,layout=(3,1),dpi=500)

#===========================================================================================
                                    Tests!
                        Null-divergence is not satisfied?
===========================================================================================#
# using Test
# #* Checking constants
# @test LGM.mu_a == LGM.phi_a - LGM.phi_h*LGM.tau_a/LGM.tau_h
# @test LGM.mu_h == LGM.phi_h*LGM.tau_a/LGM.tau_h
# #* All of these are solved quantities - So they're satisfied by construction
# @test LGM.Hh*ph ≈ -LGM.Gh*dph       #! Checking Thermal BEM System
# @test LGM.Ha*pa ≈ -LGM.Ga*dpa       #! Checking Acoustics BEM system
# @test LGM.Hv*v  ≈ -LGM.Gv*dvn       #! Checking viscous BEM system
# @test LGM.tau_a*pa ≈ -LGM.tau_h*ph  #! Checking Isothermal condition
# #* These are the extra conditions
# @test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0
#! Checking null-divergence (off by a factor??)
# @test LGM.Dr*v ≈ -LGM.Nd'*dvn
# #* We can fix this by solving a rectangular system instead
# rhs_v = -[LGM.Hv*v; LGM.Dr*v]
# V    = [LGM.Gv; LGM.Nd']
# dvn2 = lsqr(V,rhs_v;atol=1e-15,btol=1e-13,verbose=true);
# @test LGM.Dr*v ≈ -LGM.Nd'*dvn2            # Now this is satisfied
# @test LGM.Hv*v ≈ -LGM.Gv*dvn2 atol = 1e-5 # But the BEM system has a higher error

# For discontinuous linear elements something is weird...
c1 = sort(unique(mesh.physics_topology[1,:]))
c2 = sort(unique(mesh.physics_topology[2,:]))
c3 = sort(unique(mesh.physics_topology[3,:]))
scatter(ang_axis[c1],real(v_rAN_V[c1]))
scatter!(ang_axis[c1],real.(v_n0[c1]))
scatter!(ang_axis[c2],real.(v_n0[c2]))
scatter!(ang_axis[c3],real.(v_n0[c3]))

using Statistics
quantile(real.(v_rAN_V[c1]) ./ real.(v_n0[c1]))
quantile(real.(v_rAN_V[c2]) ./ real.(v_n0[c2]))
quantile(real.(v_rAN_V[c3]) ./ real.(v_n0[c3]))


c4 = sort(unique(mesh.physics_topology[4,:]))
c5 = sort(unique(mesh.physics_topology[5,:]))
c6 = sort(unique(mesh.physics_topology[6,:]))
scatter!(ang_axis[c4],real.(v_n0[c4]))
scatter!(ang_axis[c5],real.(v_n0[c5]))
scatter!(ang_axis[c6],real.(v_n0[c6]))

quantile(real.(v_rAN_V[c4]) ./ real.(v_n0[c4]))
quantile(real.(v_rAN_V[c5]) ./ real.(v_n0[c5]))
quantile(real.(v_rAN_V[c6]) ./ real.(v_n0[c6]))



# fun
nx = mesh.normals[1,:]
ny = mesh.normals[2,:]
nz = mesh.normals[3,:]
tx = mesh.tangents[1,:]
ty = mesh.tangents[2,:]
tz = mesh.tangents[3,:]
sx = mesh.sangents[1,:]
sy = mesh.sangents[2,:]
sz = mesh.sangents[3,:]
D = [sparse(Diagonal(tx)) sparse(Diagonal(ty)) sparse(Diagonal(tz));
     sparse(Diagonal(sx)) sparse(Diagonal(sy)) sparse(Diagonal(sz));
     sparse(Diagonal(nx)) sparse(Diagonal(ny)) sparse(Diagonal(nz))]

T = D'*D # Should be the identity matrix?
x = rand(size(T,1))
norm(T*x - x)./norm(x)
sum(T,dims=2)

S = D*D'
x = rand(size(S,1))
norm(S*x - x)./norm(x)
sum(S,dims=2)

# T[abs.(T .- 1.0) .< 1e-15]

D
D2 = [sparse(Diagonal(tx)) sparse(Diagonal(ty)) sparse(Diagonal(tz));
      sparse(Diagonal(sx)) sparse(Diagonal(sy)) sparse(Diagonal(sz));
      spzeros(n,3n)]

M = D'*D2
sum(M,dims=2)

N = D*LGM.Dc


data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=abs.(data_viz))
wgl.Colorbar(fig[1,2],label="|p|"); fig
