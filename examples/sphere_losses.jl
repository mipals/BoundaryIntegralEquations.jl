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
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes")
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_fine");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_finer");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_finest");
tri_mesh_file = joinpath(mesh_path,"sphere_1m_35k");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_77k");
# tri_mesh_file = joinpath(mesh_path,"examples/meshes/sphere_1m_77k")
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[2])

# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/test_binary.ply");
#==========================================================================================
                                3d Visualization
==========================================================================================#
# using Meshes, MeshViz
# import WGLMakie as wgl
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
# @time BB = IntegralEquations.LossyBlockMatrix(mesh,freq;blockoutput=true,depth=1)
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
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
# Creating the "outer"-struct representation of the system matrix
# outer = LossyOneVariableOuter(mesh,BB,freq;fmm_on=true)
outer = LossyOneVariableOuter(mesh,freq)
bt    = compute_lossy_rhs(outer,vn0,vt10,vt20)
# Solving the problem using the "outer" struct
@time pa,pa_hist = gmres(outer,bt;verbose=true,log=true)
# Generating analytical solution

n = length(vn0)
coordinates = [radius*ones(n,1) acos.(xyzb[3,:]/radius)]
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                    BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)

# Plotting
K = 10
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM-global",marker=:cross,markersize=2,color=:black)
# scatter!(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:red)
ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
#===========================================================================================
                                Checking results
===========================================================================================#
ghpainv,gh_hist = gmres(outer.Gh,outer.Hh*pa;verbose=true,log=true)
gapainv,ga_hist = gmres(outer.Ga,outer.Ha*pa;verbose=true,log=true)
ghpainv = ghpainv*ϕₕ*τₐ/τₕ
gapainv = gapainv*ϕₐ
# gapainv = gmres!(x0,outer.Ga,outer.Ha*pa;verbose=true)*ϕₐ
vx = normals[1,:] .* vn0 + tangent1[1,:] .* vt10 + tangent2[1,:] .* vt20 -
    (normals[1,:] .* (gapainv - ghpainv)         +
   (tangent1[1,:] .* (outer.Dt1*pa)              +
    tangent2[1,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vy = normals[2,:] .* vn0 + tangent1[2,:] .* vt10 + tangent2[2,:] .* vt20 -
    (normals[2,:] .* (gapainv - ghpainv)         +
   (tangent1[2,:] .* (outer.Dt1*pa)              +
    tangent2[2,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vz = normals[3,:] .* vn0 + tangent1[3,:] .* vt10 + tangent2[3,:] .* vt20 -
    (normals[3,:] .* (gapainv - ghpainv)         +
   (tangent1[3,:] .* (outer.Dt1*pa)              +
    tangent2[3,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))

# Local components of the viscous velocity on the boundary
v_n0   = -(normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz) # The normal points inwards
v_t1   =  tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   =  tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)
#===========================================================================================
                            Testing new approach
===========================================================================================#
# using Test
# Dx,Dy,Dz = BoundaryIntegralEquations.shape_function_derivatives(mesh;global_derivatives=true)
# ph = -outer.τₐ/outer.τₕ*pa
# nx = outer.nx
# ny = outer.ny
# nz = outer.nz
# dpa = gmres(outer.Ga,outer.Ha*pa;verbose=true)
# dph = gmres(outer.Gh,outer.Hh*ph;verbose=true)
# vx0 = mesh.normals[1,:]*0.0
# vy0 = mesh.normals[2,:]*0.0
# vz0 = mesh.normals[3,:]*u₀
# # Velocity coupling seems true
# @test abs.(outer.ϕₐ*(Dx*pa + Diagonal(nx)*dpa) + outer.ϕₕ*(Dx*ph + Diagonal(nx)*dph) + vx) ≈ vx0 atol=1e-5
# @test abs.(outer.ϕₐ*(Dy*pa + Diagonal(ny)*dpa) + outer.ϕₕ*(Dy*ph + Diagonal(ny)*dph) + vy) ≈ vy0 atol=1e-5
# @test abs.(outer.ϕₐ*(Dz*pa + Diagonal(nz)*dpa) + outer.ϕₕ*(Dz*ph + Diagonal(nz)*dph) + vz) ≈ vz0 atol=1e+0
# # Null-divergence?
# dvx = gmres(outer.Gv,outer.Hv*vx;verbose=true)
# dvy = gmres(outer.Gv,outer.Hv*vy;verbose=true)
# dvz = gmres(outer.Gv,outer.Hv*vz;verbose=true)
# # Also looks ok?
# @test abs.(Dx*vx + Diagonal(nx)*dvx + Dy*vy + Diagonal(ny)*dvy + Dz*vz + Diagonal(nz)*dvz) ≈ vx0 atol=1e-3
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
    BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,[radius*ones(M,1) acos.(xyzb[3,:]/radius)];S=1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)[1:K:end]

# Plotting
plt1 = scatter(ang_axis[1:K:end],abs.(pa[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|p|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
plt2 = scatter(ang_axis[1:K:end],abs.(v_n0[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2)
# ylims!((-5e-6,5e-6))
plt3 = scatter(ang_axis[1:K:end],abs.(vt_sum[1:K:end]),label="BEM",marker=:cross,markersize=2,color=:black)
# scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|Vt|")
plot(plt1,plt2,plt3,layout=(3,1))
# plt = plot(plt1,plt2,plt3,layout=(3,1),dpi=300)
