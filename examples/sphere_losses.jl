#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using SpecialFunctions, LinearAlgebra, IntegralEquations, Plots, IterativeSolvers
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/sphere_1m"
tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                physics_order=tri_physics_orders[2])
#==========================================================================================
        3d Visualization - Seems highly unstable on M1 chips. Problems with GLMakie?
==========================================================================================#
# using Meshes, MeshViz
# ##choose a Makie backend
# import GLMakie as Mke
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 500.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=-1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
@time BB = IntegralEquations.LossyBlockMatrix(mesh,freq;blockoutput=true,depth=0)
# Fa,Ba,C0 = assemble_parallel!(mesh,ka,mesh.sources;m=3,n=3)
# cond(Matrix(BB.Aᵥ))
# cond(Matrix(BB.Bᵥ))
# cond(Matrix(BB.Aₕ))
# cond(Matrix(BB.Bₕ))

xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=-1)

# Creating RHS
vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
# Creating the "outer"-struct representation of the system matrix
outer = LossyOneVariableOuter(mesh,BB,freq)
bt    = compute_lossy_rhs(outer,vn0,vt10,vt20)
# Solving the problem using the "outer" struct
@time pa,pa_hist = gmres(outer,bt;verbose=true,log=true)
#===========================================================================================
                                Checking results
===========================================================================================#
ghpainv,gh_hist = gmres(outer.Gh,outer.Hh*pa;verbose=true,log=true)
gapainv,ga_hist = gmres(outer.Ga,outer.Ha*pa;verbose=true,log=true)
ghpainv = ghpainv*ϕₕ*τₐ/τₕ
gapainv = gapainv*ϕₐ
# gapainv = gmres!(x0,outer.Ga,outer.Ha*pa;verbose=true)*ϕₐ
vx = normals[1,:] .* vn0 + tangent1[1,:] .* vt10 + tangent2[1,:] .* vt20 -
    (normals[1,:] .* (gapainv - ghpainv)       +
   (tangent1[1,:] .* (outer.Dt1*pa)   +
    tangent2[1,:] .* (outer.Dt2*pa))*(ϕₐ - ϕₕ*τₐ/τₕ))
vy = normals[2,:] .* vn0 + tangent1[2,:] .* vt10 + tangent2[2,:] .* vt20 -
    (normals[2,:] .* (gapainv - ghpainv)       +
   (tangent1[2,:] .* (outer.Dt1*pa)   +
    tangent2[2,:] .* (outer.Dt2*pa))*(ϕₐ - ϕₕ*τₐ/τₕ))
vz = normals[3,:] .* vn0 + tangent1[3,:] .* vt10 + tangent2[3,:] .* vt20 -
    (normals[3,:] .* (gapainv - ghpainv)         +
   (tangent1[3,:] .* (outer.Dt1*pa)      +
    tangent2[3,:] .* (outer.Dt2*pa))*(ϕₐ - ϕₕ*τₐ/τₕ))

# Local components of the viscous velocity on the boundary
v_n0   =  normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz
v_t1   = tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   = tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
    IntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,[radius*ones(M,1) acos.(xyzb[3,:]/radius)];S=-1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)

# Plotting
plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:yellow)
# scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:yellow)
ylabel!("|p|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",markersize=2)
title!("Frequency = $(freq)")
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:yellow)
# scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:yellow)
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",markersize=2)
# ylims!((0,1e-5))
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:yellow)
# scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:yellow)
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",markersize=2)
xlabel!("Angle"); ylabel!("|Vt|")
# plot(plt1,plt2,plt3,layout=(3,1))
plt = plot(plt1,plt2,plt3,layout=(3,1),dpi=300)
