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
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finest"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])
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
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
@time BB = IntegralEquations.LossyBlockMatrix(mesh,freq;blockoutput=true,depth=1)

# The error comes from the single-layer potential
Ga = FMMGOperator(mesh,ka;n=3,eps=1e-6,offset=0.2,nearfield=true)
xa = randn(size(Ga,1))
maximum(abs.(Ga*xa - BB.Bₐ*xa))
x1 = ones(ComplexF64,size(Ga,1))
maximum(abs.(Ga*x1 - BB.Bₐ*x1))


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
outer     = LossyOneVariableOuter(mesh,BB,freq;fmm_on=false)
outer_fmm = LossyOneVariableOuter(mesh,BB,freq;fmm_on=true,n=3)


import IntegralEquations: l_GH_r!, l_D_r!
vbn = vn0
vt1 = vt10
vt2 = vt20
l_GH_r!(outer.Gv,outer.luGv, outer.Hv, outer.nx, outer.ny, outer.nz, outer.tx, outer.ty, outer.tz,   outer.tmp1, outer.tmp2, outer.x1, vt1, true, outer.lu_on)
l_D_r!(outer.Dt1,        outer.tx, outer.ty, outer.tz, outer.tx, outer.ty, outer.tz, outer.tmp1, outer.tmp2,         outer.x1, vt1, false)
l_D_r!(outer.Dt2,        outer.sx, outer.sy, outer.sz, outer.tx, outer.ty, outer.tz, outer.tmp1, outer.tmp2,         outer.x1, vt1, false)
# Second multiplication
l_GH_r!(outer.Gv, outer.luGv, outer.Hv, outer.nx, outer.ny, outer.nz, outer.sx, outer.sy, outer.sz,   outer.tmp1, outer.tmp2, outer.x2, vt2, true, outer.lu_on)
l_D_r!(outer.Dt1,         outer.tx, outer.ty, outer.tz, outer.sx, outer.sy, outer.sz, outer.tmp1, outer.tmp2,         outer.x2, vt2, false)
l_D_r!(outer.Dt2,         outer.sx, outer.sy, outer.sz, outer.sx, outer.sy, outer.sz, outer.tmp1, outer.tmp2,         outer.x2, vt2, false)


bt        = compute_lossy_rhs(outer,vn0,vt10,vt20);
bt_fmm    = compute_lossy_rhs(outer_fmm,vn0,vt10,vt20)
bt ./ bt_fmm
# Solving the problem using the "outer" struct
@time pa     = gmres(outer,bt;verbose=true)
@time pa_fmm = gmres(outer_fmm,bt;verbose=true,maxiter=50)
pa ./ pa_fmm
#===========================================================================================
                                Checking results
===========================================================================================#
ghpainv     = gmres(outer.Gh,outer.Hh*pa;verbose=true)*ϕₕ*τₐ/τₕ
gapainv     = gmres(outer.Ga,outer.Ha*pa;verbose=true)*ϕₐ
ghpainv_fmm = gmres(outer_fmm.Gh,outer_fmm.Hh*pa_fmm;verbose=true)*ϕₕ*τₐ/τₕ
# gapainv_fmm = gmres(outer_fmm.Ga,outer_fmm.Ha*pa_fmm;verbose=true)*ϕₐ
# ghpainv ./ ghpainv_fmm
ghpainv = ghpainv_fmm
# gapainv = gapainv_fmm
# gapainv ./ gapainv_fmm
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
                                   Plotting solutions
===========================================================================================#
# Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
    IntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,[radius*ones(M,1) acos.(xyzb[3,:]/radius)];S=1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)

# Plotting
# plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|p|"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
# plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|Vn|"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
# ylims!((-1e-6,1e-6))
# plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|Vt|")
# plot(plt1,plt2,plt3,layout=(3,1))
# plt = plot(plt1,plt2,plt3,layout=(3,1),dpi=300)
