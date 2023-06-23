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
tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finest"
# tri_mesh_file = "examples/meshes/sphere_1m_35k"
# tri_mesh_file = "examples/meshes/sphere_1m_77k"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[2])
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 500.0                                   # Frequency                 [Hz]
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents

vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#
# fmm of
BB = LossyBlockMatrix(mesh,freq;blockoutput=true)
outer = LossyOneVariableOuter(mesh,BB,freq;fmm_on=false)
#cond(Matrix(outer))
rhs    = compute_lossy_rhs(outer,vn0,vt10,vt20)
pa = gmres(outer,rhs;verbose=true);
#pa = Matrix(outer)\rhs

# Generating analytical solution
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)]
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = coordinates[:,2]*180.0/pi
perm = sortperm(ang_axis)

# Plotting pressure
K = 1
gr(size=(600,500))
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 1n",marker=:x,markersize=2,color=:black,dpi=400);
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=1,color=:blue);
ylabel!("Re(Pa)");
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]")
# savefig("pa1x1_$(M)DOFs_$(freq)Hz.png")
#===========================================================================================
                                                Reconstructing unknowns
===========================================================================================#
ph  = - outer.τₐ/outer.τₕ * pa
dph,dph_hist = gmres(outer.Gh,outer.Hh*ph;verbose=true,log=true)
dpa,dpa_hist = gmres(outer.Ga,outer.Ha*pa;verbose=true,log=true)
vx = normals[1,:] .* vn0 + tangent1[1,:] .* vt10 + tangent2[1,:] .* vt20 -
    (normals[1,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)   +
   (tangent1[1,:] .* (outer.Dt1*pa)              +
    tangent2[1,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vy = normals[2,:] .* vn0 + tangent1[2,:] .* vt10 + tangent2[2,:] .* vt20 -
    (normals[2,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)    +
   (tangent1[2,:] .* (outer.Dt1*pa)              +
    tangent2[2,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vz = normals[3,:] .* vn0 + tangent1[3,:] .* vt10 + tangent2[3,:] .* vt20 -
    (normals[3,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)    +
   (tangent1[3,:] .* (outer.Dt1*pa)              +
    tangent2[3,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))

# Local components of the viscous velocity on the boundary
v_n0   =  normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz 
v_t1   =  tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   =  tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("Re(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1_Real_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("Imag(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1_Imag_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(-v_rAN_V[perm]),label="Analytical",linewidth=2); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2);
xlabel!("Angle"); ylabel!("|Vt|");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1_Abs_$(M)DOFs_$(freq)Hz.png")
#===========================================================================================
                                   Same with FMM (plus comparison)
===========================================================================================#
# solution
outer = LossyOneVariableOuter(mesh,freq)
rhs    = compute_lossy_rhs(outer,vn0,vt10,vt20)
pa_fmm = gmres(outer,rhs;verbose=true)

# FMM approximation error
eps = (sum(abs.(pa_fmm-pa).^2)/M).^(1/2)
refNorm = (sum(abs.(pa).^2)/M).^(1/2)
epsRel = eps/refNorm

# plotting pressure
K = 1
gr(size=(600,500))
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 1n",marker=:x,markersize=2,color=:black,dpi=400);
scatter!(ang_axis[1:K:end],real.(pa_fmm[1:K:end]),label="FMBEM",marker=:cross,markersize=2,color=:red,markerstrokewidth=0.5,dpi=400);
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=1,color=:blue);
ylabel!("Re(Pa)");
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]")
# savefig("pa1x1FMM_$(M)DOFs_$(freq)Hz.png")

# reconstruct unknowns
ph  = - outer.τₐ/outer.τₕ * pa_fmm
dph,dph_hist = gmres(outer.Gh,outer.Hh*ph;verbose=true,log=true)
dpa,dpa_hist = gmres(outer.Ga,outer.Ha*pa_fmm;verbose=true,log=true)
vx = normals[1,:] .* vn0 + tangent1[1,:] .* vt10 + tangent2[1,:] .* vt20 -
    (normals[1,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)   +
   (tangent1[1,:] .* (outer.Dt1*pa_fmm)              +
    tangent2[1,:] .* (outer.Dt2*pa_fmm))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vy = normals[2,:] .* vn0 + tangent1[2,:] .* vt10 + tangent2[2,:] .* vt20 -
    (normals[2,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)    +
   (tangent1[2,:] .* (outer.Dt1*pa_fmm)              +
    tangent2[2,:] .* (outer.Dt2*pa_fmm))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vz = normals[3,:] .* vn0 + tangent1[3,:] .* vt10 + tangent2[3,:] .* vt20 -
    (normals[3,:] .* (outer.ϕₐ*dpa + outer.ϕₕ*dph)    +
   (tangent1[3,:] .* (outer.Dt1*pa_fmm)              +
    tangent2[3,:] .* (outer.Dt2*pa_fmm))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))

# Local components of the viscous velocity on the boundary
v_n0_fmm  =  normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz 
v_t1  =  tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2  =  tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum_fmm = sqrt.(v_t1.^2 + v_t2.^2)


# plotting
plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,real.(pa_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,real.(v_n0_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,real.(vt_sum_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("Re(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1FMM_Real_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,imag.(pa_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,imag.(v_n0_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,imag.(vt_sum_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("Imag(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1FMM_Imag_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,abs.(pa_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,abs.(v_n0_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);# minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
scatter!(ang_axis,abs.(vt_sum_fmm),label="FMBEM",marker=:x,markersize=1,color=:red,markerstrokewidth=0.5);
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("|Vt|");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all1x1FMM_Abs_$(M)DOFs_$(freq)Hz.png")