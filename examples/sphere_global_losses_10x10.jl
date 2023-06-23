#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using Plots
using IterativeSolvers
#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/sphere_1m_coarser"
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
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
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals

v0 = [zeros(2M); u₀*ones(M)]
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3);
rhs = [zeros(7M);v0];

# dense format
LGM_dense  = BoundaryIntegralEquations._full10(LGM);
#cond(LGM_dense)
sol = LGM_dense\rhs;
pa = sol[1:M];

# Generating analytical solution
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
radius = 1.0                                     # Radius of sphere_1m       [m]
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)]
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = coordinates[:,2]*180.0/pi
perm = sortperm(ang_axis)

# Plotting pressure
K = 1
gr(size=(600,500))
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 10n",marker=:cross,markersize=2,color=:black,dpi=400);
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=1,color=:blue);
ylabel!("Re(Pa)");
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]");
#savefig("paGlobal10x10_$(M)DOFs_$(freq)Hz.png")
#===========================================================================================
                                Reconstructing unknowns
===========================================================================================#
v = sol[4M+1:7M]
# Local components of the viscous velocity on the boundary
v_n0   = LGM.Nd'*v 
v_t = v + LGM.Nd*v_n0 # Computing the tangential velocity by substracting the normal information
vt_sum = sqrt.(v_t[0M+1:1M].^2 + v_t[1M+1:2M].^2 + v_t[2M+1:3M].^2)
#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Plotting
plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("Re(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal10x10_Real_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal10x10_Imag_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("|Vt|");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
savefig("allGlobal10x10_Abs_$(M)DOFs_$(freq)Hz.png")