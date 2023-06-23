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
                                                       physics_order=tri_physics_orders[1])
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
                        BEM matrix assembly and direct solution of the 10-variable system
===========================================================================================#
BB = LossyBlockMatrix(mesh,freq;blockoutput=false)
#cond(Matrix(BB))
rhs = [zeros(5M); # Ap + Bv = 0  (for all 3 modes)
        zeros(M); # τₐpₐ + τₕpₕ = 0  (fsothermal coupling)
             vn0; # ϕₐ(∇n)pₐ  + ϕₕ(∇n)pₕ  + (vn)ᵥ  = vn (no-slip)
            vt10;
            vt20;
        zeros(M)] # vx + vy + vz = 0 (Null Divergence)
sol = Matrix(BB)\rhs # gmres doesn't work
pa = sol[0M+1:1M]

# Generating analytical solution
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)]
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
ang_axis = coordinates[:,2]*180.0/pi
perm = sortperm(ang_axis)

# Plotting
K = 1
gr(size=(600,500))
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 10n",marker=:x,markersize=2,color=:black,dpi=400);
ylabel!("Re(Pa) [Pa]"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]");
#savefig("pa10x10_$(M)DOFs_$(freq)Hz.png")
#===========================================================================================
                                                Reconstructing unknowns
===========================================================================================#
vx = sol[4*M+1:5*M]
vy = sol[6*M+1:7*M]
vz = sol[8*M+1:9*M]

# Local components of the viscous velocity on the boundary
v_n0   =  normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz 
v_t1   =  tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   =  tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2) # inconsistency doesn't matter due to summation
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
#savefig("all10x10_Real_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("Imag(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all10x10_Imag_$(M)DOFs_$(freq)Hz.png")

plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(freq) Hz");
plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(-v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue); # minus due to inconsistency with analytical solution >> v_rAN_V pointing outwards
plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle"); ylabel!("|Vt|");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))
#savefig("all10x10_Abs_$(M)DOFs_$(freq)Hz.png")