#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using Plots
using IterativeSolvers
using JLD2
#==========================================================================================
                                Loading Mesh
==========================================================================================#
# Triangular Meshes
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes")
mesh_file = "sphere_1m_coarser"
mesh_file = "sphere_1m_coarse"
# mesh_file = "sphere_1m"
# mesh_file = "sphere_1m_fine"
# mesh_file = "sphere_1m_finer"
tri_mesh_file = joinpath(mesh_path,mesh_file);
@time mesh = load3dTriangularComsolMesh(tri_mesh_file)
#==========================================================================================
                                    3d Visualization
==========================================================================================#
# using MeshViz
# import WGLMakie as wgl
# wgl.set_theme!(resolution=(800, 800))
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
ω = 2π*freq
depth = 1
#===========================================================================================
                        Computing Condition Numbers
===========================================================================================#
if Sys.isapple()
    data_directory = "/Users/mpasc/OneDrive - Danmarks Tekniske Universitet/paper1"
else if Sys.islinux()

end
data_file = joinpath(data_directory, mesh_file, "condition_numbers_depth$(depth).jld2")
if !isfile(data_file)
    # Computing matrices
    LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=depth);
    # Creating the different systems
    F10 = BoundaryIntegralEquations._full10(LGM)
    F4  = BoundaryIntegralEquations._full4(LGM)
    F1  = BoundaryIntegralEquations._full1(LGM)
    #
    @info "Computing Conditon Number of 1x1"
    condF1  = cond(F1)
    @info "Computing Conditon Number of 4x4"
    condF4  = cond(F4)
    @info "Computing Conditon Number of 10x10"
    condF10 = cond(F10)

    jldsave(data_file, condF1=condF1,condF4=condF4,condF10=condF10)

    # Plottin for a quick sanity check
    K=1
    gr(size=(600,300))
    scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM-global",marker=:cross,
            markersize=2,color=:black,dpi=600,xtickfontsize=15,ytickfontsize=15,
            legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
    ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
    title!("Frequency = $(freq)")
    xlabel!("Angle")
    savefig(joinpath(data_directory, mesh_file,"pa.png"))
else
    @warn "Loading data "
    f = jldopen(dataFile)
    condF1 = f["condF1"]
    condF4 = f["condF4"]
    condF10 = f["condF10"]
end


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
dph  = -gmres(LGM.Gh,LGM.Hh*ph)
v   = v0 - (LGM.phi_a*LGM.Dc*pa  +
            LGM.phi_a*LGM.Nd*dpa +
            LGM.phi_h*LGM.Dc*ph  +
            LGM.phi_h*LGM.Nd*dph)
dvn = -gmres(LGM.Gv,LGM.Hv*v;verbose=true)
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
