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
# mesh_files = ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer"]
# mesh_files = ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m","sphere_1m_fine"]
mesh_files = ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m"] # Just as a start

# mesh_files = ["sphere_1m_coarser"]

for mesh_file in mesh_files

    tri_mesh_file = joinpath(mesh_path,mesh_file);
    @time mesh = load3dTriangularComsolMesh(tri_mesh_file)
    #==========================================================================================
                                    Setting up constants
    ==========================================================================================#
    freq   = 100.0                                   # Frequency                 [Hz]
    rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
    k      = 2*π*freq/c                              # Wavenumber                [1/m]
    radius = 1.0                                     # Radius of sphere_1m       [m]
    ω = 2π*freq
    depth = 1

    xyzb = mesh.sources
    n = size(xyzb,2)
    u₀ = 1e-2
    normals  = mesh.normals
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
    #===========================================================================================
                            Computing Condition Numbers
    ===========================================================================================#
    if Sys.isapple()
        data_directory = "/Users/mpasc/OneDrive - Danmarks Tekniske Universitet/paper1"
    elseif Sys.islinux()
        data_directory = "/home/mpasc/OneDrive/paper1"
    end
        for depth in [1]
        data_file = joinpath(data_directory, mesh_file, "results_depth$(depth)_freq$(Int(freq)).jld2")
        if !isfile(data_file)
            @info "Computing setup: $(data_file)"
            # Computing matrices
            LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=depth);
            # Creating the different systems
            F10 = BoundaryIntegralEquations._full10(LGM)
            F4  = BoundaryIntegralEquations._full4(LGM)
            F1  = BoundaryIntegralEquations._full1(LGM)
            #
            @info "Computing Conditon Number of 1x1 (n=$(n))"
            condF1  = cond(F1)
            @info "Computing Conditon Number of 4x4 (n=$(n))"
            condF4  = cond(F4)
            @info "Computing Conditon Number of 10x10 (n=$(n))"
            condF10 = cond(F10)

            # Plottin for a quick sanity check
            v0 = [u₀*ones(n); zeros(2n)]
            coordinates = [radius*ones(n,1) acos.(xyzb[1,:]/radius)]
            rhs = LGM.Ga*gmres(LGM.inner,LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0);verbose=true);
            pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ)
            ang_axis = coordinates[:,2]*180.0/pi
            sol1  = F1\rhs
            sol4  = F4\[zeros(n); v0]
            sol10 = F10\[zeros(7n);v0]
            # Extracting pressures
            pa1  = sol1
            pa4  = sol4[1:n]
            pa10 = sol10[1:n]

            perm = sortperm(ang_axis)
            gr(size=(600,300))
            scatter(ang_axis,real.(pa1),label="BEM-1",marker=:square,
                    markersize=2,dpi=600,xtickfontsize=15,ytickfontsize=15,
                    legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            scatter!(ang_axis,real.(pa4),label="BEM-4",marker=:triangle,
                    markersize=2,xtickfontsize=15,ytickfontsize=15,
                    legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            scatter!(ang_axis,real.(pa10),label="BEM-10",marker=:cross,
                    markersize=2,xtickfontsize=15,ytickfontsize=15,
                    legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            ylabel!("Re(p)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
            title!("Frequency = $(freq)"); xlabel!("Angle")
            savefig(joinpath(data_directory, mesh_file,"all_pa_depth$(depth)_freq$(Int(freq)).png"))

            # Computing normal velocity
            Aview = @view F4[1n+1:4n,0n+1:1n] # Because I am lazy i use the bottom equation from F4 directly
            v1  = v0 - Aview*pa1    # This is the bottom equation of F4
            v4  = sol4[1n+1:4n]
            v10 = sol10[4n+1:7n]

            v_n1  = LGM.Nd'*v1
            v_n4  = LGM.Nd'*v4
            v_n10 = LGM.Nd'*v10

            scatter(ang_axis,real.(v_n1),label="BEM-1",marker=:square,
                markersize=2,dpi=600,xtickfontsize=15,ytickfontsize=15,
                legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            scatter!(ang_axis,real.(v_n4),label="BEM-4",marker=:triangle,
                markersize=2,xtickfontsize=15,ytickfontsize=15,
                legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            scatter!(ang_axis,real.(v_n10),label="BEM-10",marker=:cross,
                markersize=2,xtickfontsize=15,ytickfontsize=15,
                legendfontsize=15,titlefontsize=15,xlabelfontsize=15,ylabelfontsize=15)
            ylabel!("Re(p)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
            title!("Frequency = $(freq)"); xlabel!("Angle")
            savefig(joinpath(data_directory, mesh_file,"all_vn_depth$(depth)_freq$(Int(freq)).png"))

            # Saving data
            jldsave(data_file,
                    condF1=condF1,
                    condF4=condF4,
                    condF10=condF1,
                    pa1=pa1,
                    pa4=pa4,
                    pa10=pa10,
                    v_n1=v_n1,
                    v_n4=v_n4,
                    v_n10=v_n10,
                    ang_axis=ang_axis,
                    pasAN=pasAN
                    )
        else
            @warn "Loading data "
            f = jldopen(data_file)
            condF1 = f["condF1"]
            condF4 = f["condF4"]
            condF10 = f["condF10"]
            close(f)
        end
    end
end
