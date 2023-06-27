using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            default = 100.0
            arg_type = Float64
        "--compute_full_solution"
            help = raw"computing full solution"
            default = false
            arg_type = Bool
        "--mesh_file"
            help = "mesh file"
            default = "sphere_1m_coarser"
            arg_type = String

    end
    return parse_args(s)
end
function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return parsed_args
end

intputArguments = main()
freq = intputArguments["freq"]
compute_full_solution = intputArguments["compute_full_solution"]
mesh_file = intputArguments["mesh_file"]

#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using Plots
using IterativeSolvers
using BenchmarkTools
using JLD2
using Statistics, DataFrames


### Auxiliary functions for benchmarking
preprocess_trial(t::BenchmarkTools.Trial, id::AbstractString) = (id=id,
        minimum=minimum(t.times),
        median=median(t.times),
        maximum=maximum(t.times),
        allocations=t.allocs,
        memory_estimate=t.memory)

function denseSolve(LGM::LossyGlobalOuter,rhs::Vector{Float64})
    LGM_dense  = BoundaryIntegralEquations._full4(LGM);
    sol = LGM_dense\rhs;
end

function reconstructUnknowns(LGM::LossyGlobalOuter,sol::Vector{ComplexF64},M::Int64)
    # Local components of the viscous velocity on the boundary
    v_n0 = LGM.Nd'*sol[1M+1:4M]; 
    v_t = sol[1M+1:4M] + LGM.Nd*v_n0; # Computing the tangential velocity by substracting the normal information
    vt_sum = sqrt.(v_t[0M+1:1M].^2 + v_t[1M+1:2M].^2 + v_t[2M+1:3M].^2);
end


#=============================ß============================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
#bool_reconstruction = true;
# Triangular Meshes

# tri_mesh_file = "examples/meshes/sphere_1m_coarser"
# tri_mesh_file = "examples/meshes/sphere_1m_coarse"
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finest"
# tri_mesh_file = "examples/meshes/sphere_1m_35k"
# tri_mesh_file = "examples/meshes/sphere_1m_77k"
#mesh_files =  ["sphere_1m_coarser"]#,"sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer"];
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");

#for mesh_file in mesh_files
#mesh_file = mesh_files[1]
tri_mesh_file = joinpath(mesh_path,mesh_file);
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                    physics_order=tri_physics_orders[2]);
#==========================================================================================
                                Setting up constants
==========================================================================================#
#freq   = 1000.0;                                   # Frequency                 [Hz]
#==========================================================================================
                            Creating excitation vector
==========================================================================================#
xyzb = mesh.sources;
M  = size(xyzb,2);
u₀ = 1e-2;
v0 = [zeros(2M); u₀*ones(M)];
#===========================================================================================
                        BEM matrix assembly and (iterative) solution of the 1-variable system
===========================================================================================#

output = DataFrame()
@info "Computing LGM"
push!(output, preprocess_trial(@benchmark(LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false),evals=20), "tLGM"))
LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=1,n=3,progress=false);

@info "Computing RHS"
push!(output,preprocess_trial(@benchmark([zeros(M); v0],evals=20),"tRHS"))
rhs = [zeros(M); v0];

# dense
@info "Assembling and solving dense system"
push!(output, preprocess_trial(@benchmark(denseSolve(LGM,rhs),evals=20), "tDense"))
LGM_dense  = BoundaryIntegralEquations._full4(LGM);
sol = LGM_dense\rhs;
pa = sol[1:M];
#cond(LGM_dense)

# Generating analytical solution
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1);
radius = 1.0;                                     # Radius of sphere_1m       [m]
coordinates = [radius*ones(M,1) acos.(xyzb[3,:]/radius)];
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
                BoundaryIntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,coordinates;S=1,kv=kᵥ);
ang_axis = coordinates[:,2]*180.0/pi;
perm = sortperm(ang_axis);

# Plotting pressure
K = 1;
gr(size=(600,500));
scatter(ang_axis[1:K:end],real.(pa[1:K:end]),label="BEM - 4n Dense",marker=:cross,markersize=2,color=:black,dpi=400);
plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=1,color=:blue);
ylabel!("Re(Pa)");
title!("Frequency = $(freq) Hz");
xlabel!("Angle [deg]");
savefig("paGlobal4x4_$(M)DOFs_$(Int(freq))Hz.png")

if compute_full_solution == true
    #===========================================================================================
                                    Reconstructing unknowns
    ===========================================================================================#
    @info "Reconstructing unknowns"
    push!(output, preprocess_trial(@benchmark(reconstructUnknowns(LGM,sol,M),evals=20), "tRec"))
    v = sol[1M+1:4M]; 
    # Local components of the viscous velocity on the boundary
    v_n0   = LGM.Nd'*v; 
    v_t = v + LGM.Nd*v_n0; # Computing the tangential velocity by substracting the normal information
    vt_sum = sqrt.(v_t[0M+1:1M].^2 + v_t[1M+1:2M].^2 + v_t[2M+1:3M].^2);
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
    savefig("allGlobal4x4_Real_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,imag.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Pa)"); plot!(ang_axis[perm],imag.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,imag.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("Imag(Vn)"); plot!(ang_axis[perm],imag.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,imag.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],imag.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("Imag(Vt)");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal4x4_Imag_$(M)DOFs_$(Int(freq))Hz.png")

    plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Pa|"); plot!(ang_axis[perm],abs.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
    title!("Frequency = $(freq) Hz");
    plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
    ylabel!("|Vn|"); plot!(ang_axis[perm],abs.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
    plot!(ang_axis[perm],abs.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
    xlabel!("Angle [deg]"); ylabel!("|Vt|");
    plt4 = plot(plt1,plt2,plt3,layout=(3,1))
    savefig("allGlobal4x4_Abs_$(M)DOFs_$(Int(freq))Hz.png")

    # Saving data
    jldsave("runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    runtimes=output,
    pa=pa,
    pasAN=pasAN,
    v_rAN_V=v_rAN_V,
    v_thetaAN_V=v_thetaAN_V,
    ang_axis=ang_axis,
    v_n0=v_n0,
    vt_sum=vt_sum
    )

    print(output)

else
    # Saving data
    jldsave("runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2", 
    runtimes=output,
    pa=pa,
    pasAN=pasAN,
    ang_axis=ang_axis
    );

    print(output)
end

#end