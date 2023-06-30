using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--freq"
            help = "frequency"
            required = true
            arg_type = Float64
        "--depth"
            help = "depth of sparse matrices"
            arg_type = Int64
            default = 1
        "--fmm"
            help = "computing solution with the FMM"
            default = true
            arg_type = Bool
        "--hmatrix"
            help = "computing solution using the HMatrix"
            default = false
            arg_type = Bool
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
hmatrix_on = intputArguments["hmatrix"]
frequency = intputArguments["freq"]
fmm_on = intputArguments["fmm"]
depth = intputArguments["depth"]
#==========================================================================================
                                Loading packages
==========================================================================================#
# # Importing related packages
frequency = 100.0
depth = 1
fmm_on = true
hmatrix_on = false
using LinearAlgebra, BoundaryIntegralEquations, IterativeSolvers, Plots, JLD2
# # Setting up constants
ρ₀,c,_,_,_,kᵥ,_,_,_,_,_,_ = visco_thermal_constants(;freq=frequency,S=1);
a       = 1.0;                       # Radius of sphere        [m]
ω       = 2π*frequency;              # Angular frequency       [rad/s]
k       = 2π*frequency/c;            # Wavenumber              [1/m]
v₀      = 1e-2;                      # Initial velocity        [m/s]
# # Analytical solution
β = kᵥ*a; # As defined in (6.9.19) in Temkin
b = k*a;  # As defined in (6.9.19) in Temkin
# Now computing ``A_1``using (6.9.25)
A₁ = -v₀/(3im*k)*b^3*exp(-im*b)*(3β+3im-im*β^2)/(β^2*(b^2-2)-b^2+im*(β*b^2+2b*β^2));
# Using the two equations described earlier we can compute the expressions containing ``B_1`` as follows
h1_ka  = exp(im*b)/(b^2)*(-b - im);           # Spherical hankel function
dh1_ka = exp(im*b)*(2b + im*(2 - b^2))/(b^3); # Derivative of spherical hankel function
B1h1   = -(v₀/(3im*k) - A₁*dh1_ka)/2;         # Top equation in the system shown earlier
B1dh1  = -(v₀/(3im)*a - A₁*h1_ka);            # Bottom equation in the system shown earlier
# Using the above the analytical expression can be evaluted on the surface
θ_analytical  = collect(0:0.01:π)
p_analytical  = -3.0*ρ₀*c*k*A₁*(h1_ka).*(cos.(θ_analytical));
vn_analytical = -6im*k*B1h1*cos.(θ_analytical);  # Analytical normal velocity
vt_analytical = -3im/a*B1dh1*sin.(θ_analytical); # Analytical tangential velocity


# # Loading the mesh
data_directory = "/work3/mpasc/fmm_paper"
if Sys.isapple()
    data_directory = "/Users/mpasc/OneDrive - Danmarks Tekniske Universitet/paper_fmm"
elseif Sys.islinux()
    data_directory = "/home/mpasc/OneDrive/paper_fmm"
end
mesh_files = ["sphere_1m_fine"]
# mesh_files = ["sphere_1m_fine",   "sphere_1m_finer", "sphere_1m_extremely_fine",
#               "sphere_1m_finest", "sphere_1m_35k",   "sphere_1m_77k"]
            #   "sphere_136k", "sphere_195k", "sphere_305k", "sphere_377k", "sphere_478k"]
# # Iterative Solution through BEM
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
info_string = "depth$(depth)_freq$(Int(frequency))_fmm$(Int(fmm_on))_hmatrix$(Int(hmatrix_on))"
for mesh_file in mesh_files
    @info mesh_file
    mesh = load3dTriangularComsolMesh(joinpath(mesh_path,mesh_file))
    data_file = joinpath(data_directory, mesh_file, "results_$(info_string).jld2")
    if !isfile(data_file)
        LGO = LossyGlobalOuter(mesh,frequency;fmm_on=fmm_on,hmatrix_on=hmatrix_on,depth=depth,n=3,progress=true);
        n   = size(LGO, 2);          # Number of degrees of freedom of the 1x1 system
        v0  = [zeros(2n); v₀*ones(n) ];
        rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*v0 - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*v0),verbose=true);

        time_pa  = @elapsed pa, hist_pa = gmres(LGO,rhs;verbose=true,log=true);
        time_ph  = @elapsed ph  = -LGO.tau_a/LGO.tau_h * pa
        time_dpa = @elapsed dpa, hist_dpa = gmres(LGO.Ga,-(LGO.Ha*pa);log=true) # <- This is the bottleneck...
        time_dph = @elapsed dph, hist_dph = gmres(LGO.Gh,-(LGO.Hh*ph);log=true)
        time_v   = @elapsed v   = v0 - (LGO.phi_a*LGO.Dc*pa +
                                     LGO.phi_a*LGO.Nd*dpa   +
                                     LGO.phi_h*LGO.Dc*ph    +
                                     LGO.phi_h*LGO.Nd*dph)
        time_dvn = @elapsed dvn, hist_dvn = gmres(LGO.Gv,-(LGO.Hv*v);log=true)

        # Normal compontent of the viscous flow
        v_n0 = LGO.Nd'*v;
        # Computing the tangential velocity by substracting the normal information
        v_t = v + LGO.Nd*v_n0;
        vt_sum = sqrt.(v_t[0n+1:1n].^2 + v_t[1n+1:2n].^2 + v_t[2n+1:3n].^2);
        # Plotting real part of
        θ = acos.(mesh.coordinates[3,:]/a)
        scatter(θ,abs.(pa),label="BEM",markersize=1)
        plot!(θ_analytical,abs.(p_analytical),label="Analytical",linewidth=2);
        xlabel!("θ (rad)"); ylabel!("|p|"); title!("Frequency = $(frequency)")
        savefig(joinpath(data_directory, mesh_file,"pa_$(info_string).png"))
        # Plotting normal velocity
        scatter(θ,abs.(v_n0),label="BEM",markersize=1)
        plot!(θ_analytical,abs.(vn_analytical),label="Analytical",linewidth=2)
        xlabel!("θ (rad)"); ylabel!("|Vₙ|"); title!("Frequency = $(frequency)")
        savefig(joinpath(data_directory, mesh_file,"vn_$(info_string).png"))
        # Plotting tangential velocity
        scatter(θ,abs.(vt_sum),label="BEM",markersize=1);
        plot!(θ_analytical,abs.(vt_analytical),label="Analytical",linewidth=2);
        xlabel!("θ (rad)"); ylabel!("|Vₜ|"); title!("Frequency = $(frequency)")
        savefig(joinpath(data_directory, mesh_file,"vt_$(info_string).png"))

        # Saving data
        jldsave(data_file,
                frequency = frequency,
                vt_sum = vt_sum,
                v_n0 = v_n0,
                dpa = dpa,
                dvn = dvn,
                dph = dph,
                pa = pa,
                ph = ph,
                v = v,
                n = n,
                time_dpa = time_dpa,
                time_dph = time_dph,
                time_dvn = time_dvn,
                time_pa = time_pa,
                time_ph = time_ph,
                time_v = time_v,
                hist_pa = hist_pa,
                hist_dpa = hist_dpa,
                info = varinfo(r"LGM")
                )
    else
        @warn "Loading data "
        f = jldopen(data_file)
        frequency = f["frequency"]
        vt_sum = f["vt_sum"]
        v_n0 = f["v_n0"]
        pa  = f["pa"]
        dpa = f["dpa"]
        ph  = f["ph"]
        dph = f["dph"]
        v   = f["v"]
        dv  = f["dv"]
        n   = f["n"]
        info = f["info"]
        time_dpa = f["time_dpa"]
        time_dph = f["time_dph"]
        time_dvn = f["time_dvn"]
        time_pa = f["time_pa"]
        time_ph = f["time_ph"]
        time_v = f["time_v"]
        hist_pa = f["hist_pa"]
        hist_dpa = f["hist_dpa"]
        close(f)
    end
end
