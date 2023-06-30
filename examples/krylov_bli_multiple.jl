using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--fmin"
            help = "minimum primary frequency"
            required = true
            arg_type = Float64
        "--freq"
            help = "expansion frequency"
            required = true
            arg_type = Float64
        "--fmax"
            help = "maximum primary frequency"
            required = true
            arg_type = Float64
        "--fminsweep"
            help = "minimum response frequency"
            arg_type = Float64
            default = 200.0
        "--fmaxsweep"
            help = "maximum  response frequency"
            default = 600.0
            arg_type = Float64
        "--M"
            help = "number of Taylor terms"
            arg_type = Int64
            default = 20
        "--L"
            help = "number of primary frequencies"
            arg_type = Int64
            default = 1
        "--q"
            help = "number of Krylov vectors per primary frequency"
            required = true
            arg_type = Int64
        "--compute_full_solution"
            help = "computing full solution"
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
fmin = intputArguments["fmin"]
fmax = intputArguments["fmax"]
freq = intputArguments["freq"]
fminsweep = intputArguments["fminsweep"]
fmaxsweep = intputArguments["fmaxsweep"]
M = intputArguments["M"]
L = intputArguments["L"]
q = intputArguments["q"]
compute_full_solution = intputArguments["compute_full_solution"]
#==========================================================================================
                                Loading packages
==========================================================================================#
using LinearAlgebra,JLD2,OpenBEM,ProgressMeter
import OpenBEM: amb2prop,VTconst,krylovBLI!,newFullSingularAssembleParallel!,assembleBLI!
import OpenBEM: arnoldiBLI, krylovBasisBLI,reducedOrderModel,applyROMTaylorExpansion
import OpenBEM: compute_k1_k2
# println(versioninfo())
#==========================================================================================
                            Helper function
==========================================================================================#
masterFrequencies(fmin,fmax,Mmaster)   = collect(range(fmin,fmax,length=Mmaster))
masterWavenumbers(fmin,fmax,Mmaster,c) = 2.0*π*masterFrequencies(fmin,fmax,Mmaster)/c
absorptionCoefficient(p1,p2,k,d) = 1.0 - abs2((exp(-im*k*d) - p2/p1)/(p2/p1 - exp(im*k*d)))
#==========================================================================================
           Setting up the problem (mesh, constants and boundary conditions)
                    Note: We need a quadratic mesh for this to work
==========================================================================================#
## Defining relevant constants
rho,c,cf,gamma,nu,alfa,Cp,Cv,lambda,beta = amb2prop(freq)
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = VTconst(;freq=freq,S=-1)
k  = 2.0*π*freq/c                               # Wavenumber
u0 = 1.0/(rho*c)                                # Initial Velocity
ω  = 2.0*π*freq                                 # Angular frequency
v0 = -im*ω*rho*u0                               # Initial velocity
dV = sqrt(2.0*μ/(rho*ω))                        # Viscous BL-thickness
dT = sqrt(2*lambda/(ω*rho*Cp))                  # Thermal BL-thickness
k1 = -dV*(im - 1.0)/2.0                         # Constant 1
k2 = +dT*k^2*(im - 1.0)/2.0*(gamma - 1.0)       # Constant 2
#==========================================================================================
                                    COMSOl Test Mesh
==========================================================================================#
data_directory = "/work3/mpasc/"
boundaryConditionEntityRes = [3] .- 1
mesh_file = "multipleResonatorsFine"
# mesh_file = "resonatorFinal"
meshRes,entsRes = load3dTriangularComsolMesh("OpenBEM.jl/examples/meshes/comsol/" * mesh_file,
                                        entites=true,physicsType="geometry")
VelBC = Bool.(sum(entsRes .∈ boundaryConditionEntityRes',dims=2))[:]
BLIBC = .!VelBC
VelID = sort(unique(meshRes.physicsTopology[:,VelBC]))
BLIID = sort(unique(meshRes.physicsTopology[:,BLIBC]))
#==========================================================================================
                            Series-expansion of the Kernels
==========================================================================================#
Nk = 401
k0 = k   # Expansion frequency
file_name = "ROMBEM_M$(M)_q$(q)_L$(L)_freq$(Int(freq))_fmin$(Int(fmin))_fmax$(Int(fmax)).jld2"
data_file = data_directory * mesh_file * "/" * file_name
if isfile(data_file)
    f = jldopen(data_file)
    meshRes = f["mesh"]
    entsRes = f["ents"]
    Aout  = f["Abasis"]
    Hout  = f["Hbasis"]
    Tout  = f["Tbasis"]
    bout  = f["bbasis"]
    freq  = f["freq"]
    VelBC = f["VelBC"]
    BLIBC = f["BLIBC"]
    klist = f["klist"]
    k0 = f["k0"]
    q  = f["q"]
    M  = f["M"]
    V  = f["V"]
    # Frequency Sweep
    Nk = f["Nk"]
    x1 = f["x1"]
    x2 = f["x2"]
    d  = f["d"]
    sources = f["sources"]
    kSweep  = f["kSweep"]
    timeROMSweep    = f["timeROMSweep"]
    absorptionCoeffients = f["absorptionCoeffients"]
    close(f)
    VelID = sort(unique(meshRes.physicsTopology[:,VelBC]))
    BLIID = sort(unique(meshRes.physicsTopology[:,BLIBC]))
else
    if L == 1
        klist = [k0]
    else
        klist = masterWavenumbers(fmin,fmax,L,c) # Primary frequency range
    end
    # Computing Krylov basis
    timeKrylovAssembly = @elapsed V = krylovBasisBLI(meshRes,klist,q,v0;n=3,m=10,
                                                eps=1-4,VelBC=VelBC,BLIBC=BLIBC,c=c)
    # Computing ROM (Compressing rows)
    timeTaylorAssembly = @elapsed Aout,Tout,Hout,bout = krylovBLI!(
                                    meshRes,k0,meshRes.sources;V=V,
                                    singCheck=true,n=3,m=10,M=M,velocityBC=VelBC,BLIBC=BLIBC)
    # Computing Frequency Sweep
    x1 = [0.0; 0.2; 0.0]
    x2 = [0.0; 0.3; 0.0]
    d  = norm(x1-x2)
    sources = [x1 x2]
    kSweep  = masterWavenumbers(fminsweep,fmaxsweep,Nk,c)
    absorptionCoeffients = similar(kSweep)
    pointPressure = zeros(ComplexF64,length(kSweep),2)
    timeROMSweep = @elapsed begin
        prog = Progress(length(kSweep), 0.1, "Frequency Sweep - ROM: \t\t", 25)
        for i = 1:length(kSweep)
            ks = kSweep[i]
            As,Hs,Ts,bs = applyROMTaylorExpansion(Aout,Hout,Tout,bout,ks,k0)

            ks1,ks2 = compute_k1_k2(ks*c/(2π))
            ps = V*((As + ks1*Ts + ks2*Hs)\(-v0*bs))

            Fa,Ga,Ca = newFullSingularAssembleParallel!(meshRes,ks,sources;
                                                        n=3,m=20,velocityBC=VelBC)

            ba = -sum(Ga[:,VelID,1],dims=2)[:]*v0
            P  = (Fa[:,:,1])*ps - ba

            absorptionCoeffients[i] = absorptionCoefficient(P[1],P[2],ks,d)
            pointPressure[i,1] = P[1]
            pointPressure[i,2] = P[2]
            next!(prog)
        end
    end
    # Saving timing
    kst = kSweep[1]
    timeROMAssembly = @elapsed begin
        Ast,Hst,Tst,bst = applyROMTaylorExpansion(Aout,Hout,Tout,bout,kst,k0)
    end
    timeROMSolve = @elapsed begin
        ks1t,ks2t = compute_k1_k2(kst*c/(2π))
        pst = V*((Ast + ks1t*Tst + ks2t*Hst)\(-v0*bst))
    end
    # Saving results
    jldsave(data_file,Abasis=Aout,Hbasis=Hout,Tbasis=Tout,bbasis=bout,
                    V=V,klist=klist,q=q,M=M,freq=freq,k0=k0,
                    mesh=meshRes,ents=entsRes,VelBC=VelBC,BLIBC=BLIBC,
                    timeKrylovAssembly=timeKrylovAssembly,
                    timeTaylorAssembly=timeTaylorAssembly,
                    timeROMAssembly=timeROMAssembly,
                    timeROMSolve=timeROMSolve,
                    timeROMSweep=timeROMSweep,x1=x1,x2=x2,d=d,sources=sources,Nk=Nk,
                    kSweep=kSweep,absorptionCoeffients=absorptionCoeffients,
                    pointPressure=pointPressure)
end
#==========================================================================================
                    Computing frequency sweep using full model
==========================================================================================#
N = length(fminsweep:fmaxsweep)
if compute_full_solution
results_name = "/Results_N$(N)_fminsweep$(Int(fminsweep))_fmaxsweep$(Int(fmaxsweep)).jld2"
results_file = data_directory * mesh_file * results_name
if isfile(results_file)
    results = jldopen(results_file)
    kBLI    = results["kBLI"]
    ACBLI   = results["ACBLI"]
    sources = results["sources"]
    timeFullSweep    = results["timeFullSweep"]
    surfaceSolution  = results["surfaceSolution"]
    close(results)
else
    kBLI  = masterWavenumbers(fminsweep,fmaxsweep,N,c)
    ACBLI = similar(kBLI)
    surfaceSolution = zeros(ComplexF64,size(meshRes.normals,2),length(kBLI))
    timeFullSweep = @elapsed begin
        prog  = Progress(length(kBLI), 0.1, "Frequency Sweep - Results: \t\t", 25)
        for i = 1:length(kBLI)
            ks = kBLI[i]
            ks1,ks2 = compute_k1_k2(ks*c/(2π))

            Als,bs = krylovBLI!(meshRes,ks,meshRes.sources,Diagonal(ones(size(meshRes.sources,2))),ks1,ks2;M=1,v0=v0,
                    n=3,m=10,velocityBC=VelBC,BLIBC=BLIBC,singCheck=true,progress=false)
            ps = Als\bs               # BEM
            # Saving Surface Solution
            surfaceSolution[:,i] = ps
            # Computing field point solution
            Fa,Ga,Ca = newFullSingularAssembleParallel!(meshRes,ks,sources;
                                                    n=3,m=20,velocityBC=VelBC)
            ba = -sum(Ga[:,VelID,1],dims=2)[:]*v0
            P  = (Fa[:,:,1])*ps - ba

            ACBLI[i] = absorptionCoefficient(P[1],P[2],ks,d)
            if isnan(ACBLI[i])
                error("Absorption Coefficient computed to be a NaN")
                break
            end
            next!(prog)
        end
    end
    ks = kBLI[1]
    ks1,ks2 = compute_k1_k2(ks*c/(2π))

    timeFullAssembly = @elapsed begin
        Als,bs = krylovBLI!(meshRes,ks,meshRes.sources,Diagonal(ones(size(meshRes.sources,2))),ks1,ks2;M=1,v0=v0,
                n=3,m=10,velocityBC=VelBC,BLIBC=BLIBC,singCheck=true,progress=false)
    end
    timeFullSolve = @elapsed begin
        ps = Als\bs               # BEM
    end
    jldsave(results_file,
            kBLI=kBLI,ACBLI=ACBLI,surfaceSolution=surfaceSolution,
            x1=x1,x2=x2,d=d,sources=sources,timeFullSweep=timeFullSweep,
            timeFullAssembly=timeFullAssembly,
            timeFullSolve=timeFullSolve,
            mesh=meshRes,ents=entsRes,VelBC=VelBC,BLIBC=BLIBC)
end
end
#==========================================================================================
                    Computing frequency sweep using ROM model
==========================================================================================#
sweep_name = "/Sweep_q$(q)_N$(N)_f0$(Int(freq))_fmin$(Int(fmin))_fmax$(Int(fmax)).jld2"
sweep_file = data_directory * mesh_file * sweep_name
if isfile(sweep_file)
    sweep = jldopen(sweep_file)
    AC = sweep["AC"]
    romSurfaceSolution = sweep["romSurfaceSolution"]
    close(sweep)
else
    romSurfaceSolution  = zeros(ComplexF64,size(V,1),length(kSweep),size(Aout,3))
    AC = zeros(length(kSweep),size(Aout,3))
    prog = Progress(size(Aout,3)*length(kSweep), 0.1, "Frequency Sweep - Taylor: \t", 25)
    for i = 1:length(kSweep)
        for j = 1:size(Aout,3)
            ks = kSweep[i]
            ks1,ks2 = compute_k1_k2(ks*c/(2π))

            Ao = @view Aout[:,:,1:j]
            Ho = @view Hout[:,:,1:j]
            To = @view Tout[:,:,1:j]
            bo = @view bout[:,1:j]
            As,Hs,Ts,bs = applyROMTaylorExpansion(Ao,Ho,To,bo,ks,k0)
            ps = V*((As + ks1*Ts + ks2*Hs)\(-v0*bs))
            # Saving surface solution
            romSurfaceSolution[:,i,j] = ps
            # Computing field point solution
            Fa,Ga,Ca = newFullSingularAssembleParallel!(meshRes,ks,sources;
                                                        n=3,m=20,velocityBC=VelBC)
            ba = -sum(Ga[:,VelID,1],dims=2)[:]*v0
            P = (Fa[:,:,1])*ps - ba
            # Evaluating absorption coefficient
            AC[i,j] = absorptionCoefficient(P[1],P[2],ks,d)
            next!(prog)
        end
    end
    jldsave(sweep_file,AC=AC,romSurfaceSolution=romSurfaceSolution)
end
