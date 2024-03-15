# # Cube with vibrating sides (Interior)
# # Importing related packages
import BoundaryIntegralEquations: arnoldi_basis
using BoundaryIntegralEquations # For BIEs
using IterativeSolvers          # For gmres
using LinearAlgebra             # For Diagonal
using ProgressMeter
using PGFPlots
using Meshes                    # For 3d mesh plots
using Plots                     # For 2d plots
using JLD2
# # Setting up constants
frequency = 350.0;      # Frequency                [Hz]
c  = 343.0;             # Speed up sound           [m/s]
ρ₀ = 1.21;              # Mean density             [kg/m^3]
Z₀ = ρ₀*c;              # Characteristic impedance [Rayl]
v₀ = 1.0;               # Speed in the x-direction [m/s]
# k  = 2π*frequency/c;    # Wavenumber
k = 5.079581127681914
# # Loading and visualizing the triangular cube mesh
# First we define the path to the mesh file
# mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_extremely_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_extra_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_coarser");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_coarse");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube");
# Now we read in the mesh
physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=:linear, physics_order=physics_orders[5])
# We now define the mesh entity that contain the boundary condition. In this case it is boundary 0.
bc_ents = [0];
# wgl.save(joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","docs","src","examples","3d_cube_wave_viz.png"),fig) #hide
# ![](3d_cube_wave_viz.png)
# # Analytical Solution
# The analytical description of the interior pressure in unit cube with the side at ``x=0`` be applied a velocity of ``v_{0}``
# ```math
#  p_\text{analytical}(x) = -\frac{\mathrm{i}Z_0v_{0}\cos(k(1 - x))}{\sin(k)}
# ```
# where ``Z_0`` is the characteristic impedance and ``k`` is the wavenumber.
# We now compute the analytical expression is computed at the points
x_ana = collect(0.00:0.01:1)
p_analytical = -im*Z₀*v₀*cos.(k*(1 .- x_ana))/(sin(k));
# ## Solution using the dense BEM
# In order to apply the boundary coniditon we must find the nodes corresponding to surface entity 0
bc_elements = 0 .== mesh.entities;
bc2 = sort(unique(mesh.physics_topology[:,bc_elements]));
# We now define a vector corresponding to the velocities at `bc2`
v_bem = zeros(ComplexF64,size(mesh.sources,2));
v_bem[bc2] .= im*Z₀*k*v₀*mesh.normals[1,bc2];

function cube_krylov_basis(mesh,klist;v_bem=[],eps=1-4,n_gauss=3,verbose=true,progress=true)
    if isempty(v_bem)
        error("No velocities given")
    end
    n_sources = size(mesh.sources,2)
    V  = zeros(ComplexF64,n_sources, 0)         # Preallocating the total Krylov system
    nK = length(klist)
    solutions = zeros(ComplexF64,n_sources,nK)  # list of solutions
    qlist = zeros(Int64,length(klist))
    if progress == true
        prog = Progress(nK, 0.2, "Assembling Krylov vectors:\t", 50)
    end
    for i = 0:nK-1
        k = klist[i+1]

        Gf = FMMGOperator(mesh,k,depth=2,n_gauss=n_gauss);
        Ff = FMMFOperator(mesh,k,depth=2,n_gauss=n_gauss);
        Hf = Ff - 0.5I ;
        rhs = -Gf*(k*v_bem);                     # Computing right-hand side

        pa_fmm,history = gmres(Hf,rhs;verbose=verbose,log=true);
        solutions[:,i+1] = pa_fmm
        V = [V arnoldi_basis(Hf,rhs,history.iters)]
        qlist[i+1] = history.iters
        if progress == true
            next!(prog)
        end
    end
    U,S,_ = svd(V)                              # Extracting projection matrix
    idx = S .> eps  # Only use projection direcitons with singular values larger than `eps`
    return U[:,idx], solutions, qlist
end
function create_chebyshev_basis(mesh,k_secondary,U;n_gauss=3,v_bem=[],progress=true)
    if isempty(v_bem)
        error("No velocities given")
    end
    temp = zeros(eltype(U), size(U))
    A_output = zeros(eltype(U),size(U,2),size(U,2),length(k_secondary))
    b_output = zeros(eltype(U),size(U,2),length(k_secondary))
    if progress == true
        prog = Progress(length(k_secondary), 0.2, "Assembling Chebyshev basis:\t", 50)
    end
    for (index,k) in enumerate(k_secondary)
        Gf = FMMGOperator(mesh,k,depth=2,n_gauss=n_gauss);
        Ff = FMMFOperator(mesh,k,depth=2,n_gauss=n_gauss);
        rhs = -Gf*v_bem;
        # Hf = Ff - 0.5I;
        mul!(temp,Ff,U)
        A_output[:,:,index] = U'*(temp - U/2)
        # b_output[:,index] = k*(U'*rhs)
        b_output[:,index] = (U'*rhs)
        if progress == true
            next!(prog)
        end
    end
    return A_output, b_output
end
function chebyshev_matrix_coefficients(A_reduced)
    M = size(A_reduced,3)
    Cj = similar(A_reduced)
    for j = 0:M-1
        Cj[:,:,j+1] = 2/M*sum(i -> A_reduced[:,:,i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return Cj
end
function chebyshev_vector_coefficients(b_reduced)
    M = size(b_reduced,2)
    bj = similar(b_reduced)
    for j = 0:M-1
        bj[:,j+1] = 2/M*sum(i -> b_reduced[:,i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return bj
end
function chebyshev_eval(x,M)
    output = ones(M)
    output[2] = x
    for i = 3:M
        output[i] = 2x*output[i-1] - output[i-2]
    end
    return output
end
function eval_chebyshev_matrix(Cj,coeffs)
    output = -Cj[:,:,1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*Cj[:,:,index]
    end
    return output
end
function eval_chebyshev_vector(bj,coeffs)
    output = -bj[:,1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*bj[:,index]
    end
    return output
end
function assemble_chebyshev(Cj,bj,k)
    A = eval_chebyshev_matrix(Cj,chebyshev_eval(k,size(Cj,3)))
    b = eval_chebyshev_vector(bj,chebyshev_eval(k,size(bj,2)))
    return A,b
end

data_directory = "/Users/mpasc/OneDrive - Danmarks Tekniske Universitet/JASA"
data_file = data_directory * "/cube/cube.jld"
if isfile(data_file)
    f = jldopen(data_file)
    mesh = f["mesh"]
    v_bem_constant = f["v_bem_constant"]
    X_fieldpoints = f["X_fieldpoints"]
    p_chebyshev = f["p_chebyshev"]
    k_secondary = f["k_secondary"]
    k_priamry = f["k_primary"]
    A_reduced = f["A_reduced"]
    klist = f["klist"]
    qlist = f["qlist"]
    kmin = f["kmin"]
    kmax = f["kmax"]
    sols = f["sols"]
    Cj = f["Cj"]
    U = f["U"]
    M = f["M"]
    L = f["L"]
    c = f["c"]
    ρ₀ = f["rho0"]
    Z₀ = f["z0"]
    v₀ = f["v0"]
    # x = f["x"]
    # N = f["N"]
    # t_basis = f["t_basis"]
    # t_sweep = f["t_sweep"]
    # t_reducing = f["t_reducing"]
    # nsweep = f["nsweep"]
    close(f)
else
    v_bem_constant = zeros(ComplexF64,size(mesh.sources,2));
    v_bem_constant[bc2] .= im*Z₀*v₀*mesh.normals[1,bc2];
    # v_bem_constant = v_bem/k
    # First we define a set of primary wavenumbers for which we compute the reduced basis (`U`) and the solution (`sols`)
    L = 5;           # Number of primary wavenumbers used to compute the ROM Basis
    k_primary = 2π*(LinRange(300,400,L))/c;  # Defining the primary frequencies
    U,sols,qlist = cube_krylov_basis(mesh,k_primary;v_bem=v_bem_constant,verbose=true);
    @info "Reduced basis size: $(size(U,2)) | Reduction in DOF: $(1 - size(U,2)/size(U,1)) %";
    # Then we define the secondary wavenumbers, i.e. the wavenumbers for which we compute the matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega)\right)\mathbf{U}``
    kmin = 2π*250/c;
    kmax = 2π*450/c;
    g(k)    = 2/(kmax - kmin)*k .- (kmax + kmin )/(kmax - kmin);
    ginv(ω) = (kmax - kmin)/2*ω .+ (kmin + kmax)/2;
    M  = 25;                                # Number of terms in the Chebyshev approximation
    ωᵢ = cos.(π*(collect(0:M-1) .+ 1/2)/M); # Zeros of the Chebyshev polynomials
    k_secondary = ginv.(ωᵢ);                # Mapping [-1,1] to [kmin, kmax]
    #
    t_basis = @elapsed begin
        A_reduced,b_reduced = create_chebyshev_basis(mesh,k_secondary,U;v_bem=v_bem_constant);
    end
    # Now from the basis matrices we can compute the Chebyshev coefficient matrices (``\mathbf{C}_j``) as
    t_reducing = @elapsed begin
        Cj = chebyshev_matrix_coefficients(A_reduced);
        bj = chebyshev_vector_coefficients(b_reduced);
    end

    klist = collect(kmin:0.05:kmax);

    N = 20;
    x = collect(0.1:(0.9-0.1)/(N-1):0.9);
    y = 0.5ones(N);
    z = 0.5ones(N);
    X_fieldpoints = [x'; y'; z'];
    p_chebyshev = zeros(ComplexF64,N,length(klist));
    prog = Progress(length(klist), 0.2, "Frequency sweep:\t", 50)
    t_sweep = @elapsed begin
    for (i,k) in enumerate(klist)
        k_scaled = g(k);                       # Scaling the wavenumber to the interval [-1,1]
        A,b = assemble_chebyshev(Cj,bj,k_scaled)
        p_romi  = U*(A\(k*b));
        # Using the described field-points we can assemble the (dense) matrices for the field point evaluation
        Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=5,m=5,progress=false);
        # Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
        pf_field = Fp*p_romi + Gp*(k*v_bem_constant);
        p_chebyshev[:,i] = pf_field;
        next!(prog)
    end
    end

    jldsave(data_file,mesh=mesh,p_chebyshev=p_chebyshev,klist=klist,
                k_secondary=k_secondary,c=c,k_primary=k_primary,
                U=U,sols=sols,A_reduced=A_reduced,X_fieldpoints=X_fieldpoints,
                M=M,L=L,kmin=kmin,kmax=kmax,Cj=Cj,bj=bj,qlist=qlist,
                v_bem_constant=v_bem_constant, rho0= ρ₀,z0=Z₀,v0=v₀, N=N,x=x,
                t_basis=t_basis, t_reducing=t_reducing, t_sweep=t_sweep,nsweep=length(klist))
end

x = X_fieldpoints[1,:]
N = length(x)

for i in eachindex(klist)
    x_ana = collect(0.00:0.01:1)
    p_analytical = -im*Z₀*v₀*cos.(klist[i]*(1 .- x_ana))/(sin(klist[i]));
    Plots.plot(x_ana,abs.(p_analytical./Z₀), label="Analytical")
    # scatter!(x,imag.(pf_field./Z₀), label="FMM")
    Plots.scatter!(x,abs.(p_chebyshev[:,i]./Z₀), label="ROSEBEM")
    ylabel!("|p/Z₀|"); xlabel!("r/a"); title!("k = $(round(klist[i],digits=2))")
    Plots.savefig(data_directory * "/cube/pressure_klist/pressure_xpos$(i)" )
end

for x_pos in 1:size(X_fieldpoints,2)
    k_ana = kmin:0.01:kmax
    p_analytical = -im*Z₀*v₀*cos.(k_ana*(1 .- x[x_pos]))./(sin.(k_ana));
    Plots.plot(k_ana,abs.(p_analytical./Z₀), label="Analytical", yaxis=:log, linewdith=2)
    Plots.scatter!(klist,abs.(p_chebyshev[x_pos,:]./Z₀), label="ROSEBEM")
    ylabel!("|p/Z₀|"); xlabel!("k"); title!("x = $(round(x[x_pos],digits=2))")
    Plots.savefig(data_directory * "/cube/pressure_x_position/pressure_xpos$(x_pos)" )
end



x_pos = 3
k_ana = kmin:0.01:kmax
p_analytical = -im*Z₀*v₀*cos.(k_ana*(1 .- x[x_pos]))./(sin.(k_ana));
Plots.plot(k_ana*c/(2π),abs.(p_analytical./Z₀), label="Analytical", yaxis=:log, linewdith=2)
Plots.scatter!(klist*c/(2π),abs.(p_chebyshev[x_pos,:]./Z₀), label="ROM", linestyle=:dash,linewidth=2)
ylabel!("|p/Z₀|"); xlabel!("f"); title!("x = $(round(x[x_pos],digits=2))")
Plots.scatter!(k_secondary*c/(2π),abs.(p_secondary ./ Z₀), label="Nodes")
Plots.plot!(k_ana*c/(2π), abs.(p_scalar./Z₀), label="Interpolation")

function chebyshev_scalar_coefficients(b_reduced)
    M = length(b_reduced)
    bj = similar(b_reduced)
    for j = 0:M-1
        bj[j+1] = 2/M*sum(i -> b_reduced[i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return bj
end
function eval_chebyshev_scalar(bj,coeffs)
    output = -bj[1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*bj[index]
    end
    return output
end

g(k)    = 2/(kmax - kmin)*k .- (kmax + kmin )/(kmax - kmin);
ginv(ω) = (kmax - kmin)/2*ω .+ (kmin + kmax)/2;
M  = 25;                                # Number of terms in the Chebyshev approximation
ωᵢ = cos.(π*(collect(0:M-1) .+ 1/2)/M); # Zeros of the Chebyshev polynomials
k_secondary = ginv.(ωᵢ);                # Mapping [-1,1] to [kmin, kmax]
p_secondary = -im*Z₀*v₀*cos.(k_secondary*(1 .- x[x_pos]))./(sin.(k_secondary));
b_secondary = chebyshev_scalar_coefficients(p_secondary)
p_scalar = zeros(ComplexF64,length(k_ana))
for (i,k) in enumerate(k_ana)
    k_scaled = g(k);
    p_scalar[i] = eval_chebyshev_scalar(b_secondary,chebyshev_eval(k_scaled,M))
end




x_pos = 5
k_ana = kmin:0.001:kmax
# k_ana = klist
p_analytical = -im*Z₀*v₀*cos.(k_ana*(1 .- x[x_pos]))./(sin.(k_ana));
plot(k_ana*c/(2π),20*log10.(abs.(p_analytical)/(2e-5)), label="Analytical", linewdith=2)
scatter!(klist*c/(2π),20*log10.(abs.(p_chebyshev[x_pos,:])/(2e-5)), label="ROM", linestyle=:dash,linewidth=2)
ylabel!("|p/Z₀|"); xlabel!("f"); title!("x = $(round(x[x_pos],digits=2))")


x_pos = 5
k_ana = kmin:0.01:kmax
p_analytical = -im*Z₀*v₀*cos.(k_ana*(1 .- x[x_pos]))./(sin.(k_ana));
plot(k_ana*c/(2π),abs.(p_analytical./Z₀), label="Analytical", yaxis=:log, linewdith=2)
scatter!(klist*c/(2π),abs.(p_chebyshev[x_pos,:]./Z₀), label="ROM", linestyle=:dash,linewidth=2)
ylabel!("|p/Z₀|"); xlabel!("f"); title!("x = $(round(x[x_pos],digits=2))")

using PGFPlots
x_pos = 5
k_ana = kmin:0.01:kmax
p_analytical = -im*Z₀*v₀*cos.(k_ana*(1 .- x[x_pos]))./(sin.(k_ana));
plt = Axis(PGFPlots.Plots.Linear(k_ana*c/(2π),abs.(p_analytical./Z₀),legendentry="Analytical", mark="none",style="thick",color="black"),
            xlabel="Frequency [Hz]",ylabel="Pressure [Pa]",width="12cm", height="8cm",ymin=0.0,ymode="log")
push!(plt,PGFPlots.Plots.Linear(klist*c/(2π),abs.(p_chebyshev[x_pos,:]./Z₀),legendentry="ROSEBEM",mark="*",markSize=2,onlyMarks=true))
plt

print(tikzCode(plt))



using WriteVTK
bc_ents = sort(unique(mesh.topology[:,bc_elements]));
not_bc_ents = sort(unique(mesh.topology[:,.!bc_elements]));

dataFileVTK = "/Users/mpasc/Dropbox/Apps/ShareLaTeX/preprint_scalable_rosebem/meshes/cube_bc_ents"
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, con) for con in eachcol(mesh.topology[:,bc_elements])]
vtkfile = vtk_grid(dataFileVTK, mesh.coordinates[1,:],
                                mesh.coordinates[2,:],
                                mesh.coordinates[3,:], cells)
vtk_save(vtkfile)

point_file = "/Users/mpasc/Dropbox/Apps/ShareLaTeX/preprint_scalable_rosebem/cube"
vtk_grid(point_file, X_fieldpoints[1,:], X_fieldpoints[2,:], X_fieldpoints[3,:]) do vtk
    vtk["my_point_data"] = 1.0
end




data_directory = "/Users/mpasc/OneDrive - Danmarks Tekniske Universitet/JASA"
data_file = data_directory * "/cube/cube_full.jld"
if isfile(data_file)
    f = jldopen(data_file)
    mesh = f["mesh"]
    v_bem_constant = f["v_bem_constant"]
    X_fieldpoints = f["X_fieldpoints"]
    p_chebyshev = f["p_chebyshev"]
    k_secondary = f["k_secondary"]
    k_priamry = f["k_primary"]
    A_reduced = f["A_reduced"]
    klist = f["klist"]
    qlist = f["qlist"]
    kmin = f["kmin"]
    kmax = f["kmax"]
    sols = f["sols"]
    Cj = f["Cj"]
    U = f["U"]
    M = f["M"]
    N = f["N"]
    L = f["L"]
    c = f["c"]
    x = f["x"]
    ρ₀ = f["rho0"]
    Z₀ = f["z0"]
    v₀ = f["v0"]
    t_sweep_full = f["t_sweep"]
    nsweep = f["nsweep"]
    close(f)
else
    v_bem_constant = zeros(ComplexF64,size(mesh.sources,2));
    v_bem_constant[bc2] .= im*Z₀*v₀*mesh.normals[1,bc2];
    # Then we define the secondary wavenumbers, i.e. the wavenumbers for which we compute the matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega)\right)\mathbf{U}``
    kmin = 2π*250/c;
    kmax = 2π*450/c;

    N = 20;
    x = collect(0.1:(0.9-0.1)/(N-1):0.9);
    y = 0.5ones(N);
    z = 0.5ones(N);
    X_fieldpoints = [x'; y'; z'];
    p_full = zeros(ComplexF64,N,length(klist));
    prog = Progress(length(klist), 0.2, "Frequency sweep:\t", 50)
    t_sweep = @elapsed begin
    for (i,k) in enumerate(klist)
        k = klist[i+1]
        Gf = FMMGOperator(mesh,k,depth=2,n_gauss=n_gauss);
        Ff = FMMFOperator(mesh,k,depth=2,n_gauss=n_gauss);
        Hf = Ff - 0.5I ;
        rhs = -Gf*(k*v_bem);                     # Computing right-hand side
        pa_fmm,_ = gmres(Hf,rhs;verbose=false);

        # Using the described field-points we can assemble the (dense) matrices for the field point evaluation
        Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=5,m=5,progress=false);
        # Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
        pf_field = Fp*pa_fmm + Gp*(k*v_bem_constant);
        p_full[:,i] = pf_field;
        next!(prog)
    end
    end

    jldsave(data_file,mesh=mesh,p_full=p_full,klist=klist,t_sweep=t_sweep,nsweep=length(klist))
end


v_bem_constant = zeros(ComplexF64,size(mesh.sources,2));
v_bem_constant[bc2] .= im*Z₀*v₀*mesh.normals[1,bc2];
# Then we define the secondary wavenumbers, i.e. the wavenumbers for which we compute the matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega)\right)\mathbf{U}``
kmin = 2π*250/c;
kmax = 2π*450/c;

n_gauss = 3
N = 20;
x = collect(0.1:(0.9-0.1)/(N-1):0.9);
y = 0.5ones(N);
z = 0.5ones(N);
X_fieldpoints = [x'; y'; z'];
p_full = zeros(ComplexF64,N,length(klist));
# prog = Progress(length(klist), 0.2, "Frequency sweep:\t", 50)
k = kmin
Gf = FMMGOperator(mesh,k,depth=2,n_gauss=n_gauss);
Ff = FMMFOperator(mesh,k,depth=2,n_gauss=n_gauss);
Hf = Ff - 0.5I ;
rhs = -Gf*(k*v_bem);                     # Computing right-hand side
pa_fmm = gmres(Hf,rhs;verbose=true);

# Using the described field-points we can assemble the (dense) matrices for the field point evaluation
Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=5,m=5,progress=false);
# Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
pf_field = Fp*pa_fmm + Gp*(k*v_bem_constant);
p_full[:,i] = pf_field;
    # next!(prog)
# end
# end


F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=5,m=5,progress=true)

(G'*x) ./ (Gf'*x)

(G*x) ./ (Gf*x)

H = F - 0.5I

pa_fmm = gmres(Gf,rhs;verbose=true);
pa = gmres(G,rhs;verbose=true);
