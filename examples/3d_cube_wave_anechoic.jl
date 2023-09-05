# # Cube with anechoic condition (Interior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations, IterativeSolvers, Plots, Meshes, SpecialFunctions
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading Mesh
# Loading and visualizing the triangular cube mesh
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_extremely_coarse");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_extra_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_coarser");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube");
physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=:linear, physics_order=physics_orders[5])
bc_pre = [1] .- 1; # The .-1 is due to COMSOL 0-indexing of exported entities
bc_ane = [6] .- 1; # The .-1 is due to COMSOL 0-indexing of exported entities
# Creating simple meshes for full mesh and boundary condition
simple_mesh = create_bc_simple_mesh(mesh,[bc_pre; bc_ane],false);
simple_pre  = create_bc_simple_mesh(mesh,bc_pre);
simple_ana  = create_bc_simple_mesh(mesh,bc_ane);
# We now plot the mesh with the pressure condition shown in red and the anechoic condition shown in blue
viz(simple_pre;showfacets=true,color=:red)
viz!(simple_mesh;showfacets=true,alpha=0.1)
viz!(simple_ana;showfacets=true,color=:blue)
wgl.current_figure()
# # Setting up constants
frequency = 54.59;                              # Frequency                [Hz]
c  = 343;                                       # Speed up sound           [m/s]
ρ₀ = 1.21;                                      # Mean density             [kg/m^3]
Z₀ = ρ₀*c;                                      # Characteristic impedance [Rayl]
P₀ = 1.0;                                       # Pressure of plane wave   [m/s]
l  = 1.0;                                       # Length of cube           [m]
k  = 2π*frequency/c;                            # Computing the wavenumber
# # Analytical Solution
# The analytical description of the interior pressure of a planewave is given by
# ```math
#  p(x) = P_0\exp(-\mathrm{i}kx),
# ```
# where ``P_0`` is the magnitude of the planewave.
# We now compute the analytical expression is computed at the points
x_ana = collect(0.0:0.01:1)
p_analytical = P₀*exp.(-im*k*x_ana);
# ## Solution using the dense BEM
# We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=4,m=4,progress=false);
H  = Diagonal(C) - F; # Adding the integral free term
# Now we find the correct indicies
bc_pre = 0 .== mesh.entities;
bc_ane = 5 .== mesh.entities;
bc1 = sort(unique(mesh.physics_topology[:,bc_pre]));
bc2 = sort(unique(mesh.physics_topology[:,bc_ane]));
# First we define the ``p=1`` bc at ``x=0``
ps = zeros(ComplexF64,length(C));
ps[bc1] .= 1.0;
bs = -(H*ps); # Computing right-hand side
# Then we define the anechoic condition at ``x=1``.
Z = zeros(ComplexF64,length(C));
Z[bc1] .=  1.0;  # Unknown velocities
Z[bc2] .= -im*k; # Impedance condition
#
Hmask = ones(ComplexF64,length(C));
Hmask[bc1] .= 0.0; # Known pressures
A = H*Diagonal(Hmask) + G*Diagonal(-Z);
# Using this we can solve for the surface pressure
z_bem = gmres(A,bs);
# Extracting correct surface pressures
p_bem = copy(z_bem);
p_bem[bc1] .= 1.0;
# Extracting correct surface velocities
v_bem = zeros(ComplexF64,length(C));
v_bem[bc1] = z_bem[bc1];
v_bem[bc2] = p_bem[bc2]*(-im*k);

# ## Solution using the FMM-BEM
# First we define the two operators
Gf = FMMGOperator(mesh,k,depth=2);
Ff = FMMFOperator(mesh,k,depth=2);
Hf = 0.5I - Ff;
# Then we compute the right-hand side and the linear system similar to before. Note that the system matrix will be represented by a collection of linear maps
bf = -(Hf*ps);
Af = Hf*Diagonal(Hmask) + Gf*Diagonal(-Z)
# Similary to before we can now solve the linear system using an iterative scheme
z_fmm = gmres(Af,bf);
# Finally we extract the pressure (and assert the pressures at boundary 1)
p_fmm = copy(z_fmm);
p_fmm[bc1] .= 1.0;
# as well as the surface velocities
v_fmm = zeros(ComplexF64,length(C));
v_fmm[bc1] = z_fmm[bc1];
v_fmm[bc2] = p_fmm[bc2]*(-im*k);

# # Field evaluations
# In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of points (as columns). Thi is done by chosing ``N`` linearly spaced points defined by ``(x,0.5,0.5)`` with ``x\in[0.1,0.9]``.
N = 20;
x = collect(0.1:(0.9-0.1)/(N-1):0.9);
y = 0.5ones(N);
z = 0.5ones(N);
X_fieldpoints = [x'; y'; z'];
# Using the described field-points we can assemble the (dense) matrices for the field point evaluation
Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
# Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
p_field  = Fp*p_bem + Gp*v_bem;
pf_field = Fp*p_fmm + Gp*v_fmm;
# Plotting the field point pressure it can be seen that all 3 methods find the correct pressures
plot(x_ana,real.(p_analytical/P₀),  label="Re-Analytical")
scatter!(x,real.(p_field/P₀),   label="Re-Full", markershape=:rect)
scatter!(x,real.(pf_field/P₀),  label="Re-FMM")
plot!(x_ana,imag.(p_analytical/P₀), label="Im-Analytical")
scatter!(x,imag.(p_field/P₀),   label="Im-Full", markershape=:rect)
scatter!(x,imag.(pf_field/P₀),  label="Im-FMM")
xlabel!("r/a")
# The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=real.(data_viz/P₀))
#src wgl.scatter!(mesh.sources[1,:],mesh.sources[2,:],mesh.sources[3,:],color=abs.(p_fmm/Z₀))
wgl.Colorbar(fig[1,2],label="Re(p/p₀)"); fig
