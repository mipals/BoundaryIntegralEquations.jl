# # Cube with vibrating sides (3D - Interior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations,IterativeSolvers, Plots, MeshViz, SpecialFunctions
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading and visualizing the triangular cube mesh
# First we define the path to the mesh file
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube");
# Now we read in the mesh
mesh = load3dTriangularComsolMesh(mesh_file)
# We now define the mesh entity that contain the boundary condition. In this case it is boundary 0.
bc_ents = [0];
# Now two simple meshes is created. One for the boundary condition and one for the remaining part of the mesh
simple_bc   = create_bc_simple_mesh(mesh,bc_ents);
simple_mesh = create_bc_simple_mesh(mesh,bc_ents,false);
# Using the simple meshes we can visualize the mesh, with boundary 0 (where velocity condition will be applied) shown in red
viz(simple_mesh;showfacets=true)
viz!(simple_bc;showfacets=true,color=:red)
wgl.current_figure()
# # Setting up constants
frequency = 100.0;      # Frequency                [Hz]
c  = 343.0;             # Speed up sound           [m/s]
ρ₀ = 1.21;              # Mean density             [kg/m^3]
Z₀ = ρ₀*c;              # Characteristic impedance [Rayl]
v₀ = 1.0;               # Speed in the x-direction [m/s]
k  = 2π*frequency/c;    # Wavenumber
# # Analytical Solution
# The analytical description of the interior pressure in unit cube with the side at ``x=0`` be applied a velocity of ``v_{0}``
# ```math
#  p_\text{analytical}(x) = -\frac{\mathrm{i}Z_0v_{0}\cos(k(1 - x))}{\sin(k)}
# ```
# where ``Z_0`` is the characteristic impedance and ``k`` is the wavenumber.
# We start by defining a ``N`` linearly spaced points defined by ``(x,0.5,0.5)`` with ``x\in[0.1,0.9]``.
N = 20;
x = collect(0.1:(0.9-0.1)/(N-1):0.9);
y = 0.5ones(N);
z = 0.5ones(N);
# We now compute the analytical expression is computed at the points
p_analytical = -im*Z₀*v₀*cos.(k*(1 .- x))/(sin(k));
# ## Solution using the dense BEM
# We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,progress=false);
H  = Diagonal(C) - F; # Adding the integral free term
# In order to apply the boundary coniditon we must find the nodes corresponding to surface entity 0
bc_elements = 0 .== mesh.entities;
bc2 = sort(unique(mesh.physics_topology[:,bc_elements]));
# We now define a vector corresponding to the velocities at `bc2`
v_bem = zeros(ComplexF64,length(C));
v_bem[bc2] .= im*Z₀*k*v₀*mesh.normals[1,bc2];
# Using this we can define the right-hand side
b = G*v_bem;
# Finally the pressures can be computed by solving a linear system of equations
p_bem = gmres(H,b;verbose=true);
# ## Solution using the FMM-BEM
# Similarly we can compute the BEM solution using the Fast Multipole Operators
Gf = FMMGOperator(mesh,k);
Ff = FMMFOperator(mesh,k);
Hf = Diagonal(C) - Ff;
bf = Gf*v_bem;                     # Computing right-hand side
p_fmm = gmres(Hf,bf;verbose=true); # Solving the linear system

# # Field evaluations
# In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of points (as columns)
X_fieldpoints = [x'; y'; z'];
# Using the described field-points we can assemble the (dense) matrices for the field point evaluation
Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
# Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
p_field  = Fp*p_bem + Gp*v_bem;
pf_field = Fp*p_fmm + Gp*v_bem;
# Plotting the field point pressure it can be seen that all 3 methods find the correct pressures
plot(x,abs.(p_analytical./Z₀), label="Analytical")
scatter!(x,abs.(p_field./Z₀),  label="Full", markershape=:rect)
scatter!(x,abs.(pf_field./Z₀), label="FMM")
ylabel!("|p/Z₀|"); xlabel!("r/a")

# The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.
simple_mesh = create_simple_mesh(mesh)                          # Creating simple mesh
p_corners = p_bem[sort!(unique(mesh.physics_topology[1:3,:]))]  # Removing quadratic parts
fig, ax, hm  = viz(simple_mesh;showfacets=true, color=abs.(p_corners/Z₀))
wgl.Colorbar(fig[1,2],label="|p|"); fig
