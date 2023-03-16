# # Pulsating sphere (3D - Interior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations,IterativeSolvers, Plots, MeshViz, SpecialFunctions
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading Mesh
# Loading and visualizing the triangular (spherical) mesh
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(mesh_file)
# # Setting up constants
frequency = 100.0;      # Frequency                [Hz]
c  = 343.0;             # Speed up sound           [m/s]
ρ₀ = 1.21;              # Mean density             [kg/m^3]
Z₀ = ρ₀*c;              # Characteristic impedance [Rayl]
vr = 1.0;               # Speed in the r-direction [m/s]
a  = 1.0;               # Radius of sphere_1m      [m]
k  = 2π*frequency/c;    # Wavenumber
# # Analytical Solution
# The analytical description of the interior pressure of an z-oscilating sphere is given by
# ```math
#  p_\text{analytical}(r) = \mathrm{i}Z_0v_r\left(\frac{a}{r}\right)\frac{ka\sin(kr)}{ka\cos(ka) - \sin(ka)}
# ```
# where ``a`` is the sphere radius and ``r (\leq a)`` is the distance from origo to the point ``(x,y,z)``.
# We now generate ``N`` points in 3-dimensional space
N = 20;
x = zeros(N);
y = zeros(N);
z = collect(0.1:(0.9-0.1)/(N-1):0.9);
# From this the analytical expression is computed
p_analytical = im*Z₀*vr*(a./z).*(k*a*sin.(k*z))/(k*a*cos(k*a) - sin(k*a));
# # Solution using the BEM
# We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,progress=false);
H  = Diagonal(1.0 .- C) - F; # Adding the integral free term
# We now setup the right-hand side
vs = im*Z₀*k*vr*ones(length(C));
b  = G*vs; # Computing right-hand side
# Using this we can solve for the surface pressure
p_bem = gmres(H,b;verbose=true);
# Similarly we can compute the BEM solution using the Fast Multipole Operators
Gf = FMMGOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2);
Ff = FMMFOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2);
Hf = Diagonal(1.0 .- C) - Ff;
bf = Gf*vs;                        # Computing right-hand side
p_fmm = gmres(Hf,bf;verbose=true); # Solving the linear system
# Finally the same system is solved using the H-matrix operators
Gh = HGOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2);
Fh = HFOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2);
Hh = Diagonal(1.0 .- C) - Fh;
bh = Gh*vs;                      # Computing right-hand side
p_h = gmres(Hh,bh;verbose=true); # Solving the linear system

# # Field evaluations
# In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of points (as columns)
X_fieldpoints = [x'; y'; z'];
# Using the described field-points we can assemble the (dense) matrices for the field point evaluation
Fp,Gp,_ = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
# Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
p_field  = Fp*p_bem + Gp*vs;
pf_field = Fp*p_fmm + Gp*vs;
ph_field = Fp*p_h   + Gp*vs;
# Plotting the field point pressure it can be seen that all 3 methods find the correct pressures
plot(z,abs.(p_analytical./Z₀),label="Analytical")
scatter!(z,abs.(p_field./Z₀),  label="Full", markershape=:rect)
scatter!(z,abs.(pf_field./Z₀), label="FMM")
scatter!(z,abs.(ph_field./Z₀), label="H-Matrix", markershape=:xcross,color=:black)
ylabel!("p/Z₀"); xlabel!("r/a")

# The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.
simple_mesh = create_simple_mesh(mesh)                          # Creating simple mesh
p_corners = p_bem[sort!(unique(mesh.physics_topology[1:3,:]))]  # Removing quadratic parts
fig, ax, hm = viz(simple_mesh;showfacets=true, color=abs.(p_corners/Z₀))
wgl.Colorbar(fig[1,2],label="|p|"); fig
