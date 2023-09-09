# # Oscilating sphere (Interior)
# # Importing related packages
using BoundaryIntegralEquations # For BIEs
using IterativeSolvers          # For gmres
using SpecialFunctions          # For Bessel functions
using LinearAlgebra             # For Diagonal
using Meshes                    # For 3d mesh plots
using Plots                     # For 2d plots
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Setting up constants
frequency = 100.0;      # Frequency                [Hz]
c  = 343.0;             # Speed up sound           [m/s]
ρ₀ = 1.21;              # Mean density             [kg/m^3]
Z₀ = ρ₀*c;              # Characteristic impedance [Rayl]
vz = 1.0;               # Speed in the z-direction [m/s]
a  = 1.0;               # Radius of sphere_1m      [m]
k  = 2π*frequency/c;    # Wavenumber
# # Loading Mesh
# Loading and visualizing the triangular (spherical) mesh
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarser");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(mesh_file)
# # Analytical Solution
# The analytical description of the interior pressure of an z-oscilating sphere is given by
# ```math
#  p_\text{analytical}(x,y,z) = -3v_z z \mathrm{i}kZ_0\frac{j_1(kr)}{kr(j_0(ka) - 2j_2(ka))},
# ```
# where ``j_i`` is the spherical bessel function of the first kind, ``r`` is the distance from origo to the point ``(x,y,z)`` and ``v_z`` is the oscilating velocity in the ``z``-direction. For the points in the interior where the pressure will evaluated ``z`` will be equal to ``r`` as the remaining two coordinates are set to 0.
# From this the analytical expression is computed
z_ana = collect(0.0:0.01:1)
p_analytical = -3*vz*z_ana*im*k*Z₀.*(sphericalbesselj.(1,k*z_ana))./(k*z_ana*(sphericalbesselj.(0,k*a) - 2*sphericalbesselj.(2,k*a)));
# # Solution using the BEM
# We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=3,m=3,progress=false);
H  = Diagonal(1.0 .- C) - F; # Adding the integral free term
# We now setup the right-hand side
vs = im*Z₀*k*mesh.normals[3,:]*vz;
b = G*vs; # Computing right-hand side
# Using this we can solve for the surface pressure
p_bem = gmres(H,b;verbose=true);
# Similarly we can compute the BEM solution using the Fast Multipole Operators
Gf = FMMGOperator(mesh,k);
Ff = FMMFOperator(mesh,k);
Hf = Diagonal(1.0 .- C) - Ff;
bf = Gf*vs;                        # Computing right-hand side
p_fmm = gmres(Hf,bf;verbose=true); # Solving the linear system

# # Field evaluations
# In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of ``N`` points (as columns) in 3-dimensional space
N = 20;
x = zeros(N);
y = zeros(N);
z = collect(0.1:(0.9-0.1)/(N-1):0.9);
X_fieldpoints = [x'; y'; z'];
# Using the described field-points we can assemble the (dense) matrices for the field point evaluation
Fp,Gp,_ = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
# Using the matrices we can evalute the pressure at the field points using the previously found surface pressure
p_field  = Fp*p_bem + Gp*vs;
pf_field = Fp*p_fmm + Gp*vs;
# Plotting the field point pressure it can be seen that all 3 methods find the correct pressures
plot(z_ana,abs.(p_analytical./Z₀),  label="Analytical")
scatter!(z,abs.(p_field./Z₀),       label="Full", markershape=:rect)
scatter!(z,abs.(pf_field./Z₀),      label="FMM")
ylabel!("p/Z₀"); xlabel!("r/a")

# The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=abs.(data_viz/Z₀))
wgl.Colorbar(fig[1,2],label="|p/Z₀|"); fig
