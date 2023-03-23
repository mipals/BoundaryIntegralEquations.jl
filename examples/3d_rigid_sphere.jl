# # Rigid sphere scattering (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations
using IterativeSolvers, Plots, MeshViz, SpecialFunctions, LegendrePolynomials
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading Mesh
# Loading and visualizing the triangular (spherical) mesh
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarser");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(mesh_file;physics_order=:disctrilinear)
# # Setting up constants
frequency = 100.0;                              # Frequency                [Hz]
c  = 340;                                       # Speed up sound           [m/s]
a  = 1.0;                                       # Radius of sphere_1m      [m]
k  = 2π*frequency/c;                            # Wavenumber
P₀ = 1.0;                                       # Magnitude of planewave
# # Analytical Solution
# The analytical solution of the scattering of a sphere by plane wave can be computed as (Ihlenburg1998)
# ```math
#  p_\text{analytical}(r, \theta) = P_0\left(\exp(\mathrm{i}kr\cos(\theta)) - \sum_{n}^\infty \mathrm{i}^n(2n+1)\frac{j_n^{'}(ka)}{h_n^{'}(ka)}P_n(\cos(\theta))h_n(kr)\right),
# ```
# where ``j_n, h_n`` and ``P_n`` are respectively the spherical Bessel function of the first kind, the Hankel function of the first kind and the Legendre polynomial of degree ``n``.
# To make the implementation easier we defin the following helper functions
dsp_j(n,z) = n/z*sphericalbesselj.(n,z) - sphericalbesselj.(n+1,z); # Derivative of j
dsp_y(n,z) = n/z*sphericalbessely.(n,z) - sphericalbessely.(n+1,z); # Derivative of y
sp_h(n,z)  = sphericalbesselj.(n,z) + im*sphericalbessely.(n,z);    # Hankel function (h)
dsp_h(n,z) = dsp_j.(n,z) + im*dsp_y.(n,z);                          # Derivative of h
# Using the helper functions we can define the a function for the coefficients
c_n(n,ka)  = (im)^n*(2n + 1)*(dsp_j.(n,ka)./dsp_h.(n,ka));
# The total pressure as the sum of scattered and incident pressure
θ_analytical = collect(0:0.01:π);   # Angles where the analytical solution will be evaluated
N_truncation = 50;                  # Truncation of the sum
p_s = P₀*sum(n -> -c_n(n,k*a) .* Pl.(cos.(θ_analytical), n) .* sp_h.(n,k*a), 0:N_truncation);
p_i = P₀*exp.(im*k*cos.(θ_analytical))
p_analytical = p_s + p_i;
# # Solution using the BEM
# Computing incident pressure at the collocation points
pI = P₀*exp.(im*k*mesh.sources[3,:]);
# Computing the BEM solution using dense matrices
Fp,_,Cp = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,progress=false);
Hp = Fp + Diagonal(1.0 .- Cp);
p_bem = gmres(Hp,pI;verbose=true);
# Computing BEM solution using the FMM operator
Ff = FMMFOperator(mesh,k);
Hf = Ff + 0.5I;
p_fmm = gmres(Hf,pI;verbose=true);
# Computing BEM solution using the H operator
Fh = HFOperator(mesh,k);
Hh = Fh + 0.5I;
p_h = gmres(Hh,pI;verbose=true);

# # Plotting solution
# For plotting purposes we must order the solution w.r.t the angles
surface_angles = acos.(mesh.sources[3,:]/a);
perm = sortperm(surface_angles);
# Plotting real part of pressure
plot(θ_analytical, real.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("re(p)");
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting imaginary part of pressure
plot(θ_analytical, imag.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("im(p)");
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting absolute pressure values
plot(θ_analytical, abs.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("|p|");
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)

# Plotting the solution on the sphere. Note that the mesh is plotted as linear due to the underlying library used, not because the mesh itself is linear.
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=abs.(data_viz))
wgl.Colorbar(fig[1,2],label="|p|"); fig
