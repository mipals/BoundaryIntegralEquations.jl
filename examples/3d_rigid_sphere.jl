# # Rigid sphere scattering (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations,IterativeSolvers, Plots, MeshViz
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
c = 340;                                        # Speed up sound           [m/s]
a = 1.0;                                        # Radius of sphere_1m      [m]
θ = [π/2 0.0];                                  # Angles of incoming wave  [radians]
k = 2π*frequency/c;                             # Computing the wavenumber
# # Analytical Solution
# The analytical solution of the scattering of a sphere by plane wave can be computed as, Ihlenburg1998
# ```math
#  p_s(r, \theta) = -P_0\sum_{n}^\infty \mathrm{i}^n(2n+1)\frac{j_n^{'}(ka)}{h_n^{'}(ka)}P_n(\cos(\theta))h_n(kr)
# ```
# To use the formula we must compute the angles of the collocation points
# Now we can compute the analytical solution
θ_analtrical = collect(0:.01:π)
p_analytical, _ = BoundaryIntegralEquations.plane_wave_scattering_sphere(k,a,a,θ_analtrical,1e-6);
# # Solution using the BEM
# Computing incident pressure
pI = BoundaryIntegralEquations.incoming_wave(θ,1.0,mesh.sources,k);
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
surface_angles = acos.(mesh.sources[1,:]/a);
perm = sortperm(surface_angles);
# Plotting real part of pressure
plot(θ_analtrical, real.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("re(p)");
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting imaginary part of pressure
plot(θ_analtrical, imag.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("im(p)");
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting absolute pressure values
plot(θ_analtrical, abs.(p_analytical),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|p|");
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)

# Plotting the solution on the sphere. Note that the mesh is plotted as linear due to the underlying library used, not because the mesh itself is linear.
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=abs.(data_viz))
wgl.Colorbar(fig[1,2],label="|p|"); fig
