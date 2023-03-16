# # Rigid sphere scattering (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations,IterativeSolvers, Plots, MeshViz
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
surface_angles = acos.(mesh.sources[1,:]/a);
# Now we can compute the analytical solution
p_analytical, _ = BoundaryIntegralEquations.plane_wave_scattering_sphere(k,a,a,surface_angles,1e-6);
# # Solution using the BEM
# Computing incident pressure
pI = BoundaryIntegralEquations.incoming_wave(θ,1.0,mesh.sources,k);
# Computing the BEM solution using dense matrices
Fp,_,Cp = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,progress=false);
Hp = Diagonal(1.0 .- Cp) + Fp;
p_bem = gmres(Hp,pI;verbose=true);
# Computing BEM solution using the FMM operator
Hf = FMMHOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2)
p_fmm = gmres(Hf,pI;verbose=true);
# Computing BEM solution using the H operator
Hh = HHOperator(mesh,k;nearfield=true,n_gauss=3,offset=0.2)
p_h = gmres(Hh,pI;verbose=true);

# # Plotting solution
# For plotting purposes we must order the solution w.r.t the angles
perm = sortperm(surface_angles);
# Plotting real part of pressure
plot(surface_angles[perm], real.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("re(p)");
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting imaginary part of pressure
plot(surface_angles[perm], imag.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("im(p)");
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting absolute pressure values
plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|p|");
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)

# Plotting solution
simple_mesh = create_simple_mesh(mesh)
cols = p_bem[sort!(unique(mesh.physics_topology[1:3,:]))]
fig, ax, hm = viz(simple_mesh;showfacets=true, color=abs.(cols))
wgl.Colorbar(fig[1,2],label="|p|"); fig
