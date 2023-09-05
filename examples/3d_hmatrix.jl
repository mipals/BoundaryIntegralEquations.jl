# # Scattering of a hard sphere (H-Matrix)
# # Importing related packages
using LinearAlgebra
using BoundaryIntegralEquations
using IterativeSolvers
using Plots
using Meshes
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading Mesh
# Defining possible combinations of geometry and physics orders
geometry_orders     = [:linear,:quadratic];
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
# Loading and visualizing the triangular (spherical) mesh
#src tri_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
tri_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src tri_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src tri_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                physics_order=tri_physics_orders[2])
simple_tri_mesh = create_simple_mesh(mesh)
viz(simple_tri_mesh;showfacets=true)
# # Setting up constants
# Defining the frequency of interest
freq = 100.0;                                    # Frequency [Hz]
# Computing constants from frequency
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1);
zk = Complex(ka);                                # Making sure that the wavenumber is complex
radius = 1.0 ;                                   # Radius of sphere_1m       [m]
angles = [π/2 0.0];                              # Angles of incoming wave   [radians]
# # Solution using the BEM
# Computing incident pressure
pI = BoundaryIntegralEquations.incoming_wave(angles,1.0,mesh.sources,zk);
# Computing the BEM solution using dense matrices
Fp,_,Cp = assemble_parallel!(mesh,zk,mesh.sources,n=2,m=2,progress=false);
Ap = Diagonal(1.0 .- Cp) + Fp;
p_bem = gmres(Ap,pI;verbose=true);
# Computing BEM soltuiong using the FMM operator
Ah = BoundaryIntegralEquations.HFOperator(mesh,zk;nearfield=true,n_gauss=3,offset=0.2) + 0.5I;
p_h = gmres(Ah,pI;verbose=true);

# # Plotting solution
# First setting up analytical solution
surface_angles = acos.(mesh.sources[1,:]/radius)
perm = sortperm(surface_angles)
p_analytical, _ = BoundaryIntegralEquations.plane_wave_scattering_sphere(zk,radius,1.0,surface_angles,1e-6);
# Plotting real part of pressure
plot(surface_angles[perm], real.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("re(p)");
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting imaginary part of pressure
plot(surface_angles[perm], imag.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("im(p)");
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
# Plotting absolute pressure values
plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|p|");
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_h[perm]),label="H-matrix",linestyle=:dash,linewidth=2)
