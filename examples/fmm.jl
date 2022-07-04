using Pkg
Pkg.add(PackageSpec(url="https://github.com/flatironinstitute/FMM3D.git", subdir="julia/"))
#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using FMM3D
using LinearAlgebra
using IntegralEquations
using IterativeSolvers
using Plots
#==========================================================================================
            Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])
#==========================================================================================
                            Setting up constants
==========================================================================================#
# Defining Frequency
freq = 100.0
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
zk = Complex(kp)
radius = 1.0                                    # Radius of sphere_1m       [m]
# Computing incident pressure
angles = [π/2 0.0]                              # Angles of incoming wave   [radians]
pI = IntegralEquations.incoming_wave(angles,1.0,mesh.sources,zk)
#==========================================================================================
                            Defining LossyBlockMatrix
            Block matrix corresponding to 5 BEM systems and 5 constraints
==========================================================================================#
Fp,Gp,Cp = assemble_parallel!(mesh,zk,mesh.sources,n=2,m=2,sparse=false);
Ap = Fp + Diagonal(1.0 .- Cp);
# Ag = FMMGOperator(mesh,zk)
Af = FMMFOperator(mesh,zk;n=3)
# xg = ones(eltype(Ag),size(Ag,1))
xf = ones(eltype(Af),size(Af,1))

# Checking multiplation
@time maximum(abs.(Af*xf - Ap*xf))
# @time maximum(abs.(Ag*xg - Gp*xg))

# p = rand(eltype(xf),size(xf))
# @time IntegralEquations.nodes_to_gauss!(Af.tmp,Af.element_interpolation,Af.physics_topology,p)

# Solving the scattering problem
p_bem = gmres(Ap,pI;verbose=true);
p_fmm = gmres(Af,pI;verbose=true);

# Plotting solution
surface_angles = acos.(mesh.sources[1,:]/radius)
perm = sortperm(surface_angles)
p_analytical, _ = IntegralEquations.plane_wave_scattering_sphere(zk,radius,1.0,surface_angles,1e-6)
plot(surface_angles[perm], real.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot(surface_angles[perm], imag.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)

plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
