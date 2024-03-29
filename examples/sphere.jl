#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using BoundaryIntegralEquations
using LinearAlgebra
using Plots
using SpecialFunctions
using IterativeSolvers
import BoundaryIntegralEquations: incoming_wave, plane_wave_scattering_sphere, psca, pinc
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
quad_physics_orders = [:linear,:geometry,:discquadconstant,:discquadlinear,:discquadquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                        physics_order=tri_physics_orders[2])
# Quadrilateral Meshes
# quad_mesh_file = "examples/meshes/quad_sphere"
# quad_mesh_file = "examples/meshes/quad_sphere_1m_fine"
# quad_mesh_file = "examples/meshes/quad_sphere_1m_finer"
# mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=geometry_orders[1],
                                            # physics_order=quad_physics_orders[5])
#==========================================================================================
            3d Visualization - Seems highly unstable on M1. Problems with GLMakie?
==========================================================================================#
# using Meshes, MeshViz
# # ##choose a Makie backend
# # import GLMakie as Mke
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq  = 1000.0                                   # Frequency                 [Hz]
c     = 340.0                                   # Speed of sound            [m/s]
k     = 2*π*freq/c                              # Wavenumber                [1/m]
angles = [π/2 0.0]                              # Angles of incoming wave   [radians]
radius = 1.0                                    # Radius of sphere_1m       [m]
# Computing incident pressure
pI = BoundaryIntegralEquations.incoming_wave(angles,1.0,mesh.sources,k)
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
# @time Fp,_,Cp = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,gOn=false,sparse=false);
@time Fp,Gp,Cp = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,sparse=false);
#==========================================================================================
            Setting up a linear system and solving for the pressure
==========================================================================================#
# Assembling linear system (we )
# Ap    = Fp + Diagonal(1.0 .- Cp)
SF,SG = assemble_parallel!(mesh,k,mesh.sources;sparse=true,depth=2,progress=true,offset=0.15)
# Flu = lu(SF + Diagonal(1.0 .- Cp))

# p_bem = gmres(Ap,pI;verbose=true,Pl=Flu);
# p_bem = gmres(Ap,pI;verbose=true);

Glu = lu(SG)
p_bem = gmres(Gp,pI;verbose=true,Pr=Glu);
p_bem = gmres(Gp,pI;verbose=true);

T = Matrix(SG)\Gp
gmres(T,pI;verbose=true);

# p_bem = Ap\pI
## Plotting pressure on surface nodes
surface_angles = acos.(mesh.sources[1,:]/radius)
perm = sortperm(surface_angles)
p_analytical, _ = BoundaryIntegralEquations.plane_wave_scattering_sphere(k,radius,1.0,surface_angles,1e-6)
plot(surface_angles[perm], real.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot(surface_angles[perm], imag.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
#==========================================================================================
                            Adding CHIEF points
==========================================================================================#
# Defining CHIEF points inside of sphere
x = collect(range(-radius*0.30,0.60*radius,length=5))
y = collect(range(-radius*0.65,0.50*radius,length=5))
z = collect(range(-radius*0.11,0.55*radius,length=5))
xyzb_chief=[x y z]';
# Computing incident pressure on nodes including CHIEF points
pI_chief  = incoming_wave(angles,1.0,[mesh.sources xyzb_chief],k)
# Assembling field point matrices
Fpc,Gpc,_ = assemble_parallel!(mesh,k,xyzb_chief;n=4,m=4);
# Defining BEM systen
A_chief = [Ap;Fpc]
p_chief = A_chief\pI_chief

## Plotting pressure on surface nodes
plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
# scatter!(surface_angles,abs.(p_chief),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_chief[perm]),label="BEM",linestyle=:dash,linewidth=2)
title!("Frequency = $(freq) (Hz)")
#==========================================================================================
                                    Figure 1 + 2 from:
    The Burton and Miller Method: Unlocking Another Mystery of Its Coupling Parameter
==========================================================================================#
using LaTeXStrings
freqs = collect(1:10:1500)
klist = freqs/340*2π
angle = 0.0
pscat12 = psca(2.0*radius,angle,klist;N=60,R=R)
pinco12 = pinc(2.0*radius,angle,klist)
ptot12 = pinco12 + pscat12
plot(freqs,abs.(ptot12),label="Analytical Solution",legend=:topleft)
ylims!((0.7,1.5))
xlabel!("Frequency [Hz]")
ylabel!(L"p/p_0")
#==========================================================================================
                                    Figure 3 from:
    The Burton and Miller Method: Unlocking Another Mystery of Its Coupling Parameter
==========================================================================================#
angle  = π
pscat3 = psca(2.0*radius,angle,klist;N=60,R=radius)
pinco3 = pinc(2.0*radius,angle,klist)
ptot3 = pinco3 + pscat3
plot(freqs,abs.(ptot3),label="Analytical Solution",legend=:topleft)
ylims!((0.50,1.6))
xlabel!("Frequency [Hz]")
ylabel!(L"p/p_0")



import BoundaryIntegralEquations: get_topology,number_of_elements,get_coordinates,getpolar_gaussian
import BoundaryIntegralEquations: create_rotated_element, copy_interpolation_nodes!,set_interpolation_nodes!
physics_function = mesh.physics_function

beta = BoundaryIntegralEquations.get_beta(physics_function)
tmp  = DiscontinuousTriangularLinear(physics_function,beta)
offr = BoundaryIntegralEquations.get_offset(physics_function)
# Dealing with corners
DO = Diagonal(ones(3))
nodesX1,nodesY1,weights1 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[1.00;offr;offr],DO,1e-6)
nodesX2,nodesY2,weights2 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[offr;1.00;offr],DO,1e-6)
nodesX3,nodesY3,weights3 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[offr;offr;1.00],DO,1e-6)
nodesX4,nodesY4,weights4 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[0.50;0.50;offr],DO,1e-6)
nodesX5,nodesY5,weights5 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[offr;0.50;0.50],DO,1e-6)
nodesX6,nodesY6,weights6 = BoundaryIntegralEquations.singular_triangle_integration(tmp,3,[0.50;offr;0.50],DO,1e-6)
