#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using SpecialFunctions, LinearAlgebra, IntegralEquations, Plots
import IntegralEquations: incoming_wave, plane_wave_scattering_sphere, psca, pinc
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
# Triangular Meshes
tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# mesh  = load3dTriangularComsolMesh(tri_mesh_file)
# mesh  = load3dTriangularComsolMesh(tri_mesh_file;physics_order=:disctriconstant)
# mesh  = load3dTriangularComsolMesh(tri_mesh_file;physics_order=:disctrilinear)
# mesh  = load3dTriangularComsolMesh(tri_mesh_file;physics_order=:disctriquadratic)
# mesh  = load3dTriangularComsolMesh(tri_mesh_file;physics_order=:disctriconstant,geometry_order=:linear)

# Quadrilateral Meshes
quad_mesh_file = "examples/meshes/quad_sphere_1m_finer"
# quad_mesh_file = "examples/meshes/quad_sphere_1m_finer"
mesh = load3dQuadComsolMesh(quad_mesh_file)
# mesh = load3dQuadComsolMesh(quad_mesh_file;physics_order=:linear)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq  = 250.0                                   # Frequency                 [Hz]
c     = 340.0                                   # Speed of sound            [m/s]
k     = 2*π*freq/c                              # Wavenumber                [1/m]
el_wl =   6*freq/c                              # Nyquist sampling
angles = [π/2 0.0]                              # Angles of incoming wave   [radians]
radius = 1.0                                    # Radius of sphere_1m       [m]
# Defining CHIEF points inside of sphere
x = collect(range(-radius*0.30,0.60*radius,length=5))
y = collect(range(-radius*0.65,0.50*radius,length=5)) 
z = collect(range(-radius*0.11,0.55*radius,length=5))
xyzb_chief=[x y z]';
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
Fp,Gp,Cp  = assemble_parallel!(mesh,k,mesh.sources);
Fpc,Gpc,_ = assemble_parallel!(mesh,k,xyzb_chief;n=2,m=2);
#==========================================================================================
            Setting up a linear system and solving for the pressure
==========================================================================================#
if typeof(mesh.shape_function) <: Quadrilateral
    Ap  = (Fp + Diagonal(1.0 .- Cp))
    FAp = [Ap;Fpc]
else
    Ap  = (-Fp + Diagonal(1.0 .+ Cp))
    FAp = [Ap;-Fpc]
end
pI  = incoming_wave(angles,1.0,[mesh.sources xyzb_chief],k)
p_bem = FAp\pI

## Plotting pressure on surface nodes
surface_angles = acos.(-mesh.sources[1,:]/radius)
perm = sortperm(surface_angles)
p_analytical, _ = plane_wave_scattering_sphere(k,radius,1.0,surface_angles,1e-6)
plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
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
