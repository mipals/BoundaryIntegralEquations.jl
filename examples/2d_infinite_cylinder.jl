# # Scattering of infinite cylinder (2D)
# # Importing relevant packages
using BoundaryIntegralEquations
using Plots
using LinearAlgebra
import BoundaryIntegralEquations: mesh_circle, plane_wave_scattering_circle

# # Setting up the system
n_elements = 20;    # Number of elements
freq = 100.0;       # Frequency                 (Hz)
c    = 340.0;       # Speed of sound            (m/s)
k    = 2*π*freq/c;  # Wavenumber                (1/m)
r    = 1.0;         # Circle radius             (m)
# # Loading and plotting meshes
mesh_lin       = mesh_circle(ContinuousCurveLinear(3),n_elements;radius=r)
mesh_quad      = mesh_circle(ContinuousCurveQuadratic(3),n_elements;radius=r)
mesh_disc_con  = mesh_circle(ContinuousCurveLinear(3),DiscontinuousCurveConstant(3),n_elements;radius=r)
mesh_disc_lin  = mesh_circle(ContinuousCurveQuadratic(3),DiscontinuousCurveLinear(3),n_elements;radius=r)
mesh_disc_quad = mesh_circle(ContinuousCurveQuadratic(3),DiscontinuousCurveQuadratic(3),n_elements;radius=r)
meshes = [mesh_lin,mesh_quad,mesh_disc_con,mesh_disc_lin,mesh_disc_quad]
mesh_types = ["Linear Geometry, Continuous Linear Physics",
              "Quadratic Geometry, Continuous Quadratic Physics",
              "Linear Geometry, Discontinuous Constant Physics",
              "Quadratic Geometry, Discontinuous Linear Physics",
              "Quadratic Geometry, Discontinuous Quadratic Physics"]
begin #hide
for (mesh,plot_title) in zip(meshes,mesh_types)
    plt1=plot(mesh.coordinates[1,:],mesh.coordinates[2,:],aspect_ratio=1,label=false)
    scatter!(mesh.sources[1,:],mesh.sources[2,:],label="Collocation"); xlabel!("x"); ylabel!("y")
    scatter!(mesh.coordinates[1,:],mesh.coordinates[2,:],aspect_ratio=1,label="Geometry nodes")
    title!(plot_title)
    display(plt1)
end
end #hide

# # Solving the scattering of hard infinite cylinder for each mesh
# We first start with computing the analytical solution
src_angle = angle.(mesh_disc_quad.sources[1,:] + im*mesh_disc_quad.sources[2,:])
p = plane_wave_scattering_circle(src_angle,k*r,150)
plt_pressure = plot(src_angle,abs.(p),label="Analytical")
# We now loop over all elements
begin #hide
for mesh in meshes
    F,_,C=assemble_parallel!(mesh,k,mesh.sources;n=4,gOn=false,progress=false);
    pI = exp.(im*k*mesh.sources[1,:])
    src_angle = angle.(mesh.sources[1,:] + im*mesh.sources[2,:])
    pB = (F + Diagonal(C .- 1.0))\pI
    plot!(plt_pressure,src_angle,abs.(pB),label="BEM",linestyle=:dash)
end
end #hide
# Displaying the solution
current() #hide
