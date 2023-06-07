# # Rigid sphere scattering (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations
using IterativeSolvers, MeshViz, SpecialFunctions, LegendrePolynomials, Plots
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(600, 600))
using JSServe                           #hide
Page(exportable=true, offline=true)     #hide
# # Loading Mesh
# Loading and visualizing the triangular (spherical) mesh
# mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarser");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
# mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
# mesh = load3dTriangularComsolMesh(mesh_file;physics_order=:disctrilinear)
mesh = load3dTriangularComsolMesh(mesh_file;physics_order=:disctriquadratic)
# # Setting up constants
frequency = 250.0;                              # Frequency                [Hz]
c  = 343;                                       # Speed up sound           [m/s]
a  = 1.0;                                       # Radius of sphere_1m      [m]
k  = 2π*frequency/c;                            # Wavenumber
P₀ = 1.0;                                       # Magnitude of planewave
# # Defining radius and colatitude angles for points where the analytical solution is evaluated
r  = 1.0;
θ_analytical = collect(0:0.01:π);
# # Analytical Solution
# The analytical solution of the scattering of a sphere by plane wave can be computed as (Ihlenburg1998)
# ```math
#  p_\text{analytical}(r, \theta) = P_0\left(\exp(\mathrm{i}kr\cos(\theta)) - \sum_{n=1}^\infty \mathrm{i}^n(2n+1)\frac{j_n^{'}(ka)}{h_n^{'}(ka)}P_n(\cos(\theta))h_n(kr)\right),
# ```
# where ``j_n, h_n`` and ``P_n`` are respectively the spherical Bessel function of the first kind, the Hankel function of the first kind and the Legendre polynomial of degree ``n``.
# To make the implementation easier we defin the following helper functions
dsp_j(n,z) = n/z*sphericalbesselj(n,z) - sphericalbesselj(n+1,z); # Derivative of j
dsp_y(n,z) = n/z*sphericalbessely(n,z) - sphericalbessely(n+1,z); # Derivative of y
sp_h(n,z)  = sphericalbesselj(n,z) + im*sphericalbessely(n,z);    # Hankel function (h)
dsp_h(n,z) = dsp_j(n,z) + im*dsp_y(n,z);                          # Derivative of h
# Using the helper functions we can define the a function for the coefficients
c_n(n,ka)  = (im)^n*(2n + 1)*(dsp_j.(n,ka)./dsp_h.(n,ka));
# The analytical pressure as definde in Ihlenburg1998
function p_analytical(θ_analytical,r,k;N_trunc = 50,R=1)
    p_s = P₀*sum(n -> -c_n(n,k*R) .* Pl.(cos.(θ_analytical), n) .* sp_h.(n,k*r), 0:N_trunc)
    p_i = P₀*exp.(im*k*r*cos.(θ_analytical))
    return p_i, p_s
end
# Evaluating the incident and scattered pressure at surface points
p_i, p_s = p_analytical(θ_analytical,r,k);
p_t = p_s + p_i; # The total pressure as the sum of scattered and incident pressure
# # Solution using the ROSEBEM
# First we setup of the model
k0 = 2π*225/c;  # Expansion wavenumber (for the Taylor series)
M = 25          # Number of terms in the Taylor series
L = 3           # Number of primary frequencyies used to compute the ROM Basis
klist = 2π*(LinRange(100,300,L))/c  # Defining the primary frequencies
# We now compute the ROM basis (V) and the solution at each frequency (sols)
V,sols,_ = BoundaryIntegralEquations.scattering_krylov_basis(mesh,klist;P₀=P₀)
# Given the ROM basis we can now compute the Taylor series including M terms
Fv, Gv, Cv = BoundaryIntegralEquations.taylor_assemble!(mesh,k0,mesh.sources,mesh.shape_function;n=2,m=2,M=M,gOn=false,V=V)
# In order to test the Taylor series we evaluate and compute the solution at k
Ft = BoundaryIntegralEquations.apply_taylor_expansion(Fv,k,k0); # Taylor series of F at k
Ht = Ft + V'*Diagonal(1.0 .- Cv)*V;     # Adding the integral free term to the Taylor series
pI = P₀*exp.(im*k*mesh.sources[3,:]);   # Compute incident wave (rhs) at frequency k
p_rom = V*(Ht\(V'*pI));                 # Solving the ROM system + project back using V
# # Plotting solution
# For plotting purposes we must order the solution w.r.t the angles
surface_angles = acos.(mesh.sources[3,:]/a);
perm = sortperm(surface_angles);
# Plotting real part of pressure
plot(θ_analytical, real.(p_t),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("re(p)");
plot!(surface_angles[perm],real.(p_rom[perm]),label="ROM",linestyle=:dash,linewidth=2)
# Plotting imaginary part of pressure
plot(θ_analytical, imag.(p_t),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("im(p)");
plot!(surface_angles[perm],imag.(p_rom[perm]),label="ROM",linestyle=:dash,linewidth=2)
# Plotting absolute pressure values
plot(θ_analytical, abs.(p_t),label="Analytical",linewidth=2)
xlabel!("Angle (rad)"); ylabel!("|p|");
plot!(surface_angles[perm],abs.(p_rom[perm]),label="ROM",linestyle=:dash,linewidth=2)

# We can use `sols` to evaluate the errors as
errors = zeros(length(klist))
for (i,k) in enumerate(klist)
    pI = P₀*exp.(im*k*mesh.sources[3,:]);
    Ft = BoundaryIntegralEquations.apply_taylor_expansion(Fv,k,k0)
    Ht = Ft + V'*Diagonal(1.0 .- Cv)*V
    p_rom = V*(Ht\(V'*pI));
    errors[i] = norm(sols[:,i] - p_rom)/norm(sols[:,i])
end
println(errors)

# # Evaluating the pressure for many frequencies
# The main advantange of using the ROSEBEM is when evaluating many frequencies. As such we here defining a list of frequencies (and the corresponding wavenumber)
frequencies = collect(10:1:400)
ks = frequencies/340*2π
# In order to avoid spurious frequencies we additionally defining some so-called CHIEF points
src_chief = 0.9*rand(3,10)/sqrt(3)
# Furthermore we want to evaluate the pressure at two field points: One directly in front and one directly in the back of the sphere at double the radius. As such
X_fieldpoints = [[0.0;0.0;2*r] [0.0;0.0;-2*r]]
# We now pre-allocate the output of pressure at the two field points
p1,p2 = zeros(ComplexF64,length(ks)), zeros(ComplexF64,length(ks))
@time begin
for (i,k) in enumerate(ks)
    # Evaluating Taylor series (lhs) + new incident wave (rhs)
    Ft = BoundaryIntegralEquations.apply_taylor_expansion(Fv,k,k0)
    Ht = Ft + V'*Diagonal(1.0 .- Cv)*V
    pI = P₀*exp.(im*k*mesh.sources[3,:]);
    # Computing BEM equation for the CHIEF points
    Fc,_,_  = assemble_parallel!(mesh,k,src_chief,n=2,m=2,progress=false);
    p_chief = P₀*exp.(im*k*src_chief[3,:]);
    # Adding CHIEF points to the Taylor system
    Ht = [Ht; Fc*V]
    p_rom = V*(Ht\[(V'*pI);p_chief]);
    # Evluating the pressure at the two field points
    Fp,_,_ = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
    p_field  = -Fp*p_rom;
    # Saving results
    p1[i] = p_field[1];
    p2[i] = p_field[2];
    println("Frequency: $(frequencies[i])Hz \t | Iteration $(i)\t/$(length(ks))")
end
end

# We now compare the ROSEBEM solution with the analytical solution
p_i, p_s = p_analytical(0,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft)
plot!(frequencies,abs.(p_i + p1),label="ROSEBEM",legend=:topleft)
ylims!((0.7,1.5))
xlabel!("Frequency [Hz]")
ylabel!("|p/p0|")

p_i, p_s = p_analytical(π,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft)
plot!(frequencies,abs.(p_i + p2),label="ROSEBEM",legend=:topleft)
ylims!((0.50,1.6))
xlabel!("Frequency [Hz]")
ylabel!("|p/p0|")
