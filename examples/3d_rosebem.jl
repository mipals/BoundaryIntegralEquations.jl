# # ROSEBEM: Rigid sphere scattering
# The Reduced Order Series Expansion Boundary Element Method (ROSEBEM) is a technique to significantly increase the computational efficiency of multifrequency BEM problems [panagiotopoulos2020](@cite) [Paltorp2023](@cite).
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations
using IterativeSolvers, SpecialFunctions, LegendrePolynomials, Plots
import BoundaryIntegralEquations: scattering_krylov_basis, taylor_assemble!, apply_taylor_expansion
# # Loading Mesh
# Loading and visualizing the triangular (spherical) mesh
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarser");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(mesh_file;physics_order=:disctriquadratic)
# # Setting up constants
frequency = 250.0;                              # Frequency                [Hz]
c  = 343;                                       # Speed up sound           [m/s]
a  = 1.0;                                       # Radius of sphere_1m      [m]
k  = 2π*frequency/c;                            # Wavenumber
P₀ = 1.0;                                       # Magnitude of planewave
# Furthmore the radius and colatitude angles for points where the analytical solution is evaluated
r  = 1.0;
θ_analytical = collect(0:0.01:π);
# # Analytical Solution
# The analytical solution of the scattering of a sphere by plane wave can be computed as ([ihlenburg1998a](@cite))
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
# Function to evalauate the analytical pressure  ([ihlenburg1998a](@cite))
function p_analytical(θ_analytical,r,k;N_trunc = 50,R=1)
    p_s = P₀*sum(n -> -c_n(n,k*R) .* Pl.(cos.(θ_analytical), n) .* sp_h.(n,k*r), 0:N_trunc)
    p_i = P₀*exp.(im*k*r*cos.(θ_analytical))
    return p_i, p_s
end
# # Solution using the ROSEBEM
# The ROSEBEM is based on a Taylor series expansion of the BEM matrices (``\mathbf{F}_m(k_0) \in \mathbb{C}^{n\times n}``) and a reduced basis (``\mathbf{U} \in \mathbb{C}^{n\times \ell}``).
# Using the Taylor expansion the BEM system at a different frequency (``k``) than the expansion frequency (``k_0``) can be assembled as
# ```math
#  \left(\mathbf{U}^{\text{H}}\text{diag}(\mathbf{c})\mathbf{U} + \sum_{m=0}^{M-1}\frac{(k-k_0)^m}{m!}\mathbf{U}^{\text{H}}\mathbf{F}_m(k_0)\mathbf{U}\right)\mathbf{p}_\ell = \mathbf{U}^{\text{H}}\mathbf{p}_\text{incident}(k),
# ```
# where ``\mathbf{p}_\ell`` is the so-called reduced variable. Note that in practice the ``\mathbf{F}_m``-matrice should never be stored directly. Instead ``\mathbf{U}^{\text{H}}\mathbf{F}_m\mathbf{U} \in \mathbb{C}^{\ell\times\ell}`` should be stored as it requires significantly less memory (as ``\ell < n``).
# After ``\mathbf{p}_\ell`` is computed the full pressure can be extracted as ``\mathbf{p} = \mathbf{U}\mathbf{p}_\ell``.
#
# First we setup of the model
k0 = 2π*225/c;  # Expansion wavenumber (for the Taylor series)
M = 25          # Number of terms in the Taylor series
L = 3           # Number of primary frequencyies used to compute the ROM Basis
klist = 2π*(LinRange(100,300,L))/c;  # Defining the primary frequencies
# For the computation of the reduced basis (`V`) solution (`sols`) at each frequency is used.
U,sols,_ = scattering_krylov_basis(mesh,klist;P₀=P₀,verbose=false,progress=false);
println("Reduced basis size: $(size(U,2)) | Reduction in DOF: $(1 - size(U,2)/size(U,1)) %");
# Given the ROM basis we can now compute the Taylor series including M terms
Fm, _, Cm = taylor_assemble!(mesh,k0,mesh.sources,mesh.shape_function;n=2,m=2,M=M,gOn=false,U=U,progress=false);
# In order to test the Taylor series we evaluate and compute the solution at k
Ft = apply_taylor_expansion(Fm,k,k0); # Taylor series of F at k
Ht = Ft + U'*Diagonal(1.0 .- Cm)*U;     # Adding the integral free term to the Taylor series
pI = P₀*exp.(im*k*mesh.sources[3,:]);   # Compute incident wave (rhs) at frequency k
p_rom = U*(Ht\(U'*pI));                 # Solving the ROM system + project back using V
# # Plotting solution
# For plotting purposes the solution is ordered w.r.t the colatitude
surface_angles = acos.(mesh.sources[3,:]/a);
perm = sortperm(surface_angles);
# Additionally we compute analytical incident and scattered pressure at the colatitude angles
p_i, p_s = p_analytical(θ_analytical,r,k);
p_t = p_s + p_i; # The total pressure as the sum of scattered and incident pressure
# Plotting the real, imaginarg and absolute values of the pressure
plot(θ_analytical, real.(p_t),label="Analytical",linewidth=2,color=:black)
plot!(θ_analytical, imag.(p_t),label=false,linewidth=2,color=:black);
plot!(θ_analytical, abs.(p_t),label=false,linewidth=2,color=:black);
plot!(surface_angles[perm],real.(p_rom[perm]),label="Re(ROSEBEM)",linestyle=:dash,linewidth=2);
plot!(surface_angles[perm],imag.(p_rom[perm]),label="Im(ROSEBEM)",linestyle=:dash,linewidth=2);
plot!(surface_angles[perm],abs.(p_rom[perm]),label="abs(ROSEBEM)",linestyle=:dash,linewidth=2);
xlabel!("Angle (rad)"); ylabel!("re(p) / im(p) / abs(p)")
# We can use `sols` to evaluate the errors at the primary frequencies
errors = zeros(length(klist))
for (e,k) in enumerate(klist)
    pIe = P₀*exp.(im*k*mesh.sources[3,:]);
    Fte = BoundaryIntegralEquations.apply_taylor_expansion(Fm,k,k0)
    Hte = Fte + U'*Diagonal(1.0 .- Cm)*U
    p_rome = U*(Hte\(U'*pIe));
    errors[e] = norm(sols[:,e] - p_rome)/norm(sols[:,e])
end
println(errors)

# # Evaluating the pressure for many frequencies
# The main advantange of using the ROSEBEM is when evaluating many frequencies. As such we here defining a list of frequencies (and the corresponding wavenumber)
frequencies = collect(10:1:400);
ks = frequencies/340*2π;
# In order to avoid spurious frequencies we additionally define some so-called CHIEF points
src_chief = 0.9*rand(3,10)/sqrt(3);
# Furthermore we want to evaluate the pressure at two field points: One directly in front and one directly in the back of the sphere at double the radius. As such
X_fieldpoints = [[0.0;0.0;2*r] [0.0;0.0;-2*r]];
# We now pre-allocate the output of pressure at the two field points
p1 = zeros(ComplexF64,length(ks));
p2 = zeros(ComplexF64,length(ks));
p1_chief = zeros(ComplexF64,length(ks));
p2_chief = zeros(ComplexF64,length(ks));
for (i,k) in enumerate(ks)
    Fp,_,_  = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false); # BEM matrix at field point
    Fti = apply_taylor_expansion(Fm,k,k0); # Evaluate Taylor series at k
    Hti = Fti + U'*Diagonal(1.0 .- Cm)*U;   # Adding integral free term to lhs
    pIi = P₀*exp.(im*k*mesh.sources[3,:]);  # Evaluating incident at colloation points (rhs)
    p_romi = U*(Hti\(U'*pIi));               # Computing ROM solution without CHIEF points
    p_field = -Fp*p_romi;   # Evaluating pressure at the two field points
    p1[i] = p_field[1];     # Saving pressure 1 (no CHIEF points)
    p2[i] = p_field[2];     # Saving pressure 2 (no CHIEF points)
    F_chief,_,_  = assemble_parallel!(mesh,k,src_chief,n=2,m=2,progress=false); # CHIEF-point BEM
    p_chief = P₀*exp.(im*k*src_chief[3,:]);                                     # CHIEF-point rhs
    Hti     = [Hti; F_chief*U];             # Adding CHIEF system to ROM system
    p_romi_chief  = U*(Hti\[(U'*pIi);p_chief]);   # Adding CHIEF rhs to ROM rhs
    p_field_chief = -Fp*p_romi_chief;   # Evaluating pressure at the two field points
    p1_chief[i] = p_field_chief[1];     # Saving pressure at point 1 (with CHIEF points)
    p2_chief[i] = p_field_chief[2];     # Saving pressure at point 2 (with CHIEF points)
end

# # Comparing ROSEBEM solution with the analytical solution
# First we start with the point located directly behind the sphere
p_i, p_s = p_analytical(0,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft,linewidth=2)
plot!(frequencies,abs.(p_i + p1),label="ROSEBEM",legend=:topleft,linestyle=:dash,linewidth=2)
plot!(frequencies,abs.(p_i + p1_chief),label="ROSEBEM-CHIEF",legend=:topleft,linestyle=:dash,linewidth=2)
ylims!((0.7,1.5)); xlabel!("Frequency [Hz]"); ylabel!("|p/p0|")
# While for the point located directly in front of the sphere the results look as follows
p_i, p_s = p_analytical(π,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft,linewidth=2)
plot!(frequencies,abs.(p_i + p2),label="ROSEBEM",legend=:topleft,linestyle=:dash,linewidth=2)
plot!(frequencies,abs.(p_i + p2_chief),label="ROSEBEM-CHIEF",legend=:topleft,linestyle=:dash,linewidth=2)
ylims!((0.50,1.6)); xlabel!("Frequency [Hz]"); ylabel!("|p/p0|")

# # Bibliography
# ```@bibliography
# Pages = []
# Canonical = false

# panagiotopoulos2020
# Paltorp2023
# ihlenburg1998a
# ```
