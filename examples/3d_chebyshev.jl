# # ROSEBEM - Chebyshev
# The Reduced Order Series Expansion Boundary Element Method (ROSEBEM) is a technique to significantly increase the computational efficiency of multifrequency BEM problems [panagiotopoulos2020](@cite) [Paltorp2023](@cite) [panagiotopoulos2022a](@cite). In this example we look at the scattering of a rigid sphere and how the ROSEBEM using a Chebyshev series expansion can be used to accelerate the computational efforts.
# # Importing related packages
using BoundaryIntegralEquations # For BIEs
using LegendrePolynomials       # For Legendre Polynomials
using SpecialFunctions          # For Bessel functions
using IterativeSolvers          # For gmres
using LinearAlgebra             # For Diagonal
using Meshes                    # For using `viz`
using Plots                     # For 2d plots
import GLMakie as wgl           # For 3d plotting
wgl.set_theme!(resolution=(800, 800)) #hide
# # Setting up constants
c  = 343;                           # Speed up sound           (m/s)
a  = 1.0;                           # Radius of sphere_1m      (m)
P₀ = 1.0;                           # Magnitude of planewave   (Pa)
r  = 1.0;                           # Radius of sphere         (m)
θ_analytical = collect(0:0.01:π);   # Colaltitude angles
# # Loading Mesh
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarser");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m");
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_fine");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_finer");
#src mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_extremely_fine");
mesh = load3dTriangularComsolMesh(mesh_file;physics_order=:disctriquadratic)
simple_mesh = create_simple_mesh(mesh);
fig = viz(simple_mesh;showfacets=true);
wgl.save(joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","docs","src","examples","3d_chebyshev.png"),fig) #hide
# ![](3d_chebyshev.png)
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
# The ROSEBEM is based on a series expansion of the BEM matrices and a reduced basis (``\mathbf{U} \in \mathbb{C}^{n\times \ell}``). In [another example](3d_rosebem.md) the Taylor series expansion the BEM system is used. In this example the Chebyshev is applied (similar to the approach in [panagiotopoulos2022a](@cite)). In short this means approximating the systems as
# ```math
#   \left(\sum_{j=0}^{M-1}\mathbf{C}_jc_j\left(g(k)\right)\right)\mathbf{x}(k) = \mathbf{U}^{\text{H}}\mathbf{p}_\text{incident}(k)
# ```
# where ``c_j(\omega)`` are the Chebyshev polynomials of the first kind, ``g(k)`` is a function that maps ``k`` from the frequency range of interst onto ``[-1,1]``, and ``\mathbf{x}`` is the unknown vector. Futhermore the coefficient matrices (``\mathbf{C}_j``) is computed similarly to a standard Chebyshev approximation
# ```math
# \mathbf{C}_j = \frac{2}{M }\sum_{i=0}^{M-1}\left(\mathbf{U}^{\text{H}}\mathbf{A}\left(g^{-1}(\omega_i\right)\mathbf{U}\right)c_j(\omega_i), \quad \omega_i \in [-1,1]
# ```
# with ``\omega_i = \cos\left(\frac{\pi\left(i + \frac{1}{2}\right)}{M}\right)``, ``i = 0,\dots,M-1`` and ``g^{-1}(\omega)`` being the inverse of ``g`` i.e. the function that maps ``[-1,1]`` to the frequency range of interest e.g. mapping to ``[k_\text{min}, k_\text{max}]`` we have that
# ```math
# \begin{aligned}
# g(k)           &= \frac{2}{k_\text{max} - k_\text{min}}k      - \frac{k_\text{max} + k_\text{min}}{k_\text{max} - k_\text{min}}\\
# g^{-1}(\omega) &= \frac{k_\text{max} - k_\text{min}}{2}\omega + \frac{k_\text{max} + k_\text{min}}{2}
# \end{aligned}
# ```
# Note that the above computes the coefficients from the reduced systems, i.e. after ``\mathbf{U}`` as been applied. In practice this means that an acceleration technique can be used generate the matrices. Another important aspect is that the above is non-intrusive, making it easy to implement in existing BEM software.
#
# We now define the relevant functions required to compute the reduced basis matrices and Chebyshev coefficient matrices. Notice that this is all done using a "matrix-free" approach, meaning that the large ``\mathbf{A}``-matrices is never stored in memory. Instead the multiplication with the matrix is approximated using the fast multipole method.
function create_basis_matrices(mesh,k_secondary,U;n_gauss=3)
    temp = zeros(eltype(U), size(U))
    output = zeros(eltype(U),size(U,2),size(U,2),length(k_secondary))
    for (index,k) in enumerate(k_secondary)
        F = FMMFOperator(mesh,k;n_gauss=n_gauss); # We're actually interseted in H = F + 0.5I
        mul!(temp,F,U) # The mul! currently only works directly on the FMMFOperator
        output[:,:,index] = U'*(temp + U/2) # We add half of a U due to 0.5I in the H-operator.
    end
    return output
end
function create_chebyshev_coefficients(A_reduced)
    M = size(A_reduced,3)
    Cj = similar(A_reduced)
    for j = 0:M-1
        Cj[:,:,j+1] = 2/M*sum(i -> A_reduced[:,:,i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return Cj
end
function chebyshev_eval(x,M)
    output = ones(M)
    output[2] = x
    for i = 3:M
        output[i] = 2x*output[i-1] - output[i-2]
    end
    return output
end
function eval_chebyshev_matrix(Cj,coeffs)
    output = -Cj[:,:,1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*Cj[:,:,index]
    end
    return output
end
function assemble_chebyshev(Cj,k)
    return eval_chebyshev_matrix(Cj,chebyshev_eval(k,size(Cj,3)))
end
# First we define a set of primary wavenumbers for which we compute the reduced basis (`U`) and the solution (`sols`)
L = 3;           # Number of primary wavenumbers used to compute the ROM Basis
k_primary = 2π*(LinRange(100,300,L))/c;  # Defining the primary frequencies
U,sols,_ = scattering_krylov_basis(mesh,k_primary;P₀=P₀,verbose=false,progress=false);
@info "Reduced basis size: $(size(U,2)) | Reduction in DOF: $(1 - size(U,2)/size(U,1)) %";
# Then we define the secondary wavenumbers, i.e. the wavenumbers for which we compute the matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega)\right)\mathbf{U}``
kmin = 2π*10/c;
kmax = 2π*400/c;
g(k)    = 2/(kmax - kmin)*k .- (kmax + kmin )/(kmax - kmin);
ginv(ω) = (kmax - kmin)/2*ω .+ (kmin + kmax)/2;
M  = 25;                                # Number of terms in the Chebyshev approximation
ωᵢ = cos.(π*(collect(0:M-1) .+ 1/2)/M); # Zeros of the Chebyshev polynomials
k_secondary = ginv.(ωᵢ);                # Mapping [-1,1] to [kmin, kmax]
# Using the defined functions we can create the basis matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega_i)\right)\mathbf{U}``
A_reduced = create_basis_matrices(mesh,k_secondary,U);
# Now from the basis matrices we can compute the Chebyshev coefficient matrices (``\mathbf{C}_j``) as
Cj  = create_chebyshev_coefficients(A_reduced);
# In order to avoid spurious frequencies we additionally define some so-called CHIEF points
src_chief = 0.9*rand(3,30)/sqrt(3);
# Furthermore we want to evaluate the pressure at two field points: One directly in front and one directly in the back of the sphere at double the radius. As such
X_fieldpoints = [[0.0;0.0;2*r] [0.0;0.0;-2*r]];
# We define the frequency range of interest
frequencies = collect(10:1:400);
ks = frequencies/c*2π;
# Lastly we pre-allocate the output of pressure at the two field points
p1_chebyshev = zeros(ComplexF64,length(ks));
p2_chebyshev = zeros(ComplexF64,length(ks));
for (i,k) in enumerate(ks)
    Fp,_,_  = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false); # BEM matrix at field point
    F_chief,_,_  = assemble_parallel!(mesh,k,src_chief,n=2,m=2,progress=false); # CHIEF-point BEM
    p_chief  = P₀*exp.(im*k*src_chief[3,:]);                                    # CHIEF-point rhs
    k_scaled = g(k);                       # Scaling the wavenumber to the interval [-1,1]
    Hti = assemble_chebyshev(Cj,k_scaled); # Evaluating the series expansion at the scaled k value
    pIi = P₀*exp.(im*k*mesh.sources[3,:]); # Evaluating incident at colloation points (rhs)
    Hti     = [Hti; F_chief*U];            # Adding CHIEF system to ROM system
    p_romi  = U*(Hti\[(U'*pIi);p_chief]);  # Adding CHIEF rhs to ROM rhs
    p_field = -Fp*p_romi;                  # Evaluating pressure at the two field points
    p1_chebyshev[i] = p_field[1];          # Saving pressure at field point 1
    p2_chebyshev[i] = p_field[2];          # Saving pressure at field point 2
end

# # Comparing ROSEBEM solution with the analytical solution
# First we start with the point located directly behind the sphere
p_i, p_s = p_analytical(0,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft,linewidth=2)
plot!(frequencies,abs.(p_i + p1_chebyshev),label="Chebyshev",legend=:topleft,linestyle=:dash,linewidth=2)
ylims!((0.7,1.5)); xlabel!("Frequency [Hz]"); ylabel!("|p/p0|")
scatter!(k_secondary*340/(2π),0.7ones(length(k_secondary)),label="ωᵢ")
# While for the point located directly in front of the sphere the results look as follows
p_i, p_s = p_analytical(π,2*r,ks;N_trunc = 80)
plot(frequencies,abs.(p_i + p_s),label="Analytical Solution",legend=:topleft,linewidth=2)
plot!(frequencies,abs.(p_i + p2_chebyshev),label="Chebyshev",legend=:topleft,linestyle=:dash,linewidth=2)
ylims!((0.5,1.6)); xlabel!("Frequency [Hz]"); ylabel!("|p/p0|")
scatter!(k_secondary*340/(2π),0.5ones(length(k_secondary)),label="ωᵢ")
# # Bibliography
# ```@bibliography
# Pages = []
# Canonical = false

# panagiotopoulos2020
# Paltorp2023
# ihlenburg1998a
# panagiotopoulos2022a
# ```
