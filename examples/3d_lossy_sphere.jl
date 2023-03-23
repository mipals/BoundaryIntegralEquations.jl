# # Lossy sphere (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations, IterativeSolvers, Plots, MeshViz
# # Loading the mesh
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
#src mesh_file = joinpath(mesh_path,"sphere_1m_extra_coarse");
mesh_file = joinpath(mesh_path,"sphere_1m_coarser");
#src mesh_file = joinpath(mesh_path,"sphere_1m_coarse");
#src mesh_file = joinpath(mesh_path,"sphere_1m");
#src mesh_file = joinpath(mesh_path,"sphere_1m_fine");
mesh = load3dTriangularComsolMesh(mesh_file)
# # Setting up constants
frequency = 100.0;                    # Frequency               [Hz]
ρ₀,c,_,_,_,kᵥ,_,_,_,_,_,_ = visco_thermal_constants(;freq=frequency,S=1);
a       = 1.0;                       # Radius of sphere        [m]
ω       = 2π*frequency;              # Angular frequency       [rad/s]
k       = 2π*frequency/c;            # Wavenumber              [1/m]
v₀      = 1e-2;                      # Initial velocity        [m/s]
normals = mesh.normals;              # Normal coordinates      [m]
targets = mesh.sources;              # Source coordinates      [m]
n       = size(targets, 2);          # Number of target points
# # Analytical solution
# The following is the analytical solution of a sphere oscillating in a viscous fluid described in section 6.9 of S. Temkin, 2001.
# The solution is truncated to only include the first order terms (which are the dominating terms). Here the analytical solution for the acoustical pressure, normal velocity and tangential velocity described by
# ```math
# \begin{aligned}
#   p_\text{analytical}(r,\theta) &\approx -3\rho c A_1h_1^{(1)}(kr)\cos(\theta)\\
#   v^n_\text{analytical}(r,\theta) &\approx -\frac{6\mathrm{i}kB_1}{k_vr}h_1(k_vr)\cos(\theta)\\
#   v^t_\text{analytical}(r,\theta) &\approx -\frac{3\mathrm{i}B_1\left(k_vrh_1'(k_vr) + h_1(k_vr)\right)}{r}\sin(\theta)
# \end{aligned}
# ```
# where ``A_1`` and ``B_1`` is found by solving the following system of equations
# ```math
#   \begin{bmatrix} v_{0}/(3\mathrm{i}k) \\ av_{0}/(3\mathrm{i})\end{bmatrix} =
#   \begin{bmatrix}
#       h_1'(ka)  & -2h_1(k_va)/(ka)\\
#       h_1(ka)   & -\left(k_vah_1'(k_va) + h_1(k_va)\right)
#   \end{bmatrix}
#   \begin{bmatrix} A_1 \\ B_1\end{bmatrix}
# ```
# Solving the above naively will be numerically unstable, due the large values of ``k_v``. Instead one can find ``A_1`` using (6.9.25) and then solve for the two expressions including ``B_1`` using the two equations above. Following the notation of Temkin we introduce the following
β = kᵥ*a; # As defined in (6.9.19) in Temkin
b = k*a;  # As defined in (6.9.19) in Temkin
# Now computing ``A_1``using (6.9.25)
A₁ = -v₀/(3im*k)*b^3*exp(-im*b)*(3β+3im-im*β^2)/(β^2*(b^2-2)-b^2+im*(β*b^2+2b*β^2));
# Using the two equations described earlier we can compute the expressions containing ``B_1`` as follows
h1_ka  = exp(im*b)/(b^2)*(-b - im);           # Spherical hankel function
dh1_ka = exp(im*b)*(2b + im*(2 - b^2))/(b^3); # Derivative of spherical hankel function
B1h1   = -(v₀/(3im*k) - A₁*dh1_ka)/2;         # Top equation in the system shown earlier
B1dh1  = -(v₀/(3im)*a - A₁*h1_ka);            # Bottom equation in the system shown earlier
# Using the above the analytical expression can be evaluted on the surface
θ_analytical  = collect(0:0.01:π)
p_analytical  = -3.0*ρ₀*c*k*A₁*(h1_ka).*(cos.(θ_analytical));
vn_analytical = -6im*k*B1h1*cos.(θ_analytical);  # Analytical normal velocity
vt_analytical = -3im/a*B1dh1*sin.(θ_analytical); # Analytical tangential velocity
# # Iterative Solution through BEM
LGO = LossyGlobalOuter(mesh,frequency;fmm_on=true,depth=1,n=3,progress=false);
v0  = [v₀*ones(n); zeros(2n)];
rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*v0 - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*v0));
pa  = gmres(LGO,rhs;verbose=true);
# Finally we plot the results
θ = acos.(targets[1,:]/a);    # Angles to target points [rad]
scatter(θ,real.(pa),label="BEM-global",markersize=3);
plot!(θ_analytical,real.(p_analytical),label="Analytical",linewidth=2);
ylabel!("Re(p)"); title!("Frequency = $(frequency)"); xlabel!("θ (rad)")
# # Computing the remaining variables
# First we compute the remaining variables
ph  = -LGO.tau_a/LGO.tau_h*pa;
dpa = -gmres(LGO.Ga,LGO.Ha*pa);
dph = -gmres(LGO.Gh,LGO.Hh*ph);
v   = v0 - (LGO.phi_a*LGO.Dc*pa  +
            LGO.phi_a*LGO.Nd*dpa +
            LGO.phi_h*LGO.Dc*ph  +
            LGO.phi_h*LGO.Nd*dph);
dvn = -gmres(LGO.Gv,LGO.Hv*v);
# Normal compontent of the viscous flow
v_n0 = LGO.Nd'*v;
# Computing the tangential velocity by substracting the normal information
v_t = v + LGO.Nd*v_n0;
vt_sum = sqrt.(v_t[0n+1:1n].^2 + v_t[1n+1:2n].^2 + v_t[2n+1:3n].^2);
# Plotting real part of
scatter(θ,abs.(pa),label="BEM",markersize=3);
plot!(θ_analytical,abs.(p_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|p|"); title!("Frequency = $(frequency)")
# Plotting normal velocity
scatter(θ,abs.(v_n0),label="BEM",markersize=3);
plot!(θ_analytical,abs.(vn_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|Vₙ|"); title!("Frequency = $(frequency)")
# Plotting tangential velocity
scatter(θ,abs.(vt_sum),label="BEM",markersize=3);
plot!(θ_analytical,abs.(vt_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|Vₜ|"); title!("Frequency = $(frequency)")
