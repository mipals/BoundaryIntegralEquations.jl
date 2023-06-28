# # Lossy sphere (3D - Exterior)
# # Importing related packages
using LinearAlgebra, BoundaryIntegralEquations, IterativeSolvers, Plots
# # Loading the mesh
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
#src mesh_file = joinpath(mesh_path,"sphere_1m_extra_coarse");
# mesh_file = joinpath(mesh_path,"sphere_1m_coarser");
# mesh_file = joinpath(mesh_path,"sphere_1m_coarse");
# mesh_file = joinpath(mesh_path,"sphere_1m");
mesh_file = joinpath(mesh_path,"sphere_1m_fine");
# mesh_file = joinpath(mesh_path,"sphere_1m_finer");
mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
# mesh_file = joinpath(mesh_path,"sphere_1m_finest");
# mesh_file = joinpath(mesh_path,"sphere_1m_35k");
# mesh_file = joinpath(mesh_path,"sphere_1m_77k");
mesh = load3dTriangularComsolMesh(mesh_file)
# Mesh with 323.7k DOFs - 647.3k elements
# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/test_binary.ply")
# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/F1.stl")
# mesh = load3dTriangularMesh("/Users/mpasc/Documents/testfiles/quad_sphere_text.stl")

BoundaryIntegralEquations.get_hmin(mesh;n=10)
BoundaryIntegralEquations.get_hmax(mesh;n=10)
# # Setting up constants
frequency = 1000.0;                   # Frequency               [Hz]
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
thres = 1e-4
@time LGO_FMM = LossyGlobalOuter(mesh,frequency;thres=thres,fmm_on=true,depth=1,n=3,progress=true);
@time LGO_H   = LossyGlobalOuter(mesh,frequency;thres=thres,hmatrix_on=true,depth=1,n=3,progress=true);
println(prod(size(LGO_H))*2*8/(2^30))
varinfo(r"LGO")

v0  = [zeros(2n); v₀*ones(n) ];
rhs = LGO_FMM.Ga*gmres(LGO_FMM.inner,LGO_FMM.Dr*v0 - LGO_FMM.Nd'*gmres(LGO_FMM.Gv,LGO_FMM.Hv*v0));

@time pa_fmm = gmres(LGO_FMM,rhs;verbose=true);
@time pa_h   = gmres(LGO_H,rhs;verbose=true);

time_h = @elapsed pa_h,hist_h = gmres(LGO_H,rhs;verbose=true,log=true);

# Finally we plot the results
θ = acos.(targets[3,:]/a);    # Angles to target points [rad]
scatter(θ,real.(pa_h),label="BEM-HMatrix",markersize=3);
scatter!(θ,real.(pa_fmm),label="BEM-FMM",markersize=3);
plot!(θ_analytical,real.(p_analytical),label="Analytical",linewidth=2);
ylabel!("Re(p)"); title!("Frequency = $(frequency)"); xlabel!("θ (rad)")
# # Computing the remaining variables
# First we compute the remaining variables
@time begin
    ph_fmm  = -LGO_FMM.tau_a/LGO_FMM.tau_h*pa_fmm;
    dpa_fmm = -gmres(LGO_FMM.Ga,LGO_FMM.Ha*pa_fmm);
    dph_fmm = -gmres(LGO_FMM.Gh,LGO_FMM.Hh*ph_fmm);
    v_fmm   = v0 - (LGO_FMM.phi_a*LGO_FMM.Dc*pa_fmm  +
                LGO_FMM.phi_a*LGO_FMM.Nd*dpa_fmm +
                LGO_FMM.phi_h*LGO_FMM.Dc*ph_fmm  +
                LGO_FMM.phi_h*LGO_FMM.Nd*dph_fmm);
    dvn_fmm = -gmres(LGO_FMM.Gv,LGO_FMM.Hv*v_fmm);
    # Normal compontent of the viscous flow
    v_n0_fmm = LGO_FMM.Nd'*v_fmm;
    # Computing the tangential velocity by substracting the normal information
    v_t_fmm = v_fmm + LGO_FMM.Nd*v_n0_fmm;
    vt_sum_fmm = sqrt.(v_t_fmm[0n+1:1n].^2 + v_t_fmm[1n+1:2n].^2 + v_t_fmm[2n+1:3n].^2);
end;
# H-matrix
@time begin
    ph_h  = -LGO_H.tau_a/LGO_H.tau_h*pa_h;
    dpa_h = -gmres(LGO_H.Ga,LGO_H.Ha*pa_h);
    dph_h = -gmres(LGO_H.Gh,LGO_H.Hh*ph_h);
    v_h   = v0 - (LGO_H.phi_a*LGO_H.Dc*pa_h  +
                LGO_H.phi_a*LGO_H.Nd*dpa_h +
                LGO_H.phi_h*LGO_H.Dc*ph_h  +
                LGO_H.phi_h*LGO_H.Nd*dph_h);
    dvn_h = -gmres(LGO_H.Gv,LGO_H.Hv*v_h);
    # Normal compontent of the viscous flow
    v_n0_h = LGO_H.Nd'*v_h;
    # Computing the tangential velocity by substracting the normal information
    v_t_h = v_h + LGO_FMM.Nd*v_n0_h;
    vt_sum_h = sqrt.(v_t_h[0n+1:1n].^2 + v_t_h[1n+1:2n].^2 + v_t_h[2n+1:3n].^2);
end;
# Plotting real part of
scatter(θ,abs.(pa_fmm),label="BEM-FMM",markersize=3);
scatter!(θ,abs.(pa_h),label="BEM-HMatrix",markersize=3);
plot!(θ_analytical,abs.(p_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|p|"); title!("Frequency = $(frequency)")
# Plotting normal velocity
scatter(θ,abs.(v_n0_fmm),label="BEM-FMM",markersize=3)
scatter!(θ,abs.(v_n0_h),label="BEM-HMatrix",markersize=3)
plot!(θ_analytical,abs.(vn_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|Vₙ|"); title!("Frequency = $(frequency)")
# Plotting tangential velocity
scatter(θ,abs.(vt_sum_fmm),label="BEM-FMM",markersize=3);
scatter!(θ,abs.(vt_sum_fmm),label="BEM-HMatrix",markersize=3);
plot!(θ_analytical,abs.(vt_analytical),label="Analytical",linewidth=2);
xlabel!("θ (rad)"); ylabel!("|Vₜ|"); title!("Frequency = $(frequency)")
