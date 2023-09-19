# # Lossy Oscilating sphere (Exterior)
# # Importing related packages
using BoundaryIntegralEquations # For BIEs
using IterativeSolvers          # For gmres
using Plots                     # For plottings
# # Setting up constants
frequency = 100.0;                   # Frequency               (Hz)
ρ₀,c,_,_,_,kᵥ,_,_,_,_,_,_ = visco_thermal_constants(;freq=frequency,S=1);
a       = 1.0;                       # Radius of sphere        (m)
k       = 2π*frequency/c;            # Wavenumber              (1/m)
v₀      = 1e-2;                      # Initial velocity        (m/s)
# # Loading the mesh
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
#src mesh_file = joinpath(mesh_path,"sphere_1m_extra_coarse");
mesh_file = joinpath(mesh_path,"sphere_1m_coarser");
#src mesh_file = joinpath(mesh_path,"sphere_1m_coarse");
#src mesh_file = joinpath(mesh_path,"sphere_1m");
#src mesh_file = joinpath(mesh_path,"sphere_1m_fine");
mesh = load3dTriangularComsolMesh(mesh_file)
# For use later we extract the polar and azimuth angles from the collocation nodes
targets = mesh.sources;                         # Target coordinates      [m]
n = size(targets, 2);                           # Number of target points
x = targets[1,:];                               # x-coordinates
y = targets[2,:];                               # y-coordinates
z = targets[3,:];                               # z-coordinates
θ = acos.(z./a);                                # Polar/Colatitude angle
ϕ = acos.(x./sqrt.(x.^2 .+ y.^2)).*sign.(y);    # Azimuth angle
# # Analytical solution
# The following is the analytical solution of a sphere oscillating in a viscous fluid described in section 6.9 of Temkin ([temkin2001elements](@cite)).
# The solution is truncated to only include the first order terms (which are the dominating terms). As such the analytical expression for the acoustical pressure (``p^a``), viscous velocity in the radial direction (``v^r``), and viscous velocity in the polar direction (``v^\theta``) on the surface are given as functions of the polar angle
# ```math
# \begin{aligned}
#   p^a_\text{analytical}(\theta) &\approx -3\rho_0 c A_1h_1^{(1)}(ka)\cos(\theta)\\
#   v^r_\text{analytical}(\theta) &\approx \left(v_0 - 3\mathrm{i}k A_1h_1'(ka)\right)\cos(\theta)\\
#   v^\theta_\text{analytical}(\theta) &\approx \left(-v_0 + \frac{3\mathrm{i}}{a}A_1h_1(ka)\right)\sin(\theta)
# \end{aligned}
# ```
# where ``\rho_0`` is the ambient density of air, ``c`` is the speed of sound, ``k`` is the wavenuumber, ``a`` is the radius of the sphere, ``h_1^{(1)}`` is the spherical Hankel function of order 1, and ``v_0`` is the oscilating velocity in the ``z``-direction. Furthermore ``A_1`` is a constant computed from the radius of the sphere, the wavenumber, and the viscous wavenumber. In order to compute ``A_1`` we start by defining two constants similar to Temkin
β = a*kᵥ; # As defined in (6.9.19) in Temkin (kᵥ is the so-called viscous wavenumber)
b = a*k;  # As defined in (6.9.19) in Temkin
# Then using (6.9.25) from Temkin it is now possible to compute ``A_1`` as
A₁ = -v₀/(3im*k)*b^3*exp(-im*b)*(3β+3im-im*β^2)/(β^2*(b^2-2)-b^2+im*(β*b^2+2b*β^2));
# In order to evaluate the three analytical expressions we need to define the sphereical Hankel function as well as its derivative at order 1
sp_h(z)  = -exp(im*z)*(z + im)/(z^2);          # Spherical Hankel function of order 1
dsp_h(z) = exp(im*z)*(2z + im*(2-z^2))/(z^3);  # Derivative of Spherical Hankel function of order 1
# Using these it is now possible to evaluate the analytical expressions
θ_analytical  = collect(0:0.01:π)
pa_analytical = -3.0*ρ₀*c*k*A₁*(sp_h(k*a))*cos.(θ_analytical);   # Analytical acoustic pressure
vr_analytical = ( v₀ - 3*im*k*A₁*dsp_h(k*a))*cos.(θ_analytical); # Analytical viscous velocity in the radial direction
vθ_analytical = (-v₀ + 3*im/a*A₁*sp_h(k*a))*sin.(θ_analytical);  # Analytical viscous velocity in the polar direction
# # Solution through BEM
# The following is based on the ideas presented in [Preuss2023](@cite). A summary of theory can be found [in documentation.](../theory_lossy.md) In simple terms the problem of including the viscous and thermal losses into the BEM system boils down to solving the following linear system of equations
# ```math
# \begin{equation}
#     \underbrace{\left[\mathbf{G}_a\left(\mu_a\left(\mathbf{R}\mathbf{N}\right)^{-1}\mathbf{R}\mathbf{D}_c + \mu_h\mathbf{G}_h^{-1}\mathbf{H}_h\right) - \phi_a\mathbf{H}_a\right]}_{\text{LGO}}\mathbf{p}_a = \mathbf{G}_a\left(\underbrace{\mathbf{R}\mathbf{N}}_{\text{inner}}\right)^{-1}\mathbf{R} \mathbf{v}_s,
# \end{equation}
# ```
# where ``\mathbf{R} = \mathbf{D}_r - \mathbf{N}^\top\mathbf{G}_v^{-1}\mathbf{H}_v``. The definition of the remaining matrices can be found [here.](../theory_lossy.md) In the code LGO is implemented as a lazy linear operator that contain all matrices that goes into the linear system and evaluates the multiplication as in the algorithm below (note that (4.50) refers to the equation above).
# ![Algorithm](../figures/algorithm_nested.png)
# In the iterative scheme it is possible to apply acceleration techniques such as the fast multipole method (FMM) or hierarchical matrices instead of the dense versions of ``\mathbf{G_a}`` and ``\mathbf{H_a}``. In the following code snippet it was chosen to use the FMM by setting `fmm_on=true`.
LGO = LossyGlobalOuter(mesh,frequency;fmm_on=true,depth=1,n=3,progress=false);
vs  = [zeros(2n); v₀*ones(n)];  # Define surface velocity in the
rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*vs - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*vs));
pa_bem = gmres(LGO,rhs;verbose=true);
# To verify the acoustical pressure we plot the real part against the analytical solution
scatter(θ,real.(pa_bem),label="BEM",markersize=2);
plot!(θ_analytical,real.(pa_analytical),label="Analytical",linewidth=2);
ylabel!("Re(pᵃ)"); title!("Acoustical Pressure"); xlabel!("θ (rad)")
# From the acoustical pressure it possible to reconstruct the viscous velocity as
ph  = -LGO.tau_a/LGO.tau_h*pa_bem;
dpa = -gmres(LGO.Ga,LGO.Ha*pa_bem);
dph = -gmres(LGO.Gh,LGO.Hh*ph);
v   = vs - (LGO.phi_a*LGO.Dc*pa_bem +
            LGO.phi_a*LGO.Nd*dpa    +
            LGO.phi_h*LGO.Dc*ph     +
            LGO.phi_h*LGO.Nd*dph);
# From this we can extract the radial component of the viscous flow as
vr_bem = LGO.Nd'*v;
# Then subtracting the radial part from the viscous flow we the tangential viscous velocity.
vt = v + LGO.Nd*vr_bem;
# Finally, projecting the tangential flow onto the polar direction it follows that
vθ_bem = cos.(θ).*cos.(ϕ).*vt[1:1n] + cos.(θ).*sin.(ϕ).*vt[n+1:2n] - sin.(θ).*vt[2n+1:3n];
# Plotting real part of the radial velocity and polar velocity
p1 = scatter(θ,real.(vr_bem),label="BEM",markersize=2);
plot!(p1,θ_analytical, real.(vr_analytical),label="Analytical",linewidth=2);
ylabel!(p1,"Re(vʳ)"); title!(p1,"Viscous velocity in the radial direction")
p2 = scatter(θ,real.(vθ_bem),markersize=2,legend=false);
plot!(p2,θ_analytical,real.(vθ_analytical),linewidth=2);
xlabel!(p2,"θ (rad)"); ylabel!(p2,"Re(v⁽⁾)"); title!(p2,"Viscous velocity in the polar direction")
plot(p1,p2,layout=(2,1))

# # Bibliography
# ```@bibliography
# Pages = []
# Canonical = false

# temkin2001elements
# Preuss2023
# ```
