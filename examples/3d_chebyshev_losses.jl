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
# mesh_file = joinpath(mesh_path,"sphere_1m_coarser");
#src mesh_file = joinpath(mesh_path,"sphere_1m_coarse");
#src mesh_file = joinpath(mesh_path,"sphere_1m");
# mesh_file = joinpath(mesh_path,"sphere_1m_finer");
mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
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

using ProgressMeter
import BoundaryIntegralEquations: arnoldi_basis
function lossy_krylov_basis(mesh,freqlist;eps=1-4,n_gauss=3,verbose=true,v₀=1e-2,progress=true)
    n_sources = size(mesh.sources,2)
    V  = zeros(ComplexF64,n_sources, 0)         # Preallocating the total Krylov system
    nK = length(freqlist)
    solutions = zeros(ComplexF64,n_sources,nK)  # list of solutions
    qlist = zeros(Int64,length(freqlist))
    if progress == true
        prog = Progress(nK, 0.2, "Assembling Krylov vectors:\t", 50)
    end
    for i = 0:nK-1
        frequency = freqlist[i+1]

        vs  = [zeros(2n); v₀*ones(n)];  # Define surface velocity in the
        LGO = LossyGlobalOuter(mesh,frequency;fmm_on=true,depth=1,n=n_gauss);
        rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*vs - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*vs));

        pa_fmm,history = gmres(LGO,rhs;verbose=verbose,log=true);
        solutions[:,i+1] = pa_fmm
        V = [V arnoldi_basis(LGO,rhs,history.iters)]
        qlist[i+1] = history.iters
        if progress == true
            next!(prog)
        end
    end
    U,S,_ = svd(V)                              # Extracting projection matrix
    idx = S .> eps  # Only use projection direcitons with singular values larger than `eps`
    return U[:,idx], solutions, qlist
end
function create_chebyshev_basis(mesh,freq_secondary,U;n_gauss=3,v₀=1e-2,progress=true)
    n = size(U,1)
    vs  = [zeros(2n); v₀*ones(n)];  # Define surface velocity in the
    temp = zeros(eltype(U), size(U))
    A_output = zeros(eltype(U),size(U,2),size(U,2),length(freq_secondary))
    b_output = zeros(eltype(U),size(U,2),length(freq_secondary))
    if progress == true
        prog = Progress(length(freq_secondary), 0.2, "Assembling Chebyshev basis:\t", 50)
    end
    for (index,frequency) in enumerate(freq_secondary)
        LGO = LossyGlobalOuter(mesh,frequency;fmm_on=true,depth=1,n=n_gauss);
        for (idx,col) in enumerate(eachcol(U))
            temp[:,idx] = LGO*col
        end
        A_output[:,:,index] = U'*temp
        rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*vs - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*vs));
        b_output[:,index] = U'*rhs
        if progress == true
            next!(prog)
        end
    end
    return A_output, b_output
end
function chebyshev_matrix_coefficients(A_reduced)
    M = size(A_reduced,3)
    Cj = similar(A_reduced)
    for j = 0:M-1
        Cj[:,:,j+1] = 2/M*sum(i -> A_reduced[:,:,i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return Cj
end
function chebyshev_vector_coefficients(b_reduced)
    M = size(b_reduced,2)
    bj = similar(b_reduced)
    for j = 0:M-1
        bj[:,j+1] = 2/M*sum(i -> b_reduced[:,i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return bj
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
function eval_chebyshev_vector(bj,coeffs)
    output = -bj[:,1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*bj[:,index]
    end
    return output
end
function assemble_chebyshev(Cj,bj,k)
    A = eval_chebyshev_matrix(Cj,chebyshev_eval(k,size(Cj,3)))
    b = eval_chebyshev_vector(bj,chebyshev_eval(k,size(bj,2)))
    return A,b
end

freqlist = [150.0 250.0 350.0]
U, sols, qlist = lossy_krylov_basis(mesh,freqlist)
@info "Reduced basis size: $(size(U,2)) | Reduction in DOF: $(1 - size(U,2)/size(U,1)) %";
# Then we define the secondary wavenumbers, i.e. the wavenumbers for which we compute the matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega)\right)\mathbf{U}``
fmin = 100;
fmax = 400;
g(k)    = 2/(fmax - fmin)*k .- (fmax + fmin )/(fmax - fmin);
ginv(ω) = (fmax - fmin)/2*ω .+ (fmin + fmax)/2;
M  = 20;                                # Number of terms in the Chebyshev approximation
ωᵢ = cos.(π*(collect(0:M-1) .+ 1/2)/M); # Zeros of the Chebyshev polynomials
f_secondary = ginv.(ωᵢ);                # Mapping [-1,1] to [kmin, kmax]
# Using the defined functions we can create the basis matrices ``\mathbf{U}^\text{H}\mathbf{A}\left(g^{-1}(\omega_i)\right)\mathbf{U}``
A_reduced,b_reduced = create_chebyshev_basis(mesh,f_secondary,U);
# Now from the basis matrices we can compute the Chebyshev coefficient matrices (``\mathbf{C}_j``) as
Cj = chebyshev_matrix_coefficients(A_reduced);
bj = chebyshev_vector_coefficients(b_reduced);


p1 = 708
p2 = 100
θ_analytical = [θ[p1]; θ[p2]]

frequencies = collect(fmin:10:fmax);
p1_chebyshev = zeros(ComplexF64,length(frequencies));
p2_chebyshev = zeros(ComplexF64,length(frequencies));
prog = Progress(length(frequencies), 0.2, "Frequency sweep:\t", 50)
for (i,k) in enumerate(frequencies)
    k_scaled = g(k);                       # Scaling the wavenumber to the interval [-1,1]
    A,b = assemble_chebyshev(Cj,bj,k_scaled)
    p_romi  = U*(A\b);
    p1_chebyshev[i] = p_romi[p1];          # Saving pressure at field point 1
    p2_chebyshev[i] = p_romi[p2];          # Saving pressure at field point 2
    next!(prog)
end

fanalytical = collect(fmin:1:fmax);
pa_analytical = zeros(ComplexF64,length(fanalytical),length(θ_analytical));
for (idx,frequency) in enumerate(fanalytical)
    ρ₀,c,_,_,_,kᵥ,_,_,_,_,_,_ = visco_thermal_constants(;freq=frequency,S=1);
    a       = 1.0;                       # Radius of sphere        (m)
    k       = 2π*frequency/c;            # Wavenumber              (1/m)
    β = a*kᵥ; # As defined in (6.9.19) in Temkin (kᵥ is the so-called viscous wavenumber)
    b = a*k;  # As defined in (6.9.19) in Temkin
    # Then using (6.9.25) from Temkin it is now possible to compute ``A_1`` as
    A₁ = -v₀/(3im*k)*b^3*exp(-im*b)*(3β+3im-im*β^2)/(β^2*(b^2-2)-b^2+im*(β*b^2+2b*β^2));
    # In order to evaluate the three analytical expressions we need to define the sphereical Hankel function as well as its derivative at order 1
    sp_h(z)  = -exp(im*z)*(z + im)/(z^2);          # Spherical Hankel function of order 1
    dsp_h(z) = exp(im*z)*(2z + im*(2-z^2))/(z^3);  # Derivative of Spherical Hankel function of order 1
    # Using these it is now possible to evaluate the analytical expressions
    pa_analytical[idx,:] = -3.0*ρ₀*c*k*A₁*(sp_h(k*a))*cos.(θ_analytical);   # Analytical acoustic pressure
end

# To verify the acoustical pressure we plot the real part against the analytical solution
plot(fanalytical,real.(pa_analytical[:,1]),label="Analytical",markersize=2)
plot!(fanalytical,real.(pa_analytical[:,2]),label="Analytical",markersize=2)
scatter!(frequencies,real.(p1_chebyshev),label="BEM")
scatter!(frequencies,real.(p2_chebyshev),label="BEM")




function create_chebyshev_vectors(mesh,freq_secondary;n_gauss=3,v₀=1e-2,progress=true)
    n = size(mesh.sources,2)
    vs  = [zeros(2n); v₀*ones(n)];  # Define surface velocity in the
    p_output = zeros(eltype(U),n,length(freq_secondary))
    if progress == true
        prog = Progress(length(freq_secondary), 0.2, "Assembling Chebyshev basis:\t", 50)
    end
    for (index,frequency) in enumerate(freq_secondary)
        LGO = LossyGlobalOuter(mesh,frequency;fmm_on=true,depth=1,n=n_gauss);
        rhs = LGO.Ga*gmres(LGO.inner,LGO.Dr*vs - LGO.Nd'*gmres(LGO.Gv,LGO.Hv*vs));
        p_output[:,index] = gmres(LGO,rhs)
        if progress == true
            next!(prog)
        end
    end
    return p_output
end

p_output = create_chebyshev_vectors(mesh,f_secondary)
p_coeff = chebyshev_vector_coefficients(p_output);

frequencies_p = collect(fmin:1:fmax);
p1_chebyshev_p = zeros(ComplexF64,length(frequencies_p));
p2_chebyshev_p = zeros(ComplexF64,length(frequencies_p));
prog = Progress(length(frequencies_p), 0.2, "Frequency sweep:\t", 50)
for (i,k) in enumerate(frequencies_p)
    k_scaled = g(k);                       # Scaling the wavenumber to the interval [-1,1]
    p_romi = eval_chebyshev_vector(p_coeff,chebyshev_eval(k_scaled,size(p_coeff,2)))
    p1_chebyshev_p[i] = p_romi[p1];          # Saving pressure at field point 1
    p2_chebyshev_p[i] = p_romi[p2];          # Saving pressure at field point 2
    next!(prog)
end

plot(fanalytical,real.(pa_analytical[:,1]),label="Analytical",markersize=2)
plot!(fanalytical,real.(pa_analytical[:,2]),label="Analytical",markersize=2)
scatter!(frequencies,real.(p1_chebyshev),label="ROM")
scatter!(frequencies,real.(p2_chebyshev),label="ROM")
scatter!(frequencies_p,real.(p1_chebyshev_p),label="Interpolation")
plot!(frequencies_p,real.(p1_chebyshev_p),label="Interpolation")



function chebyshev_scalar_coefficients(b_reduced)
    M = length(b_reduced)
    bj = similar(b_reduced)
    for j = 0:M-1
        bj[j+1] = 2/M*sum(i -> b_reduced[i+1]*cos(π*j*(i+1/2)/M),0:M-1)
    end
    return bj
end
function eval_chebyshev_scalar(bj,coeffs)
    output = -bj[1]/2
    for index = eachindex(coeffs)
        output += coeffs[index]*bj[index]
    end
    return output
end

# function f(x)
#     mid = (fmax-fmin)/2 + fmin
#     if x < mid
#         return x/((fmax-fmin)/2 + fmin - x)
#     else
#         return -x/((fmax-fmin)/2 + fmin - x) + x/(fmax - x)
#     end
# end

f(x) = sin(x/5)

plot(frequencies_p, f.(frequencies_p))



M  = 40;                                # Number of terms in the Chebyshev approximation
omg = cos.(π*(collect(0:M-1) .+ 1/2)/M); # Zeros of the Chebyshev polynomials
fs = ginv.(omg);                # Mapping [-1,1] to [kmin, kmax]
y = f.(fs)
yt = chebyshev_scalar_coefficients(y)
coeff = chebyshev_eval(-1,length(yt))
eval_chebyshev_scalar(yt,coeff)
p_scalar = zeros(ComplexF64,length(frequencies_p))
for (i,k) in enumerate(frequencies_p)
    k_scaled = g(k);
    p_scalar[i] = eval_chebyshev_scalar(yt,chebyshev_eval(k_scaled,length(coeff)))
end

f(100)
plot(frequencies_p,f.(frequencies_p))
plot!(frequencies_p,real.(p_scalar))

plot(frequencies_p,abs.(f.(frequencies_p) - p_scalar)./f.(frequencies_p) ,yaxis=:log)

# # Bibliography
# ```@bibliography
# Pages = []
# Canonical = false

# temkin2001elements
# Preuss2023
# ```
