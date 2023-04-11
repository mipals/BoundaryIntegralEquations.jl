#==========================================================================================
                                Adapted from OpenBEM
==========================================================================================#
function incoming_wave(angles,P0,xyzb,k)
    # pI = zeros(size((xyzb)))
    A  = sin.(angles[1]) .* cos.(angles[2])
    B  = sin.(angles[1]) .* sin.(angles[2])
    C  = cos.(angles[1])

    dps = (A .* xyzb[1,:] + B .* xyzb[2,:] + C .* xyzb[3,:]) ./ sqrt.(A.^2 + B.^2 + C.^2)
    return P0 .* exp.(-im .* (k .* dps))
end
"""
    plane_wave_scattering_sphere(k,a,r,θ,Ieps)

```math
  p_s(r, \\theta) = -P_0\\sum_{n}^\\infty \\mathrm{i}^n(2n+1)\\frac{j_n^{'}(ka)}{h_n^{'}(ka)}P_n(\\cos(\\theta))h_n(kr)
```
"""
function plane_wave_scattering_sphere(k,a,r,θ,Ieps)

    ka = k*a
    kr = k*r
    n  = 0
    ptot = zeros(length(θ),length(r))
    accuracy = []
    rel_accuracy = Ieps*10.0
    while rel_accuracy > Ieps
        apri = (n*besselj(n - 0.5,ka) - (n + 1)*besselj(n + 1.5,ka)) ./
               (n*besselh(n - 0.5,ka) - (n + 1)*besselh(n + 1.5,ka));
        # Pmatrix = legendre(n,cos(theta));
        # Pmatrix[1,:]; # vector of legendre_n function of different arguments
        P = Pl.(cos.(θ),n)
        factor = (-im)^n*(2*n+1)*sqrt(π/2.0./kr);
        jn = besselj.(n + 0.5,kr);
        hn = besselh.(n + 0.5,kr);
        iter = (factor.*(jn - apri .* hn))*P;
        ptot = ptot+iter;
        rel_accuracy = minimum(maximum(abs.(iter ./ ptot)));
        accuracy = [accuracy; rel_accuracy];
        n=n+1;
    end
    return ptot[:], accuracy
end
"""
```math
  p_\\text{int}(r) = -\\mathrm{i}Z_0v_r(\\frac{a}{r})\\frac{ka\\sin(kr)}{ka\\cos(ka) - \\sin(ka)}
 ```
"""
function p_interior_pulsating_sphere()

end
"""
```math
    p_\\text{int}(x,y,z) = -3v_z z \\mathrm{i}kZ_0\\frac{j_1(kr)}{kr(j_0(ka) - 2j_2(ka))},
```
"""
function p_interior_oscilating_sphere()

end

pinc(r,phi,k;p0=1.0) = p0*exp.(im*k*r.*cos.(phi))
# Short-cuts for spherical bessel j and y functions
sp_j(n,z)  = sphericalbesselj.(n,z)
dsp_j(n,z) = (sphericalbesselj.(n-1,z) - sphericalbesselj.(n+1,z))/2 - sphericalbesselj.(n,z)./(2z)
sp_y(n,z)  = sphericalbessely.(n,z)
dsp_y(n,z) = (sphericalbessely.(n-1,z) - sphericalbessely.(n+1,z))/2 - sphericalbessely.(n,z)./(2z)
# Defining spherical hankel function of the first kind
sp_h(n,z)  = sp_j.(n,z) + im*sp_y.(n,z)
dsp_h(n,z) = dsp_j.(n,z) + im*dsp_y.(n,z)
# Scattered wave - (Equation 130 in Marburg)
function psca(r,phi,k;N=50,R=1.0,p0=1.0)
    # Args:
    #    r: Distance at evaluation point
    #    phi: Angle to evluation point
    #    k: Wavenumber
    # Kwargs:
    #    N: Number of terms in sum (Standard )
    #    R: Radius of sphere
    #    p0: Sound pressure amplitude of incident wave
    z = k*r
    Z = k*R
    psca = (-sp_j.(0+1,Z)./(-sp_j.(0+1,Z)-im*sp_y.(0+1,Z))) .*Pl.(cos.(phi), 0) .* sp_h.(0,z)
    for n = 1:N-1
        psca += im^n*(2.0*n + 1.0)*(dsp_j.(n,Z)./dsp_h.(n,Z)) .*Pl.(cos.(phi), n) .* sp_h.(n,z)
    end
    # F. Ihlenburg, FEA of Acoustic Scattering, Applied Mathematical Sciences, (2.1.18)
    return -p0*psca
end
