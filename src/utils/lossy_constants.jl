function visco_thermal_constants(;freq=1000,S=-1,pa=101325.0,t=20.0,Hr=50.0)

    rho,c,cf,cpcv,mu,alfa,Cp,Cv,landa,beta = ambient_to_properties(freq,pa,t,Hr)
    eta = 0.6*mu

    # Characteristic Lengths
    lh   = landa/(rho*Cp*c)
    lv   = (eta + 4.0*mu/3.0)/(rho*c)
    lvp  = mu/(rho*c)
    lvh  = lv + (cpcv-1.0)*lh
    lvhp = (cpcv-1.0)*(lh - lv)

    # Wavenumbers
    kp  = 2*pi*freq/c
    ka2 = kp.^2 ./ (1.0 .- S*im*kp*lvh - kp.^2*lh*lvhp)
    kh2 = S*im*kp ./ (lh * (1.0 .+ S*im*kp*lvhp))
    kv2 = S*im*kp/lvp

    # The result of the square 
    ka = sqrt(ka2)
    kh = sqrt(kh2)
    kv = sqrt(kv2)

    # Thermal boundary condition constants
    tau_a = (cpcv - 1.0) ./ (beta*cpcv*(1.0 + S*im*lh*ka2/kp))
    tau_h = (cpcv - 1.0) ./ (beta*cpcv*(1.0 + S*im*lh*kh2/kp))

    # Velocity boundary condition constants
    phi_a = -S*im/(rho*2.0*pi*freq * (1.0 + S*im*lv*ka2/kp))
    phi_h = -S*im/(rho*2.0*pi*freq * (1.0 + S*im*lv*kh2/kp))

    return rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu

end
function compute_k1_k2(freq,S=-1)
    rho,c,cf,gamma,nu,alfa,Cp,Cv,lambda,beta = ambient_to_properties(freq)
    ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=S)
    k  = 2.0*π*freq/c                               # Wavenumber
    ω  = 2.0*π*freq                                 # Angular frequency
    dV = sqrt(2.0*μ/(rho*ω))                        # Viscous BL-thickness
    dT = sqrt(2*lambda/(ω*rho*Cp))                  # Thermal BL-thickness
    k1 = -dV*(im - 1.0)/2.0                         # Constant 1
    k2 = +dT*k^2*(im - 1.0)/2.0*(gamma - 1.0)       # Constant 2
    return k1,k2
end
