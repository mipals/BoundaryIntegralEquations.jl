#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using SpecialFunctions, LinearAlgebra, IntegralEquations, Plots
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq  = 50.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=-1)
k     = 2*π*freq/c                              # Wavenumber                [1/m]
#==========================================================================================
                                Setting up constants
==========================================================================================#
alog = -8
blog = -1
rlog = collect(10.0.^range(alog,blog,length=100))

f(x) = exp(-im*kh*x)/(4π*x)
plot(rlog,real.(f.(rlog)),xaxis=:log)
plot(rlog,imag.(f.(rlog)),xaxis=:log)
