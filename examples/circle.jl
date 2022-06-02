using IntegralEquations,Plots,LinearAlgebra

n_elements = 1000
freq  = 1000.0                                  # Frequency                 [Hz]
c     = 340.0                                   # Speed of sound            [m/s]
k     = 2*π*freq/c                              # Wavenumber                [1/m]
r = 1.0

# mesh = IntegralEquations.mesh_circle(ContinuousCurveLinear(3),n_elements;radius=r)
mesh = IntegralEquations.mesh_circle(ContinuousCurveQuadratic(3),n_elements;radius=r)
# mesh = IntegralEquations.mesh_circle(ContinuousCurveLinear(3),DiscontinuousCurveConstant(3),n_elements;radius=r)
# mesh = IntegralEquations.mesh_circle(ContinuousCurveQuadratic(3),DiscontinuousCurveLinear(3),n_elements;radius=r)
# mesh = IntegralEquations.mesh_circle(ContinuousCurveQuadratic(3),DiscontinuousCurveQuadratic(3),n_elements;radius=r)
# @time F,G,C=assemble_parallel!(mesh,k,mesh.sources;n=4,fOn=false,gOn=false,cOn=false);
@time F,G,C=assemble_parallel!(mesh,k,mesh.sources;n=4,gOn=false);

src = mesh.sources
pI = exp.(im*k*src[1,:])
ϕ = angle.(src[1,:] + im*src[2,:])
r = 1.0
p = IntegralEquations.plane_wave_scattering_circle(ϕ,k*r,150)

pB = (F + Diagonal(C .- 1.0))\pI

plot(ϕ,abs.(p),label="Analytical")
plot!(ϕ,abs.(pB),label="BEM")
