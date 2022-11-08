#==========================================================================================
                                    Packages
==========================================================================================#
using FastGaussQuadrature, LinearAlgebra, SpecialFunctions, ForwardDiff, Plots
#==========================================================================================
                                Helper Functions
==========================================================================================#
JacMul!(j,w)    = @inbounds for i = eachindex(j) j[i]   = j[i]*w[i]                      end
dist!(x,y,r)    = @inbounds for i = eachindex(r) r[i]   = hypot(x[1,i]-y[1],x[2,i]-y[2]) end
jacobian!(j,n)  = @inbounds for i = eachindex(j) j[i]   = hypot(n[1,i],      n[2,i])     end
normalize!(n,j) = @inbounds for i = eachindex(j) n[1,i] = n[1,i]/j[i];n[2,i]=n[2,i]/j[i] end
normals!(n,dX)  = @inbounds for i = axes(n,2)    n[1,i] = -dX[2,i];   n[2,i]=dX[1,i]     end
#==========================================================================================
                                Basis Functions
==========================================================================================#
linear(ξ)    = [(1 .- ξ)/2; (1 .+ ξ)/2]
quadratic(ξ) = [ξ .* (ξ .- 1)/2; 1 .- ξ .^2; ξ .* (ξ .+ 1)/2]
constant(ξ)  = ones(eltype(ξ),1,length(ξ))
β = gausslegendre(2)[1][end]; discLinear(ξ)    = linear(ξ/β)
δ = gausslegendre(3)[1][end]; discQuadratic(ξ) = quadratic(ξ/δ)
# Derivatives this is fine for smaller problems. Inefficient for larger computations.
basisFunctionDerivative(func,ξ) = hcat(ForwardDiff.derivative.(func,ξ)...)
#==========================================================================================
                                    Kernels
==========================================================================================#
# Green's function
G!(k,r,int) = @inbounds for i = eachindex(r) int[i] = im/4*hankelh1(0,k*r[i]) end
# Normal derivative of Green's function
function F!(x,y,k,n,r,int)
    @inbounds for i = eachindex(r)
        int[i] = -k*im*hankelh1(1,k*r[i])*(n[1,i]*(x[1,i] - y[1]) +
                                           n[2,i]*(x[2,i] - y[2]))/(4r[i])
    end
    return int
end
# C-constant
function C!(x,y,n,r,int)
    @inbounds for i = eachindex(r)
        int[i] = (n[1,i]*(x[1,i] - y[1]) + n[2,i]*(x[2,i] - y[2]))/(2π*r[i]^2)
    end
    return int
end
#==========================================================================================
                            Mesh interpolation routine
==========================================================================================#
function meshInterpolation!(Fslice,Gslice,Cslice,elementInterpolation,elementDerivatives,
                            physicsInterpolation,weights,
                            coordinates,source,k,r,interpolation,dX,normals,jacobian,integrand)
    # Element/Geometric interpolation
    mul!(interpolation,coordinates,elementInterpolation)    # Interpolating on element
    mul!(dX,coordinates,elementDerivatives)                 # Computing tangents
    dist!(interpolation,source,r)                           # Computing distances
    normals!(normals,dX)                                    # Computing normal direction
    jacobian!(jacobian,normals)                             # Computing Jacobian
    normalize!(normals,jacobian)                            # Normalizing normal
    JacMul!(jacobian,weights)                               # NOW JACOBIAN IS JAC*WEIGHTS
    # Physics Interpolation
    G!(k,r,integrand)                                       # Evaluate Greens function
    JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
    mul!(Gslice,integrand,physicsInterpolation',true,true)  # Integration with basis func
    F!(interpolation,source,k,normals,r,integrand)          # Evaluate ∂ₙGreens function
    JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
    mul!(Fslice,integrand,physicsInterpolation',true,true)  # Integration with basis func
    C!(interpolation,source,normals,r,integrand)            # Evaluate G₀
    Cslice[1] += dot(integrand,jacobian)                    # Integration of G₀
end
#==========================================================================================
                                   Assembly
==========================================================================================#
function assemble(coordinates,topology,k,n;interior=false,order=-1,fieldPoints=[])
    # Defining correct interpolation of geometry and physics (this is a mess. Sorry)
    physicsTopology,sources  = (order ∉ [0,1,2] ? (topology,coordinates) : createPhysics(order,topology,coordinates))
    fieldPoints          = (isempty(fieldPoints) ? sources : fieldPoints)
    basisFunction(ξ)     = (size(topology,1) == 2 ? linear(ξ) : quadratic(ξ))
    physicsFunction(ξ)   = (order ∉ [0,1,2] ? basisFunction(ξ) : (order == 0 ? discConstant(ξ) : (order == 1 ? discLinear(ξ) : discQuadratic(ξ))))
    # Quadrature points
    nodes,weights        = gausslegendre(n)
    # Element/geometry interpolation
    elementInterpolation = basisFunction(nodes')
    elementDerivatives   = basisFunctionDerivative(basisFunction,nodes')
    # Physics interpolation
    physicsInterpolation = physicsFunction(nodes')
    # Pre-allocation
    F,G = (zeros(ComplexF64, size(fieldPoints,2), size(sources,2)) for _ in 1:2)
    C   =  zeros(ComplexF64, size(fieldPoints,2))
    @inbounds @views Threads.@threads for fieldNode = axes(fieldPoints,2)
        # Creating views
        Cslice = C[fieldNode:fieldNode]
        source = fieldPoints[:,fieldNode]
        # Pre-allocation
        r, jacobian              = (zeros(n)   for _ in 1:2)
        interpolation,dX,normals = (zeros(2,n) for _ in 1:3)
        integrand                = zeros(ComplexF64,1,n)
        @inbounds @views for element = axes(topology,2)
            # Extracting element information
            elementNodes        = topology[:,element]
            elementCoordinates  = coordinates[:,elementNodes]
            physicsNodes        = physicsTopology[:,element]
            # views to correct part of matrices
            Fslice = F[fieldNode:fieldNode,physicsNodes]
            Gslice = G[fieldNode:fieldNode,physicsNodes]
            # Computing integral(s) over element
            meshInterpolation!(Fslice,Gslice,Cslice,elementInterpolation,elementDerivatives,
                                physicsInterpolation,weights,elementCoordinates,source,k,r,
                                interpolation,dX,normals,jacobian,integrand)
        end
    end
    return F,G,(interior ? C : 1.0 .+ C)
end
#==========================================================================================
                                        Meshing
==========================================================================================#
function createPhysics(ord,topology,coordinates)
    bF(ξ) = (size(topology,1) == 2 ? linear(ξ) : quadratic(ξ))
    physicsTopology = reshape(collect(1:(ord+1)*size(topology,2)), ord+1, size(topology,2))
    sources         = zeros(2,(ord+1)*size(topology,2))
    geometryElement = (ord == 0 ? bF(0) : (ord == 1 ? [bF(-β) bF(β)] : [bF(-δ) bF(0) bF(δ)]))
    for element = axes(topology,2)
        elementCoordinates  = @view coordinates[:,topology[:,element]]
        elementSources      = @view physicsTopology[:,element]
        sources[:,elementSources] = elementCoordinates * geometryElement
    end
    return physicsTopology,sources
end
function createTopology(ord,nElements)
    top = ones(Int64, ord+1, nElements)
    for i = 1:ord+1; top[i,1:length(i:ord:ord*nElements)] = i:ord:ord*nElements end
    return top
end
#==========================================================================================
                            Infinite Cylinder Helper Functions
                This for the analytical solution of an infinite cylinder
==========================================================================================#
pressure(mn,ka) = (mn == 0 ? atan(-besselj(1,ka)/bessely(1,ka)) :
            atan((besselj(mn-1,ka)-besselj(mn+1,ka))/(bessely(mn+1,ka)-bessely(mn-1,ka))))
function createCircle(topology,radius=1.0)
    angles = collect(range(0,2π,length=1+length(unique(topology))))
    return radius*[cos.(angles[1:end-1])'; sin.(angles[1:end-1])']
end
function cylscat(ϕ,ka,nterm=10)
    IntI,IntS = (ones(length(ϕ))*besselj(0,ka), zeros(length(ϕ)))
    for m=1:nterm-1; IntI = IntI + 2*im .^(m) * besselj(m,ka) .* cos.(m*(ϕ)) end
    for m=0:nterm-1;
        IntS = IntS - (m == 0 ? 1.0 : 2.0)*im .^(m+1.0) * exp(-im*pressure(m,ka)) *
                    sin.(pressure(m,ka)) * (besselj(m,ka)+im*bessely(m,ka)) .* cos.(m*(ϕ))
    end
    return IntI + IntS
end
#==========================================================================================
        Test case: Scattering of an infinite cylinder with hard boundaries (v=0 on ∂Ω)
==========================================================================================#
# Physical properties
freq = 600.0    # Frequency of interest [Hertz]
c = 340.0       # Speed up sound [m/s]
k = 2*π*freq/c  # Wavenumber [m^(-1)]

# geometryOrder: Denotes the order of the geometry
#  1 = Continuous Linear Elements
#  2 = Continuous Quadratic Elements
geometryOrder = 2

# physicsOrder: Denotes the order of the physics
# -1 = Same as geometry
#  0 = Discontinuous Constant Elements
#  1 = Discontinuous Linear Elements
#  2 = Discontinuous Quadratic Elements
physicsOrder = 2

# Setting the number of elements
# A rule of thumb says 6 elements per wavelength
el_wl = 6*freq/c
# The circumference of a circle is 2π. We can compute the number of elements required as
nElements = Int(ceil(el_wl*2π))
# nElements = 150
# NOTE: You will see problems for linear geoemtry elements.
# Theese errors are simply because the approximation of the cylinder is not good enough,
# and not because the number of elements can not describe the underlying physics.

### Actual BEM computations
topology    = createTopology(geometryOrder,nElements)   # Create a topology matrix
coordinates = createCircle(topology)                    # Create the coordinates of the mesh
F,G,C       = assemble(coordinates,topology,k,5;order=physicsOrder) # Assemble BEM matrices
# t1=createPhysics(-physicsOrder,topology,coordinates)[2]
# Collocation points depend on the physics
src = (physicsOrder ∉ [0,1,2] ? coordinates : createPhysics(physicsOrder,topology,coordinates)[2])
pIncident   = exp.(im*k*src[1,:])
# diag(C)p = G*v - F*p + pIncident
# With boundary condition v=0 we have that: p = (diag(C) + F)\pIncident
ps = (F + Diagonal(C))\pIncident

# Compute angles of the collocation points
θ  = angle.(src[1,:] + im*src[2,:]); θ[θ .< 0] .+= 2π
# Compute analytical solution
pA = cylscat(θ,k,150)
# Plot Results
plot(θ,abs.(pA),label="Analytical")
scatter!(θ,abs.(ps),label="BEM")
title!("Frequency = $(freq)")
xlabel!("Angle")
ylabel!("|p|")

# notes
nodes,weights  = gausslegendre(20)
m1 = [1;0.0]
m2 = [1.5;0.5]
m3 = [1;1.0]
coordinates = [m1 m2 m3]
elementInterpolation = quadratic(nodes')
interpolation = similar(coordinates*elementInterpolation)
mul!(interpolation,coordinates,elementInterpolation)
scatter(coordinates[1,:],coordinates[2,:])
plot!(interpolation[1,:],interpolation[2,:])
scatter!(interpolation[1,:],interpolation[2,:])
