#==========================================================================================
                                    README
==========================================================================================#
# This file contains a simple implementation of 2D-BEM for acoustics.
# It supports constant, linear and quadratic discontinuous and continuous elements.
#==========================================================================================
                                    Packages
==========================================================================================#
using FastGaussQuadrature   # For Numerical Integration
using LoopVectorization     # For fast small matrix-matrix multiplications
using LinearAlgebra         # For performing standard Linear Algebra
using SpecialFunctions      # For Hankel functions
using Plots                 # For visualizing the results
#==========================================================================================
                                Helper Functions
==========================================================================================#
JacMul!(j,w)    = @inbounds for i = axes(j,1) j[i]   = j[i]*w[i]                      end
dist!(x,y,r)    = @inbounds for i = axes(r,1) r[i]   = hypot(x[1,i]-y[1],x[2,i]-y[2]) end
jacobian!(j,n)  = @inbounds for i = axes(j,1) j[i]   = hypot(n[1,i],      n[2,i])     end
normalize!(n,j) = @inbounds for i = axes(j,1) n[1,i] = n[1,i]/j[i];n[2,i]=n[2,i]/j[i] end
normals!(n,dX)  = @inbounds for i = axes(n,2) n[1,i] = -dX[2,i];   n[2,i]=dX[1,i]     end
#==========================================================================================
                                Basis Functions
==========================================================================================#
linear(ξ)    = [(1 .- ξ)/2; (1 .+ ξ)/2]
quadratic(ξ) = [ξ .* (ξ .- 1)/2; 1 .- ξ .^2; ξ .* (ξ .+ 1)/2]
constant(ξ)  = ones(eltype(ξ),1,length(ξ))
β = gausslegendre(2)[1][end]; discLinear(ξ)    = linear(ξ/β)
δ = gausslegendre(3)[1][end]; discQuadratic(ξ) = quadratic(ξ/δ)
#==========================================================================================
                                    Kernels
==========================================================================================#
G!(k,r,int) = @inbounds for i = eachindex(r) int[i] = im/4*hankelh1(0,k*r[i]) end
function F!(x,y,k,n,r,int)
    @inbounds for i = eachindex(r)
        int[i] = -k*im*hankelh1(1,k*r[i])*(n[1,i]*(x[1,i] - y[1]) +
                                           n[2,i]*(x[2,i] - y[2]))/(4r[i])
    end
    return int
end
function C!(x,y,n,r,int)
    @inbounds for i = eachindex(r)
        int[i] = (n[1,i]*(x[1,i] - y[1]) +
                  n[2,i]*(x[2,i] - y[2]))/(2π*r[i]^2)
    end
    return int
end
#==========================================================================================
                                        Meshing
==========================================================================================#
function create_topology(ord,nElements)
    top = ones(Int64, ord+1, nElements)
    for i = 1:ord+1; top[i,1:length(i:ord:ord*nElements)] = i:ord:ord*nElements end
    return top
end
function mygemm_vec!(C, A, B)
    @inbounds @fastmath for n ∈ eachindex(C)
        Cmn = zero(eltype(C))
        for k ∈ eachindex(A)
            Cmn += A[k] * B[k,n]
        end
        C[n] += Cmn
    end
end
#==========================================================================================
                            Mesh interpolation routine
==========================================================================================#
function compute_integrands!(Fslice,Gslice,Cslice,interpolation,physics_interpolation,
                            source,k,r,normals,jacobian,integrand,gOn,fOn,cOn)
    dist!(interpolation,source,r)                           # Computing distances
    if gOn # Only compute G if you need it
        G!(k,r,integrand)                                       # Evaluate Greens function
        JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
        # mul!(Gslice,integrand,physics_interpolation!',true,true)  # Integration with basis func
        mygemm_vec!(Gslice,integrand,physics_interpolation')
    end
    if fOn # Only compute F if you need it
        F!(interpolation,source,k,normals,r,integrand)          # Evaluate ∂ₙGreens function
        JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
        # mul!(Fslice,integrand,physics_interpolation!',true,true)  # Integration with basis func
        mygemm_vec!(Fslice,integrand,physics_interpolation')
    end
    if cOn # Only compute c if you need it
        C!(interpolation,source,normals,r,integrand)            # Evaluate G₀
        Cslice[1] += dot(integrand,jacobian)                    # Integration of G₀
    end
end
#==========================================================================================
                    Creating continuous mesh topology and a circle geometry
==========================================================================================#
function createPhysicsTopology(N,physicsOrder)
    return physicsOrder < 0 ? create_topology(-physicsOrder,N) : reshape(collect(1:(physicsOrder+1)*N), physicsOrder+1, N)
end
function create_circle(topology,radius=1.0)
    angles = collect(range(0,2π,length=1+length(unique(topology))))
    return radius*[cos.(angles[1:end-1])'; sin.(angles[1:end-1])']
end

#==========================================================================================
                            Assembly - Spline
==========================================================================================#
function assemble(spline,k,n;interior=false,gOn=true,fOn=true,cOn=true,
                physicsOrder=-1,fieldPoints=[],nElements=10,N=10)
    # Setting up element division of Spline
    physics_function(ξ) = (physicsOrder == -2 ? quadratic(ξ)  :
                          (physicsOrder == -1 ? linear(ξ)     :
                          (physicsOrder ==  0 ? constant(ξ)   :
                          (physicsOrder ==  1 ? discLinear(ξ) : discQuadratic(ξ)))))
    nInterpolations = max(abs(physicsOrder),1)
    tInterp         = range(1,N,length=nInterpolations*nElements+1)
    nodes,weights   = scaledgausslegendre(n,0.0,1.0)
    physicsTopology = createPhysicsTopology(nElements,physicsOrder)
    # Compute stuff for every Gaussian node
    physics_interpolation! = physics_function(gausslegendre(n)[1]')
    tGauss = kron(ones(nElements),nodes) + kron(tInterp[1:nInterpolations:end-1],ones(n))
    Interpolation = spline(tGauss)
    DerivInterp   = derivative(spline,tGauss)
    Weights  = kron(ones(nElements),weights)/(nElements*sum(weights))*(N-1)
    Normals  = similar(DerivInterp);   normals!(Normals,DerivInterp) # Allocate & compute
    Jacobian = zeros(size(Normals,2)); jacobian!(Jacobian,Normals)   # Allocate & compute
    normalize!(Normals,Jacobian) # Normalizing using the Jacobian (length of normal)
    JacMul!(Jacobian,Weights)    # NB! Now jacobian = jacobian*weights (abuse of notation)
    # Computing the physics-nodes by interpolation of the spline
    tInterp = (physicsOrder == 0 ?  tInterp[1:end-1]    .+ 0.5*(tInterp[2]-tInterp[1])       :
              (physicsOrder == 1 ? [tInterp[1:end-1]'   .+ 0.5*(tInterp[2]-tInterp[1])*β     ;
                                    tInterp[2:end]'     .- 0.5*(tInterp[2]-tInterp[1])*β][:] :
              (physicsOrder == 2 ? [tInterp[1:2:end-1]' .+ 0.5*(tInterp[3]-tInterp[1])*β     ;
                                    tInterp[2:2:end]'                                        ;
                                    tInterp[3:2:end]'   .- 0.5*(tInterp[3]-tInterp[1])*β][:]
                                    : tInterp[1:end-1])))
    fieldPoints = (isempty(fieldPoints) ? spline(tInterp) : fieldPoints)
    # Preallocation
    F,G = (zeros(ComplexF64, size(fieldPoints,2), length(tInterp)) for _ in 1:2)
    C   =  zeros(ComplexF64, size(fieldPoints,2))
    @inbounds @views Threads.@threads for fieldNode = 1:size(fieldPoints,2)
        Cslice    = C[fieldNode:fieldNode]
        source    = fieldPoints[:,fieldNode]
        r         = zeros(n)
        integrand = zeros(ComplexF64,1,n)
        @inbounds @views for element = 1:nElements
            physicsNodes = physicsTopology[:,element]
            Fslice = F[fieldNode:fieldNode,physicsNodes]
            Gslice = G[fieldNode:fieldNode,physicsNodes]
            interpolation = Interpolation[:,n*(element-1)+1:n*element]
            jacobian = Jacobian[n*(element-1)+1:n*element]
            normals  = Normals[:,n*(element-1)+1:n*element]
            compute_integrands!(Fslice,Gslice,Cslice,interpolation,physics_interpolation!,
                                source,k,r,normals,jacobian,integrand,gOn,fOn,cOn)
        end
    end
    return F,G,(interior ? C : 1.0 .+ C),fieldPoints
end
#==========================================================================================
                        Infinite Cylinder Helper Functions
==========================================================================================#
pressure(mn,ka) = (mn == 0 ? atan(-besselj(1,ka)/bessely(1,ka)) :
            atan((besselj(mn-1,ka)-besselj(mn+1,ka))/(bessely(mn+1,ka)-bessely(mn-1,ka))))
function cylscat(ϕ,ka,nterm=10)
    IntI,IntS = (ones(length(ϕ))*besselj(0,ka), zeros(length(ϕ)))
    for m=1:nterm-1; IntI = IntI + 2*im .^(m) * besselj(m,ka) .* cos.(m*(ϕ)) end
    for m=0:nterm-1;
        IntS = IntS - (m == 0 ? 1.0 : 2.0)*im .^(m+1.0) * exp(-im*pressure(m,ka)) *
                    sin.(pressure(m,ka)) * (besselj(m,ka)+im*bessely(m,ka)) .* cos.(m*(ϕ))
    end
    return IntI + IntS
end
