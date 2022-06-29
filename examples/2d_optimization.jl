#==========================================================================================
                                    Packages
==========================================================================================#
using FastGaussQuadrature, LinearAlgebra, SpecialFunctions, FiniteDifferences, Dierckx, Plots
#==========================================================================================
                                Helper Functions
==========================================================================================#
JacMul!(j,w)    = @inbounds for i = 1:length(j) j[i]   = j[i]*w[i]                      end
dist!(x,y,r)    = @inbounds for i = 1:length(r) r[i]   = hypot(x[1,i]-y[1],x[2,i]-y[2]) end
jacobian!(j,n)  = @inbounds for i = 1:length(j) j[i]   = hypot(n[1,i],      n[2,i])     end
normalize!(n,j) = @inbounds for i = 1:length(j) n[1,i] = n[1,i]/j[i];n[2,i]=n[2,i]/j[i] end
normals!(n,dX)  = @inbounds for i = 1:size(n,2) n[1,i] = -dX[2,i];   n[2,i]=dX[1,i]     end
#==========================================================================================
                                Basis Functions
==========================================================================================#
linear(ξ)    = [0.5*(1.0 .- ξ); 0.5*(1.0 .+ ξ)]
quadratic(ξ) = [0.5 * ξ .* (ξ .- 1.0); 1.0 .- ξ .^2; 0.5 * ξ .* (ξ .+ 1.0)]
constant(ξ)  = ones(1,length(ξ))
β = gausslegendre(2)[1][end]; discLinear(ξ)    = linear(ξ/β)
δ = gausslegendre(3)[1][end]; discQuadratic(ξ) = quadratic(ξ/δ)
#==========================================================================================
                                    Kernels
==========================================================================================#
G!(k,r,int)       = @inbounds for i = 1:length(r) int[i] =  im*0.25*hankelh1.(0,k*r[i]) end
F!(x,y,k,n,r,int) = @inbounds for i = 1:length(r) int[i] = -im*0.25*hankelh1.(1,k*r[i])*k.*(n[1,i]*(x[1,i]-y[1])+n[2,i]*(x[2,i]-y[2]))/r[i] end
C!(x,y,n,r,int)   = @inbounds for i = 1:length(r) int[i] = (n[1,i]*(x[1,i]-y[1])+n[2,i]*(x[2,i]-y[2]))/(2π*r[i]^2) end
#==========================================================================================
                                        Meshing
==========================================================================================#
function create_topology(ord,nElements)
    top = ones(Int64, ord+1, nElements)
    for i = 1:ord+1; top[i,1:length(i:ord:ord*nElements)] = i:ord:ord*nElements end
    return top
end
# Test
using LoopVectorization
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
function compute_integrands!(Fslice,Gslice,Cslice,interpolation,physics_interpolation!,
                            source,k,r,normals,jacobian,integrand,gOn,fOn,cOn)
    dist!(interpolation,source,r)                           # Computing distances
    if gOn # Only compute G if you need it
        G!(k,r,integrand)                                       # Evaluate Greens function
        JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
        # mul!(Gslice,integrand,physics_interpolation!',true,true)  # Integration with basis func
        mygemm_vec!(Gslice,integrand,physics_interpolation!')
    end
    if fOn # Only compute F if you need it
        F!(interpolation,source,k,normals,r,integrand)          # Evaluate ∂ₙGreens function
        JacMul!(integrand,jacobian)                             # Multiply by jacobian*weights
        # mul!(Fslice,integrand,physics_interpolation!',true,true)  # Integration with basis func
        mygemm_vec!(Fslice,integrand,physics_interpolation!')
    end
    if cOn # Only compute c if you need it
        C!(interpolation,source,normals,r,integrand)            # Evaluate G₀
        Cslice[1] += dot(integrand,jacobian)                    # Integration of G₀
    end
end
#==========================================================================================
                            Assembly - Spline
==========================================================================================#
function scaledgausslegendre(n,a=0.0,b=1.0)
    node, w = gausslegendre(n)
    return 0.5*node*(b-a) .+ 0.5*(b+a) , 0.5*(b-a)*w
end
function createPhysicsTopology(N,physicsOrder)
    return physicsOrder < 0 ? create_topology(-physicsOrder,N) : reshape(collect(1:(physicsOrder+1)*N), physicsOrder+1, N)
end
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
                                Optimization
==========================================================================================#
kappa(spline,t) = abs(det([derivative(spline,t) derivative(spline,t;nu=2)]))/hypot(derivative(spline,t)...)^3
function objective(x,spl,k,n,fieldpoints;interior=false,gOn=true,fOn=true,cOn=true,
                    physicsOrder=-1,nElements=10,N=10)
    spline = deepcopy(spl)
    Nspline = size(spline.c,2)
    spline.c[:,2:Nspline-2] = reshape(x,2,Nspline-3)
    spline.c[:,[1 Nspline-1 Nspline]] = spline.c[:,[Nspline-2 2 3]]
    F,_,C,src    = assemble(spline,k,n;interior=interior,gOn=gOn,fOn=fOn,cOn=cOn,
                                        physicsOrder=physicsOrder,nElements=nElements,N=N)
    Fp,Gp,Cp,_   = assemble(spline,k,n;interior=interior,gOn=gOn,fOn=fOn,cOn=cOn,
                physicsOrder=physicsOrder,fieldPoints=fieldpoints,nElements=nElements,N=N)
    return sum(abs.(Diagonal(Cp)\((Fp - Gp)*((F + Diagonal(C))\exp.(im*k*src[1,:])))))
end
obj(x)     = objective(x,spline,k,n,fieldpoints;gOn=false,physicsOrder=physicsOrder,nElements=nElements,N=N)
gradObj(x) = grad(central_fdm(2,1),obj,x)[1]


freq = 200.0; c=340.0;k = 2.0*pi*freq/c; n = 6; nElements = 20; physicsOrder = 0;
N = 19+1; coordinates = createCircle(collect(1:N-1)); coordinates = [coordinates coordinates[:,1]]
spline = ParametricSpline(1:size(coordinates,2),coordinates;periodic=true)
Nspline = size(spline.c,2)
x = range(-2.2,-2.0,length=5)
y = range(-0.1, 0.1,length=5)
fieldpoints = [kron(x,ones(length(y)))'; kron(ones(length(x)),y)']
x0 = (spline.c[:,2:size(spline.c,2)-2])[:]
# anim = Plots.Animation()
for i = 1:10
    x0 = x0 + 0.002*gradObj(x0)
    spline.c[:,2:Nspline-2] = reshape(x0,2,Nspline-3)
    spline.c[:,[1 Nspline-1 Nspline]] = spline.c[:,[Nspline-2 2 3]]
    println(obj(x0))
    # scatter(fieldpoints[1,:],fieldpoints[2,:],aspect_ratio=1,legend=:topleft,label="Fieldpoints",markerstyle=:dot)
    plt = scatter(fieldpoints[1,:],fieldpoints[2,:],aspect_ratio=1,legend=:topleft,label="Fieldpoints",markerstyle=:dot)
    scatter!(plt,eachcol(spline.c')...,label="Control Points")
    plot!(plt,eachcol(spline(1:0.1:get_knots(spline)[end])')...,label="Geometry")
    # display(plt)
    # Plots.frame(anim)
end
# gif(anim, "testAnim.gif", fps = 4)

gradObj(x)  = grad(central_fdm(2,1),obj,x)[1]
gradbobj(x) = grad(backward_fdm(2,1),obj,x)[1]
gradfobj(x) = grad(forward_fdm(2,1),obj,x)[1]
gradcobj(x) = grad(central_fdm(2,1),obj,x)[1]

@time gradbobj(x0);
@time gradfobj(x0);
@time gradcobj(x0);
@time gradObj(x0);

#==========================================================================================
                                Remeshing?
==========================================================================================#
function compute_segment_lengths(spl,N)
    t_knots = get_knots(spl)
    segment_lenths = zeros(length(t_knots)-1)
    jacobian = zeros(N)
    for idx = 1:length(segment_lenths)
        nodes,weights=scaledgausslegendre(N,t_knots[idx],t_knots[idx+1])
        jacobian!(jacobian,derivative(spl,nodes))
        segment_lenths[idx] = dot(jacobian,weights)
    end
    return segment_lenths
end
compute_arc_length(spl,N) = sum(compute_segment_lengths(spl,N))
function re_meshing(spl,nControl,N=4000)
    knots = get_knots(spl)
    maxl = compute_arc_length(spl,10)/nControl

    nodes,weights = scaledgausslegendre(N,knots[1],knots[end])
    tmp = 0.0
    nNodes = length(nodes)
    interps = 1
    i = 1
    while i < nNodes
        while tmp < maxl && i < nNodes
            tmp += hypot(derivative(spl,nodes[i])...)*weights[i]
            i = i+1
        end
        interps = [interps; i]
        tmp = 0.0
    end
    coords = spl(nodes[interps]);
    coords[:,end] = coords[:,1]
    return coords
end
coo = re_meshing(spline,19)
splT = ParametricSpline(1:size(coo,2),coo;periodic=true)
tmp = compute_segment_lengths(splT,10)


# plot(eachcol(splT(1:0.1:get_knots(spline)[end])')...,aspect_ratio=1,label="Geometry")
arc_length = compute_arc_length(splT,10)
elements_per_wavelength = 6.0*freq/c
element_spacing = 1.0/elements_per_wavelength
nElements = ceil(Int,arc_length/element_spacing)

nKnots = size(coo,2)
_,_,_,srcNew = assemble(splT,k,n;physicsOrder=physicsOrder,nElements=nElements,N=nKnots,gOn=false,fOn=false)
_,_,_,srcOld = assemble(spline,k,n;physicsOrder=physicsOrder,nElements=nElements,N=nKnots,gOn=false,fOn=false)
plot(eachcol(splT(1:0.01:get_knots(splT)[end])')...,aspect_ratio=1,label="NewGeometry")
plot!(eachcol(spline(1:0.01:get_knots(spline)[end])')...,aspect_ratio=1,label="OldGeometry")
scatter!(srcOld[1,:],srcOld[2,:],label="Old")
scatter!(srcNew[1,:],srcNew[2,:],label="New")
