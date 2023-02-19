#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using HMatrices
using Plots
using IterativeSolvers
#==========================================================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes")
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_fine");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_finer");
tri_mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_finest");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_35k");
# tri_mesh_file = joinpath(mesh_path,"sphere_1m_77k");
@time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[2])
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
ω = 2π*freq
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
xyzb = mesh.sources
n = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=1)
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
Gf = FMMGOperator(mesh,k)
Gh = BoundaryIntegralEquations.HGOperator(mesh,k)
x = ones(n)
varinfo(r"G")

yf = Gf*x
yh = Gh*x

(yh - yf)./yf

@time Gf*x;
@time Gh*x;


@time xf = gmres(Gf,yf;verbose=true,maxiter=100)
@time xh = gmres(Gh,yh;verbose=true,maxiter=100)

norm(xf - x)/norm(x)
norm(xh - x)/norm(x)
norm(xh - xf)/norm(xf)

# The IFGF is currently problematic as it only accepts AbstractTrees = "0.3" (we need 0.4)
using StaticArrays
const Point3D = SVector{3,Float64}
# sample some points on a sphere
X = [Point3D(x) for x in eachcol(mesh.sources)]
Y = [Point3D(y) for y in eachcol(Gf.sources)]

## H-matrices
struct HelmholtzDoubleLayer <: AbstractMatrix{ComplexF64}
    X::Vector{Point3D}
    Y::Vector{Point3D}
    NY::Vector{Point3D} # normals at Y coordinate
    k::Float64
end
# \\frac{e^{ikr_j}}{4\\pi r_j}(ikr_j - 1) (x_j - y)\\cdot n
function Base.getindex(K::HelmholtzDoubleLayer,i::Int,j::Int)
    # r = K.X[i] - K.Y[j]
    r = K.Y[j] - K.X[i]
    d = norm(r)
    return exp(im*K.k*d)/(4π*d^3) * (im*K.k*d - 1) * dot(r, K.NY[j])
end
Base.size(K::HelmholtzDoubleLayer) = length(K.X), length(K.Y)

Hf = FMMHOperator(mesh,k)
NY = [Point3D(n) for n in eachcol(Hf.normals)]
K = HelmholtzDoubleLayer(X,Y,NY,k)
Xclt = ClusterTree(X)
Yclt = ClusterTree(Y)
Hh = assemble_hmat(K,Xclt,Yclt;comp=PartialACA(;rtol=1e-6))

@time zf = Hf*x
@time zh = -(Hh*(Hf.C*x)) + Hf.nearfield_correction*x
(zh - zf)./zf

# Just to check if the FMM is correct (it is)
# @time Fp,_,_ = assemble_parallel!(mesh,k,mesh.sources,n=2,m=2,sparse=false);
# Hp = Fp + 0.5I
# @time zp = Hp*x
# (zf - zp)./zp


# The package makes it possible to plot the rank-structure!
using Plots
plot(Hh;aspect_ratio=1)



## Alternative IFGF implementation. The developer said it was experimental.
# using IFGF
# struct HelmholtzMatrix <: AbstractMatrix{ComplexF64}
#     X::Vector{Point3D}
#     Y::Vector{Point3D}
#     k::Float64
# end
# IFGF.wavenumber(A::HelmholtzMatrix) = A.k
# # abstract matrix interface
# Base.size(A::HelmholtzMatrix) = length(X), length(Y)
# Base.getindex(A::HelmholtzMatrix,i::Int,j::Int) = A(A.X[i],A.Y[j])
# # functor interface
# function (K::HelmholtzMatrix)(x,y)
#     k = IFGF.wavenumber(K)
#     d = norm(x-y)
#     return exp(im*k*d)/(4π*d)
# end
# A = HelmholtzMatrix(X,Y,k)
# @time Gi = assemble_ifgf(A,X,Y;tol=1e-3);
# coeffs = Gf.C*x
# yi = Gi*coeffs
# yi += Gf.nearfield_correction*x

# # Relatvive error between Hmatrix with fmm
# (yi - yf)./yf

# # The memory vs.
# @time Hh*(Ga.C*x); # 1
# @time Hi*(Ga.C*x); # 2
# @time Ga*x;        # 2

# varinfo(r"Hh") # 3
# varinfo(r"Hi") # 2
# varinfo(r"Ga") # 1
