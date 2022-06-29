using FMM3D
using Random
using LinearAlgebra
using StatsBase
using IntegralEquations

# Defining Frequency
freq = 500.0
_,_,kp,ka,kh,kv,_,_,_,_,_,_ = visco_thermal_constants(;freq=freq,S=1)
zk = Complex(kp)
#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, IntegralEquations, Plots
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                physics_order=tri_physics_orders[3])

function unroll_interpolations(interpolations)
    N = length(interpolations)
    Nweights = length(interpolations[1].jacobian_mul_weights)

    weights = zeros(N*Nweights)
    interps = zeros(3,N*Nweights)
    normals = zeros(3,N*Nweights)

    for i = 1:N
        weights[(i-1)*Nweights + 1:i*Nweights]   = interpolations[i].jacobian_mul_weights
        interps[:,(i-1)*Nweights + 1:i*Nweights] = interpolations[i].interpolation
        normals[:,(i-1)*Nweights + 1:i*Nweights] = interpolations[i].normals
    end

    return interps, weights, normals
end
@time interpolations = IntegralEquations.interpolate_elements(mesh;n=2,m=2);
@time interp,weights,normals = unroll_interpolations(interpolations);
sum(weights)/4π # Total weights should add up to surface area

# Creating targets for FMM
targets = mesh.sources
w = ComplexF64.(weights)
acc = 1e-6
@time vals = hfmm3d(acc,zk,interp,charges=w,targets=targets,pgt=1);
@time Fp,Gp,Cp = assemble_parallel!(mesh,kp,mesh.sources,n=2,m=2,sparse=false);


# Checking if FMM and BEM gives the same results
N = length(Cp)
@time vals = hfmm3d(acc,zk,interp,charges=w,targets=targets,pgt=1);
(vals.pottarg/4π - Gp*ones(N))./(Gp*ones(N))

#==========================================================================================
                                Setting up constants
==========================================================================================#
radius = 1.0                                    # Radius of sphere_1m       [m]
# Computing incident pressure
pI = IntegralEquations.incoming_wave(angles,1.0,mesh.sources,zk)

rle = 4ones(Int64,N)
inverse_rle(pI,rle)

#==========================================================================================
                            Defining LossyBlockMatrix
            Block matrix corresponding to 5 BEM systems and 5 constraints
==========================================================================================#
using LinearMaps
struct FMMGOperator{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    m::Int64                            # Number of gauss nodes
    # Physical Quantities
    k::T                                # Wavenumber
    eps::Float64                        # Precision
    # Acoustic Matrices
    targets::AbstractMatrix{Float64}    # FMM-targets (size = 3,n)
    sources::AbstractMatrix{Float64}    # FMM-sources (size = 3,m)
    weights::AbstractVecOrMat{T}        # Nodal weights (size = m)
    #
    rle::AbstractVecOrMat
end
function FMMGOperator(eps,k,targets,sources,weights,N=4)
    n = size(targets,2)
    m = size(sources,2)
    rle = N*ones(Int64,n)
    return FMMGOperator(n,m,k,eps,targets,sources,weights,rle)
end
Base.size(A::FMMGOperator) = (A.n, A.n)
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMGOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    # A.tmp .= inverse(rle)
    tmp  = inverse_rle(x,A.rle)
    weights = A.weights .* tmp
    vals = hfmm3d(A.eps,zk,A.sources,charges=weights,targets=A.targets,pgt=1)
    y .= vals.pottarg/4π
end
# FMMF
struct FMMFOperator{T} <: LinearMaps.LinearMap{T}
    n::Int64                            # Number of nodes
    m::Int64                            # Number of gauss nodes
    # Physical Quantities
    k::T                                # Wavenumber
    eps::Float64                        # Precision
    # Acoustic Matrices
    targets::AbstractMatrix{Float64}    # FMM-targets (size = 3,n)
    sources::AbstractMatrix{Float64}    # FMM-sources (size = 3,m)
    normals::AbstractMatrix{T}        # Nodal weights (size = m)
    #
    rle::AbstractVecOrMat
end
function FMMFOperator(eps,k,targets,sources,normals,weights,N=4)
    n = size(targets,2)
    m = size(sources,2)
    rle = N*ones(Int64,n)
    normals = normals .* weights'
    return FMMFOperator(n,m,k,eps,targets,sources,normals,rle)
end
Base.size(A::FMMFOperator) = (A.n, A.n)
# Standard Multiplication
function LinearAlgebra.mul!(y::AbstractVecOrMat{T},
                            A::FMMFOperator{T},
                            x::AbstractVector) where {T <: ComplexF64}
    LinearMaps.check_dim_mul(y, A, x)
    # A.tmp .= inverse(rle)
    tmp  = inverse_rle(x,A.rle)
    weights = A.normals .* Transpose(tmp)
    vals = hfmm3d(A.eps,zk,A.sources,targets=A.targets,dipvecs=weights,pgt=2)
    y .= vals.pottarg/4π + 0.5*x
end

Ag = FMMGOperator(acc,zk,targets,interp,w)
Af = FMMFOperator(acc,zk,targets,interp,normals,w)

x = ones(eltype(Ag),size(Ag,1))
Ag*x
Af*x


using IterativeSolvers
# Solving scattering
Ap    = Fp + Diagonal(1.0 .- Cp);
maximum(abs.(Af*x - Ap*x))
p_bem = gmres(Ap,pI;verbose=true);
p_fmm = gmres(Af,pI;verbose=true);
p_bem./p_fmm

surface_angles = acos.(mesh.sources[1,:]/radius)
perm = sortperm(surface_angles)
p_analytical, _ = IntegralEquations.plane_wave_scattering_sphere(zk,radius,1.0,surface_angles,1e-6)
plot(surface_angles[perm], real.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],real.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],real.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)
plot(surface_angles[perm], imag.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],imag.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],imag.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)

plot(surface_angles[perm], abs.(p_analytical[perm]),label="Analytical",linewidth=2)
plot!(surface_angles[perm],abs.(p_bem[perm]),label="BEM",linestyle=:dash,linewidth=2)
plot!(surface_angles[perm],abs.(p_fmm[perm]),label="FMM",linestyle=:dash,linewidth=2)


# Solving radiating sphere - Takes long time to converge
# p_bem = gmres(Gp,pI;verbose=true);
# p_fmm = gmres(A,pI;verbose=true);
# p_bem./p_fmm
