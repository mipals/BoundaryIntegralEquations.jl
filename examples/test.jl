#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra, IntegralEquations, Plots, IterativeSolvers, FMM3D
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/sphere_1m"
tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
                                                physics_order=tri_physics_orders[2])
#==========================================================================================
        3d Visualization - Seems highly unstable on M1 chips. Problems with GLMakie?
==========================================================================================#
# using Meshes, MeshViz
# ##choose a Makie backend
# import GLMakie as Mke
# simple_mesh = create_simple_mesh(mesh)
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
freq   = 100.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
@time BB = IntegralEquations.LossyBlockMatrix(mesh,freq;blockoutput=true,depth=1)

# The error comes from the single-layer potential
A = FMMGOperator(mesh,ka;n=3,eps=1e-6,offset=0.2,nearfield=true)
x = rand(ComplexF64,size(A,1))
maximum(abs.(A*x - BB.Bₐ*x))

A = FMMFOperator(mesh,ka;n=3,eps=1e-6,offset=0.2,nearfield=true)
x = rand(ComplexF64,size(A,1));
maximum(abs.(A*x - BB.Aₐ*x))

IntegralEquations.nodes_to_gauss!(A.tmp_weights,A.element_interpolation,A.physics_topology,x)
IntegralEquations.integrand_mul!(A.tmp_weights,A.weights)
vals = hfmm3d(A.eps,A.k,A.sources,charges=A.weights,targets=A.targets,pgt=1)
y = vals.pottarg + A.nearfield_correction*x

A.targets
A.sources

function create_ga(targets,sources,k)
    Ga = zeros(ComplexF64,size(targets,2),size(sources,2))
    for i = 1:size(Ga,1), j = 1:size(Ga,2)
        r = norm(targets[:,i] - sources[:,j])
        if r > 1e-12
            Ga[i,j] = exp(im*k*r)/(4π*r)
        end
    end
    return Ga
end
Ga = create_ga(A.targets,A.sources,A.k)
maximum(abs.(Ga*A.tmp_weights - A*x))
maximum(abs.(Ga*A.tmp_weights - BB.Bₐ*x))
maximum(abs.(A*x - BB.Bₐ*x))

maximum(abs.(A*x)./abs.(BB.Bₐ*x))
minimum(abs.(A*x)./abs.(BB.Bₐ*x))
