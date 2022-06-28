#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using SpecialFunctions, LinearAlgebra, IntegralEquations, Plots, IterativeSolvers
#==========================================================================================
                Loading Mesh + Visualization (viz seems broken on M1 chips :/ )
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
# tri_mesh_file = "examples/meshes/sphere_1m"
# tri_mesh_file = "examples/meshes/sphere_1m_fine"
# tri_mesh_file = "examples/meshes/sphere_1m_finer"
tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],
                                                physics_order=tri_physics_orders[1])
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
freq   = 1000.0                                   # Frequency                 [Hz]
rho,c,kp,ka,kh,kv,ta,th,phi_a,phi_h,eta,mu = visco_thermal_constants(;freq=freq,S=-1)
k      = 2*π*freq/c                              # Wavenumber                [1/m]
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Assembling BEM matrices
==========================================================================================#
@time BB = IntegralEquations.LossyBlockMatrix(mesh,freq;blockoutput=true,depth=1)
# Fa,Ba,C0 = assemble_parallel!(mesh,ka,mesh.sources;m=3,n=3)
# cond(Matrix(BB.Aᵥ))
# cond(Matrix(BB.Bᵥ))
# cond(Matrix(BB.Aₕ))
# cond(Matrix(BB.Bₕ))

xyzb = mesh.sources
M  = size(xyzb,2)
u₀ = 1e-2
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents
ρ,c,kₚ,kₐ,kₕ,kᵥ,τₐ,τₕ,ϕₐ,ϕₕ,η,μ = visco_thermal_constants(;freq=freq,S=-1)

# Creating RHS
vn0  = u₀*normals[3,:]
vt10 = u₀*tangent1[3,:]
vt20 = u₀*tangent2[3,:]
#===========================================================================================
                        Iterative Solution of the 1-variable system
===========================================================================================#
# Creating the "outer"-struct representation of the system matrix
outer = LossyOneVariableOuter(mesh,BB,freq)
bt    = compute_lossy_rhs(outer,vn0,vt10,vt20)
# Solving the problem using the "outer" struct
@time pa,pa_hist = gmres(outer,bt;verbose=true,log=true)
#===========================================================================================
                                Checking results
===========================================================================================#
ghpainv,gh_hist = gmres(outer.Gh,outer.Hh*pa;verbose=true,log=true)
gapainv,ga_hist = gmres(outer.Ga,outer.Ha*pa;verbose=true,log=true)
ghpainv = ghpainv*ϕₕ*τₐ/τₕ
gapainv = gapainv*ϕₐ
# gapainv = gmres!(x0,outer.Ga,outer.Ha*pa;verbose=true)*ϕₐ
vx = normals[1,:] .* vn0 + tangent1[1,:] .* vt10 + tangent2[1,:] .* vt20 -
    (normals[1,:] .* (gapainv - ghpainv)         +
   (tangent1[1,:] .* (outer.Dt1*pa)              +
    tangent2[1,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vy = normals[2,:] .* vn0 + tangent1[2,:] .* vt10 + tangent2[2,:] .* vt20 -
    (normals[2,:] .* (gapainv - ghpainv)         +
   (tangent1[2,:] .* (outer.Dt1*pa)              +
    tangent2[2,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))
vz = normals[3,:] .* vn0 + tangent1[3,:] .* vt10 + tangent2[3,:] .* vt20 -
    (normals[3,:] .* (gapainv - ghpainv)         +
   (tangent1[3,:] .* (outer.Dt1*pa)              +
    tangent2[3,:] .* (outer.Dt2*pa))*(outer.ϕₐ - outer.ϕₕ*outer.τₐ/outer.τₕ))

# Local components of the viscous velocity on the boundary
v_n0   = -(normals[1,:].*vx +  normals[2,:].*vy +  normals[3,:].*vz) # The normal points inwards
v_t1   =  tangent1[1,:].*vx + tangent1[2,:].*vy + tangent1[3,:].*vz
v_t2   =  tangent2[1,:].*vx + tangent2[2,:].*vy + tangent2[3,:].*vz
vt_sum = sqrt.(v_t1.^2 + v_t2.^2)

#===========================================================================================
                                   Plotting solutions
===========================================================================================#
# Generating analytical solution
pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V =
    IntegralEquations.sphere_first_order(kₐ,c,ρ,radius,u₀,[radius*ones(M,1) acos.(xyzb[3,:]/radius)];S=-1,kv=kᵥ)
ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
perm = sortperm(ang_axis)

# Plotting
# plt1 = scatter(ang_axis,abs.(pa),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|p|"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2)
title!("Frequency = $(freq)")
# plt2 = scatter(ang_axis,abs.(v_n0),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black)
ylabel!("|Vn|"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2)
ylims!((-1e-6,1e-6))
# plt3 = scatter(ang_axis,abs.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:yellow)
scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black)
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2)
xlabel!("Angle"); ylabel!("|Vt|")
# plot(plt1,plt2,plt3,layout=(3,1))
# plt = plot(plt1,plt2,plt3,layout=(3,1),dpi=300)


# count = IntegralEquations.count_elements_pr_node(mesh)
# physics_function = mesh.physics_function
# IntegralEquations.set_nodal_interpolation!(physics_function)
# physics_function.interpolation
# d = 1e-9
# maximum(abs.((physics_function(physics_function.gauss_u' .+ d, physics_function.gauss_v') - I)/d - physics_function.derivatives_u))
# maximum(abs.((physics_function(physics_function.gauss_u', physics_function.gauss_v' .+ d) - I)/d - physics_function.derivatives_v))

linmesh  = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[1],physics_order=tri_physics_orders[1])
quadmesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2], physics_order=tri_physics_orders[2])

lin_nr  = IntegralEquations.count_elements_pr_node(linmesh)
quad_nr = IntegralEquations.count_elements_pr_node(quadmesh)

edges = sort(unique(quadmesh.topology[4:6,:]))
sum(2 .== uad_nr[edges]) == length(edges)

topology    = mesh.topology
coordinates = mesh.coordinates
nElements   = IntegralEquations.number_of_elements(mesh)
n_sources   = size(mesh.sources,2)

# Making a copy of the element type
shape_function   = deepcopy(mesh.shape_function)
physics_function = deepcopy(mesh.physics_function)

# Setting interpolations to be at the nodal positions
# (since this is where we want to compute the derivatives)
IntegralEquations.set_nodal_interpolation!(physics_function)
IntegralEquations.copy_interpolation_nodes!(shape_function,physics_function)
physics_function.derivatives_u - physics_function.derivatives_u
physics_function.derivatives_v - physics_function.derivatives_v
physics_function.interpolation - physics_function.interpolation
physics_function.weights       - physics_function.weights
physics_function.gauss_u       - physics_function.gauss_u
physics_function.gauss_v       - physics_function.gauss_v

n_shape = IntegralEquations.number_of_shape_functions(physics_function)
n_elements_pr_node = IntegralEquations.count_elements_pr_node(mesh)
# Pre-Allocation
global_gradients    = zeros(3,n_shape)   # ∇ₓ
local_gradients     = zeros(3,n_shape)   # ∇ξ
change_of_variables = zeros(3,3)         # COV-matrix
dX  = zeros(3,n_shape) # X-Derivatives
dY  = zeros(3,n_shape) # Y-Derivatives
dZ  = zeros(3,n_shape) # Z-Derivatives
using SparseArrays
Dx  = spzeros(n_sources,n_sources)  # ∂x
Dy  = spzeros(n_sources,n_sources)  # ∂y
Dz  = spzeros(n_sources,n_sources)  # ∂z

# for element = 1:nElements
    element = 1
    # Extract element and compute
    element_nodes        = @view topology[:,element]
    element_coordinates  = coordinates[:,element_nodes]
    IntegralEquations.my_mul!(dX,element_coordinates,shape_function.derivatives_u)
    IntegralEquations.my_mul!(dY,element_coordinates,shape_function.derivatives_v)
    IntegralEquations.cross_product!(dZ,dX,dY)

    physics_nodes = @view mesh.physics_topology[:,element]
    # Add nodal gradients to the shape function derivatives
    # for node = 1:n_shape
        node = 2
        change_of_variables[1,:] .= dX[:,node]
        change_of_variables[2,:] .= dY[:,node]
        change_of_variables[3,:] .= dZ[:,node]
        local_gradients[1,:]     .= physics_function.derivatives_u[:,node]
        local_gradients[2,:]     .= physics_function.derivatives_v[:,node]
        global_gradients         .= change_of_variables\local_gradients

        Dx[physics_nodes[node],physics_nodes] += global_gradients[1,:] ./
                                                 n_elements_pr_node[physics_nodes[node]]
        Dy[physics_nodes[node],physics_nodes] += global_gradients[2,:] ./
                                                 n_elements_pr_node[physics_nodes[node]]
        Dz[physics_nodes[node],physics_nodes] += global_gradients[3,:] ./
                                                 n_elements_pr_node[physics_nodes[node]]
    # end
