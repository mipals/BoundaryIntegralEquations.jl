using IntegralEquations
using LinearAlgebra
using Plots

meshFile = "examples/meshes/cylinder"
mesh = load3dTriangularComsolMesh(meshFile)
mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="linear")

mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctriconstant")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctrilinear")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctriquadratic")


element_interpolations = IntegralEquations.interpolate_elements(mesh)

S = 0.0
for i = 1:length(element_interpolations)
    S += dot(element_interpolations[i].jacobians,mesh.shape_function.weights)
end
S

@time IntegralEquations.interpolate_elements(mesh);


# quad_mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")
quad_mesh = "examples/meshes/quad_cylinder"
coords,topology,entities = read_comsol_mesh(quad_mesh,QuadrilateralQuadratic9(3,3))
quadMesh = load3dQuadComsolMesh(quad_mesh)
quadMeshgl = load3dQuadComsolMesh(quad_mesh;geometryType="quadlinear")
quadMeshpl = load3dQuadComsolMesh(quad_mesh;physicsType="linear")

quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadconstant")
quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadlinear")
quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadquadratic")

@time IntegralEquations.interpolate_elements(quadMesh);
TMP = IntegralEquations.interpolate_elements(quadMesh)

Ql = QuadrilateralLinear(3,3)
IntegralEquations.set_nodal_interpolation!(Ql)
Ql4 = QuadrilateralLinear4(3,3)
IntegralEquations.set_nodal_interpolation!(Ql4)

Qq = QuadrilateralQuadratic(3,3)
IntegralEquations.set_nodal_interpolation!(Qq)
Qq9 = QuadrilateralQuadratic9(3,3)
IntegralEquations.set_nodal_interpolation!(Qq9)


normals = mesh.normals
tangentX = similar(normals)
tangentY = similar(normals)
IntegralEquations.tangents!(normals,tangentX,tangentY)



#==========================================================================================
                                        Plotting
==========================================================================================#

# using IntegralEquations
# using Meshes
# using MeshViz
# import Makie as Mke

# mesh_file = "examples/meshes/resonatorFinal"
# tri_mesh = load3dTriangularComsolMesh(mesh_file)
# # plot_mesh(tri_mesh)

# quad_mesh_file = "examples/meshes/quad_cylinder"
# quad_mesh = load3dQuadComsolMesh(quad_mesh_file)

# SP_lin  = create_simple_mesh(tri_mesh)
# viz(SP_lin;showfacets=true)
# # SP_quad = create_simple_mesh(quad_mesh)
# bc_ents = [5]
# SP_quad = create_bc_simple_mesh(quad_mesh_file,quad_mesh,quad_mesh.shape_function,bc_ents,false)
# SP_bc   = create_bc_simple_mesh(quad_mesh_file,quad_mesh,quad_mesh.shape_function,bc_ents)
# fig1, axis1, plot1 = viz(SP_quad;showfacets=true)
# Mke.scatter!(axis1, [0.525], [0.00], [0.00],
#     color = :red,markersize = 10.0,
#     )
# viz(SP_quad;showfacets=true)
# viz!(SP_bc;showfacets=true,color=:red)

# _,topo,ents = read_comsol_mesh(quad_mesh_file,quad_mesh.shape_function)
# # .!Bool.(sum(bc_ents .∈ ents,dims=1))[:]

# verts = vertices(SP_lin)
# connec = SP_lin.topology.connec
# d1 = rand(length(verts))
# d2 = tri_mesh.coordinates[tri_mesh.topology[2,:]]
# md = meshdata(SP_lin, 
#     Dict(0 => (temperature = d1,),
#          2 => (pressure = d2,))
# )
# viz(md)


#==========================================================================================
                                Assembly
==========================================================================================#
using IntegralEquations
using LinearAlgebra
import IntegralEquations: SurfaceElement, jacobian!, interpolate_on_nodes!, interpolate_elements
import IntegralEquations: assemble_parallel!, create_rotated_element, computing_integrals!
import IntegralEquations: freens3d!, greens3d!,freens3dk0!, computing_integrand!


mesh_file = "examples/meshes/resonatorFinal"
tri_mesh = load3dTriangularComsolMesh(mesh_file)
# Fl,Gl,Cl = assemble_parallel!(tri_mesh,1.0,tri_mesh.sources);

quad_mesh_file = "examples/meshes/quad_cylinder"
quad_mesh = load3dQuadComsolMesh(quad_mesh_file)
# quad_mesh = load3dQuadComsolMesh(quad_mesh_file;physicsType="linear")
# Fq,Gq,Cq = assemble_parallel!(quad_mesh,1.0,quad_mesh.sources;n=4,m=4);


import IntegralEquations: assemble_parallel_galerkin!
Fg,Gg,Cg = assemble_parallel_galerkin!(tri_mesh,1.0,tri_mesh.sources)



import IntegralEquations: number_of_elements, create_shape_function, interpolate_elements,
                            copy_interpolation_nodes!, computing_galerkin_integrals!
mesh = tri_mesh
insources = tri_mesh.sources
shape_function = tri_mesh.shape_function
ntest = 3
mtest = 3
nbasis = 4
mbasis = 4


nElements   = number_of_elements(mesh)
# nThreads    = nthreads()
# sources     = convert.(eltype(shape_function),insources)
nSources    = size(insources,2)
nNodes      = size(mesh.sources,2)
physicsTopology  = mesh.physics_topology
physics_function = mesh.physics_function
#======================================================================================

======================================================================================#
test_function       = create_shape_function(shape_function;n=ntest,m=mtest)
basis_function      = create_shape_function(shape_function;n=nbasis,m=mbasis)
test_interpolation  = interpolate_elements(mesh,test_function)
basis_interpolation = interpolate_elements(mesh,basis_function)
test_physics        = deepcopy(physics_function)
basis_physics       = deepcopy(physics_function)
copy_interpolation_nodes!(test_physics,test_function)
copy_interpolation_nodes!(basis_physics,basis_function)
# Avoiding to have gauss-node on a singularity. Should be handled differently,
# but, this is good for now.
# if typeof(shape_function) <: QuadrilateralQuadratic9
#     for (x,y) in zip(test_function.gauss_u,shape_function1.gauss_v)
#         if isapprox.(x, 0.0, atol=1e-15) && isapprox.(y, 0.0, atol=1e-15)
#             error("Gauss Node Singularity.")
#         end
#     end
# end
#======================================================================================
                    Preallocation of return values & Intermediate values
======================================================================================#
# 
F = zeros(ComplexF64, nSources, nNodes)
G = zeros(ComplexF64, nSources, nNodes)
C = zeros(ComplexF64, nSources)

# Preallocation according to the number of threads
#======================================================================================
                                Assembly
======================================================================================#
@inbounds for test_element = 1:nElements
    # Access source
    # Every thread has access to parts of the pre-allocated matrices
    integrand = zeros(ComplexF64, ntest*mtest, nbasis*nbasis)
    r         = zeros(            ntest*mtest, nbasis*nbasis)
    test_nodes = @view physicsTopology[:,test_element]
    @inbounds for basis_element = 1:nElements
        if basis_element == test_element
            continue
        end
        basis_nodes = @view physicsTopology[:,basis_element]
        # Acces submatrix of the BEM matrix
        submatrixF = @view F[test_nodes,basis_nodes]
        submatrixG = @view G[test_nodes,basis_nodes]
        subvectorC = @view C[test_nodes]
        # Interpolating on the mesh. 
        # test_physics,test_interpolation,
        #  basis_physics,basis_interpolation,
        #  submatrixF,submatrixG,subvectorC,k,integrand,r)
        computing_galerkin_integrals!(test_physics,test_interpolation[test_element],
                                        basis_physics,basis_interpolation[basis_element],
                                        submatrixF,submatrixG,subvectorC,k,integrand,r)
        @assert all(.!isnan.(r))
        @assert all(.!isnan.(submatrixG))
        @assert all(.!isnan.(submatrixF))
    end
    println(test_element)
end

1.0 + 1.0

# gx = 16
# gy = 16
# nx = 9
# ny = 9
# wx = ones(gx)
# Nx = ones(gx,nx)'
# wy = ones(gy)
# Ny = ones(gy,ny)'
# Gxy = ones(gy,gx)
# @time sum(Ny[:,i]*(wx'*(Gxy[i,:] .* Nx')) for i = 1:gy);
# function summing!(Nx,wx,Ny,wy,Gxy,tmp)
#     for i = 1:length(wy)
#         tmp .+= Ny[:,i]*(wx'*(Gxy[i,:] .* Nx'))
#     end
# end
# tmp = zeros(ny,nx)
# @time summing!(Nx,wx,Ny,wy,Gxy,tmp)
# tmp - sum(Ny[:,i]*(wx'*(Gxy[i,:] .* Nx')) for i = 1:gy)
# @time sum(Ny[:,i]*(wx'*(Gxy[i,:] .* Nx')) for i = 1:gy);

# @inbounds for i = 1:size(Integrand,1),j = 1:size(Integrand,2)
#     println((i,j))
# end


# import IntegralEquations: compute_distances!
# nsources = 3
# ngauss = 16
# r = zeros(nsources,ngauss)
# interpolation = ones(3,ngauss) .* collect(range(0.0,1.0,length=ngauss))'
# sources = zeros(3,nsources)
# compute_distances!(r,interpolation,sources)
# r

# r = ones(1,size(integrand,2))
# greens3d!(integrand,r,1.0)
# R = ones(size(Integrand))
# greens3d!(Integrand,R,1.0)

# normals = ones(3,length(r))
# sources = ones(3,size(Integrand,1))
# interpolation = ones(3,length(r))
# freens3d!(Integrand,R,interpolation,sources,normals,1.0)
# Integrand

