using IntegralEquations
using LinearAlgebra
using Plots
using Meshes
using MeshViz
using GLMakie
set_theme!(resolution=(1200, 1200))

mesh_file = "examples/meshes/cylinder"
mesh = load3dTriangularComsolMesh(mesh_file)
simple_mesh = create_simple_mesh(mesh)
viz(simple_mesh;showfacets=true)


mesh = load3dTriangularComsolMesh(mesh_file;geometryType="TriLinear")
mesh = load3dTriangularComsolMesh(mesh_file;physicsType="linear")

mesh = load3dTriangularComsolMesh(mesh_file;physicsType="disctriconstant")
mesh = load3dTriangularComsolMesh(mesh_file;physicsType="disctrilinear")
mesh = load3dTriangularComsolMesh(mesh_file;physicsType="disctriquadratic")


simple_mesh = create_simple_mesh(mesh)
viz(simple_mesh;showfacets=true)

element_interpolations = IntegralEquations.interpolate_elements(mesh)

S = 0.0
for i = 1:length(element_interpolations)
    S += dot(element_interpolations[i].jacobians,mesh.shape_function.weights)
end
S

@time IntegralEquations.interpolate_elements(mesh);


# quad_mesh = load3dTriangularComsolMesh(mesh_file;geometryType="TriLinear")
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

using IntegralEquations
using Meshes
using MeshViz
using WGLMakie
set_theme!(resolution=(1200, 1200))

mesh_file = "examples/meshes/sphere_1m"
tri_mesh = load3dTriangularComsolMesh(mesh_file)

quad_mesh_file = "examples/meshes/quad_sphere"
quad_mesh = load3dQuadComsolMesh(quad_mesh_file)

SP_lin  = create_simple_mesh(tri_mesh)
viz(SP_lin;showfacets=true)
SP_quad = create_simple_mesh(quad_mesh)
viz(SP_quad;showfacets=true)
# SP_quad = create_simple_mesh(quad_mesh)
bc_ents = [0,1]
SP_quad = create_bc_simple_mesh(quad_mesh_file,quad_mesh,quad_mesh.shape_function,bc_ents,false)
SP_bc   = create_bc_simple_mesh(quad_mesh_file,quad_mesh,quad_mesh.shape_function,bc_ents)
viz(SP_quad;showfacets=true)
viz!(SP_bc;showfacets=true,color=:red)
current_figure()

x = quad_mesh.sources[1,:]
y = quad_mesh.sources[2,:]
z = quad_mesh.sources[3,:]
u = quad_mesh.normals[1,:]
v = quad_mesh.normals[2,:]
w = quad_mesh.normals[3,:]
arrows(x, y, z, u, v, w)


# _,topo,ents = read_comsol_mesh(quad_mesh_file,quad_mesh.shape_function)
# fig1, axis1, plot1 = viz(SP_quad;showfacets=true)
# scatter!(axis1, [0.525], [0.00], [0.00],
#     color = :red,markersize = 10.0,
#     )
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
import IntegralEquations: assemble_parallel_galerkin!

mesh_file = "examples/meshes/resonatorFinal"
tri_mesh = load3dTriangularComsolMesh(mesh_file)
# Fl,Gl,Cl = assemble_parallel!(tri_mesh,1.0,tri_mesh.sources);
# Fl,Gl,Cl = assemble_parallel_galerkin!(tri_mesh,1.0,tri_mesh.sources)

quad_mesh_file = "examples/meshes/quad_cylinder"
# quad_mesh = load3dQuadComsolMesh(quad_mesh_file)
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;physicsType="linear")
# Fq,Gq,Cq = assemble_parallel!(quad_mesh,1.0,quad_mesh.sources;n=4,m=4);
Fq,Gq,Cq = assemble_parallel_galerkin!(quad_mesh,1.0);

import IntegralEquations: compute_distances



gx = 10
gy = 16
nx = 6
ny = 9
wx = ones(gx)
Nx = rand(gx,nx)'
wy = ones(gy)
Ny = rand(gy,ny)'
Gxy = ones(gy,gx)
integrand = ones(ny,nx)
@time computing_integrand!(integrand,Gxy,Ny,wy,Nx,wx);


i = 1
T1 = wx'*(Gxy[i,:] .* Nx')
wx'*(Gxy[i,:] .* Nx')
(wx' .* Gxy[i,:]')*Nx'
T1 = wx' .* Gxy[i,:]'
T2 = wx'*(Gxy[i,:] .* Nx')

T1 = wx .* Gxy[i,:]
T2 = Nx*T1

function jacmul!(jac,g,w)
    @inbounds for i = 1:length(jac)
        jac[i] = g[i]*w[i]
    end
end
function computing_integrand_test!(integrand,Gxy,Ny,wy,Nx,wx,temporary1,temporary2)
    for i = 1:length(wy)
        # temporary1 .= Gxy[i,:] .* wx
        jacmul!(temporary1,Gxy[i,:],wx)
        mul!(temporary2,Nx,temporary1)
        mul!(integrand,Ny[:,i],temporary2',wy[i],true)
    end
end
# @time Ny[:,i]*(wx'*(Gxy[i,:] .* Nx'))*wy[i]
@time computing_integrand!(integrand,Gxy,Ny,wy,Nx,wx);
@time computing_integrand_test!(integrand,Gxy,Ny,wy,Nx,wx,T1,T2);

tm = deepcopy(Gxy[i,:])
@time jacmul!(T1,Gxy[i,:],wx)
@time mul!(T2,Nx,T1);
i = 1
@time mul!(integrand,Ny[:,i],T2',wy[i],true);


function computing_c_integrand!(integrand,Ny,wy)
    @inbounds for i = 1:length(wy)
        mul!(integrand,Ny[:,i],Ny[:,i]',wy[i],true)
    end
end
inte = Ny[:,i]*Ny[:,i]'
@time computing_c_integrand!(inte,Ny,wy);


tmp = zeros(ny,nx)
tmp1 = wx'* Nx'
tmp2 = Ny*wy
function computing_integrand2!(integrand,Ny,wy,Nx,wx,tmp1)
    tmp1 .= wx'* Nx'
    for i = 1:length(wy)
        # integrand .+= Ny[:,i]*wy[i]*(tmp1)
        mul!(integrand,Ny[:,i],tmp1,wy[i],true)
    end
end
function computing_integrand3!(integrand,Ny,wy,Nx,wx,tmp1,tmp2)
    mul!(tmp1,wx',Nx')
    mul!(tmp2,Ny,wy)
    mul!(integrand,tmp2,tmp1,true,true)
end
@time tmp .= Ny*wy * (wx'*Nx');
@time computing_integrand3!(tmp,Ny,wy,Nx,wx,tmp1,tmp2);
@time computing_integrand2!(tmp,Ny,wy,Nx,wx,tmp1);

# Gxy = ones(gy,gx)
# @time sum(Ny[:,i]*(wx'*(Gxy[i,:] .* Nx')) for i = 1:gy);
# function summing!(Nx,wx,Ny,wy,Gxy,tmp)
#     for i = 1:length(wy)
#         tmp .+= Ny[:,i]*(wx'*(Gxy[i,:] .* Nx'))
#     end
# end
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
