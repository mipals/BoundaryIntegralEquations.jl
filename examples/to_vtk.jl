#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using Plots
using WriteVTK
#==========================================================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
mesh_path = "/Users/mpasc/Documents/bigmeshes";
# mesh_file = joinpath(mesh_path,"sphere_136k");
# mesh_file = joinpath(mesh_path,"sphere_195k");
# mesh_file = joinpath(mesh_path,"sphere_305k");
# mesh_file = joinpath(mesh_path,"sphere_377k");
# mesh_file = joinpath(mesh_path,"sphere_478k");
# mesh_file = joinpath(mesh_path,"sphere_5m_119k");
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
mesh_file = joinpath(mesh_path,"sphere_1m_fine");
mesh_file = joinpath(mesh_path,"sphere_1m_finer");
mesh_file = joinpath(mesh_path,"sphere_1m_extremely_fine");
mesh_file = joinpath(mesh_path,"sphere_1m_finest");
mesh_file = joinpath(mesh_path,"sphere_1m_35k");
mesh_file = joinpath(mesh_path,"sphere_1m_77k");

@time mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=geometry_orders[2],
                                                       physics_order=tri_physics_orders[2])

#==========================================================================================
                                    3d Visualization
==========================================================================================#
# using MeshViz
# import WGLMakie as wgl
# simple_mesh = create_simple_mesh(mesh)
# wgl.set_theme!(resolution=(1200, 1200))
# viz(simple_mesh, showfacets = true)
#==========================================================================================
                                Setting up constants
==========================================================================================#
radius = 1.0                                     # Radius of sphere_1m       [m]
#==========================================================================================
                            Write to vtk
==========================================================================================#
xyzb = mesh.sources
M  = size(xyzb,2)
normals  = mesh.normals
tangent1 = mesh.tangents
tangent2 = mesh.sangents

# Export geometry mesh and nodal normals/tangentials to paraview format (https://juliahub.com/ui/Packages/WriteVTK/D2v2J/1.12.0)
numElem = size(mesh.topology,2)
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, mesh.topology[:,i]) for i=1:numElem] # array that contains all the cells of the mesh
vtkfile = vtk_grid("/Users/mpasc/Documents/bigmeshes/sphere_$(M)DOFs", mesh.coordinates, cells)
vtkfile["nVec",VTKPointData()] = normals
vtkfile["tVec",VTKPointData()] = tangent1
vtkfile["sVec",VTKPointData()] = tangent2
outfiles = vtk_save(vtkfile)


DOFS = [13.062 ; 39.486 ; 69.038 ; 153.978 ; 271.778 ; 390.890 ; 610.058 ; 753.778 ; 955.706]*1000
2*(DOFS.^2)*8*2/(2^30)

# Export geometry error
# r = sqrt.(xyzb[1,:].^2+xyzb[2,:].^2+xyzb[3,:].^2)
# vtkfile["RadiusError", VTKPointData()] = abs.(r.-radius)
# outfiles = vtk_save(vtkfile)


# #==========================================================================================
#                             Calculate and plot geometry error
# ==========================================================================================#
# function errorCalc(calc,ref,M)

#     eps = (sum(abs.(calc.-ref).^2)/M).^(1/2)
#     refNorm = (sum(abs.(ref).^2)/M).^(1/2);
#     epsRel = eps/refNorm;
#     return epsRel
# end

# ang_axis = acos.(xyzb[3,:]/radius)*180.0/pi
# perm = sortperm(ang_axis)
# epsRel = errorCalc(r,radius,2)
# scatter(ang_axis[perm],(radius.-r)[perm],label="p2 $(M)DOFs epsRel=$(round(epsRel; digits = 7))",marker=:cross,markersize=2,color=:black,dpi=400); # for quadratic physics discretization
# title!("Geometry error");
# xlabel!("Angle [deg]");
# ylabel!("r_ref - r_num [m]")
# savefig("geometryError.png")
# #==========================================================================================
#                             Compare geometry error to other mesh
# ==========================================================================================#
# # Triangular Meshes
# # tri_mesh_file = "examples/meshes/sphere_1m_coarser"
# # tri_mesh_file = "examples/meshes/sphere_1m_coarse"
# # tri_mesh_file = "examples/meshes/sphere_1m"
# # tri_mesh_file = "examples/meshes/sphere_1m_fine"
# # tri_mesh_file = "examples/meshes/sphere_1m_finer"
# tri_mesh_file = "examples/meshes/sphere_1m_extremely_fine"
# # tri_mesh_file = "examples/meshes/sphere_1m_4p5k"
# # tri_mesh_file = "examples/meshes/sphere_1m_35k"
# # tri_mesh_file = "examples/meshes/sphere_1m_77k"
# @time mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=geometry_orders[2],
#                                                        physics_order=tri_physics_orders[1])
# xyzb = mesh.sources
# M2  = size(xyzb,2)

# r2 = sqrt.(xyzb[1,:].^2+xyzb[2,:].^2+xyzb[3,:].^2)
# ang_axis2 = acos.(xyzb[3,:]/radius)*180.0/pi
# perm2 = sortperm(ang_axis2)
# epsRel2 = errorCalc(r2,radius,2)

# scatter(ang_axis[perm],(radius.-r)[perm],label="p2 $(M)DOFs epsRel=$(round(epsRel; digits = 7))",marker=:cross,markersize=2,color=:black,dpi=400);
# scatter!(ang_axis2[perm2],(radius.-r2)[perm2],label="p1 $(M2)DOFs epsRel=$(round(epsRel2; digits = 7))",marker=:cross,markersize=2,color=:red);
# title!("Geometry error");
# xlabel!("Angle [deg]");
# ylabel!("r_ref - r_num [m]")
# savefig("geometryErrorComp.png")
