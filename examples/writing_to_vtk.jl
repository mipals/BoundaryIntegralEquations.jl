#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using IntegralEquations
#==========================================================================================
                                Loading Mesh
==========================================================================================#
geometry_orders     = [:linear,:quadratic]
tri_physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
# Triangular Meshes
mesh_file = "multipleResonatorsFine"
# mesh_file = "resonatorFinalFine"
@time mesh = load3dTriangularComsolMesh("examples/meshes/" * mesh_file;
                                            geometry_order=geometry_orders[2],
                                            physics_order=tri_physics_orders[2])
#==========================================================================================
                                    3d Visualization
==========================================================================================#
using MeshViz
import WGLMakie as wgl
simple_mesh = create_simple_mesh(mesh)
wgl.set_theme!(resolution=(800, 800))
viz(simple_mesh, showfacets = true)
#==========================================================================================
                                    Write to VTK
==========================================================================================#
using WriteVTK
points = mesh.coordinates
cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE, con) for con in eachcol(mesh.topology)]

vtk_grid("examples/meshes/vtk_outputs/" * mesh_file, points, cells) do vtk
    # add datasets...
end


x1 = [0.0; 0.2; 0.0]
x2 = [0.0; 0.3; 0.0]
x = [x1[1]; x2[1]]
y = [x1[2]; x2[2]]
z = [x1[3]; x2[3]]
vtk_grid("examples/meshes/vtk_outputs/points",x,y,z) do vtk
end
