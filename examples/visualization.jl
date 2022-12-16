#==========================================================================================
                                Visualization of full mesh
==========================================================================================#
using IntegralEquations
using MeshViz
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(1000, 1000))

# Triangular Mesh
tri_mesh_file = "examples/meshes/calibrationplate"
# tri_mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=:linear)
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=:quadratic)
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)

# Quadrilateral Mesh
quad_mesh_file = "examples/meshes/quad_cylinder"
# quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:linear)
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:quadratic)
simple_quad_mesh = create_simple_mesh(quad_mesh)
viz(simple_quad_mesh;showfacets=true)

#==========================================================================================
                                Visualization of boundary conditions
==========================================================================================#
using IntegralEquations
using MeshViz
import WGLMakie as wgl

# Triangles
tri_bc_ents = [6,7,10] .- 1
# Loading triangular mesh
tri_mesh_file = "examples/meshes/calibrationplate"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file)
# Plotting entities in `bc_ents` red
tri_simple_mesh = create_bc_simple_mesh(tri_mesh,tri_bc_ents,false)
viz(tri_simple_mesh;showfacets=true)
tri_simple_bc   = create_bc_simple_mesh(tri_mesh,tri_bc_ents)
viz(tri_simple_bc;showfacets=true,color=:red)

# Quadrilaterals
quad_bc_ents = [6,7]
# Loading quadrilateral mesh
quad_mesh_file = "examples/meshes/quad_sphere"
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:quadratic,physics_order=:linear)
# Plotting entities in `bc_ents` red
quad_simple_mesh = create_bc_simple_mesh(quad_mesh,quad_bc_ents,false)
viz(quad_simple_mesh;showfacets=true)
quad_simple_bc   = create_bc_simple_mesh(quad_mesh,quad_bc_ents)
viz!(quad_simple_bc;showfacets=true,color=:red)
# current_figure()

#==========================================================================================
                                Visualization of boundary conditions
==========================================================================================#
using IntegralEquations
using MeshViz
import WGLMakie as wgl

# Loading triangular mesh
tri_mesh_file = "examples/meshes/meta_boundary"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=:quadratic)
simple_tri_mesh = create_simple_mesh(tri_mesh)
wgl.set_theme!(resolution=(600, 600))
viz(simple_tri_mesh;showfacets=true)
