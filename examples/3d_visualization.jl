# # Visualization (3D)
# # Importing packages
using BoundaryIntegralEquations
using Meshes
import WGLMakie as wgl # WGLMakie integrates into VSCode. Other backends can also be used.
wgl.set_theme!(resolution=(800, 800))

using JSServe                       #hide
Page(exportable=true, offline=true) #hide

# # Visualizing full mesh
# First we start with a triangular Mesh
tri_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),
                            "..","examples","meshes","cylinder");
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=:linear)
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)

# A similar approach can be done for quadrilateral meshes
quad_mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),
                            "..","examples","meshes","quad_cylinder")
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:linear)
simple_quad_mesh = create_simple_mesh(quad_mesh)
viz(simple_quad_mesh;showfacets=true)

# # Visualization of mesh with boundary conditions
# Again we start with a triangluar mesh. First the Boundary entities is defined
tri_bc_ents = [6,7,10] .- 1; # The .-1 is due to COMSOL 0-indexing of exported entities
# Creating simple meshes for full mesh and boundary condition
tri_simple_mesh = create_bc_simple_mesh(tri_mesh,tri_bc_ents,false);
tri_simple_bc   = create_bc_simple_mesh(tri_mesh,tri_bc_ents);
# Plotting entities in `bc_ents` red
viz(tri_simple_mesh;showfacets=true)
viz!(tri_simple_bc;showfacets=true,color=:red)
wgl.current_figure()

# And similarly for a quadrilateral mesh. First we define the entites for boundaries
quad_bc_ents = [1,2];
# Creating simple meshes for full mesh and boundary condition
quad_simple_mesh = create_bc_simple_mesh(quad_mesh,quad_bc_ents,false);
quad_simple_bc   = create_bc_simple_mesh(quad_mesh,quad_bc_ents);
# Plotting entities in `bc_ents` red
viz(quad_simple_mesh;showfacets=true)
viz!(quad_simple_bc;showfacets=true,color=:red)
wgl.current_figure()
