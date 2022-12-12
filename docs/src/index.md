# IntegralEquations.jl Documentation

`IntegralEquations.jl` provides the basis functionalities of boundary element methods.

The package is focused on solving the Helmholtz equation for acoustical problems. 

```@example introduction
using JSServe # hide
Page(exportable=true, offline=true) # hide
```

## Quick start guide
```@example introduction
import WGLMakie as Mke # hide
Mke.set_theme!(resolution=(800, 800)) # hide
```

```@example introduction
using IntegralEquations, MeshViz
examples_path = normpath(joinpath(@__DIR__, "..", "..", "examples"));
tri_mesh_file = examples_path * "/meshes/tri_sphere"
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=:linear)
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)
```

```@example introduction
quad_mesh_file = examples_path * "/meshes/quad_sphere"
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:linear)
simple_quad_mesh = create_simple_mesh(quad_mesh)
viz(simple_quad_mesh;showfacets=true)
```



