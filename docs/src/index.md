# BoundaryIntegralEquations.jl Documentation

`BoundaryIntegralEquations.jl` provides the basis functionalities for implementing Boundary Element Methods.

Currently, the package is focused on solving the Helmholtz equation for acoustical problems through the Kirchoff-Helmholtz integral Equation

```math
c(\mathbf{y})p(\mathbf{y}) + \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial \mathbf{n} }p(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x} = \mathrm{i}\rho ck\int_\Gamma G(\mathbf{x},\mathbf{y})v_s(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x},
```
where ``G(\mathbf{x},\mathbf{y})`` is the Green's function of the Helmholtz operator.

```@example introduction
using JSServe                       # hide
Page(exportable=true, offline=true) # hide
```

## Quick start guide
```@example introduction
import WGLMakie as Mke # hide
Mke.set_theme!(resolution=(800, 800)) # hide
```

A mesh can be loaded using the following
```@example introduction
using BoundaryIntegralEquations, MeshViz
examples_path = normpath(joinpath(@__DIR__, "..", "..", "examples"));
tri_mesh_file = joinpath(examples_path, "meshes","tri_sphere")
tri_mesh = load3dTriangularComsolMesh(tri_mesh_file;geometry_order=:linear)
```
Furthermore, mesh can be visualized by creating `SimpleMesh` that is compatible with the `MeshViz` library.
```@example introduction
simple_tri_mesh = create_simple_mesh(tri_mesh)
viz(simple_tri_mesh;showfacets=true)
```
A similar approach can be used to import and visualize a mesh of quadrilaterals.
```@example introduction
quad_mesh_file = joinpath(examples_path, "meshes","quad_sphere")
quad_mesh = load3dQuadComsolMesh(quad_mesh_file;geometry_order=:linear)
```
```@example introduction
simple_quad_mesh = create_simple_mesh(quad_mesh)
viz(simple_quad_mesh;showfacets=true)
```
