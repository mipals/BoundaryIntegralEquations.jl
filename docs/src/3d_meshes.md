# Meshes in 3D

## Loading meshes
```@docs
BoundaryIntegralEquations.load_mesh_file
```

```@docs
BoundaryIntegralEquations.read_comsol_mesh
```

```@docs
BoundaryIntegralEquations.load3dTriangularMesh
```

```@docs
BoundaryIntegralEquations.load3dTriangularComsolMesh
```

```@docs
BoundaryIntegralEquations.load3dQuadComsolMesh
```

## Mesh manipulation

```@docs
BoundaryIntegralEquations.Mesh3d
```

```@docs
BoundaryIntegralEquations.remove_unused_nodes
```

```@docs
BoundaryIntegralEquations.normalize_vector!
```

```@docs
BoundaryIntegralEquations.tangents!
```

```@docs
BoundaryIntegralEquations.get_element_normals
```

```@docs
BoundaryIntegralEquations.compute_sources
```

```@docs
BoundaryIntegralEquations.set_physics_element
```

```@docs
BoundaryIntegralEquations.cross!
```

```@docs
BoundaryIntegralEquations.cross_product!
```

```@docs
BoundaryIntegralEquations.column_norms!
```

```@docs
BoundaryIntegralEquations.normalize!
```

```@docs
BoundaryIntegralEquations.jacobian!(basisElement::SurfaceFunction,coordinates,normals,tangent,sangent,jacobian)
```
