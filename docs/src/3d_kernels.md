# Kernels
```@docs
BoundaryIntegralEquations.greens3d!
```

```@docs
BoundaryIntegralEquations.freens3d!
```

```@docs
BoundaryIntegralEquations.freens3dk0!
```

```@docs
BoundaryIntegralEquations.onefunction!
```

## Derivatives for ROSEBEM
```@docs
BoundaryIntegralEquations.taylor_greens3d!(integrand,r,k0,m)
```
```@docs
BoundaryIntegralEquations.taylor_freens3d!(integrand,r,sources,collocation,normals,k0,m)
```
```@docs
BoundaryIntegralEquations.taylor_greens_gradient3d!(gradient,r,collocation,sources,k0,m)
```
```@docs
BoundaryIntegralEquations.taylor_greens_tangential_gradient3d!(gradient,r,collocation,sources,normals,k0,m)
```
