# ROSEBEM
(Reduced Order Series Expansion Boundary Element Method)

```@docs
BoundaryIntegralEquations.compute_taylor_integrals!
```

```@docs
BoundaryIntegralEquations.taylor_assemble!
```

```@docs
BoundaryIntegralEquations.apply_taylor_expansion(Abasis,k,k0)
```

```@docs
BoundaryIntegralEquations.apply_taylor_expansion(Abasis,bbasis,k,k0)
```

```@docs
BoundaryIntegralEquations.arnoldi_basis(A,b,q)
```

```@docs
BoundaryIntegralEquations.scattering_krylov_basis(mesh,klist;eps=1-4,n_gauss=3,verbose=true,P₀=1,progress=true)
```