# Fast Operators

## Helper Functions
```@docs
BoundaryIntegralEquations.partial_assemble_parallel!
```

```@docs
BoundaryIntegralEquations.unroll_interpolations
```

```@docs
BoundaryIntegralEquations.create_coefficient_map
```

```@docs
BoundaryIntegralEquations.scale_columns!
```

```@docs
BoundaryIntegralEquations.setup_fast_operator
```

```@docs
BoundaryIntegralEquations.create_single_layer_matrix
```

```@docs
BoundaryIntegralEquations.evaluate_targets(A::FMMFOperator,x,targets)
```

```@docs
BoundaryIntegralEquations.evaluate_targets(A::FMMFOperator,k,x,targets)
```

```@docs
BoundaryIntegralEquations.evaluate_targets(A::FMMGOperator,x,targets)
```

```@docs
BoundaryIntegralEquations.evaluate_targets(A::FMMGOperator,k,x,targets)
```

## Fast Multipole Operators

```@docs
BoundaryIntegralEquations.FMMGOperator
```

```@docs
BoundaryIntegralEquations.FMMHOperator
```

## H-Matrix Operators
```@docs
BoundaryIntegralEquations.HGOperator
```

```@docs
BoundaryIntegralEquations.HFOperator
```
