# BoundaryIntegralEquations.jl

`BoundaryIntegralEquations.jl` provides the basis functionalities for implementing Boundary Element Methods. Currently, the package is focused on solving the Helmholtz equation for acoustical problems through the Kirchoff-Helmholtz integral Equation

```math
\zeta(\mathbf{x})p(\mathbf{x}) + \int_\Gamma \frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial \mathbf{n}(\mathbf{y})}p(\mathbf{y})\ \mathrm{d}S_\mathbf{y} -
    \mathrm{i} \rho_0 c k \int_\Gamma G(\mathbf{x},\mathbf{y})v_\mathbf{n}(\mathbf{y})\ \mathrm{d}S_\mathbf{y} = 0.
```
where 
```math 
    G(\mathbf{x},\mathbf{y}) = \frac{\exp\left(\mathrm{i}k\|\mathbf{x} - \mathbf{y}\|_2\right)}{4\pi\|\mathbf{x} - \mathbf{y}\|_2},
``` 
is the Green's function of the Helmholtz operator, ``k`` is the wavenumber, ``\rho_0`` is the ambient density and ``c`` is the speed of sound.

## Mesh formats
* COMSOLs *.mphtxt* files (best for applying boundary conditions)
* .obj, .ply, .stl, .off, .2DM through [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl).

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
