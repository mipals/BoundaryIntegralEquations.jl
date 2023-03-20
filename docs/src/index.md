# BoundaryIntegralEquations.jl

`BoundaryIntegralEquations.jl` provides the basis functionalities for implementing Boundary Element Methods.

Currently, the package is focused on solving the Helmholtz equation for acoustical problems through the Kirchoff-Helmholtz integral Equation

```math
c(\mathbf{y})p(\mathbf{y}) + \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial \mathbf{n} }p(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x} = \mathrm{i}\rho ck\int_\Gamma G(\mathbf{x},\mathbf{y})v_s(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x},
```
where ``G(\mathbf{x},\mathbf{y})`` is the Green's function of the Helmholtz operator.

For the Fast Multipole Method this package utilizes the Julia interfaces for the Flatiron Institute Fast Multipole Libraries: [2D](https://github.com/mipals/FMM2D.jl), [3D](https://github.com/flatironinstitute/FMM3D/tree/master/julia).

For H-matrices the package utilizes the [HMatrices.jl](https://github.com/WaveProp/HMatrices.jl).

**N.B. The package is still under heavy development.**

## Installation
The package can be downloaded directly from GitHub 

```julia
using Pkg
Pkg.add(url="https://github.com/mipals/BoundaryIntegralEquations.jl")
```

## Element types
* (Dis)continuous (Constant, Linear and Quadratic) Line elements
* (Dis)continuous (Constant, Linear and Quadratic) Triangular Elements
* (Dis)continuous (Constant, Linear and Quadratic) Quadrilateral Elements

## Mesh formats
* COMSOLs *.mphtxt* files (best for applying boundary conditions)
* .obj, .ply, .stl, .off, .2DM through [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl).

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
