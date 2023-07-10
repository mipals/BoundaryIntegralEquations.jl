# BoundaryIntegralEquations.jl

| **Documentation** | **Tests** | **CodeCov** | **Lifecycle** | **Aqua** |
|:-----------------:|:---------:|:-----------:|:-------------:|:--------:|
|[![](https://img.shields.io/badge/docs-online-blue.svg)](https://mipals.github.io/BoundaryIntegralEquations.jl/dev/)| [![Build Status](https://github.com/mipals/BoundaryIntegralEquations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/BoundaryIntegralEquations.jl/actions/workflows/CI.yml?query=branch%3Amain) | [![Coverage](https://codecov.io/gh/mipals/BoundaryIntegralEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/BoundaryIntegralEquations.jl)| ![](https://img.shields.io/badge/Lifecycle-Unstable-yellow)| [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) |


BoundaryIntegralEquations.jl provides the basic building blocks required for implementing Boundary Element Methods (BEMs). Currently, it supplies the discretization of the Kirchhoff–Helmholtz integral equation found in acoustical applications.

The package supplies support for acceleration methods such as e.g. the Fast Multipole Method (FMM) and Hierarchical matrices. For the FMM the package utilizes the Flatiron Institute Fast Multipole Libraries ([2D](https://github.com/mipals/FMM2D.jl), [3D](https://github.com/flatironinstitute/FMM3D/tree/master/julia)) whereas for Hierarchical matrices it uses the [HMatrices.jl](https://github.com/WaveProp/HMatrices.jl) package.

**N.B. The package is still under heavy development.**

## Installation
The package can be installed directly from GitHub 

```julia
using Pkg
Pkg.add(url="https://github.com/mipals/BoundaryIntegralEquations.jl")
```

## Mesh formats
* COMSOLs *.mphtxt* files (including entity description)
* `obj`, `stl`, `ply`, `off` and `2DM` through [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl).

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
