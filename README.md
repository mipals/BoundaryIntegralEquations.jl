# IntegralEquations.jl

[![Build Status](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mipals/IntegralEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/IntegralEquations.jl)

This is a package provides the building blocks required for discretizing and solving (surface) integral equations using the techniques of the Boundary Element Method (BEM). The focus so far has been on the solution of the Kirchhoff–Helmholtz integral equation found in acoustical applications. However, as the package provides an interface to various element is should be easy to extended to handle all kinds of 

## Supported Element Types
* Continuous (Linear and Quadratic) Triangular Elements
* Continuous (Linear and Quadratic) Quadrilateral Elements
* Discontinuous (Constant, Linear and Quadratic) Triangular Elements
* Discontinuous (Constant, Linear and Quadratic) Quadrilateral Elements

## Supported meshes
* COMSOLs *.mphtxt* files. 

## Roadmap
* More basis functions
    - [ ] Higher order basis functions. 
    - [ ] Legendre basis (also high order).
* Support for more mesh files. 
    - [ ] Gmesh files.
* Better support for Boundary Conditions 
    - [ ] (Normal) Velocity Conditions for continuous elements.
    - [ ] Impedance Conditions.
    - [ ] Boundary Layer Impedance Conditions.
* Better post-processing 
    - [ ] Create data structure to save results.
    - [ ] Automatically plot surface pressure/velocities.

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
