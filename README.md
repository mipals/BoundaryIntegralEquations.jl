# <img src="integral_equations.png" alt="drawing" width="400"/>

[![Build Status](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mipals/IntegralEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/IntegralEquations.jl)

The IntegralEquations.jl can be used to solve integral equations through the Boundary Element Methods (FEM). Currently, it only supplies the discretization of the Kirchhoff–Helmholtz integral equation found in acoustical applications,

$$
c(\mathbf{y})p(\mathbf{y}) + \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial\mathbf{n}_\mathbf{x}}p(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x} = \mathrm{i}\rho ck\int_\Gamma G(\mathbf{x},\mathbf{y})v_s(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x}.
$$

However, the package provides an easy-to-use API to the building blocks required for discretizing and solving general (surface) integral equations using the techniques of the Boundary Element Method (BEM). All the user has to define is a different integrand. 

## Supported Element Types
* Continuous (Linear and Quadratic) Triangular Elements
* Continuous (Linear and Quadratic) Quadrilateral Elements
* Discontinuous (Constant, Linear and Quadratic) Triangular Elements
* Discontinuous (Constant, Linear and Quadratic) Quadrilateral Elements

## Supported meshes
* COMSOLs *.mphtxt* files. 

## Roadmap
* More basis functions
    - [ ] Higher order basis functions. Fairly low-hanging fruit, but could be problematic due to on-element-singularities.
    - [ ] Non-interpolatory basis such as [Higher-order Legendre basis functions.](https://ieeexplore.ieee.org/document/1353496) (seems to only be easily implemented for quadrilaterals).
* Support for more mesh files. 
    - [ ] Gmesh files.
* Precision Reduction. Could be troublesome because of the singularities.
* Better singularity handling.
    - [ ] SauterSchwab integration (Required for the full implementation of Galerkin assembly).
    - [ ] Taylor series expansion.
    - [ ] Handling hypersingular integrals correctly (efficiently).
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
