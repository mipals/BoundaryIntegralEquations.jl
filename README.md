# <img src="integral_equations.png" alt="drawing" width="400"/>

[![Build Status](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/IntegralEquations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mipals/IntegralEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/IntegralEquations.jl)

The IntegralEquations.jl provides the basic building blocks for the Boundary Element Method (BEM). Currently, it only supplies the discretization of the Kirchhoff–Helmholtz integral equation found in acoustical applications

$$
c(\mathbf{y})p(\mathbf{y}) + \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{y})}{\partial\mathbf{n}_\mathbf{x}}p(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x} = \mathrm{i}\rho ck\int_\Gamma G(\mathbf{x},\mathbf{y})v_s(\mathbf{x})\ \mathrm{d}\Gamma_\mathbf{x}.
$$

For the Fast Multipole Method this package utilizes the Julia interfaces for the Flatiron Institute Fast Multipole Libraries: [2D](https://github.com/mipals/FMM2D.jl), [3D](https://github.com/flatironinstitute/FMM3D/tree/master/julia).

## Supported Element Types
* (Dis)continuous (Constant, Linear and Quadratic) Line elements
* (Dis)continuous (Constant, Linear and Quadratic) Triangular Elements
* (Dis)continuous (Constant, Linear and Quadratic) Quadrilateral Elements

## Supported meshes
* COMSOLs *.mphtxt* files. 

## Roadmap / Notes
* Systematic way of applying boundary conditions
* Better support for quadrilaterals
    - [ ] Adaptive integration for lossy formulation.
* A performant reading of COMSOL mesh files.
* Support for more mesh files. 
    - [ ] Gmsh files.
    - [ ] Ply files through [PlyIO.jl](https://github.com/JuliaGeometry/PlyIO.jl).
    - [ ] Any many more...
* Support for single precision numbers. ([Inspired by Bempp-cl](https://www.mscroggs.co.uk/papers/2021-cise.pdf), generally needed for GPU)
* More basis functions
    - [ ] Higher order basis functions. Fairly low-hanging fruit, but could be problematic due to on-element-singularities.
    - [ ] Non-interpolatory basis such as [Higher-order Legendre basis functions.](https://ieeexplore.ieee.org/document/1353496) (seems to only be easily implemented for quadrilaterals).
* Better singularity handling.
* Better post-processing 
    - [ ] Create data structure to save results.
    - [ ] Automatically plot surface pressure/velocities.
* Galerkin support
    - [ ] Requires SauterSchwab integration. The package [SauterSchwabQuadrature.jl](https://github.com/ga96tik/SauterSchwabQuadrature.jl) is used by BEAST.jl. It might be easy to integrate here.
* Other integration techniques
    - [ ] [SparseGrids.jl](https://github.com/robertdj/SparseGrids.jl) / [DistributedSparseGrids.jl](https://github.com/baxmittens/DistributedSparseGrids.jl)
    - [ ] Adaptive Hierarchical Approaches

## Similar Packages
* [BEAST.jl](https://github.com/krcools/BEAST.jl): Boundary Element Analysis and Simulation Toolkit. A general toolkit, with a focus on electromagnetics. Limitations with respect to element orders and only supplies Galerkin assembly. 
* [NESSie.jl](https://github.com/tkemmer/NESSie.jl): Nonlocal Electrostatics in Structured Solvents. A specialized package written specifically for Nonlocal protein eletrostatics. 
