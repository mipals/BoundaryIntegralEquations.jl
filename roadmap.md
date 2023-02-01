# Roadmap / Ideas / Issues
* Systematic way of applying boundary conditions
* Better support for quadrilaterals
    - [ ] Adaptive integration for lossy formulation.
* A performant reading of COMSOL mesh files.
* Support for more mesh files. 
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