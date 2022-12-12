#==========================================================================================
                            Defining Various Element (Sub)Types
==========================================================================================#
abstract type SurfaceFunction            <: ShapeFunction               end
abstract type Triangular                 <: SurfaceFunction             end
abstract type Quadrilateral              <: SurfaceFunction             end
abstract type ContinuousTriangular       <: Triangular                  end
abstract type DiscontinuousTriangular    <: Triangular                  end
abstract type ContinuousQuadrilateral    <: Quadrilateral               end
abstract type DiscontinuousQuadrilateral <: Quadrilateral               end
abstract type QuadrilateralLagrange      <: ContinuousQuadrilateral     end
abstract type QuadrilateralSerendipity   <: ContinuousQuadrilateral     end
#==========================================================================================
                                Triangular Elements
==========================================================================================#
"""
Defines a Triangle with node placements         \n
                                          3     \n
                                          1 2   \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
1 - u - v   \\newline
    u       \\newline
    v
\\end{bmatrix},\\quad u \\in [0,1], v\\in [0,1-u].
```
"""
mutable struct TriangularLinear{T<:AbstractFloat} <: ContinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Triangle with node placements             \n
                                          3         \n
                                          6  5      \n
                                          1  4  2   \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
          (1 - v - u)(1 - 2v - 2u)   \\newline
            u(2u - 1)                \\newline
            v(2v - 1)                \\newline
            4u(1 - v - u)            \\newline
            4uv                      \\newline
            4v(1 - v - u)
\\end{bmatrix},\\quad u \\in [0,1], v\\in [0,1-u].
```
"""
mutable struct TriangularQuadratic{T<:AbstractFloat} <: ContinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct TriangularCubic{T<:AbstractFloat} <: ContinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Triangle with node placement on the middle of the triangle.

```math
\\mathbf{N}(u,v) = [1]
```
"""
mutable struct DiscontinuousTriangularConstant{T<:AbstractFloat} <: DiscontinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct DiscontinuousTriangularLinear{T<:AbstractFloat} <: DiscontinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    beta::T
end
mutable struct DiscontinuousTriangularQuadratic{T<:AbstractFloat} <: DiscontinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    beta::T
end
#==========================================================================================
                                Quadrilateral Elements
==========================================================================================#
"""
Defines a Quadrilateral with node placements         \n
                                            4 -- 3   \n
                                            |    |   \n
                                            1 -- 2   \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
          (1 - u)(1 - v)/4 \\newline
          (1 + u)(1 - v)/4 \\newline
          (1 + u)(1 + v)/4 \\newline
          (1 - u)(1 + v)/4 \\newline
\\end{bmatrix}, \\quad u, v \\in [-1, 1].
```
"""
mutable struct QuadrilateralLinearSerendipity{T<:AbstractFloat} <: QuadrilateralSerendipity
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Quadrilateral with node placements        \n
                                            4  7  3 \n
                                            8     6 \n
                                            1  5  2 \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
            [(1 - u)(v - 1)(u + v + 1)/4   \\newline
            (1 + u)(v - 1)(v - u + 1)/4 \\newline
            (1 + u)(v + 1)(u + v - 1)/4 \\newline
            (u - 1)(v + 1)(u - v + 1)/4 \\newline
            (1 - v)(1 - u^2)/2          \\newline
            (1 + u)(1 - v^2)/2          \\newline
            (1 + v)(1 - u^2)/2          \\newline
            (1 - u)(1 - v^2)/2
\\end{bmatrix}, \\quad u, v \\in [-1, 1].
```
"""
mutable struct QuadrilateralQuadratic{T<:AbstractFloat} <: QuadrilateralSerendipity
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Quadrilateral with node placements         \n
                                            3 -- 4   \n
                                            |    |   \n
                                            1 -- 2   \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
          (1 - u)(1 - v)/4 \\newline
          (1 + u)(1 - v)/4 \\newline
          (1 - u)(1 + v)/4 \\newline
          (1 + u)(1 + v)/4 \\newline
\\end{bmatrix}, \\quad u, v \\in [-1, 1].
```
"""
mutable struct QuadrilateralLinear4{T<:AbstractFloat} <: ContinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Quadrilateral with node placements         \n
                                            3  9  4  \n
                                            6  7  8  \n
                                            1  5  2  \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =
\\begin{bmatrix}
             u(1 - u)v(1 - v)/4          \\newline
            -u(1 + u)v(1 - v)/4        \\newline
            -u(1 - u)v(1 + v)/4        \\newline
            u(1 + u)v(1 + v)/4         \\newline
            -(1 + u)(1 - u)v(1 - v)/2  \\newline
            -u(1 - u)(1 + v)(1 - v)/2  \\newline
            (1 - u^2)(1 - v^2)         \\newline
            u(1 + u)(1 + v)(1 - v)/2   \\newline
            (1 + u)(1 - u)(1 + v)v/2
\\end{bmatrix}, \\quad u, v \\in [-1, 1]
```
"""
mutable struct QuadrilateralQuadraticLagrange{T<:AbstractFloat} <: QuadrilateralLagrange
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct QuadrilateralCubicLagrange{T<:AbstractFloat} <: QuadrilateralLagrange
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Quadrilateral with node placements         \n
                                              ---    \n
                                            |  1  |  \n
                                              ---    \n
and basis functions defined as

```math
\\mathbf{N}(u,v) =  [1], \\quad u,v \\in [-1,1 ]
```
"""
mutable struct DiscontinuousQuadrilateralConstant{T<:AbstractFloat} <: DiscontinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
"""
Defines a Quadrilateral with node placements         \n
                                              ---    \n
                                            | 3 4 |  \n
                                            | 1 2 |  \n
                                              ---    \n
and basis functions defined as

```math
\\mathbf{N}(u,v) = quadrilateralLinear4(u/(1-\\alpha),v/(1-\\alpha)), \\quad u,v \\in [-1,1 ]
```
"""
mutable struct DiscontinuousQuadrilateralLinear4{T<:AbstractFloat} <: DiscontinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
"""
Defines a Quadrilateral with node placements            \n
                                              -------   \n
                                            | 3  9  4 | \n
                                            | 6  7  8 | \n
                                            | 1  5  2 | \n
                                              -------   \n
and basis functions defined as

```math
\\mathbf{N}(u,v) = quadrilateralQuadraticLagrange(u/(1-\\alpha),v/(1-\\alpha)), \\quad u,v \\in [-1,1 ]
```
"""
mutable struct DiscontinuousQuadrilateralQuadraticLagrange{T<:AbstractFloat} <: DiscontinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
#==========================================================================================
                                    Legendre Elements
==========================================================================================#
mutable struct QuadrilateralLegendre{T<:AbstractFloat} <: ContinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    M::Int64
    N::Int64
end
#==========================================================================================
            Generalized Elements: https://academic.csuohio.edu/duffy_s/CVE_512_11.pdf
                This is just an idea. Requires that we find an higher-order mesher
==========================================================================================#
mutable struct QuadrilateralLagrangeM{T<:AbstractFloat} <: QuadrilateralLagrange
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    M::Int64
end
mutable struct QuadrilateralSerendipityM{T<:AbstractFloat} <: QuadrilateralSerendipity
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    M::Int64
end
#==========================================================================================
                                Show
==========================================================================================#
function Base.show(io::IO, ::MIME"text/plain", surface_function::SurfaceFunction)
    println(io, "SurfaceFunction Defined by:    \t $(typeof(surface_function))")
    println(io, "Number of gauss nodes:         \t $(length(surface_function.gauss_u))")
end
#==========================================================================================
                                    Utility functions
==========================================================================================#
number_of_shape_functions(::TriangularLinear)                            = 3
number_of_shape_functions(::TriangularQuadratic)                         = 6
number_of_shape_functions(::DiscontinuousTriangularConstant)             = 1
number_of_shape_functions(::DiscontinuousTriangularLinear)               = 3
number_of_shape_functions(::DiscontinuousTriangularQuadratic)            = 6
number_of_shape_functions(::QuadrilateralLinearSerendipity)              = 4
number_of_shape_functions(::QuadrilateralQuadratic)                      = 8
number_of_shape_functions(::QuadrilateralLinear4)                        = 4
number_of_shape_functions(::QuadrilateralQuadraticLagrange)              = 9
number_of_shape_functions(::DiscontinuousQuadrilateralConstant)          = 1
number_of_shape_functions(::DiscontinuousQuadrilateralLinear4)           = 4
number_of_shape_functions(::DiscontinuousQuadrilateralQuadraticLagrange) = 9
# Overloading base
Base.eltype(::Type{QuadrilateralQuadratic{T}})                      where {T} = T
Base.eltype(::Type{QuadrilateralLinearSerendipity{T}})              where {T} = T
Base.eltype(::Type{QuadrilateralQuadraticLagrange{T}})              where {T} = T
Base.eltype(::Type{QuadrilateralLinear4{T}})                        where {T} = T
Base.eltype(::Type{TriangularQuadratic{T}})                         where {T} = T
Base.eltype(::Type{TriangularLinear{T}})                            where {T} = T
Base.eltype(::Type{DiscontinuousTriangularConstant{T}})             where {T} = T
Base.eltype(::Type{DiscontinuousTriangularLinear{T}})               where {T} = T
Base.eltype(::Type{DiscontinuousTriangularQuadratic{T}})            where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralConstant{T}})          where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralLinear4{T}})           where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralQuadraticLagrange{T}}) where {T} = T
#==========================================================================================
                         Introduction of some niceness features
==========================================================================================#
function (surface_function::SurfaceFunction)(ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunction(surface_function,ξ,η)
end
function (surface_function::SurfaceFunction)(x)
    return basisFunction(surface_function,x...)
end
function (surface_function_derivative::ShapeFunctionDerivative{T})(ξ,η) where {T <: SurfaceFunction}
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunctionDerivative(surface_function_derivative.shape_function,ξ,η)
end
function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: SurfaceFunction}
    return coordinates * K.surface_function.derivatives_u, coordinates * K.surface_function.derivatives_v
end
#==========================================================================================
    Computing derivatives using ForwardDiff: Eases the implementation of new element
    Note: This is only used as a fallback if `basisFunctionDerivative` is not defined
==========================================================================================#
function basisFunctionDerivative(surface_function::SurfaceFunction,ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    if length(ξ) == 1
        return basisFunctionDerivative(surface_function,[[ξ η]])
    else
        return basisFunctionDerivative(surface_function,[[x y] for (x,y) in zip(ξ,η)])
    end
end
function basisFunctionDerivative(surface_function::SurfaceFunction,nodes)
    jacobians = hcat(ForwardDiff.jacobian.(Ref(surface_function),nodes)...)
    return jacobians[:,1:2:end], jacobians[:,2:2:end] # Split into dξ and dη
end
# Work-in-progress
function basisFunctionSecondOrderDerivative(surface_function::SurfaceFunction,ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    if length(ξ) == 1
        return basisFunctionSecondOrderDerivative(surface_function,[[ξ η]])
    else
        return basisFunctionSecondOrderDerivative(surface_function,[[x y] for (x,y) in zip(ξ,η)])
    end
end
function basisFunctionSecondOrderDerivative(surface_function::SurfaceFunction,nodes)
    n_shape = number_of_shape_functions(surface_function)
    d2dx    = zeros(n_shape,length(nodes))
    d2dy    = zeros(n_shape,length(nodes))
    d2dxdy  = zeros(n_shape,length(nodes))
    for (i,node) in enumerate(nodes)
        H = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(surface_function,x)',node)
        d2dx[:,i]   = H[1:2:end,1]
        d2dy[:,i]   = H[2:2:end,2]
        d2dxdy[:,i] = H[2:2:end,1]
    end
    return d2dx,d2dy,d2dxdy
end
#==========================================================================================
                        Setting interpolation nodes equal to the
==========================================================================================#
# Triangular Elements
get_nodal_nodes_u(::TriangularLinear)    = [0.0; 1.0; 0.0]
get_nodal_nodes_v(::TriangularLinear)    = [0.0; 0.0; 1.0]
get_nodal_nodes_u(::TriangularQuadratic) = [0.0; 1.0; 0.0; 0.5; 0.5; 0.0]
get_nodal_nodes_v(::TriangularQuadratic) = [0.0; 0.0; 1.0; 0.0; 0.5; 0.5]
get_nodal_nodes_u(::DiscontinuousTriangularConstant) = [1/3]
get_nodal_nodes_v(::DiscontinuousTriangularConstant) = [1/3]
get_nodal_nodes_u(sf::DiscontinuousTriangularLinear) = [sf.beta; 1-2*sf.beta; sf.beta]
get_nodal_nodes_v(sf::DiscontinuousTriangularLinear) = [sf.beta; sf.beta; 1-2*sf.beta]
function get_nodal_nodes_u(sf::DiscontinuousTriangularQuadratic)
    return [sf.beta; 1-2*sf.beta; sf.beta; (1-sf.beta)/2; (1-sf.beta)/2; sf.beta]
end
function get_nodal_nodes_v(sf::DiscontinuousTriangularQuadratic)
    return [sf.beta; sf.beta; 1-2*sf.beta; sf.beta; (1-sf.beta)/2; (1-sf.beta)/2]
end
# Quadrilateral Elements
get_nodal_nodes_u(::QuadrilateralLinearSerendipity) = [-1.0; 1.0; 1.0;-1.0]
get_nodal_nodes_v(::QuadrilateralLinearSerendipity) = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(::QuadrilateralLinear4)           = [-1.0; 1.0;-1.0; 1.0]
get_nodal_nodes_v(::QuadrilateralLinear4)           = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(::QuadrilateralQuadratic)  = [-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0;-1.0]
get_nodal_nodes_v(::QuadrilateralQuadratic)  = [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0]
get_nodal_nodes_u(::DiscontinuousQuadrilateralConstant)  = [0.0]
get_nodal_nodes_v(::DiscontinuousQuadrilateralConstant)  = [0.0]
function get_nodal_nodes_u(::QuadrilateralQuadraticLagrange)
    return [-1.0; 1.0;-1.0; 1.0; 0.0;-1.0; 0.0; 1.0; 0.0]
end
function get_nodal_nodes_v(::QuadrilateralQuadraticLagrange)
    return [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 0.0; 0.0; 1.0]
end
function get_nodal_nodes_u(sf::DiscontinuousQuadrilateralLinear4)
    return [sf.alpha-1; 1-sf.alpha; sf.alpha-1; 1-sf.alpha]
end
function get_nodal_nodes_v(sf::DiscontinuousQuadrilateralLinear4)
    return [sf.alpha-1; sf.alpha-1; 1-sf.alpha; 1-sf.alpha]
end
function get_nodal_nodes_u(sf::DiscontinuousQuadrilateralQuadraticLagrange)
    return [sf.alpha-1; 1-sf.alpha; sf.alpha-1; 1-sf.alpha; 0; sf.alpha-1; 0; 1-sf.alpha; 0]
end
function get_nodal_nodes_v(sf::DiscontinuousQuadrilateralQuadraticLagrange)
    return [sf.alpha-1; sf.alpha-1; 1-sf.alpha; 1-sf.alpha; sf.alpha-1; 0; 0; 0; 1-sf.alpha]
end
"""
    setInterpolationNodal!(surface_function::SurfaceFunction)
Sets the interpolation nodes `shape_function` to be the nodal positions.
"""
function set_nodal_interpolation!(surface_function::SurfaceFunction)
    nodes_u = get_nodal_nodes_u(surface_function)
    nodes_v = get_nodal_nodes_v(surface_function)
    set_interpolation_nodes!(surface_function,nodes_u,nodes_v)
end
function set_interpolation_nodes!(surface_function,nodes_u,nodes_v)
    surface_function.gauss_u = nodes_u
    surface_function.gauss_v = nodes_v
    surface_function.derivatives_u,surface_function.derivatives_v = basisFunctionDerivative(
                                                                        surface_function,
                                                                        nodes_u',
                                                                        nodes_v')
    surface_function.interpolation = basisFunction(surface_function,nodes_u',nodes_v')
end
function set_interpolation_nodes!(surface_function::SurfaceFunction,physics_function::SurfaceFunction)
    surface_function.weights = physics_function.weights
    nodes_u = get_nodal_nodes_u(physics_function)
    nodes_v = get_nodal_nodes_v(physics_function)
    set_interpolation_nodes!(surface_function,nodes_u,nodes_v)
end
function set_interpolation_nodes!(surface_function::SurfaceFunction,gauss_u,gauss_v,weights)
    surface_function.weights = weights
    set_interpolation_nodes!(surface_function,gauss_u,gauss_v)
end
"""
    copy_interpolation_nodes!(physics_function::SurfaceFunction,surfaceFunction::Triangular)

Sets interpolating nodes of `physics_function` to be the same as the `surfaceFunction`.
"""
function copy_interpolation_nodes!(physics_function::SurfaceFunction,surface_function::SurfaceFunction)
    physics_function.gauss_u = surface_function.gauss_u
    physics_function.gauss_v = surface_function.gauss_v
    physics_function.weights = surface_function.weights
    interpolate_on_nodes!(physics_function)
end
function interpolate_on_nodes!(physics_function::SurfaceFunction)
    physics_function.interpolation = physics_function(physics_function.gauss_u',
                                                      physics_function.gauss_v')
    physics_function.derivatives_u, physics_function.derivatives_v = physics_function'(
                                                                physics_function.gauss_u',
                                                                physics_function.gauss_v'
                                                                )
end
