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
mutable struct TriangularLinear{T<:AbstractFloat} <: ContinuousTriangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
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
mutable struct QuadrilateralLinearSerendipity{T<:AbstractFloat} <: QuadrilateralSerendipity
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct QuadrilateralQuadratic{T<:AbstractFloat} <: QuadrilateralSerendipity
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct QuadrilateralLinear4{T<:AbstractFloat} <: ContinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
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
mutable struct DiscontinuousQuadrilateralConstant{T<:AbstractFloat} <: DiscontinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct DiscontinuousQuadrilateralLinear4{T<:AbstractFloat} <: DiscontinuousQuadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
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
    println(io, "Number of basis functions:     \t $(length(surface_function))")
end
#==========================================================================================
                                    Utility functions
==========================================================================================#
get_nodes(surface_function::SurfaceFunction)        = surface_function.gauss_u, 
                                                      surface_function.gauss_v
get_derivatives(surface_function::SurfaceFunction)  = surface_function.derivatives_u, 
                                                      surface_function.derivatives_v
number_of_shape_functions(surface_function::TriangularLinear)                            = 3
number_of_shape_functions(surface_function::TriangularQuadratic)                         = 6
number_of_shape_functions(surface_function::DiscontinuousTriangularConstant)             = 1
number_of_shape_functions(surface_function::DiscontinuousTriangularLinear)               = 3
number_of_shape_functions(surface_function::DiscontinuousTriangularQuadratic)            = 6
number_of_shape_functions(surface_function::QuadrilateralLinearSerendipity)              = 4
number_of_shape_functions(surface_function::QuadrilateralQuadratic)                      = 8
number_of_shape_functions(surface_function::QuadrilateralLinear4)                        = 4
number_of_shape_functions(surface_function::QuadrilateralQuadraticLagrange)              = 9
number_of_shape_functions(surface_function::DiscontinuousQuadrilateralConstant)          = 1
number_of_shape_functions(surface_function::DiscontinuousQuadrilateralLinear4)           = 4
number_of_shape_functions(surface_function::DiscontinuousQuadrilateralQuadraticLagrange) = 9

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
get_nodal_nodes_u(sf::TriangularLinear)     = [0.0; 1.0; 0.0]
get_nodal_nodes_v(sf::TriangularLinear)     = [0.0; 0.0; 1.0]
get_nodal_nodes_u(sf::TriangularQuadratic)  = [0.0; 1.0; 0.0; 0.5; 0.5; 0.0]
get_nodal_nodes_v(sf::TriangularQuadratic)  = [0.0; 0.0; 1.0; 0.0; 0.5; 0.5]
get_nodal_nodes_u(sf::DiscontinuousTriangularConstant)  = [1.0/3.0]
get_nodal_nodes_v(sf::DiscontinuousTriangularConstant)  = [1.0/3.0]
get_nodal_nodes_u(sf::DiscontinuousTriangularLinear)    = [sf.beta; 1.0-2.0*sf.beta; sf.beta]
get_nodal_nodes_v(sf::DiscontinuousTriangularLinear)    = [sf.beta; sf.beta; 1.0-2.0*sf.beta]
get_nodal_nodes_u(sf::DiscontinuousTriangularQuadratic) = [sf.beta; 1.0-2.0*sf.beta; sf.beta; (1.0-sf.beta)/2.0; (1.0-sf.beta)/2.0; sf.beta]
get_nodal_nodes_v(sf::DiscontinuousTriangularQuadratic) = [sf.beta; sf.beta; 1.0-2.0*sf.beta; sf.beta; (1.0-sf.beta)/2.0; (1.0-sf.beta)/2.0]
# Quadrilateral Elements
get_nodal_nodes_u(sf::QuadrilateralLinearSerendipity)       = [-1.0; 1.0; 1.0;-1.0]
get_nodal_nodes_v(sf::QuadrilateralLinearSerendipity)       = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(sf::QuadrilateralLinear4)                 = [-1.0; 1.0;-1.0; 1.0]
get_nodal_nodes_v(sf::QuadrilateralLinear4)                 = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(sf::QuadrilateralQuadratic)  = [-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0;-1.0]
get_nodal_nodes_v(sf::QuadrilateralQuadratic)  = [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0]
get_nodal_nodes_u(sf::QuadrilateralQuadraticLagrange) = [-1.0; 1.0;-1.0; 1.0; 0.0;-1.0; 0.0; 1.0; 0.0]
get_nodal_nodes_v(sf::QuadrilateralQuadraticLagrange) = [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 0.0; 0.0; 1.0]
get_nodal_nodes_u(sf::DiscontinuousQuadrilateralConstant)   = [0.0]
get_nodal_nodes_v(sf::DiscontinuousQuadrilateralConstant)   = [0.0]
get_nodal_nodes_u(sf::DiscontinuousQuadrilateralLinear4)    = [sf.alpha-1.0; 1.0-sf.alpha; sf.alpha-1.0; 1.0-sf.alpha]
get_nodal_nodes_v(sf::DiscontinuousQuadrilateralLinear4)    = [sf.alpha-1.0; sf.alpha-1.0; 1.0-sf.alpha; 1.0-sf.alpha]
get_nodal_nodes_u(sf::DiscontinuousQuadrilateralQuadraticLagrange) = [sf.alpha-1.0; 1.0-sf.alpha; sf.alpha-1.0; 1.0-sf.alpha;          0.0; sf.alpha-1.0; 0.0; 1.0-sf.alpha; 0.0]
get_nodal_nodes_v(sf::DiscontinuousQuadrilateralQuadraticLagrange) = [sf.alpha-1.0; sf.alpha-1.0; 1.0-sf.alpha; 1.0-sf.alpha; sf.alpha-1.0;          0.0; 0.0;          0.0; 1.0-sf.alpha]
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
    physics_function.interpolation = physics_function(physics_function.gauss_u',physics_function.gauss_v')
    physics_function.derivatives_u, physics_function.derivatives_v = physics_function'(physics_function.gauss_u',physics_function.gauss_v')
end
