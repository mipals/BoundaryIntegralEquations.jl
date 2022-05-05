#==========================================================================================
                            Defining Various Element (Sub)Types
==========================================================================================#
abstract type ShapeFunction                                 end
abstract type CurveFunction              <: ShapeFunction   end
abstract type SurfaceFunction            <: ShapeFunction   end
abstract type Triangular                 <: SurfaceFunction end
abstract type Quadrilateral              <: SurfaceFunction end
abstract type DiscontinuousTriangular    <: Triangular      end
abstract type DiscontinuousQuadrilateral <: Quadrilateral   end
abstract type QuadrilateralLagrange      <: Quadrilateral   end
abstract type QuadrilateralSerendipity   <: Quadrilateral   end
#==========================================================================================
                            A general definition of a derivative
==========================================================================================#
struct ShapeFunctionDerivative{T <: ShapeFunction} <: ShapeFunction
    shape_function::T
    function ShapeFunctionDerivative(shape_function::T) where {T<:ShapeFunction}
		new{T}(shape_function)
	end
end
#==========================================================================================
                                Curve elements for 2D
==========================================================================================#
mutable struct CurveLinear{T<:AbstractFloat} <: CurveFunction 
    weights::AbstractArray{T,1}
    gauss::AbstractArray{T,1}
    derivatives::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct CurveQuadratic{T<:AbstractFloat} <: CurveFunction
    weights::AbstractArray{T,1}
    gauss::AbstractArray{T,1}
    derivatives::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
#==========================================================================================
                                Triangular Elements
==========================================================================================#
mutable struct TriangularLinear{T<:AbstractFloat} <: Triangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct TriangularQuadratic{T<:AbstractFloat} <: Triangular
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct TriangularCubic{T<:AbstractFloat} <: Triangular
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
mutable struct QuadrilateralLinear{T<:AbstractFloat} <: Quadrilateral
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
mutable struct QuadrilateralLinear4{T<:AbstractFloat} <: Quadrilateral
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
mutable struct DiscontinuousQuadrilateralConstant{T<:AbstractFloat} <: Quadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct DiscontinuousQuadrilateralLinear4{T<:AbstractFloat} <: Quadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
mutable struct DiscontinuousQuadrilateralQuadraticLagrange{T<:AbstractFloat} <: Quadrilateral
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
mutable struct QuadrilateralLegendre{T<:AbstractFloat} <: Quadrilateral
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
get_nodes(shape_function::CurveFunction)          = shape_function.gauss
get_nodes(shape_function::SurfaceFunction)        = shape_function.gauss_u, 
                                                    shape_function.gauss_v
get_derivatives(shape_function::CurveFunction)    = shape_function.derivatives
get_derivatives(shape_function::SurfaceFunction)  = shape_function.derivatives_u, 
                                                    shape_function.derivatives_v
get_weights(shape_function::ShapeFunction)::Array{Float64,1}        = shape_function.weights
get_interpolation(shape_function::ShapeFunction)::Array{Float64,2}  = shape_function.interpolation
number_of_shape_functions(shape_function::TriangularLinear)                     = 3
number_of_shape_functions(shape_function::TriangularQuadratic)                  = 6
number_of_shape_functions(shape_function::DiscontinuousTriangularConstant)      = 1
number_of_shape_functions(shape_function::DiscontinuousTriangularLinear)        = 3
number_of_shape_functions(shape_function::DiscontinuousTriangularQuadratic)     = 6
number_of_shape_functions(shape_function::QuadrilateralLinear)                  = 4
number_of_shape_functions(shape_function::QuadrilateralQuadratic)               = 8
number_of_shape_functions(shape_function::QuadrilateralLinear4)                 = 4
number_of_shape_functions(shape_function::QuadrilateralQuadraticLagrange)              = 9
number_of_shape_functions(shape_function::DiscontinuousQuadrilateralConstant)   = 1
number_of_shape_functions(shape_function::DiscontinuousQuadrilateralLinear4)    = 4
number_of_shape_functions(shape_function::DiscontinuousQuadrilateralQuadraticLagrange) = 9

Base.eltype(::Type{QuadrilateralQuadratic{T}})                  where {T} = T
Base.eltype(::Type{QuadrilateralLinear{T}})                     where {T} = T
Base.eltype(::Type{QuadrilateralQuadraticLagrange{T}})                 where {T} = T
Base.eltype(::Type{QuadrilateralLinear4{T}})                    where {T} = T
Base.eltype(::Type{TriangularQuadratic{T}})                     where {T} = T
Base.eltype(::Type{TriangularLinear{T}})                        where {T} = T
Base.eltype(::Type{DiscontinuousTriangularConstant{T}})         where {T} = T
Base.eltype(::Type{DiscontinuousTriangularLinear{T}})           where {T} = T
Base.eltype(::Type{DiscontinuousTriangularQuadratic{T}})        where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralConstant{T}})      where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralLinear4{T}})       where {T} = T
Base.eltype(::Type{DiscontinuousQuadrilateralQuadraticLagrange{T}})    where {T} = T

Base.length(basisFunction::SurfaceFunction) = number_of_shape_functions(basisFunction)
#==========================================================================================
                         Introduction of some niceness features 
==========================================================================================#
# Making the structs "callable"
function (shape_function::CurveFunction)(ξ)
    return basisFunction(shape_function,ξ)
end
function (shape_function::SurfaceFunction)(ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunction(shape_function,ξ,η)
end
function (shape_function::SurfaceFunction)(x)
    return basisFunction(shape_function,x...)
end
function (shape_function_derivative::ShapeFunctionDerivative{T})(ξ)   where {T <: ShapeFunction}
    return basisFunctionDerivative(shape_function_derivative.shape_function,ξ)
end
function (shape_function_derivative::ShapeFunctionDerivative{T})(ξ,η) where {T <: SurfaceFunction}
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunctionDerivative(shape_function_derivative.shape_function,ξ,η)
end
# Create a jacobian function from a ShapeFunction
jacobian(shape_function::ShapeFunction) = coordinates -> jacobian(shape_function,coordinates)::Array{Float64,2}
# Defining multiplication for certain types
(*)(coordinates::AbstractArray,K::ShapeFunction)::Array{Float64,2} = coordinates * K.interpolation
function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: CurveFunction}
    return coordinates * K.shape_function.derivatives
end
function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: SurfaceFunction}
    return coordinates * K.shape_function.derivatives_u, coordinates * K.shape_function.derivatives_v
end
# This is slight "Abuse-of-notation" Using the transpose operator for derivatives
adjoint(shape_function::T) where {T<:ShapeFunction} = ShapeFunctionDerivative(shape_function)
#==========================================================================================
    Computing derivatives using ForwardDiff: Eases the implementation of new element
    Note: This is only used as a fallback if `basisFunctionDerivative` is not defined
==========================================================================================#
function basisFunctionDerivative(shape_function::SurfaceFunction,ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    if length(ξ) == 1
        return basisFunctionDerivative(shape_function,[[ξ η]])
    else
        return basisFunctionDerivative(shape_function,[[x y] for (x,y) in zip(ξ,η)])
    end
end
function basisFunctionDerivative(shape_function::SurfaceFunction,nodes)
    jacobians = hcat(ForwardDiff.jacobian.(Ref(shape_function),nodes)...)
    return jacobians[:,1:2:end], jacobians[:,2:2:end] # Split into dξ and dη
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
get_nodal_nodes_u(sf::QuadrilateralLinear)     = [-1.0; 1.0; 1.0;-1.0]
get_nodal_nodes_v(sf::QuadrilateralLinear)     = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(sf::QuadrilateralLinear4)    = [-1.0; 1.0;-1.0; 1.0]
get_nodal_nodes_v(sf::QuadrilateralLinear4)    = [-1.0;-1.0; 1.0; 1.0]
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
    setInterpolationNodal!(shape_function::SurfaceFunction)
Sets the interpolation nodes `shape_function` to be the nodal positions.
"""
function set_nodal_interpolation!(shape_function::SurfaceFunction)
    nodes_u = get_nodal_nodes_u(shape_function)
    nodes_v = get_nodal_nodes_v(shape_function)
    set_interpolation_nodes!(shape_function,nodes_u,nodes_v)
end
function set_interpolation_nodes!(shape_function,nodes_u,nodes_v)
    shape_function.gauss_u = nodes_u
    shape_function.gauss_v = nodes_v
    shape_function.derivatives_u,shape_function.derivatives_v = basisFunctionDerivative(
                                                                        shape_function,
                                                                        nodes_u',
                                                                        nodes_v')
    shape_function.interpolation = basisFunction(shape_function,nodes_u',nodes_v')
end
function set_interpolation_nodes!(shape_function::SurfaceFunction,physics_function::SurfaceFunction)
    nodes_u = get_nodal_nodes_u(physics_function)
    nodes_v = get_nodal_nodes_v(physics_function)
    set_interpolation_nodes!(shape_function,nodes_u,nodes_v)
end
"""
    copy_interpolation_nodes!(physics_function::SurfaceFunction,surfaceFunction::Triangular)

Sets interpolating nodes of `physics_function` to be the same as the `surfaceFunction`.
"""
function copy_interpolation_nodes!(physics_function::SurfaceFunction,shape_function::SurfaceFunction)
    physics_function.gauss_u = shape_function.gauss_u
    physics_function.gauss_v = shape_function.gauss_v
    physics_function.weights = shape_function.weights
    interpolate_on_nodes!(physics_function)
end
function interpolate_on_nodes!(physics_function::SurfaceFunction)
    physics_function.interpolation = physics_function(physics_function.gauss_u',physics_function.gauss_v')
    physics_function.derivatives_u, physics_function.derivatives_v = physics_function'(physics_function.gauss_u',physics_function.gauss_v')
end
