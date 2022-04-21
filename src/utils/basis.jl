abstract type ShapeFunction end
struct ShapeFunctionDerivative{T <: ShapeFunction} <: ShapeFunction
    basis::T
    function ShapeFunctionDerivative(basis::T) where {T<:ShapeFunction}
		new{T}(basis)
	end
end
abstract type CurveFunction              <: ShapeFunction   end
abstract type SurfaceFunction            <: ShapeFunction   end
abstract type Triangular                 <: SurfaceFunction end
abstract type Quadrilateral              <: SurfaceFunction end
abstract type DiscontinuousTriangular    <: Triangular      end
abstract type DiscontinuousQuadrilateral <: Quadrilateral   end
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
mutable struct QuadrilateralQuadratic{T<:AbstractFloat} <: Quadrilateral
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
mutable struct QuadrilateralQuadratic9{T<:AbstractFloat} <: Quadrilateral
    weights::AbstractArray{T,1}
    gauss_u::AbstractArray{T,1}
    gauss_v::AbstractArray{T,1}
    derivatives_u::AbstractArray{T,2}
    derivatives_v::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
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
                                    Utility functions
==========================================================================================#
get_nodes(basis::CurveFunction)          = basis.gauss
get_nodes(basis::SurfaceFunction)        = basis.gauss_u, basis.gauss_v
get_derivatives(basis::CurveFunction)    = basis.derivatives
get_derivatives(basis::SurfaceFunction)  = basis.derivatives_u, basis.derivatives_v
get_weights(basis::ShapeFunction)::Array{Float64,1}        = basis.weights
get_interpolation(basis::ShapeFunction)::Array{Float64,2}  = basis.interpolation
number_of_shape_functions(basisFunction::TriangularLinear)                 = 3
number_of_shape_functions(basisFunction::TriangularQuadratic)              = 6
number_of_shape_functions(basisFunction::DiscontinuousTriangularConstant)  = 1
number_of_shape_functions(basisFunction::DiscontinuousTriangularLinear)    = 3
number_of_shape_functions(basisFunction::DiscontinuousTriangularQuadratic) = 6
number_of_shape_functions(basisFunction::QuadrilateralLinear)              = 4
number_of_shape_functions(basisFunction::QuadrilateralQuadratic)           = 8
number_of_shape_functions(basisFunction::QuadrilateralLinear4)             = 4
number_of_shape_functions(basisFunction::QuadrilateralQuadratic9)          = 9

Base.eltype(::Type{QuadrilateralQuadratic{T}})              where {T} = T
Base.eltype(::Type{QuadrilateralLinear{T}})                 where {T} = T
Base.eltype(::Type{QuadrilateralQuadratic9{T}})             where {T} = T
Base.eltype(::Type{QuadrilateralLinear4{T}})                where {T} = T
Base.eltype(::Type{TriangularQuadratic{T}})                 where {T} = T
Base.eltype(::Type{TriangularLinear{T}})                    where {T} = T
Base.eltype(::Type{DiscontinuousTriangularConstant{T}})     where {T} = T
Base.eltype(::Type{DiscontinuousTriangularLinear{T}})       where {T} = T
Base.eltype(::Type{DiscontinuousTriangularQuadratic{T}})    where {T} = T

Base.length(basisFunction::SurfaceFunction) = number_of_shape_functions(basisFunction)
#==========================================================================================
                         Introduction of some niceness features 
==========================================================================================#
# Making the structs "callable"
function (basis::CurveFunction)(ξ)
    return basisFunction(basis,ξ)
end
function (basis::SurfaceFunction)(ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunction(basis,ξ,η)
end
function (basis::SurfaceFunction)(x)
    return basisFunction(basis,x...)
end
function (derivative::ShapeFunctionDerivative{T})(ξ)   where {T <: ShapeFunction}
    return basisFunctionDerivative(derivative.basis,ξ)
end
function (derivative::ShapeFunctionDerivative{T})(ξ,η) where {T <: SurfaceFunction}
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    return basisFunctionDerivative(derivative.basis,ξ,η)
end
# Create a jacobian function from a ShapeFunction
jacobian(basis::ShapeFunction) = coordinates -> jacobian(basis,coordinates)::Array{Float64,2}
# Defining multiplication for certain types
(*)(coordinates::AbstractArray,K::ShapeFunction)::Array{Float64,2} = coordinates * K.interpolation
function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: CurveFunction}
    return coordinates * K.basis.derivatives
end
function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: SurfaceFunction}
    return coordinates * K.basis.derivatives_u, coordinates * K.basis.derivatives_v
end
# This is slight "Abuse-of-notation" Using the transpose operator for derivatives
adjoint(basis::T) where {T<:ShapeFunction} = ShapeFunctionDerivative(basis)
#==========================================================================================
    Computing derivatives using ForwardDiff: Eases the implementation of new element
    Note: This is only used as a fallback if `basisFunctionDerivative` is not defined
==========================================================================================#
function basisFunctionDerivative(surfaceFunction::SurfaceFunction,ξ,η)
    if length(ξ) != length(η)
        throw(ArgumentError("Input does not have equal length."))
    end
    if length(ξ) == 1
        return basisFunctionDerivative(surfaceFunction,[[ξ η]])
    else
        return basisFunctionDerivative(surfaceFunction,[[x y] for (x,y) in zip(ξ,η)])
    end
end
function basisFunctionDerivative(surfaceFunction::SurfaceFunction,nodes)
    jacobians = hcat(ForwardDiff.jacobian.(Ref(surfaceFunction),nodes)...)
    return jacobians[:,1:2:end], jacobians[:,2:2:end] # Split into dξ and dη
end

#==========================================================================================
                        Setting interpolation nodes equal to the  
==========================================================================================#
"""
    setInterpolationNodal!(surfaceFunction::SurfaceFunction)

Sets the interpolation nodes `surfaceFunction` to be the nodal positions.
"""
function set_nodal_interpolation!(surfaceFunction::SurfaceFunction)
    nodes_u = get_nodal_nodes_u(surfaceFunction)
    nodes_v = get_nodal_nodes_v(surfaceFunction)
    set_interpolation_nodes!(surfaceFunction,nodes_u,nodes_v)
end
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

get_nodal_nodes_u(sf::QuadrilateralLinear)     = [-1.0; 1.0; 1.0;-1.0]
get_nodal_nodes_v(sf::QuadrilateralLinear)     = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(sf::QuadrilateralLinear4)    = [-1.0; 1.0;-1.0; 1.0]
get_nodal_nodes_v(sf::QuadrilateralLinear4)    = [-1.0;-1.0; 1.0; 1.0]
get_nodal_nodes_u(sf::QuadrilateralQuadratic)  = [-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0;-1.0]
get_nodal_nodes_v(sf::QuadrilateralQuadratic)  = [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 1.0; 0.0]
get_nodal_nodes_u(sf::QuadrilateralQuadratic9) = [-1.0; 1.0;-1.0; 1.0; 0.0;-1.0; 0.0; 1.0; 0.0]
get_nodal_nodes_v(sf::QuadrilateralQuadratic9) = [-1.0;-1.0; 1.0; 1.0;-1.0; 0.0; 0.0; 0.0; 1.0]


function set_interpolation_nodes!(surfaceFunction,nodes_u,nodes_v)
    surfaceFunction.gauss_u = nodes_u
    surfaceFunction.gauss_v = nodes_v
    surfaceFunction.derivatives_u,surfaceFunction.derivatives_v = basisFunctionDerivative(
                                                                        surfaceFunction,
                                                                        nodes_u',
                                                                        nodes_v')
    surfaceFunction.interpolation = basisFunction(surfaceFunction,nodes_u',nodes_v')
end
