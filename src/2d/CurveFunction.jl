abstract type CurveFunction              <: ShapeFunction   end
abstract type DiscontinuousCurveFunction <: CurveFunction   end

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
mutable struct DiscontinuousCurveConstant{T<:AbstractFloat} <: DiscontinuousCurveFunction 
    weights::AbstractArray{T,1}
    gauss::AbstractArray{T,1}
    derivatives::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
end
mutable struct DiscontinuousCurveLinear{T<:AbstractFloat} <: DiscontinuousCurveFunction 
    weights::AbstractArray{T,1}
    gauss::AbstractArray{T,1}
    derivatives::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
mutable struct DiscontinuousCurveQuadratic{T<:AbstractFloat} <: DiscontinuousCurveFunction
    weights::AbstractArray{T,1}
    gauss::AbstractArray{T,1}
    derivatives::AbstractArray{T,2}
    interpolation::AbstractArray{T,2}
    alpha::T
end
#==========================================================================================
                                Show
==========================================================================================#
function Base.show(io::IO, ::MIME"text/plain", curve_function::CurveFunction)
    println(io, "SurfaceFunction Defined by:    \t $(typeof(curve_function))")
    println(io, "Number of gauss nodes:         \t $(length(curve_function.gauss))")
    println(io, "Number of basis functions:     \t $(length(curve_function))")
end
#==========================================================================================
                         Introduction of some niceness features 
==========================================================================================#
# Making the structs "callable"
function (curve_function::CurveFunction)(ξ)
    return basisFunction(shape_function,ξ)
end
function (shape_function_derivative::ShapeFunctionDerivative{T})(ξ)   where {T <: CurveFunction}
    return basisFunctionDerivative(shape_function_derivative.shape_function,ξ)
end

function (*)(coordinates::AbstractArray,K::ShapeFunctionDerivative{T}) where {T <: CurveFunction}
    return coordinates * K.curve_function.derivatives
end

get_derivatives(curve_function::CurveFunction)    = curve_function.derivatives

number_of_shape_functions(curve_function::CurveLinear)                  = 2
number_of_shape_functions(curve_function::CurveQuadratic)               = 3
number_of_shape_functions(curve_function::DiscontinuousCurveConstant)   = 1
number_of_shape_functions(curve_function::DiscontinuousCurveLinear)     = 2
number_of_shape_functions(curve_function::DiscontinuousCurveQuadratic)  = 3

#==========================================================================================
                        Overloading Base
==========================================================================================#
Base.eltype(::Type{CurveLinear{T}})                 where {T} = T
Base.eltype(::Type{CurveQuadratic{T}})              where {T} = T
Base.eltype(::Type{DiscontinuousCurveConstant{T}})  where {T} = T
Base.eltype(::Type{DiscontinuousCurveLinear{T}})    where {T} = T
Base.eltype(::Type{DiscontinuousCurveQuadratic{T}}) where {T} = T

#==========================================================================================
                        Modifications of CurveFunction structs 
==========================================================================================#
function set_interpolation_nodes!(curve_function::CurveFunction,nodes)
    curve_function.nodes         = nodes
    curve_function.derivatives   = basisFunctionDerivative(shape_function,nodes')
    curve_function.interpolation = basisFunction(shape_function,nodes')
end
function set_interpolation_nodes!(curve_function::CurveFunction,physics_function::CurveFunction)
    curve_function.weights = physics_function.weights
    nodes = get_nodal_nodes(physics_function)
    set_interpolation_nodes!(shape_function,nodes)
end
function interpolate_on_nodes!(curve_function::CurveFunction)
    curve_function.interpolation = shape_function(curve_function.gauss)
    curve_function.derivatives   = shape_function'(curve_function.gauss')
end

get_nodes(curve_function::CurveFunction)                     = curve_function.gauss
get_nodal_nodes(curve_function::CurveLinear)                 = [-1.0; 1.0]
get_nodal_nodes(curve_function::CurveQuadratic)              = [-1.0; 0.0; 1.0]
get_nodal_nodes(curve_function::DiscontinuousCurveConstant)  = [0.0]
get_nodal_nodes(curve_function::DiscontinuousCurveLinear)    = [-curve_function.alpha; 
                                                                 curve_function.alpha]
get_nodal_nodes(curve_function::DiscontinuousCurveQuadratic) = [-curve_function.alpha; 0.0; 
                                                                 curve_function.alpha]
