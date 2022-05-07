#==========================================================================================
                Abstract Type for a shape function in 2d and 3d
==========================================================================================#
abstract type ShapeFunction end
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
                Properties for shape functions in 2d and 3d
==========================================================================================#
get_weights(shape_function::ShapeFunction)::Array{Float64,1}        = shape_function.weights
get_interpolation(shape_function::ShapeFunction)::Array{Float64,2}  = shape_function.interpolation
Base.length(basisFunction::ShapeFunction) = number_of_shape_functions(basisFunction)

# Create a jacobian function from a ShapeFunction
jacobian(shape_function::ShapeFunction) = coordinates -> jacobian(shape_function,coordinates)::Array{Float64,2}
# Defining multiplication for certain types
(*)(coordinates::AbstractArray,K::ShapeFunction)::Array{Float64,2} = coordinates * K.interpolation
# This is slight "Abuse-of-notation" Using the transpose operator for derivatives
adjoint(shape_function::T) where {T<:ShapeFunction} = ShapeFunctionDerivative(shape_function)
