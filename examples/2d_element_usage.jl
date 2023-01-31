# # Surface Functions (2D)
# For 2D boundary integral equations the elements of interest are defined as `CurveFunction`s. These are defined as a vector of basis functions defined by a single local coordinate ($u$)
# ```math
#   \mathbf{N}(u)
# ```
# The following is s brief explanation of how the `CurveFunction`s are implemented in `BoundaryIntegralEquations.jl`. In short the `SurfaceFunction` is contain the following information of the surface element in local coordinates
# * `gauss`: Vector of $u$-values of the Gaussian points.
# * `weights`: Vector of Gaussian weights.
# * `interpolation`: Matrix with columns equal to the basis functions evauluated at the Gaussian nodes, i.e. $\mathbf{N}(u_i)$.
# * `derivatives`: Matrix with columns equal to the derivative of the basis functions with respect to $u$ evaluated at the Gaussian nodes, i.e. $\mathbf{N}_u(u_i)$.


# First we import relevant packages
using Plots
using LinearAlgebra
using BoundaryIntegralEquations
using FastGaussQuadrature
# # Linear Curve
curve_linear = ContinuousCurveLinear(10)
m1 = [1;0.0]
m2 = [1.5;0.5]
coordinates = [m1 m2]
interpolation = coordinates*curve_linear.interpolation
scatter(coordinates[1,:],coordinates[2,:], label="Interpolation Nodes", aspect_ratio=1)
scatter!(interpolation[1,:],interpolation[2,:], label = "Gaussian Nodes")
plot!(interpolation[1,:],interpolation[2,:], label="Element")


A_exact = norm(m1 - m2)
dU = coordinates * curve_linear.derivatives
A_approx = dot(sqrt.(sum(dU.^2,dims=1))', curve_linear.weights)

A_exact - A_approx

# # Quadratic Curves
m3 = [1;1.0]
coordinates = [m1 m2 m3]
curve_quadratic = ContinuousCurveQuadratic(10)
interpolation = coordinates*curve_quadratic.interpolation
scatter(coordinates[1,:],coordinates[2,:], label="Interpolation Nodes", aspect_ratio=1)
scatter!(interpolation[1,:],interpolation[2,:], label = "Gaussian Nodes")
plot!(interpolation[1,:],interpolation[2,:], label="Element")

# # Discontinuous Curves
