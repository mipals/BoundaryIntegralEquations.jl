# # Surface Functions
# For 3D boundary integral equations the elements of interest are defined as `SurfaceFunction`s. These are defined as a vector of basis functions defined by two local coordinates, $u$ and $v$
# ```math
#   \mathbf{N}(u,v)
# ```
# The `SurfaceFunction` contain the following information of the surface element in local coordinates
# * `gauss_u`: Vector of $u$-values of the Gaussian points
# * `gauss_v`: Vector $v$-values of the Gaussian points
# * `weights`: Vector of Gaussian weights
# * `interpolation`: Matrix with columns equal to the basis functions evauluated at the Gaussian nodes, i.e. $\mathbf{N}(u_i,v_i)$
# * `derivatives_u`: Matrix with columns equal to the derivative of the basis functions with respect to $u$ evaluated at the Gaussian nodes, i.e. $\mathbf{N}_u(u_i,v_i)$
# * `derivatives_v`: Matrix with columns equal to the derivative of the basis functions with respect to $v$ evaluated at the Gaussian nodes, i.e. $\mathbf{N}_v(u_i,v_i)$
# All this information is required when computing the underyling boundary integrals.

# # Importing relevant packages
using Plots
gr()
using LinearAlgebra
using BoundaryIntegralEquations
#hide using DisplayAs
# # Triangles
# The simplest `SurfaceFunction` is `TriangularLinear` which defines a linear interpolation of a triangle. A `TriangularLinear` with 3 gaussian points is created as follows
linear_triangular = TriangularLinear(3)
# The Gaussian points is defined in the local coordinates as
scatter(linear_triangular.gauss_u, linear_triangular.gauss_v,
        label = "Gauss Points", aspect_ratio = :equal, xlabel="u",ylabel="v")
# While the Gaussian points are important for the integral of the local surfaces the nodal points of the reference triangle are important when computing normals and derivatives at the nodal positions, i.e. the positions where only one basis function is equal 1 and the rest are 0. This can be done by using the function `set_nodal_interpolation!` that is defined on all `SurfaceFunction`
nodal_triangle = TriangularLinear(3)
BoundaryIntegralEquations.set_nodal_interpolation!(nodal_triangle) # Should be the identity matrix
# Using this we can plot the full reference triangle
scatter!(nodal_triangle.gauss_u, nodal_triangle.gauss_v,
        label="Nodal Points")
plot!(Shape(nodal_triangle.gauss_u, nodal_triangle.gauss_v),
        fillalpha=0.2, label="Reference Triangle")
# When discritizing boundary integrals over surfaces in 3d an important aspect is the mapping from the local $(u,v)$-coordinates to global $(x,y)$ coordinates. This transformation is linear with respect to the basis functions, but not nessecrarily the local coordinates, and is given by
# ```math
#   \mathbf{x}^e = \mathbf{X}^e\mathbf{N}(u,v).
# ```
# As an example we take the triangle with corners equal to the standard basis vectors
x1 = [1.0;0.0;0.0]
x2 = [0.0;1.0;0.0]
x3 = [0.0;0.0;1.0]
X = [x1 x2 x3]
# Using the mapping from local to global coordinates the Gaussian points can be mapped onto the $(x,y)$-space as
interp = X*linear_triangular.interpolation
scatter3d(interp[1,:],interp[2,:],interp[3,:],label="Gauss Points")
scatter3d!(X[1,:],X[2,:],X[3,:], label="Nodal points")
plot3d!([X[1,:];x1[1]],[X[2,:];x1[2]],[X[3,:];x1[3]],label="Global Triangle")

# The main purpose of this package is the discretization of surface integrals using elements. The simplest of such integration is where the integrand is the constant 1, as it equals the area of the surface. For a single element this can be seen as
# ```math
#   \int_{\Gamma_e} 1 \ \mathrm{d}\Gamma_\mathbf{y} = \int_0^{1-u}\int_0^1 \text{jacobian}(u,v)\ \mathrm{d}u\mathrm{d}v \approx \sum_{i=1}^{Q} \text{jacobian}(u_i,v_i) w_i,
# ```
# where $n_g$ is the number of Gaussian points $(u_i, v_i)$, $w_i$ is the corresponding weight and the jacobian is defined as
# ```math
# \text{jacobian} = \left\|\left(\frac{\mathrm{d}\mathbf{x}^e}{\mathrm{d}u}\right) \times \left(\frac{\mathrm{d}\mathbf{x}^e}{\mathrm{d}v}\right)\right\|_2
# ```
# which can be thought of as the area deformation stemming from the mapping from local to global coordinates.
# The exact area of a flat triangle is equal half of the area of the paralllelogram spanned by two sides of the triangle. This is computed similarly to the jacobian as
# ```math
#   \text{Area}(\Delta) = \frac{1}{2}\|(\mathbf{x}^e_3 - \mathbf{x}^e_1)\times(\mathbf{x}^e_2-\mathbf{x}^e_1)\|_2.
# ```
#
# We now start by computing the exact area of the flat triangle as
A_exact = norm(cross(x3-x1,x2-x1))/2
# To compute the are using the integral we first need to compute the gradient
dU = X*linear_triangular.derivatives_u
dV = X*linear_triangular.derivatives_v
jacobian = similar(linear_triangular.weights)
for i = eachindex(jacobian)
    jacobian[i] = norm(cross(dU[:,i],dV[:,i]))
end
# As a final step the sum over the jacobian multiplied with the weights is
A_approx = dot(jacobian,linear_triangular.weights)
# In this case, since the jacobian is constant over the element, the exact area of the triangle is equal to the approximated integral
A_exact - A_approx
# # Quadrilaterals
# We can perform the same computations for a quadrilateral element
linear_quadrilateral = QuadrilateralLinear4(3,3)
# The Gaussian points is defined in the local coordinates as
scatter(linear_quadrilateral.gauss_u, linear_quadrilateral.gauss_v,
    legend=:outertop,label = "Gauss Points", aspect_ratio = :equal, xlabel="u",ylabel="v")
# While the Gaussian points are important for the integral of the local surfaces the nodes of the reference triangle are important when computing normals and derivatives at the nodal positions. As such the function `set_nodal_interpolation!` is defined on all `SurfaceFunction`s
nodal_quadrilateral = QuadrilateralLinear4(3,3)
BoundaryIntegralEquations.set_nodal_interpolation!(nodal_quadrilateral);
# Using this we can plot the full reference quadrilateral
scatter!(nodal_quadrilateral.gauss_u, nodal_quadrilateral.gauss_v,
        label="Nodes")
idx = [1;2;4;3;1] # Index set
plot!(Shape(nodal_quadrilateral.gauss_u[idx], nodal_quadrilateral.gauss_v[idx]),
        fillalpha=0.2, label="Reference Quadrilateral")

#
linear_quadrilateral = QuadrilateralLinear4(2,2)
Y = [x3 x1 x2 [1.0;1.0;-1.0]]
quad_interp = Y*linear_quadrilateral.interpolation
scatter3d(quad_interp[1,:],quad_interp[2,:],quad_interp[3,:],label="Gauss Points")
scatter3d!(Y[1,:],Y[2,:],Y[3,:],label="Nodes")
plot3d!(Y[1,idx],Y[2,idx],Y[3,idx],label="Global Quadrilateral")

#
# ```math
#   \int_{\Gamma_e} 1 \ \mathrm{d}\Gamma_\mathbf{y} = \int_{-1}^1\int_{-1}^1 \text{jacobian}(u,v)\ \mathrm{d}u\mathrm{d}v \approx \sum_{i=1}^{Q} \text{jacobian}(u_i,v_i) w_i,
# ```
# where $n_g$ is the number of Gaussian points $(u_i, v_i)$ and $w_i$ is the corresponding weights. The Jacobian is computed the same as for the triangle.
# The exact area of a flat quadrilateral is to the sum of the areas of the two flat triangles that make up the quadrilateral.
# ```math
#   \text{Area}(\square) = \frac{1}{2}\|(\mathbf{x}^e_4 - \mathbf{x}^e_1)\times(\mathbf{x}^e_2-\mathbf{x}^e_1)\|_2 + \frac{1}{2}\|(\mathbf{x}^e_4 - \mathbf{x}^e_3)\times(\mathbf{x}^e_2-\mathbf{x}^e_3)\|_2.
# ```
A_exact = (norm(cross(Y[:,4]-Y[:,1],Y[:,2]-Y[:,1])) + norm(cross(Y[:,4]-Y[:,3],Y[:,2]-Y[:,3])))/2
# To compute the are using the integral we first need to compute the gradient
dU = Y*linear_quadrilateral.derivatives_u
dV = Y*linear_quadrilateral.derivatives_v
jacobian = similar(linear_quadrilateral.weights)
for i = eachindex(jacobian)
    jacobian[i] = norm(cross(dU[:,i],dV[:,i]))
end
A_approx = dot(jacobian,linear_quadrilateral.weights)
# Again, since we're dealing with a flat triangle with a constant jacobian the approximated integral is equal to the exact solution
A_exact - A_approx

# # Discontinuous Basis Elements
# These are only relevant for interpolating a physical quantity - and not the geometry. Additonally to the same information that the continuous elements have they also have an additional parameter `beta` defining the position of the nodes on the element.
beta = 0.1
disclinear_triangular = DiscontinuousTriangularLinear(linear_triangular,beta)
scatter(disclinear_triangular.gauss_u, disclinear_triangular.gauss_v,
        label = "Gauss Points", aspect_ratio = :equal, xlabel="u",ylabel="v")
# Again the most interesting part is the nodal positions - which in this case is now inside of the element
disc_nodal_triangle = DiscontinuousTriangularLinear(disclinear_triangular,beta)
BoundaryIntegralEquations.set_nodal_interpolation!(disc_nodal_triangle) # Should be the identity matrix
# Using this we can plot the full reference triangle
scatter!(disc_nodal_triangle.gauss_u, disc_nodal_triangle.gauss_v,
        label="Nodal Points")
plot!(Shape(nodal_triangle.gauss_u, nodal_triangle.gauss_v),
        fillalpha=0.2, label="Reference Triangle")
# When discritizing boundary integrals over surfaces in 3d an important aspect is the mapping from the local $(u,v)$-coordinates to global $(x,y)$ coordinates. This transformation is linear with respect to the basis functions, but not nessecrarily the local coordinates, and is given by
# ```math
#   \mathbf{x}^e = \mathbf{X}^e\mathbf{N}(u,v).
# ```
# As an example we take the triangle with corners equal to the standard basis vectors
x1 = [1.0;0.0;0.0]
x2 = [0.0;1.0;0.0]
x3 = [0.0;0.0;1.0]
X = [x1 x2 x3]
# Using the mapping from local to global coordinates the Gaussian points can be mapped onto the $(x,y)$-space as
interp = X*disclinear_triangular.interpolation
G = X*disc_nodal_triangle([0.0;1.0;0.0]',[0.0;0.0;1.0]')
idx = [1;2;3;1]
scatter3d(interp[1,:],interp[2,:],interp[3,:],label="Gauss Points")
scatter3d!(X[1,:],X[2,:],X[3,:], label="Nodal Points")
plot3d!(G[1,idx],G[2,idx],G[3,idx],label="Global Triangle")
