# # 3D Element Usage
# For 3D boundary integral equations the elements of interest are `SurfaceFunctions`.
# These are defined by two local coordinates, $u$ and $v$
# ```math
#   \mathbf{N}(u,v)
# ```

import WGLMakie as Mke
# Set the default resolution to something that fits the Documenter theme
Mke.set_theme!(resolution=(800, 800))

using IntegralEquations
import IntegralEquations: set_nodal_interpolation!
linear_triangular = TriangularLinear(3)
# Plotting Gaussian Integration Points
Mke.scatter(linear_triangular.gauss_u, linear_triangular.gauss_v)
# We can instead
reference_triangle = TriangularLinear(3)
set_nodal_interpolation!(reference_triangle)
Mke.scatter!(reference_triangle.gauss_u, reference_triangle.gauss_v)
Mke.lines!([reference_triangle.gauss_u; 0.0],[reference_triangle.gauss_v; 0.0])
Mke.current_figure()
# Weights should add to area of triangle
sum(linear_triangular.weights) ≈ 0.5
