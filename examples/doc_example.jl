# # Integration on an annulus
# In this example, we explore integration of the function:
# ```math
#   f(x,y) = \frac{x^3}{x^2+y^2-\frac{1}{4}},
# ```
# over the annulus defined by $\{(r,\theta) : \rho < r < 1, 0 < \theta < 2\pi\}$
# with parameter $\rho = \frac{2}{3}$. We will calculate the integral:
# ```math
#   \int_0^{2\pi}\int_{\frac{2}{3}}^1 f(r\cos\theta,r\sin\theta)^2r{\rm\,d}r{\rm\,d}\theta,
# ```
# by analyzing the function in an annulus polynomial series.
# We analyze the function on an $N\times M$ tensor product grid defined by:
# ```math
# \begin{aligned}
# r_n & = \sqrt{\cos^2\left[(n+\tfrac{1}{2})\pi/2N\right] + \rho^2 \sin^2\left[(n+\tfrac{1}{2})\pi/2N\right]},\quad{\rm for}\quad 0\le n < N,\quad{\rm and}\\
# \theta_m & = 2\pi m/M,\quad{\rm for}\quad 0\le m < M;
# \end{aligned}
# ```


using JSServe # hide
Page(exportable=true, offline=true) # hide
# Our function $f$ on the annulus:


import WGLMakie as Mke
# Set the default resolution to something that fits the Documenter theme
Mke.set_theme!(resolution=(800, 400))
Mke.scatter(1:4, color=1:4)
