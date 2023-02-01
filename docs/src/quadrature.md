# Quadrature
The underlying mathematics of the BEM is that of surface integrals. In most cases it is not possible to compute these analytically. Instead, a quadrature scheme is used to approximate the integrals. To do so the element in global coordinates is transformed into local coordinates using the Jacobian function describing the local deformation from the transformation. For 2D and 3D this deformation looks as follows

```math
\begin{alignat*}{3}
    \text{2D:} & \quad \text{jacobian}(u)   &&= \left\|\mathbf{X}^e \frac{\mathrm{d}\mathbf{N}(u)}{\mathrm{d}u}\right\|_2 \quad && \text{(length deformation)}\\
    \text{3D:} & \quad \text{jacobian}(u,v) &&= \left\|\left(\mathbf{X}^e \frac{\mathrm{d}\mathbf{N}(u,v)}{\mathrm{d}u}\right) \times \left(\mathbf{X}^e \frac{\mathrm{d}\mathbf{N}(u,v)}{\mathrm{d}v}\right) \right\|_2 \quad && \text{(area deformation)}
\end{alignat*}
```

Applying this the surface integral can be transformed into local coordinates and computed using a suitable quadrature scheme as

```math
    \int_{\Gamma^e} f(\mathbf{x})\ \mathrm{d}S_\mathbf{x} = \int_\mathcal{L} \text{jacobian}(\mathbf{u})f(\mathbf{u})\ \mathrm{d}\mathbf{u} \approx \sum_{i=1}^{Q} \text{jacobian}(\mathbf{u}_i)w_if(\mathbf{u}_i),
```

where $\mathbf{u}_i$ is the quadrature point with corresponding weights $w_i$.