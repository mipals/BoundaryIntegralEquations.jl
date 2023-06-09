#==========================================================================================
        Acoustical Green's Function and Derivatives - Vectors (used for collocation)
==========================================================================================#
"""
    greens3d!(integrand,r,k)

Green's function for the Helmholtz Equation in 3d:

```math
\\frac{e^{ikr_j}}{4\\pi r_j}
```
"""
function greens3d!(integrand,r,k)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = exp(im*k*r[i,j])/(4π*r[i,j])
    end
    return integrand
end

"""
    freens3d!(integrand,r,interpolation,sources,normals,k)

Normal derivative of the 3D Helmholtz Green's function with respect to interpolation nodes:

```math
\\frac{e^{ikr_j}}{4\\pi r_j^3}(ikr_j - 1) (x_j - y)\\cdot n
```
"""
function freens3d!(integrand,r,interpolation,sources,normals,k)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = exp(im*k*r[i,j])*(im*k*r[i,j] - 1)*
                        (normals[1,i]*(interpolation[1,i] - sources[1,j]) +
                         normals[2,i]*(interpolation[2,i] - sources[2,j]) +
                         normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4π*r[i,j]^3)
    end
    return integrand
end

"""
    freens3dk0!(integrand,r,interpolation,sources,normals)

freens3d! with k=0.

```math
-\\frac{(x_j - y)\\cdot n}{4\\pi r_j}
```
"""
function freens3dk0!(integrand,r,interpolation,sources,normals)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] =  -(normals[1,i]*(interpolation[1,i] - sources[1,j]) +
                            normals[2,i]*(interpolation[2,i] - sources[2,j]) +
                            normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4π*r[i,j]^3)
    end
    return integrand
end

"""
    freens3dk0_to_freens3d!(integrand,r,k)

Multiplies each term in `integrand` from `freens3dk0!` such that they equal `freens3d!`.
"""
function freens3dk0_to_freens3d!(integrand,r,k)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = integrand[i,j]*exp(im*k*r[i,j])*(1 - im*k*r[i,j])
    end
    return integrand
end
#==========================================================================================
                                    Kernels for debugging.
==========================================================================================#
"""
    onefunction!(integrand,r,k)
    onefunction!(integrand,r,interpolation,source,normals)
    onefunction!(integrand,r,interpolation,source,normals,k)

Fills `integrand` with ones. Other inputs disregarded. Can be used to compute surface areas.
"""
onefunction!(integrand,r,k)                              = fill!(integrand,1)
onefunction!(integrand,r,interpolation,source,normals)   = fill!(integrand,1)
onefunction!(integrand,r,interpolation,source,normals,k) = fill!(integrand,1)

#==========================================================================================
                                Taylor Expansions of Kernels
==========================================================================================#
"""
    taylor_greens3d!(integrand,r,k,m)

``m``th derivative of the Green's function w.r.t. the wavenumber.

```math
G^{(m)}(r_j,k_0) = \\frac{(\\mathrm{i}r)^m\\exp(\\mathrm{i}k_0r)}{4\\pi r_j}
```
"""
function taylor_greens3d!(integrand,r,k0,m)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = (im*r[i,j])^m*exp(im*k0*r[i,j])/(4π*r[i,j])
    end
    return integrand
end

"""
    taylor_freens3d!(integrand,r,interpolation,sources,normals,k0,m)

``m``th derivative of the normal derivative of the Green's function w.r.t. the wavenumber.

```math
n\\cdot \\nabla G^{(m)}(r_j,k_0)
```
"""
function taylor_freens3d!(integrand,r,sources,collocation,normals,k0,m)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = exp(im*k0*r[i,j])*(im*r[i,j])^m*(1 - im*k0*r[i,j] - m)*
                        (normals[1,i]*(collocation[1,j] - sources[1,i]) +
                         normals[2,i]*(collocation[2,j] - sources[2,i]) +
                         normals[3,i]*(collocation[3,j] - sources[3,i]))/(4π*r[i,j]^3)
    end
    return integrand
end

"""
    taylor_greens_gradient3d!(integrand,r,interpolation,sources,normals,k0,m)

Gradient of the ``m``th derivative of the Green's function w.r.t. the wavenumber.

```math
\\nabla G^{(m)}(r_j,k_0)
```
"""
function taylor_greens_gradient3d!(gradient,r,collocation,sources,k0,m)
    @fastmath @inbounds for i = 1:size(gradient,2)
        gradient[1,j] = collocation[1,j] - sources[1,i]
        gradient[2,j] = collocation[2,j] - sources[2,i]
        gradient[3,j] = collocation[3,j] - sources[3,i]
        gradient[:,j] .*= exp(im*k0*r[i,j])*(im*r[i,j])^m*(1 - im*k0*r[i,j] - m)/(4π*r[i,j]^3)
    end
    return gradient
end

"""
    taylor_greens_tangential_gradient3d!(integrand,r,interpolation,sources,normals,k0,m)

Tangential gradient of the ``m``th derivative of the Green's function w.r.t. the wavenumber.

```math
(\\mathbf{I} - \\mathbf{n}\\mathbf{n}^\\top)\\nabla G^{(m)}(r_j,k_0)
```
"""
function taylor_greens_tangential_gradient3d!(gradient,r,collocation,sources,normals,k0,m)
    taylor_greens_gradient3d!(gradient,r,collocation,sources,k0,m)
    @fastmath @inbounds for j = 1:size(gradient,2)
        gradient[:,j] .-= normals[:,j]*(normals[:,j]*gradient[:,j])
    end
    return gradient
end
