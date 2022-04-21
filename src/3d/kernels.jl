#==========================================================================================
                    Acoustical Green's Function and Derivatives
==========================================================================================#
"""
    greens3d!(integrand,r,k)

Green's function for the Helmholtz Equation in 3d.
"""
function greens3d!(integrand,r,k)
    @inbounds for i = 1:length(r)
        integrand[i] = exp(-im*k*r[i])/(4.0*π*r[i])
    end
end
"""
    freens3d!(integrand,r,interpolation,source,normals,k)

Normal derivative of the 3D Helmholtz Green's function with respect to interpolation nodes.
"""
function freens3d!(integrand,r,interpolation,source,normals,k)
    @inbounds for i = 1:length(integrand)
        integrand[i] = -exp(-im*k*r[i])*(1.0 + im*k*r[i])*
                        (normals[1,i]*(interpolation[1,i] - source[1]) + 
                         normals[2,i]*(interpolation[2,i] - source[2]) + 
                         normals[3,i]*(interpolation[3,i] - source[3]))/(4.0*π*r[i]^3)
    end
end
"""
    freens3dk0!(integrand,r,interpolation,source,normals)

freens3d! with k=0.
"""
function freens3dk0!(integrand,r,interpolation,source,normals)
    @inbounds for i = 1:length(r)
        integrand[i] = -(normals[1,i]*(interpolation[1,i] - source[1]) + 
                         normals[2,i]*(interpolation[2,i] - source[2]) + 
                         normals[3,i]*(interpolation[3,i] - source[3]))/(4.0*π*r[i]^3)
    end
end
#==========================================================================================
                                    Kernels for debugging. 
==========================================================================================#
"""
    onefunction!

Fills `integrand` with ones. Used for computing surface areas. 
"""
onefunction!(integrand,r) = fill!(integrand,1.0)