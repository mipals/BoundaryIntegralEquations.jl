#==========================================================================================
        Acoustical Green's Function and Derivatives - Vectors (used for collocation)
==========================================================================================#
"""
    greens3d!(integrand,r,k)

Green's function for the Helmholtz Equation in 3d.
"""
function greens3d!(integrand,r,k)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = exp(im*k*r[i,j])/(4.0*π*r[i,j])
    end
end
"""
    freens3d!(integrand,r,interpolation,source,normals,k)

Normal derivative of the 3D Helmholtz Green's function with respect to interpolation nodes.
"""
function freens3d!(integrand,r,interpolation,sources,normals,k)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = exp(im*k*r[i,j])*(1 - im*k*r[i,j])*
                        (normals[1,i]*(interpolation[1,i] - sources[1,j]) +
                         normals[2,i]*(interpolation[2,i] - sources[2,j]) +
                         normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4π*r[i,j]^3)
    end
end
"""
    freens3dk0!(integrand,r,interpolation,source,normals)

freens3d! with k=0.
"""
function freens3dk0!(integrand,r,interpolation,sources,normals)
    @fastmath @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = -(normals[1,i]*(interpolation[1,i] - sources[1,j]) +
                           normals[2,i]*(interpolation[2,i] - sources[2,j]) +
                           normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4π*r[i,j]^3)
    end
end
"""
    freens3dk0_to_freens3d!(integrand,r,k)

Multiplies each term in `integrand` from `freens3dk0!` such that they equal `freens3d!`.
"""
function freens3dk0_to_freens3d!(integrand,r,k)
    @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
        integrand[i,j] = integrand[i,j]*exp(im*k*r[i,j])*(1 - im*k*r[i,j])
    end
end

#==========================================================================================
                                    Kernels for debugging.
==========================================================================================#
"""
    onefunction!

Fills `integrand` with ones. Can be used to compute surface areas.
"""
onefunction!(integrand,r,k)                              = fill!(integrand,1.0)
onefunction!(integrand,r,interpolation,source,normals)   = fill!(integrand,1.0)
onefunction!(integrand,r,interpolation,source,normals,k) = fill!(integrand,1.0)



#==========================================================================================
        Acoustical Green's Function and Derivatives - Matrices (used for Galerkin)
==========================================================================================#
# """
#     matrix_greens3d!(integrand,r,k)

# Green's function for the Helmholtz Equation in 3d.
# """
# function matrix_greens3d!(integrand,r,k)
#     @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
#         integrand[i,j] = exp(-im*k*r[i,j])/(4.0*π*r[i,j])
#     end
# end
# """
#     matrix_freens3d!(integrand,r,interpolation,source,normals,k)

# Normal derivative of the 3D Helmholtz Green's function with respect to interpolation nodes.
# """
# function matrix_freens3d!(integrand,r,interpolation,sources,normals,k)
#     @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
#         integrand[i,j] = -exp(-im*k*r[i])*(1.0 + im*k*r[i])*
#                         (normals[1,i]*(interpolation[1,i] - sources[1,j]) +
#                          normals[2,i]*(interpolation[2,i] - sources[2,j]) +
#                          normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4.0*π*r[i,j]^3)
#     end
# end
# """
#     matrix_freens3dk0!(integrand,r,interpolation,source,normals)

# matrix_freens3d! with k=0.
# """
# function matrix_freens3dk0!(integrand,r,interpolation,sources,normals)
#     @inbounds for i = 1:size(integrand,1), j = 1:size(integrand,2)
#         integrand[i,j] = -(normals[1,i]*(interpolation[1,i] - sources[1,j]) +
#                            normals[2,i]*(interpolation[2,i] - sources[2,j]) +
#                            normals[3,i]*(interpolation[3,i] - sources[3,j]))/(4.0*π*r[i,j]^3)
#     end
# end
