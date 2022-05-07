#==========================================================================================
                             Kernels in 2D                                    
==========================================================================================#
# Green's function
function G_kernel2d(x,y,k,normals)
    r = sqrt.((sum((x .- y).^2, dims=1)))
    return im*0.25*hankelh1.(0,k*r)
end
# Directional derivative of the Green's function
function F_kernel2d(x,y,k,normals)
    r = sqrt.((sum((x .- y).^2, dims=1)))
    return -im*0.25*k*hankelh1.(1,k*r) .* sum(normals .* (x .- y),dims=1)./r
end

#==========================================================================================
                            Mutating kernels in 2D                          
==========================================================================================#
# Green's function
function G_kernel2d!(x,y,k,normals,r)
    r .= sqrt.((sum((x .- y).^2, dims=1))) 
    r .= im*0.25*hankelh1.(0,k*r)
end
# Directional derivative of the Green's function
function F_kernel2d!(x,y,k,normals,r)
    r .= sqrt.((sum((x .- y).^2, dims=1)))
    r .= -im*0.25*k*hankelh1.(1,k*r) .* sum(normals .* (x .- y),dims=1)./r
end
