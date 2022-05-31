#==========================================================================================
                             Kernels in 2D
==========================================================================================#
function G!(k,r,int)
    @inbounds for i = 1:length(r)
        int[i] = im*0.25*hankelh1(0,k*r[i])
    end
end

function F!(x,y,k,n,r,int)
     @inbounds for i = 1:length(r)
        int[i] = -im*0.25*hankelh1(1,k*r[i])*k*(n[1,i]*(x[1,i] - y[1]) +
                                                n[2,i]*(x[2,i] - y[2]))/r[i]
    end
end
function C!(x,y,n,r,int)
    @inbounds for i = 1:length(r)
        int[i] = (n[1,i]*(x[1,i] - y[1]) +
                  n[2,i]*(x[2,i] - y[2]))/(2π*r[i]^2)
    end
end
