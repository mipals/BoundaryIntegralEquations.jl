#==========================================================================================
                             Kernels in 2D
==========================================================================================#
function G!(int,k,r)
    @fastmath @inbounds for i = 1:length(r)
        int[i] = im*0.25*hankelh1(0,k*r[i])
    end
end

function C!(int,x,y,n,r)
    @fastmath @inbounds for i = 1:length(r)
        int[i] = (n[1,i]*(x[1,i] - y[1]) +
                  n[2,i]*(x[2,i] - y[2]))/(2π*r[i]^2)
    end
end

function F!(int,x,y,k,n,r)
    @fastmath @inbounds for i = 1:length(r)
        int[i] = -im*hankelh1(1,k*r[i])*k*(n[1,i]*(x[1,i] - y[1]) +
                                        n[2,i]*(x[2,i] - y[2]))/(4r[i])
    end
end

function C_to_F!(int,k,r)
    @fastmath @inbounds for i = 1:length(r)
        int[i] = -π*k*r[i]*im*int[i]*hankelh1(1,k*r[i])/2
    end
end
