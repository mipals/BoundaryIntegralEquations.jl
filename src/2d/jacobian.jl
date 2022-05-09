function jacobian_mul_weights!(j,w)
    @inbounds for i = 1:length(j) 
        j[i] = j[i]*w[i]
    end
end
function compute_distance!(x,y,r)
    @inbounds for i = 1:length(r) 
        r[i] = hypot(x[1,i] - y[1], x[2,i] - y[2]) 
    end
end
function jacobian!(j,n)
     @inbounds for i = 1:length(j) 
        j[i] = hypot(n[1,i], n[2,i])     
    end 
end
function normalize!(n,j)
    @inbounds for i = 1:length(j) 
        n[1,i] = n[1,i]/j[i]
        n[2,i] = n[2,i]/j[i] 
    end
end
function compute_normals!(n,dX)
    @inbounds for i = 1:size(n,2) 
        n[1,i] = -dX[2,i]
        n[2,i] =  dX[1,i]     
    end
end
