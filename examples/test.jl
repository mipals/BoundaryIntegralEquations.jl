using SparseArrays

function fill_matrix!(M)
    for i = 1:size(M,1), j = 1:size(M,2)
        M[i,j] = 1.0
    end
end
function fill_matrix_partially!(M,p=0.01)
    for i = 1:size(M,1), j = 1:Int(ceil(size(M,2)*p))
        M[i,j] = 1.0
    end
end

A = zeros(3000,3000)
@time fill_matrix!(A)
@time fill_matrix_partially!(A)

B = spzeros(size(A))
@time fill_matrix_partially!(B,0.05)
# We definately can not just start with a sparse matrix!
# We need to compute Coordinates I,J,V and build sparse(I,J,V)...
