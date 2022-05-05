#==========================================================================================
                Generating arbitrary orders of Lagrange Quadrilateral Elements
==========================================================================================#
using Symbolics, LinearAlgebra
# Number polynomial contributions
n = 4
# Number of interpolation nodes
N = n^2
# Order of the Lagrange polynomials
order = n - 1
# Defining a general Lagrange Polynomial
func_lagrange(x,y) = kron(x.^(0:order),y.^(0:order))

# Defining Nodal Nodes
x_lagrange = kron(ones(n),collect(range(-1,1,length=n)))
y_lagrange = kron(collect(range(-1,1,length=n)),ones(n))
# Interpolating the general Lagrange Polynomial on the Nodal nodes
X_lagrange = zeros(N,N)
for i = 1:N
    X_lagrange[:,i] = func_lagrange(x_lagrange[i],y_lagrange[i])
end
# Computing coefficients for the Lagrange Polynomials
C_lagrange = X_lagrange\I
# Checking if the interpolation is correct (testing)
Y_lagrange = similar(X_lagrange)
for i = 1:N
    Y_lagrange[:,i] = C_lagrange*func_lagrange(x_lagrange[i],y_lagrange[i])
end
isapprox(Y_lagrange, I)

# Generating code for the generated lagrange polynomial using the Symbolics package
@variables x y
element_lagrange = C_lagrange*func_lagrange(x,y)
f_lagrange = build_function(element_lagrange,x,y)
# write("test_functions.jl", string(f_lagrange[2]))
fun_lagrange! = eval(f_lagrange[2])
# Creating an allocating alternative
function fun_lagrange(u,v,N=N)
    @assert(length(u) == length(v))
    output = zeros(N,N)
    for i = 1:length(u)
        out = @view output[:,i]
        fun_lagrange!(out, u[i],v[i])
    end
    return output
end
isapprox(fun_lagrange(x_lagrange,y_lagrange),I)

#==========================================================================================
                Generating arbitrary orders of Serendipity Quadrilateral Elements
==========================================================================================#
# Number polynomial contributions
m = 3
# Number of interpolation nodes
M = 4*m - 4
# Order of the serendipity polynomials
order = m - 1
# Defining a general serendipity polynomial
func_serendipity(x,y) = [x.^(0:order); y.*x.^(0:order); y.^(2:order); y.^(2:order).*x]

# Defining Nodal Nodes
side1 = collect(range(-1, 1,length=m))
side2 = collect(range( 1, 1,length=m))
side3 = collect(range( 1,-1,length=m))
side4 = collect(range(-1,-1,length=m))
x_serendipity = [side1; side2[2:end]; side3[2:end]; side4[2:end-1]]
y_serendipity = [side4; side1[2:end]; side2[2:end]; side3[2:end-1]]
# Interpolating the general Serendipity Polynomial on the Nodal nodes
X_serendipity = zeros(M,M)
for i = 1:M
    X_serendipity[:,i] = func_serendipity(x_serendipity[i],y_serendipity[i])
end
# Computing coefficients for the serendipity polynomials
C_serendipity = X_serendipity\I
# Checking if interpolation is correct (testing)
Y_serendipity = similar(X_serendipity)
for i = 1:M
    Y_serendipity[:,i] = C_serendipity*func_serendipity(x_serendipity[i],y_serendipity[i])
end
isapprox(Y_serendipity, I)

# Generating code for the generated lagrange polynomial using the Symbolics package
@variables s t
element_serendipity = C_serendipity*func_serendipity(s,t)
f_serendipity = build_function(element_serendipity,s,t)
fun_serendipity! = eval(f_serendipity[2])
# write("test_functions.jl", string(f_lagrange[2]))
# Creating an allocating alternative
function fun_serendipity(u,v,M=M)
    @assert(length(u) == length(v))
    output = zeros(M,M)
    for i = 1:length(u)
        out = @view output[:,i]
        fun_serendipity!(out, u[i],v[i])
    end
    return output
end
fun_serendipity(x_serendipity,y_serendipity)
