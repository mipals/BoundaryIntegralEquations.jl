# Getting the computer to understand functions
In the most basic terms computers are only capable of understanding numbers, meaning that they are inherently unable to solve equations where the unknowns are functions. This is a problem when trying to solve BIEs. To resolve this issue the functions are instead approximated using parametrizations for which the coefficients (numbers) are unknown. Intuitively these numbers are exactly what the computer is then tasked to find. For element methods this parametrization is chosen to be the simplest possible: A linear combination of functions

```math
p(\mathbf{x}) \approx \mathbf{T}(\mathbf{x})\mathbf{p} 
    = 
    \begin{bmatrix} T_1(\mathbf{x}) &  T_2(\mathbf{x}) &  \dots  &  T_n(\mathbf{x}) \end{bmatrix}
    \begin{bmatrix} p_1             \\ p_2             \\ \vdots \\ p_n             \end{bmatrix}. 
```

Note that linearity is with respect to the unknown parameters ``\mathbf{p}``, but not necessarily in the known basis functions ``\mathbf{T}(\mathbf{x})``. A useful, but not necessary, property is called the Cardinal property stating that 

```math
T_i(\mathbf{x}_j) = \delta_{ij} = 
    \begin{cases}
        1, \quad i = j\\
        0, \quad i \neq j    
    \end{cases}.
```

The usefulness is due the property making the ``j``th value of ``\mathbf{p}`` interpretable as the value of the function at point ``\mathbf{x}_j`` as

```math
\mathbf{T}(\mathbf{x}_j)\mathbf{p} = p_j.
```

One might ask: *How does the above relate to The Boundary Element Method?* The answer is that the functions ``\mathbf{T}_i`` are chosen to be simpler functions with support equal to only a few subdomains of the original domain. These subdomains are commonly referred to as elements. 