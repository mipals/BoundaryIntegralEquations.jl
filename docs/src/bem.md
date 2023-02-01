# The Boundary Element Method
In simple terms the Boundary Element Method (BEM) is a method for solving Boundary Integral Equations through a discretization of both the unknown function and domain. In this thesis the main focus is solving the BIE through the so-called collocation approach. This approach is a particular case of the Galerkin approach where the test function is equal to the sum of Dirac-delta functions

```math
\phi(\mathbf{y}) =  \mathbf{a}^\top
                    \begin{bmatrix}
                        \delta\left(\mathbf{y} - \mathbf{z}_1\right)\\ 
                        \delta\left(\mathbf{y} - \mathbf{z}_2\right)\\ 
                        \vdots\\
                        \delta\left(\mathbf{y} - \mathbf{z}_n\right)\\ 
                    \end{bmatrix},
```

where ``\mathbf{a} \in \mathbb{C}^n`` is a vector of arbitrary constants and ``\mathbf{z}_i \in \mathbb{R}^3`` are the nodal mesh positions. 

```math
\mathbf{a}^\top
\begin{bmatrix}
    c(\mathbf{z}_1)p(\mathbf{z}_1)\\ 
    c(\mathbf{z}_2)p(\mathbf{z}_2)\\ 
    \vdots \\ 
    c(\mathbf{z}_n)p(\mathbf{z}_n)
\end{bmatrix} 
+ 
\mathbf{a}^\top
\begin{bmatrix}
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_1)}{\partial\mathbf{n}(\mathbf{x})}p(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_2)}{\partial\mathbf{n}(\mathbf{x})}p(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_n)}{\partial\mathbf{n}(\mathbf{x})} p(\mathbf{x})\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}
=
\mathbf{a}^\top
\begin{bmatrix}
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_1)v_f(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_2)v_f(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_n)v_f(\mathbf{x})\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}.
```

In its current form the solution to the above equation is a \text{function} which a computer can not be tasked to find. Instead, the problem is made computationally tractable by parametrizing the functions ``p`` and ``v_s`` as shown previously. As such the above reduces to

```math
\mathbf{a}^\top
\begin{bmatrix}
    c(\mathbf{z}_1)\mathbf{T}(\mathbf{z}_1)\mathbf{p}\\ 
    c(\mathbf{z}_2)\mathbf{T}(\mathbf{z}_2)\mathbf{p}\\ 
    \vdots \\ 
    c(\mathbf{z}_n)\mathbf{T}(\mathbf{z}_n)\mathbf{p}
\end{bmatrix}
+ 
\mathbf{a}^\top
\begin{bmatrix}
    \int_\Gamma \frac{\partial G(\mathbf{x}, \mathbf{z}_1)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\mathbf{p}\ \mathrm{d}S_\mathbf{x}\\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_2)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\mathbf{p}\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_n)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\mathbf{p}\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}
\approx
\mathbf{a}^\top
\begin{bmatrix}
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_1)\mathbf{T}(\mathbf{x})\mathbf{v}_f\ \mathrm{d}S_\mathbf{x}\\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_2)\mathbf{T}(\mathbf{x})\mathbf{v}_f\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_n)\mathbf{T}(\mathbf{x})\mathbf{v}_f\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}.
```

By the cardinal property of the basis functions it follows that

```math
\mathbf{a}^\top
\left(
\text{diag}\left(\begin{bmatrix}
    c(\mathbf{z}_1)\\ 
    c(\mathbf{z}_2)\\ 
    \vdots \\ 
    c(\mathbf{z}_n)
\end{bmatrix}\right)
+ 
\begin{bmatrix}
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_1)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_2)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    \int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_n)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}\right)\mathbf{p}
\approx
\mathbf{a}^\top
\begin{bmatrix}
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_1)\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_2)\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}\\
    \vdots \\
    sk\int_\Gamma G(\mathbf{x},\mathbf{z}_n)\mathbf{T}(\mathbf{x})\ \mathrm{d}S_\mathbf{x}
\end{bmatrix}\mathbf{v}_f.
```

Since the above has to hold for all ``\mathbf{a} \in \mathbb{C}^n`` it follows that
```math
\left(\text{diag}\left(\begin{bmatrix}
    c(\mathbf{z}_1)\\ 
    c(\mathbf{z}_2)\\ 
    \vdots \\ 
    c(\mathbf{z}_n)
\end{bmatrix}\right)
+ 
\mathbf{F}\right)\mathbf{p}
=
sk\mathbf{G}\mathbf{v}_s.
```

For simplification purposes the above is sometimes written as

```math
\mathbf{H}\mathbf{p} + \mathbf{G}\partial_\mathbf{n}\mathbf{p}.
```

We now briefly explain how to compute the ``j``th row of ``\mathbf{H}``. To do so element description of the unknown function used differently as ``\mathbf{p}^e`` is not part of the equation. Instead, we use the split as follows 

```math
\begin{aligned}
\int_\Gamma\frac{\partial G(\mathbf{x}, \mathbf{z}_j)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})\mathbf{p}\ \mathrm{d}S_\mathbf{x} 
&\approx \sum_{e=1}^{N}\left(\int_{\Gamma^e}\frac{\partial G(\mathbf{x}, \mathbf{z}_j)}{\partial\mathbf{n}(\mathbf{x})}\mathbf{T}(\mathbf{x})(\mathbf{L}^e)^\top\mathbf{L}^e\ \mathrm{d}S_\mathbf{x}\right) \mathbf{p}\\
&\approx \left(\underbrace{\sum_{e=1}^{N}\left(\sum_{i=1}^{Q(\mathbf{z}_j,e)}\frac{\partial G(\mathbf{x}^e(\mathbf{u}_i), \mathbf{z}_j)}{\partial\mathbf{n}(\mathbf{x}^e(\mathbf{u}_i))}\text{jacobian}(\mathbf{u}_i)w_i\mathbf{T}^e(\mathbf{u}_i)\right)\mathbf{L}^e}_{j\text{th row of }\mathbf{H}}\right)\mathbf{p},
\end{aligned}
```

where the number of quadrature points, ``Q(\mathbf{z}_1,e)``, is depended on element location with respect to the collocation point ``\mathbf{z}_1``. Note that approximation in the first line of the above is due to approximation of the geometry using elements while the second approximation is due to the integration error. A similar approach can be applied for the rows of ``\mathbf{H}``.

A downside of the BEM is that the resulting matrices ``\mathbf{H}, \mathbf{G} \in \mathbb{C}^{n\times n}`` are dense, meaning that storing the two matrices scales as ``\mathcal{O}(n^2)`` rendering the direct application of the method unusable large ``n``.