# Reduced Order Series Expansion Boundary Element Method (ROSEBEM)
This page includes (some of) the theoretical foundations speeding up multifrequency problems using the Reduced Order Series Expansion Boundary Element Method (ROSEBEM) [panagiotopoulos2020](@cite) [Paltorp2023](@cite).

For an example [look here.](examples/3d_rosebem.md)

## The ROSEBEM
In [panagiotopoulos2020](@cite) the boundary integrals is made frequency independent by utilizing a Taylor expansion of w.r.t the wavenumber for the Green's function and its normal derivative
```math
\begin{equation*}
    G(\mathbf{x},\mathbf{y}) \approx \sum_{m=0}^{N_v - 1} \frac{(k - k_0)^m}{m!}\left[G^{(m)}(\mathbf{x},\mathbf{y},k_0)\right],\quad \frac{\partial G(\mathbf{x},\mathbf{y})}{\partial\mathbf{n}(\mathbf{y})} \approx \sum_{m=0}^{N_p-1} \frac{(k - k_0)^m}{m!}\left[\mathbf{n}(\mathbf{y})^\top\nabla G^{(m)}(\mathbf{x},\mathbf{y},k_0)\right],
\end{equation*}
```
where $N_v$ and $N_p$ is the number of terms included in the Taylor expansion of respectively the Green's function and its normal derivative while $k_0$ is the expansion wavenumber. For simplification purposes, it is in the following chosen to set $N_p=N_v=M$. Inserting the above into the Kirchhoff-Helmholtz integral equation the following system of equations will appear
```math
\begin{equation}
    \left(\text{diag}(\mathbf{\zeta}) + \sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{F}_m(k_0)\right)\mathbf{p} \approx    
    \mathrm{i} \rho_0 c k\left(\sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{G}_m(k_0) \right)\mathbf{v}_n,
\end{equation}
```
which it sometimes referred to as the series expansion boundary element method or SEBEM for short. The important aspect of above equation is that the two series of matrices, $\mathbf{F}_m$ and $\mathbf{G}_m$, depend only on the expansion wavenumber $k_0$ and can therefore be computed in a so-called offline stage. After applying boundary conditions the system will have the form of 
```math
\begin{equation}
    \left(\sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{A}_m(k_0)\right)\mathbf{z} \approx \sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{b}_m(k_0),
\end{equation}
```
where $\mathbf{z}$ is a vector containing nodal values of either $\mathbf{p}$ or $\mathbf{v}_n$. Unfortunately, the equation above does not resolve all the issues, as it only reduces the assembly from a discretization of boundary integrals to a sum of matrices that still scales as $O(n^2)$ in both computation and storage. Additionally, for large systems, the time spent solving the linear system of equations is overshadowed by the time spent assembling the linear system, making the gain of using the equation above insignificant with respect to the total computational time. 

In order to reduce both the memory footprint and the assembly time of SEBEM a Galerkin projection can be applied. The main assumption behind this projection is that the solution, as well as the right-hand side of SEBEM is spanned by a lower dimensional subspace. This means that
```math
\begin{equation}
    \mathbf{z} \approx \mathbf{U}_\ell\mathbf{z}_\ell, \quad \mathbf{b}_m \approx \mathbf{U}_\ell \mathbf{b}_{\ell m}, \quad \ \mathbf{U}_\ell, \in \mathbb{C}^{N\times\ell}
\end{equation}
```
will be good approximations. For simplification purposes, it is often chosen to let $\mathbf{U}_\ell$ be unitary. An intuitive explanation as to why this must be a good approximation is that most BE systems are well-conditioned and iterative solvers such as the generalized minimal residual method (GMRes) converge quickly. In addition, this motivates a good place to search for both $\mathbf{U}_\ell$ is the $\ell$-Krylov subspace defined by [hdgm15](@cite)
```math
\begin{equation}
    \mathcal{K}_\ell\left(\mathbf{A}(k), \mathbf{b}(k)\right) = \text{span}\left\{\mathbf{b}(k),\ \mathbf{A}(k)\mathbf{b}(k),\ \mathbf{A}(k)^2\mathbf{b}(k),\ \dots,\ \mathbf{A}(k)^{\ell-1}\mathbf{b}(k) \right\}.
\end{equation}
```
Setting the columns of $\mathbf{U}_\ell$ equal to the Krylov vectors would result in the lower dimensional subspace approximation being a good approximation for wavenumbers close to $k$. The subspace can be made to handle a wider range of wavenumbers by combining the Krylov subspace for the wavenumbers $k_1,k_2,\dots k_L$. As such the final projection matrix, $\mathbf{U}_\ell$, is computed as the singular value decomposition (SVD) of the concatenation of the columns of all the $q$-Krylov vector subspaces 
```math
\begin{equation}
    \mathbf{U}_\ell\Sigma_\ell \mathbf{V}_\ell^{\mathsf{H}} = \text{svd}\left( \begin{bmatrix}\mathbf{K}_{k_1}& \mathbf{K}_{k_2}& \dots & \mathbf{K}_{k_L}\end{bmatrix} \right),
\end{equation}
```
where $\mathbf{K}_{k_i}$ denotes a matrix with columns equal to some Krylov vectors at wavenumber $k_i$. Now inserting the lower dimensional subspace approximatiion into SEBEM equation while also multiplying both sides with the hermitian transpose of $\mathbf{U}_\ell$ from the left it follows that
```math
\begin{equation}
    \left(\sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{U}_\ell^\mathsf{H}\mathbf{A}_m(k_0)\mathbf{U}_\ell\right)\mathbf{z}_\ell \approx \sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{U}_\ell^\mathsf{H}\mathbf{U}_\ell\mathbf{b}_{\ell m}(k_0).
\end{equation}
```
Using that $\mathbf{U}_\ell$ is unitary while defining $\mathbf{A}_{\ell m} = \mathbf{U}_\ell^\mathsf{H}\mathbf{A}_m(k_0)\mathbf{U}_\ell$ the above becomes
```math
\begin{equation}
    \left(\sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{A}_{\ell m}(k_0)\right)\mathbf{z}_\ell \approx \sum_{m=0}^{M-1}\frac{(k - k_0)^m}{m!}\mathbf{b}_{\ell m}(k_0).
\end{equation}
```
The important part of the above is that $\mathbf{A}_{\ell m}(k_0) \in \mathbb{C}^{\ell\times\ell}$ and $\mathbf{b}_{\ell m}(k_0)\in\mathbb{C}^\ell$, which means that if $\ell \ll n$ this representation requires significantly less memory than the original system. Note, however, that $\mathbf{A}_m(k_0)$ should never explicitly be stored. Instead, the rows should be projected during the assembly phase.

## Visualization
In many cases visualizations of the mathematics gives rise to a deeper understanding. Therefore, a short visual introduction to the ROSEBEM incorporating the BLI boundary condition is given. We start by introducing the following visual representation of the vector spaces
![Step1](figures/rosebem_1.png)
Using this we can visualize the SEBEM as
![Step2](figures/rosebem_2.png)
Similarly, the reduced basis assumption, for both $\mathbf{p}$ and $\partial_\mathbf{n}\mathbf{p}$, is visualized as
![Step3](figures/rosebem_3.png)
Inserting assumption from above into the SEBEM visualization while also multiplying with the Hermitian transpose of the reduced basis 
![Step4](figures/rosebem_4.png)
Finally by performing the matrix-matrix, the ROSEBEM appears as
![Step5](figures/rosebem_5.png)

## Bibliography
```@bibliography
Pages = []
Canonical = false

panagiotopoulos2020
Paltorp2023
hdgm15
```
