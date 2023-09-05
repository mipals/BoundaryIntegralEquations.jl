"""
    compute_taylor_integrals!

Approximates the integral of `freens3d!`, `greens3d!` and `freens3dk0!` multiplied by the
shapefunction over single `shape_function`, defined by `coordinates`.
"""
function compute_taylor_integrals!(physics_interpolation,interpolation_element,
                                    fOn,gOn,cOn,
                                    submatrixF,submatrixG,subvectorC,k,
                                    source,integrand,r,M)
    interpolation        = interpolation_element.interpolation
    normals              = interpolation_element.normals
    jacobian_mul_weights = interpolation_element.jacobian_mul_weights
    # Note that everything here has been written to avoid re-allocating things.
    # As such we re-use the memory of "integrand" for all computations
    compute_distances!(r,interpolation,source)
    if cOn
        ### Evaluating the G0-kernel (used for computing the C-constant)
        freens3dk0!(integrand,r,interpolation,source,interpolation_element.normals)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating the integral and adding the value to the C-vector
        sum_to_c!(subvectorC,integrand)
    end
    for m = 0:M
        if gOn
            Gview = @view submatrixG[1,:,m+1]
            ### Evaluating the G-kernel (single-layer kernel) at the global nodes
            taylor_greens3d!(integrand,r,k,m)
            # Multiplying integrand with jacobian and weights
            integrand_mul!(integrand,jacobian_mul_weights)
            # Approximating the integral and the adding the values to the BEM matrix
            mygemm_vec!(Gview,integrand,physics_interpolation)
        end
        if fOn
            Fview = @view submatrixF[1,:,m+1]
            ### Evaluating the F-kernel (double-layer kernel) at the global nodes
            taylor_freens3d!(integrand,r,interpolation,source,normals,k,m)
            # Multiplying integrand with jacobian and weights
            integrand_mul!(integrand,jacobian_mul_weights)
            # Approximating integral and adding the value to the BEM matrix
            mygemm_vec!(Fview,integrand,physics_interpolation)
        end
    end
end

"""
    taylor_assemble!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                            M=0,fOn=true,gOn=true,cOn=true,m=3,n=3,progress=true,U=I)

Return:
 * `From`: Contains the derivatives of the `F`-matrix.
 * `Grom`: Contains the derivatives of the `G`-matrix.
 * `C`: Contains the integral free term.
"""
function taylor_assemble!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                            M=0,fOn=true,gOn=true,cOn=true,m=3,n=3,progress=true,U=I)
    n_elements  = number_of_elements(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    physics_topology = mesh.physics_topology
    physics_function = mesh.physics_function
    #======================================================================================
        Introducing three elements: (The hope here is to compute singular integrals)
        The numbers corresponds to the corner for which the GP-points are clustered
        As such choose the clustering that are closest to the source point
     ————————————————————————————————————  Grid  ————————————————————————————————————————
                                           3
                                           | \
                                           1 - 2
    ======================================================================================#
    shape_function1 = create_rotated_element(shape_function,n,m,1)
    shape_function2 = create_rotated_element(shape_function,n,m,2)
    shape_function3 = create_rotated_element(shape_function,n,m,3)
    physics_function1 = deepcopy(physics_function)
    physics_function2 = deepcopy(physics_function)
    physics_function3 = deepcopy(physics_function)
    copy_interpolation_nodes!(physics_function1,shape_function1)
    copy_interpolation_nodes!(physics_function2,shape_function2)
    copy_interpolation_nodes!(physics_function3,shape_function3)

    # Computing interpolation on each element
    interpolation_list1 = interpolate_elements(mesh,shape_function1)
    interpolation_list2 = interpolate_elements(mesh,shape_function2)
    interpolation_list3 = interpolate_elements(mesh,shape_function3)

    # Copying interpolation of physics functions1
    physics_interpolation1 = copy(physics_function1.interpolation')
    physics_interpolation2 = copy(physics_function2.interpolation')
    physics_interpolation3 = copy(physics_function3.interpolation')

    # Preallocation of return values
    if typeof(U) <: UniformScaling
        n_nodes = size(mesh.sources,2)
    else
        n_nodes = size(U,2)
    end
    F = zeros(ComplexF64, n_sources, n_nodes, M+1)
    G = zeros(ComplexF64, n_sources, n_nodes, M+1)
    C = zeros(ComplexF64, n_sources)

    # Assembly loop
    if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds @threads for source_node = 1:n_sources
        # Access source
        source = sources[:,source_node]
        # Every thread has access to parts of the pre-allocated matrices
        integrand  = zeros(ComplexF64,n*m)
        r          = zeros(Float64,n*m)
        subvectorC = @view C[source_node]
        # element_coordinates = zeros(3,n_physics_functions)
        Fview = zeros(eltype(F),1,n_sources,M+1)
        Gview = zeros(eltype(G),1,n_sources,M+1)
        @inbounds for element = 1:n_elements
            # Access element topology and coordinates
            physics_nodes        = @view physics_topology[:,element]
            physics_coordinates  = @view sources[:,physics_nodes]
            # Acces submatrix of the BEM matrix
            submatrixF = @view Fview[1:1,physics_nodes,:]
            submatrixG = @view Gview[1:1,physics_nodes,:]
            # Use quadrature point clustered around the closest vertex
            close_corner = find_closest_corner(source,physics_coordinates)
            if close_corner == 1
                compute_taylor_integrals!(physics_interpolation1,interpolation_list1[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r,M)
            elseif close_corner == 2
                compute_taylor_integrals!(physics_interpolation2,interpolation_list2[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r,M)
            else
                compute_taylor_integrals!(physics_interpolation3,interpolation_list3[element],
                                        fOn,gOn,cOn,
                                        submatrixF,submatrixG,subvectorC,k,
                                        source,integrand,r,M)
            end
        end
        for m = 1:M+1
            F[source_node:source_node,:,m] = Fview[:,:,m]*U
            G[source_node:source_node,:,m] = Gview[:,:,m]*U
        end
        if progress; next!(prog); end # For the progress meter
    end

    From = zeros(eltype(F),n_nodes,n_nodes,M+1)
    Grom = zeros(eltype(F),n_nodes,n_nodes,M+1)

    for m = 1:M+1
        From[:,:,m] = U'*F[:,:,m]
        Grom[:,:,m] = U'*G[:,:,m]
    end

    return From, Grom, C

end

"""
    apply_taylor_expansion(Abasis,k,k0)

Assembling the Taylor series expansion with expansion wavenumber `k0` at a new wavenumber `k` and derivative matrices defined in `Abasis`.
"""
function apply_taylor_expansion(Abasis,k,k0)
    Mmax = size(Abasis,3) - 1
    Arom = Abasis[:,:,1]
    # Adding the Mmax Taylor terms
    @inbounds for i = 1:Mmax
        if i > 20 # Numerical issues of factorial otherwise
            Arom += (Float64((k - k0)^(i)/factorial(Int128(i))))*Abasis[:,:,i+1]
        else
            Arom += ((k - k0)^(i)/factorial(i))*Abasis[:,:,i+1]
        end
    end
    return Arom
end

"""
    apply_taylor_expansion(Abasis,bbasis,k,k0)

Assembling the Taylor series expansion with expansion wavenumber `k0` at a new wavenumber `k` and derivative matrices defined in `Abasis` and right-hand side defiend by `bbasis`. \\newline
Returns the system matrix and right-hand side, `Arom*x = brom`
"""
function apply_taylor_expansion(Abasis,bbasis,k,k0)
    Mmax = size(Abasis,3) - 1
    Arom = Abasis[:,:,1]
    brom = bbasis[:,1]
    # Adding the Mmax Taylor terms
    @inbounds for i = 1:Mmax
        if i > 20 # Numerical issues of factorial otherwise
            Arom += (Float64((k - k0)^(i)/factorial(Int128(i))))*Abasis[:,:,i+1]
            brom += (Float64((k - k0)^(i)/factorial(Int128(i))))*bbasis[:,i+1]
        else
            Arom .+= ((k - k0)^(i)/factorial(i))*Abasis[:,:,i+1]
            brom .+= ((k - k0)^(i)/factorial(i))*bbasis[:,i+1]
        end
    end
    return Arom,brom
end

"""
    arnoldi_basis(A,b,q)

Computes the first `q` Krylov vectors of the linear system `A*x=b` using the Arnodli method.
Returns a matrix `V` whose columsn are equal to the Krylov vectors.
"""
function arnoldi_basis(A,b,q)
    vp      = b/norm(b)                         # Extracting first Krylov vector
    V       = zeros(eltype(A),length(b),q)      # Preallocation of Krylov vector basis
    V[:,1]  = vp                                # Saving first Krylov vector
    Av      = similar(vp)                       # Preallocation of Matrix-vector product
    for p = 2:q
        Av .= A*vp
        vp .= Av
        for l = 1:p-1
            vp -= (V[:,l]'*Av)*V[:,l]
        end
        vp /= norm(vp)
        V[:,p] = vp
    end
    return V
end

"""
    scattering_krylov_basis(mesh,klist;eps=1-4,n_gauss=3,verbose=true,P₀=1,progress=true)

Computes a reduced basis for the scattering of an incident wave.
The basis is defined by computing the solution at the wavenumbers defined in `klist`.
The number of basis vectors at each wavenumber is chosen to be equal number of iterations of the `gmres` algorithm.
Returns:
 * `U`: The reduced basis matrix.
 * `solutions`: Colums equal to the solution at each wavenumbers in `klist`.
 * `qlist`: The number of Krylov vectors used at each wavenumber in `klist`.
"""
function scattering_krylov_basis(mesh,klist;eps=1-4,n_gauss=3,verbose=true,P₀=1,progress=true)
    n_sources = size(mesh.sources,2)
    nK       = length(klist)
    V = zeros(ComplexF64,n_sources, 0)      # Preallocating the total Krylov system
    solutions = zeros(ComplexF64,n_sources,nK) # list of solutions
    qlist = zeros(Int64,length(klist))
    if progress == true
        prog = Progress(nK, 0.2, "Assembling Krylov vectors:\t", 50)
    end
    for i = 0:nK-1
        k = klist[i+1]
        pI = P₀*exp.(im*k*mesh.sources[3,:]);
        Ff = FMMFOperator(mesh,k;n_gauss=n_gauss);
        Hf = Ff + 0.5I;
        p_fmm,history = gmres(Hf,pI;verbose=verbose,log=true);
        solutions[:,i+1] = p_fmm
        # F0, _, C0 = taylor_assemble!(mesh,k0,mesh.sources,mesh.shape_function;n=2,m=2,M=20,gOn=false)
        # V[:,i*q+1:(i+1)*q] = arnoldiBLI(A,b,q)
        V = [V arnoldi_basis(Hf,pI,history.iters)]
        qlist[i+1] = history.iters
        if progress == true
            next!(prog)
        end
    end
    U,S,_ = svd(V)                              # Extracting projection matrix
    idx = S .> eps  # Only use projection direcitons with singular values larger than `eps`
    return U[:,idx], solutions, qlist
end
