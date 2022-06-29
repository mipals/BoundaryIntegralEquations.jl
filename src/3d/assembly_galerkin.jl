"""
    computing_integrand!(integrand,Nx,wx,Ny,wy,Gxy)

Computes the double-integral ∫N(s)∫G(x,y)N(x)ᵀ ∂Γₓ∂Γₛ .
"""
function computing_integrand!(integrand,Gxy,Ny,wy,Nx,wx)
    @inbounds for i = 1:length(wy)
        mul!(integrand,Ny[:,i],wx'*(Gxy[i,:] .* Nx'),wy[i],true)
    end
end

function computing_integrand_test!(integrand,Gxy,Ny,wy,Nx,wx,temporary1,temporary2)
    @inbounds for i = 1:length(wy)
        temporary1 .= Gxy[i,:] .* wx
        mul!(temporary2,Nx,temporary1)
        mul!(integrand,Ny[:,i],temporary2',wy[i],true)
    end
end
function computing_c_integrand!(integrand,Ny,wy)
    @inbounds for i = 1:length(wy)
        mul!(integrand,Ny[:,i],Ny[:,i]',wy[i],true)
    end
end

function computing_galerkin_integrals!(test_physics,test_interpolation,
                                        basis_physics,basis_interpolation,
                                        submatrixF,submatrixG,submatrixC,k,integrand,r,
                                        Nywy,wxNx)
    # Note that everything here has been written to avoid re-allocating things.
    compute_distances!(r,test_interpolation.interpolation,basis_interpolation.interpolation)

    ### Evaluating the F-kernel (double-layer kernel) at the global nodes
    freens3d!(integrand,r,basis_interpolation.interpolation,test_interpolation.interpolation,basis_interpolation.normals,k)
    # Computing the integrand
    computing_integrand!(submatrixF,integrand,
                        test_physics.interpolation, test_interpolation.jacobian_mul_weights,
                        basis_physics.interpolation,basis_interpolation.jacobian_mul_weights)

    ### Evaluating the G-kernel (single-layer kernel) at the global nodes
    greens3d!(integrand,r,k)
    # Computing the integrand
    computing_integrand!(submatrixG,integrand,
                        test_physics.interpolation, test_interpolation.jacobian_mul_weights,
                        basis_physics.interpolation,basis_interpolation.jacobian_mul_weights)

    ### Evaluating the G0-kernel (used for computing the C-constant)
    computing_c_integrand!(submatrixC,test_physics.interpolate_elements, test_interpolation.jacobian_mul_weights)
end


"""
    assemble_parallel_galerkin!(mesh::Mesh3d,k,sources;
                        m=5,n=5,fkernel=F3d!,gkernel=G3d!,ckernel=G03d!,progress=true)

Assembles the BEM matrices for F, G and G0 kernels over the elements on the mesh.
"""
function assemble_parallel_galerkin!(mesh::Mesh3d,k;mtest=3,ntest=3,mbasis=4,nbasis=4,progress=true)
    return assemble_parallel_galerkin!(mesh::Mesh3d,k,mesh.shape_function; progress=progress,
                                    mtest=mtest,ntest=ntest,mbasis=mbasis,nbasis=nbasis)
end
function assemble_parallel_galerkin!(mesh::Mesh3d,k,shape_function::SurfaceFunction;
                            mtest=3,ntest=3,mbasis=4,nbasis=4,progress=true)
    nElements   = number_of_elements(mesh)
    nSources    = size(mesh.sources,2)
    physicsTopology  = mesh.physics_topology
    physics_function = mesh.physics_function
    #======================================================================================

    ======================================================================================#
    test_function       = create_shape_function(shape_function;n=ntest,m=mtest)
    basis_function      = create_shape_function(shape_function;n=nbasis,m=mbasis)
    test_interpolation  = interpolate_elements(mesh,test_function)
    basis_interpolation = interpolate_elements(mesh,basis_function)
    test_physics        = deepcopy(physics_function)
    basis_physics       = deepcopy(physics_function)
    copy_interpolation_nodes!(test_physics,test_function)
    copy_interpolation_nodes!(basis_physics,basis_function)
    #======================================================================================
                        Preallocation of return values & Intermediate values
    ======================================================================================#
    F = zeros(ComplexF64, nSources, nSources)
    G = zeros(ComplexF64, nSources, nSources)
    C = zeros(ComplexF64, nSources, nSources)
    #======================================================================================
                                    Assembly
    ======================================================================================#
    n_shape_test  = number_of_shape_functions(test_physics)
    n_shape_basis = number_of_shape_functions(basis_physics)
    if progress; prog = Progress(nSources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds for test_element = 1:nElements
        # Access source
        # Every thread has access to parts of the pre-allocated matrices
        integrand = zeros(ComplexF64, ntest*mtest, nbasis*nbasis)
        r         = zeros(            ntest*mtest, nbasis*nbasis)
        Nywy      = zeros(n_shape_test)
        wxNx      = zeros(1,n_shape_basis)
        test_nodes = @view physicsTopology[:,test_element]
        @inbounds for basis_element = 1:nElements
            if basis_element == test_element
                continue
            end
            basis_nodes = @view physicsTopology[:,basis_element]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[test_nodes,basis_nodes]
            submatrixG = @view G[test_nodes,basis_nodes]
            submatrixC = @view C[test_nodes,basis_nodes]
            # Interpolating on the mesh.
            # test_physics,test_interpolation,
            #  basis_physics,basis_interpolation,
            #  submatrixF,submatrixG,subvectorC,k,integrand,r)
            computing_galerkin_integrals!(test_physics,test_interpolation[test_element],
                                          basis_physics,basis_interpolation[basis_element],
                                          submatrixF,submatrixG,submatrixC,k,integrand,r,
                                          Nywy,wxNx)
            # @assert all(.!isnan.(r))
            # @assert all(.!isnan.(submatrixG))
            # @assert all(.!isnan.(submatrixF))
        end
        # println(test_element)
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C

end
