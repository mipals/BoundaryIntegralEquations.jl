"""
    computing_integrand!(integrand,Nx,wx,Ny,wy,Gxy)

Computes the double-integral ∫N(s)∫G(x,y)N(x)ᵀ ∂Γₓ∂Γₛ .
"""
function computing_integrand!(integrand,Gxy,Ny,wy,Nx,wx)
    for i = 1:length(wy)
        integrand .+= Ny[:,i]*(wx'*(Gxy[i,:] .* Nx'))*wy[i]
    end
end


function computing_galerkin_integrals!(test_physics,test_interpolation,
                                        basis_physics,basis_interpolation,
                                        submatrixF,submatrixG,subvectorC,k,integrand,r)
    # Note that everything here has been written to avoid re-allocating things. 
    # As such we re-use the memory of "integrand" for all computations
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
    # integrand_mul!(integrand,test_interpolation.jacobian_mul_weights)
    computing_integrand!(submatrixG,integrand,
                        test_physics.interpolation, test_interpolation.jacobian_mul_weights,
                        basis_physics.interpolation,basis_interpolation.jacobian_mul_weights)
    # Approximating the integral and the adding the values to the BEM matrix
    # mul!(submatrixG,integrand,Transpose(physics_function.interpolation),true,true)

    ### Evaluating the G0-kernel (used for computing the C-constant)
    # Recomputing integrand
    # freens3dk0!(integrand,r,interpolation_element.interpolation,source,interpolation_element.normals)
    # Approximating the integral and adding the value to the C-vector
    # subvectorC .+= dot(integrand, interpolation_element.jacobian_mul_weights)
end


"""
    assemble_parallel_galerkin!(mesh::Mesh3d,k,sources;
                        m=5,n=5,fkernel=F3d!,gkernel=G3d!,ckernel=G03d!,progress=true)

Assembles the BEM matrices for F, G and G0 kernels over the elements on the mesh.
"""
function assemble_parallel_galerkin!(mesh::Mesh3d,k,insources;mtest=3,ntest=3,mbasis=4,nbasis=4,progress=true)
    return assemble_parallel_galerkin!(mesh::Mesh3d,k,insources,mesh.shape_function; progress=progress,
                                    mtest=mtest,ntest=ntest,mbasis=mbasis,nbasis=nbasis)
end
function assemble_parallel_galerkin!(mesh::Mesh3d,k,insources,shape_function::SurfaceFunction;
                            mtest=3,ntest=3,mbasis=4,nbasis=4,progress=true)
    nElements   = number_of_elements(mesh)
    # nThreads    = nthreads()
    # sources     = convert.(eltype(shape_function),insources)
    nSources    = size(insources,2)
    nNodes      = size(mesh.sources,2)
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
    # Avoiding to have gauss-node on a singularity. Should be handled differently,
    # but, this is good for now.
    # if typeof(shape_function) <: QuadrilateralQuadratic9
    #     for (x,y) in zip(test_function.gauss_u,shape_function1.gauss_v)
    #         if isapprox.(x, 0.0, atol=1e-15) && isapprox.(y, 0.0, atol=1e-15)
    #             error("Gauss Node Singularity.")
    #         end
    #     end
    # end
    #======================================================================================
                        Preallocation of return values & Intermediate values
    ======================================================================================#
    # 
    F = zeros(ComplexF64, nSources, nNodes)
    G = zeros(ComplexF64, nSources, nNodes)
    C = zeros(ComplexF64, nSources)

    # Preallocation according to the number of threads
    #======================================================================================
                                    Assembly
    ======================================================================================#
    if progress; prog = Progress(nSources, 0.2, "Assembling BEM matrices: \t", 50); end
    @inbounds for test_element = 1:nElements
        # Access source
        # Every thread has access to parts of the pre-allocated matrices
        integrand = zeros(ComplexF64, ntest*mtest, nbasis*nbasis)
        r         = zeros(            ntest*mtest, nbasis*nbasis)
        test_nodes = @view physicsTopology[:,test_element]
        @inbounds for basis_element = 1:nElements
            if basis_element == test_element
                continue
            end
            basis_nodes = @view physicsTopology[:,basis_element]
            # Acces submatrix of the BEM matrix
            submatrixF = @view F[test_nodes,basis_nodes]
            submatrixG = @view G[test_nodes,basis_nodes]
            subvectorC = @view C[test_nodes]
            # Interpolating on the mesh. 
            # test_physics,test_interpolation,
            #  basis_physics,basis_interpolation,
            #  submatrixF,submatrixG,subvectorC,k,integrand,r)
            computing_galerkin_integrals!(test_physics,test_interpolation[test_element],
                                          basis_physics,basis_interpolation[basis_element],
                                          submatrixF,submatrixG,subvectorC,k,integrand,r)
            @assert all(.!isnan.(r))
            @assert all(.!isnan.(submatrixG))
            @assert all(.!isnan.(submatrixF))
        end
        println(test_element)
        if progress; next!(prog); end # For the progress meter
    end

    return F, G, C 

end
