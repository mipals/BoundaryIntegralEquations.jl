
import BoundaryIntegralEquations: set_interpolation_nodes!, gauss_points_triangle, number_of_elements, copy_interpolation_nodes!, interpolate_elements, computing_integrals!, connected_sources, sparse_assemble_parallel!,greens3d!,mygemm_vec!,integrand_mul!,compute_distances!, freens3d!
using Base.Threads
using ProgressMeter
# using Polyester
using StaticArrays
using LoopVectorization

function new_computing_integrals!(physics_interpolation,interpolation_element,
                                fOn,gOn,
                                submatrixF,submatrixG,k,
                                source,integrand,r)
    interpolation        = interpolation_element.interpolation
    normals              = interpolation_element.normals
    jacobian_mul_weights = interpolation_element.jacobian_mul_weights
    # Note that everything here has been written to avoid re-allocating things.
    # As such we re-use the memory of "integrand" for all computations
    compute_distances!(r,interpolation,source)
    if gOn
        ### Evaluating the G-kernel (single-layer kernel) at the global nodes
        greens3d!(integrand,r,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating the integral and the adding the values to the BEM matrix
        mygemm_vec!(submatrixG,integrand,physics_interpolation)
    end
    if fOn
        ### Evaluating the F-kernel (double-layer kernel) at the global nodes
        freens3d!(integrand,r,interpolation,source,normals,k)
        # Multiplying integrand with jacobian and weights
        integrand_mul!(integrand,jacobian_mul_weights)
        # Approximating integral and adding the value to the BEM matrix
        mygemm_vec!(submatrixF,integrand,physics_interpolation)
    end
end

function new_assemble_parallel!(mesh::Mesh3d,k,in_sources,shape_function::Triangular;
                                fOn=true,gOn=true,n_gauss=3,progress=true,depth=1,Ngauss=10,offset=0.5)
    n_elements  = number_of_elements(mesh)
    sources     = convert.(eltype(shape_function),in_sources)
    n_sources   = size(in_sources,2)
    n_nodes     = size(mesh.sources,2)
    physics_topology = mesh.physics_topology
    physics_function = deepcopy(mesh.physics_function)
    #======================================================================================

    ======================================================================================#
    shape_function1 = deepcopy(shape_function)
    set_interpolation_nodes!(shape_function1,gauss_points_triangle(n_gauss)...)
    set_interpolation_nodes!(physics_function,gauss_points_triangle(n_gauss)...)
    copy_interpolation_nodes!(physics_function,shape_function1)

    # Computing interpolation on each element
    interpolation_list = interpolate_elements(mesh,shape_function1)

    # Copying interpolation of physics functions1
    # physics_interpolation = SMatrix{size(physics_function.interpolation')...}(copy(physics_function.interpolation'))
    physics_interpolation = copy(physics_function.interpolation')

    # Preallocation of return values
    if (n_sources*n_nodes*2*8)/(1e9) > 10
        @warn "Each matrix is requires more than 10GB of storage"
    end
    F = zeros(ComplexF64, n_sources, n_nodes)
    G = zeros(ComplexF64, n_sources, n_nodes)
    C = zeros(ComplexF64, n_sources)

    element_connections, source_connections = connected_sources(mesh,depth)

    # Assembly loop
    # if progress; prog = Progress(n_sources, 0.2, "Assembling BEM matrices: \t", 50); end
    @threads for source_node = 1:n_sources
        # Every thread has access to parts of the pre-allocated matrices
        integrand  = zeros(ComplexF64,n_gauss)
        r          = zeros(n_gauss)
        # Access source
        source = sources[:,source_node]
        for element = 1:n_elements
            if element ∈ element_connections[source_node]
                continue
            else
                # Access element topology and coordinates
                physics_nodes = @view physics_topology[:,element]
                # Acces submatrix of the BEM matrix
                submatrixF = @view F[source_node,physics_nodes]
                submatrixG = @view G[source_node,physics_nodes]
                # Use quadrature point clustered around the closest vertex
                # tmp = interpolation_list[element]
                new_computing_integrals!(physics_interpolation,interpolation_list[element],
                                        fOn,gOn,
                                        submatrixF,submatrixG,k,
                                        source,integrand,r)
            end
        end
        # if progress; next!(prog); end # For the progress meter
    end

    # Singular corrections
    Fs,Gs = sparse_assemble_parallel!(mesh,k,in_sources,mesh.physics_function;
                fOn=fOn,gOn=gOn,depth=depth,progress=progress,offset=offset, Ngauss = Ngauss)

    F .+= Fs
    G .+= Gs

    return F, G, C

end
