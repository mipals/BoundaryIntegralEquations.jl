function compute_integrands!(Fslice,Gslice,Cslice,interpolation,physics_interpolation!,
                            source,k,r,normals,jacobian,integrand,gOn,fOn,cOn)
    compute_distances(interpolation,source,r)                    # Computing distances
    if gOn # Only compute G if you need it
        G!(k,r,integrand)                                        # Evaluate Greens function
        jacobian_mul_weights!(integrand,jacobian)                # Multiply by jacobian*weights
        mul!(Gslice,integrand,physics_interpolation!',true,true) # Integration with basis func
    end
    if fOn # Only compute F if you need it
        F!(interpolation,source,k,normals,r,integrand)           # Evaluate ∂ₙGreens function
        jacobian_mul_weights!(integrand,jacobian)                # Multiply by jacobian*weights
        mul!(Fslice,integrand,physics_interpolation!',true,true) # Integration with basis func
    end
    if cOn # Only compute c if you need it
        C!(interpolation,source,normals,r,integrand)             # Evaluate G₀
        Cslice[1] += dot(integrand,jacobian)                     # Integration of G₀ 
    end
end


function assemble_parallel!(mesh::Mesh2d,k,insources;n=4,progress=true)
    return assemble_parallel!(mesh::Mesh2d,k,insources,mesh.shape_function;
                                                    n=n,progress=progress)
end

function assemble_parallel!(mesh::Mesh2d,k,insources,shape_function::CurveFunction;
                                                    n=4,progress=true)
error("2d assembly not implemented yet")
end
