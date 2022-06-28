# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/flatironinstitute/FMM3D.git", subdir="julia"))

mutable struct HelmholtzVals
    pot
    grad

    pottarg
    gradtarg

    ier

    zk

    pgt
end

function HelmholtzVals()
    return HelmholtzVals(nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end


function h3ddir!(vals::HelmholtzVals,sources::Array{Float64},
                targets::Array{Float64};
                charges::TCN=nothing,dipvecs::TCN=nothing,
                nd::Integer=1,
                thresh::Float64=1e-16)
    # allocate memory for return values
    pottarg  = vals.pottarg
    gradtarg = vals.gradtarg
    pgt = vals.pgft
    zk = complex(vals.zk)

    # default values

    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)

    pottarg  = zero
    gradtarg = zero

    # check inputs

    if (pgt == 3)
        @warn "Hessian not implemented for Helmholtz fmm, only computing potential and gradients at targets"
        pgt = 2
    end


    pg = 0
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero end
    if (ifdipole == 0); dipvecs = zero end
    if (nt == 0); targets = 0.0 end


    if pgt == 1
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdp!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,thresh)
            else
                h3ddirectcp!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             thresh)
            end
        else
            if ifdipole == 1
                h3ddirectdp!(nd,zk,sources,
                             dipvecs,n,targets,nt,
                             pottarg,thresh)
                vals.pottarg = pottarg
            end
        end
    elseif pgt == 2
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdg!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,thresh)
            else
                h3ddirectcg!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,thresh)
            end
        else
            if ifdipole == 1
                h3ddirectdg!(nd,zk,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,thresh)
            end
        end
    end

    return vals

end
