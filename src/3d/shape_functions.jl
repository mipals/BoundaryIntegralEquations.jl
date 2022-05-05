#==========================================================================================
                                QuadrilateralLinear Surface                            
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        4 -- 3                                      
                                        |    |                                      
                                        1 -- 2                                       
==========================================================================================#
quadrilateralLinearSerendipity(u,v) = [0.25*(1.0 .- u).*(1.0 .- v);
                            0.25*(1.0 .+ u).*(1.0 .- v);
                            0.25*(1.0 .+ u).*(1.0 .+ v);
                            0.25*(1.0 .- u).*(1.0 .+ v)]
function basisFunction(surfaceFunction::QuadrilateralLinearSerendipity,u,v)
    return quadrilateralLinearSerendipity(u,v)
end
#==========================================================================================
                                QuadrilateralQuadratic Surface                          
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        4  7  3                                      
                                        8     6                                      
                                        1  5  2                                      
==========================================================================================#
quadrilateralQuadratic(u,v) = [0.25*(1.0 .- u).*(v .- 1.0).*(u + v .+ 1.0);
                               0.25*(1.0 .+ u).*(v .- 1.0).*(v - u .+ 1.0); 
                               0.25*(1.0 .+ u).*(v .+ 1.0).*(u + v .- 1.0);
                               0.25*(u .- 1.0).*(v .+ 1.0).*(u - v .+ 1.0);
                               0.50*(1.0 .- v).*(1.0 .- u.^2);
                               0.50*(1.0 .+ u).*(1.0 .- v.^2);
                               0.50*(1.0 .+ v).*(1.0 .- u.^2);
                               0.50*(1.0 .- u).*(1.0 .- v.^2)]
function basisFunction(surfaceFunction::QuadrilateralQuadratic,u,v)
    return quadrilateralQuadratic(u,v)
end
#==========================================================================================
                            QuadrilateralQuadraticLagrange (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3  9  4                                      
                                        6  7  8                                      
                                        1  5  2                                      
==========================================================================================#
quadrilateralQuadraticLagrange(u,v) = [0.25*u.*(1.0 .- u).*v.*(1.0 .- v);
                               -0.25*u.*(1.0 .+ u).*v.*(1.0 .- v);
                               -0.25*u.*(1.0 .- u).*v.*(1.0 .+ v);
                                0.25*u.*(1.0 .+ u).*v.*(1.0 .+ v);
                               -0.50*(1.0 .+ u).*(1.0 .- u).*v.*(1.0 .- v);
                               -0.50*u.*(1.0 .- u).*(1.0 .+ v).*(1.0 .- v);
                               (1.0 .- u.^2).*(1.0 .- v.^2);
                                0.50*u.*(1.0 .+ u).*(1.0 .+ v).*(1.0 .- v);
                                0.50*(1.0 .+ u).*(1.0 .- u).*(1.0 .+ v).*v]
function basisFunction(surfaceFunction::QuadrilateralQuadraticLagrange,u,v)
    return quadrilateralQuadraticLagrange(u,v)
end
#==========================================================================================
                            QuadrilateralQuadratic4 (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3 --- 4                                      
                                        |     |                                      
                                        1 --- 2                                      
==========================================================================================#
quadrilateralLinear4(u,v) = [0.25*(1.0 .- u).*(1.0 .- v);
                             0.25*(1.0 .+ u).*(1.0 .- v);
                             0.25*(1.0 .- u).*(1.0 .+ v);
                             0.25*(1.0 .+ u).*(1.0 .+ v)]
function basisFunction(surfaceFunction::QuadrilateralLinear4,u,v)
    return quadrilateralLinear4(u,v)
end
#==========================================================================================
                        DiscontinuousQuadrilateralConstant
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          ---                                       
                                        |  1  |                                      
                                          ---                                       
==========================================================================================#
discontinuousQuadrilateralConstant(u,v) = ones(1,length(u))
function basisFunction(surfaceFunction::DiscontinuousQuadrilateralConstant,u,v) 
    return discontinuousQuadrilateralConstant(u,v)
end
function discontinuousQuadrilateralConstantDerivatives(u,v)
    dNu = zeros(1,length(u))
    dNv = zeros(1,length(v))
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousQuadrilateralConstant,u,v)
    return discontinuousQuadrilateralConstantDerivatives(u,v)
end
#==========================================================================================
                            DiscontinuousQuadrilateralLinear4
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          ---
                                        | 3 4 |                                      
                                        | 1 2 |                                      
                                          ---                                       
==========================================================================================#
discontinuousQuadrilateralLinear4(u,v,alpha) = quadrilateralLinear4(u./(1.0-alpha),v./(1.0-alpha))
function basisFunction(surfaceFunction::DiscontinuousQuadrilateralLinear4,u,v) 
    return discontinuousQuadrilateralLinear4(u,v,surfaceFunction.alpha)
end
#==========================================================================================
                QuadrilateralLegendre: See https://ieeexplore.ieee.org/document/1353496
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          ---
                                        | ... |                                      
                                        | ... |                                      
                                          ---                                       
==========================================================================================#
function quadrilateralLegendre(N,M,u,v)
    @assert N >= 2
    if length(u) == 1
        Pm = vcat(collectPl.(u,lmax=M)...)
        Pn = vcat(collectPl.(v,lmax=N)...)
    else
        Pm = hcat(collectPl.(u,lmax=M)...)
        Pn = hcat(collectPl.(v,lmax=N)...)
    end
    Pm_tilde = [0.5 .- 0.5*u; 0.5 .+ 0.5*u; Pm[3:end,:] - Pm[1:end-2,:]]
    return hcat([kron(Pm_tilde[:,i],Pn[:,i]) for i = 1:size(Pn,2)]...)
end
function basisFunction(surfaceFunction::QuadrilateralLegendre,ξ,η)
    return quadrilateralLegendre(surfaceFunction.N,surfaceFunction.M,ξ,η)
end
#==========================================================================================
                            QuadrilateralQuadraticLagrange (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3  9  4                                      
                                        6  7  8                                      
                                        1  5  2                                      
==========================================================================================#
discontinuousQuadrilateralQuadraticLagrange(u,v,alpha) = quadrilateralQuadraticLagrange(u./(1.0-alpha),v./(1.0-alpha))
function basisFunction(surfaceFunction::DiscontinuousQuadrilateralQuadraticLagrange,u,v)
    return discontinuousQuadrilateralQuadraticLagrange(u,v,surfaceFunction.alpha)
end
#==========================================================================================
                    DiscontinuousQuadrilateralLagrange (Kronecker Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3    4     7  8  9     
                                1   ->         ->  4  5  6     ....
                                        1    2     1  2  3 
==========================================================================================#
quadrilateralLagrange(u,v,M) = error("QuadrilateralLagrange elements not implemented yet.")
function basisFunction(surfaceFunction::QuadrilateralLagrange,u,v)
    return discontinuousQuadrilateralLagrange(u,v,surfaceFunction.M)
end
#==========================================================================================
                                   TriangularLinear Surface                             
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3                                             
                                          .  .                                          
                                          1  .  2                                        
==========================================================================================#
triangularLinear(u,v) = [1.0 .- u .- v; u; v];
function basisFunction(surfaceFunction::TriangularLinear,u,v)
    return triangularLinear(u,v)
end
#==========================================================================================
                                TriangularQuadratic Surface                             
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3                                             
                                          6  5                                          
                                          1  4  2                                       
==========================================================================================#
triangularQuadratic(u,v) = [(1.0 .- v .- u).*(1.0 .- 2.0*v .- 2.0*u);
                            u.*(2.0*u .- 1.0);
                            v.*(2.0*v .- 1.0);
                            4.0*u.*(1.0 .- v .- u);
                            4.0*u.*v;
                            4.0*v.*(1.0 .- v .- u)]
function basisFunction(surfaceFunction::TriangularQuadratic,u,v)
    return triangularQuadratic(u,v)
end
#==========================================================================================
                            DiscontinuousConstant Surface                             
 ——————————————————————————————————————  Grid  ——————————————————————————————————————————— 
                                          .
                                          . .                                            
                                          . 1 .
                                          . .  .
==========================================================================================#
discontinuousTriangularConstant(u,v) = ones(1,length(u))
function basisFunction(surfaceFunction::DiscontinuousTriangularConstant,u,v) 
    return discontinuousTriangularConstant(u,v)
end
function discontinuousTriangularConstantDerivatives(u,v)
    dNu = zeros(1,length(u))
    dNv = zeros(1,length(v))
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousTriangularConstant,u,v)
    return discontinuousTriangularConstantDerivatives(u,v)
end
#==========================================================================================
                            DiscontinuousTriangularLinear Surface                             
 ——————————————————————————————————————  Grid  ——————————————————————————————————————————— 
                                        3     
                                        | \                                            
                                        1 - 2                                          
==========================================================================================#
gamma3(u,v)  = 1.0 .- u .- v
psik(s,β)    = (s .- β)/(1.0 - 3.0*β)
psi3(u,v,β)  =  psik(gamma3(u,v),β)
psikm(s,β)   =  ones(1,length(s))/(1.0 - 3.0*β)
psi3m(u,v,β) = -psikm(gamma3(u,v),β)

discontinuousTriangularLinear(u,v,β) = [psi3(u,v,β);
                                        psik(u,β);
                                        psik(v,β)]
function basisFunction(surfaceFunction::DiscontinuousTriangularLinear,u,v)
    return discontinuousTriangularLinear(u,v,surfaceFunction.beta)
end
function discontinuousTriangularLinearDerivatives(u,v,β)
    dNu = [psi3m(u,v,β);
            psikm(u,β);
            zeros(1,length(u))]
    dNv = [psi3m(u,v,β);
            zeros(1,length(v));
            psikm(v,β)]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousTriangularLinear,u,v)
    return discontinuousTriangularLinearDerivatives(u,v,surfaceFunction.beta)
end

#==========================================================================================
                        DiscontinuousTriangularQuadratic Surface                             
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3                                             
                                          6  5                                          
                                          1  4  2                                       
==========================================================================================#
discontinuousTriangularQuadratic(u,v,β) = [psi3(u,v,β).*(2.0*psi3(u,v,β) .- 1.0);
                                            psik(u,β).*(2.0*psik(u,β) .- 1.0);
                                            psik(v,β).*(2.0*psik(v,β) .- 1.0);
                                        4.0*psik(u,β).*psi3(u,v,β);
                                        4.0*psik(u,β).*psik(v,β);
                                        4.0*psik(v,β).*psi3(u,v,β)]
function basisFunction(surfaceFunction::DiscontinuousTriangularQuadratic,u,v)
    return discontinuousTriangularQuadratic(u,v,surfaceFunction.beta)
end
function discontinuousTriangularQuadraticDerivatives(u,v,β)
    dNu = [(4.0*psi3(u,v,β) .- 1.0) .* psi3m(u,v,β);
            (4.0*psik(u,β)    .- 1.0) .* psikm(u,β);
            zeros(1,length(u));
            4.0*psikm(u,β).*psi3(u,v,β) + 4.0*psik(u,β).*psi3m(u,v,β);
            4.0*psikm(u,β).*psik(v,β);
            4.0*psik(v,β).*psi3m(u,v,β)]
    dNv = [(4.0*psi3(u,v,β) .- 1.0) .* psi3m(u,v,β);
            zeros(1,length(v));
            (4.0*psik(v,β)    .- 1.0) .* psikm(v,β);
            4.0*psik(u,β).*psi3m(u,v,β);
            4.0*psik(u,β).*psikm(v,β);
            4.0*psikm(v,β).*psi3(u,v,β) + 4.0*psik(v,β).*psi3m(u,v,β)]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousTriangularQuadratic,u,v)
    return discontinuousTriangularQuadraticDerivatives(u,v,surfaceFunction.beta)
end
#==========================================================================================
                                Triangular Constructors                                 
==========================================================================================#
function TriangularLinear(n::Real,m::Real)
    nodes_u, nodes_v, weights = triangularQuadpoints(n,m)
    TMP = TriangularLinear(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularLinear(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function TriangularLinear(SF::Triangular)
    TMP = TriangularLinear(3,3)
    weights = SF.weights
    gauss_u = SF.gauss_u
    gauss_v = SF.gauss_v
    derivatives_u,derivatives_v = TMP'(gauss_u',gauss_v')
    interp = TMP(gauss_u',gauss_v')
    return TriangularLinear(weights,gauss_u,gauss_v,derivatives_u,derivatives_v,interp)
end
function TriangularQuadratic(n::Real,m::Real)
    nodes_u, nodes_v, weights = triangularQuadpoints(n,m)
    TMP = TriangularQuadratic(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularQuadratic(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function TriangularQuadratic(SF::Triangular)
    weights = SF.weights
    gauss_u = SF.gauss_u
    gauss_v = SF.gauss_v
    TMP     = TriangularQuadratic(3,3)
    interp  = TMP(gauss_u',gauss_v')
    derivatives_u,derivatives_v = TMP'(gauss_u',gauss_v')
    return TriangularQuadratic(weights,gauss_u,gauss_v,derivatives_u,derivatives_v,interp)
end
function DiscontinuousTriangularConstant(SF::Triangular)
    weights = SF.weights
    gauss_u = SF.gauss_u
    gauss_v = SF.gauss_v
    derivatives_u = zeros(1,length(gauss_u))
    derivatives_v = zeros(1,length(gauss_v))
    interpolation = zeros(1,length(gauss_v))
    return DiscontinuousTriangularConstant(weights,gauss_u,gauss_v,derivatives_u,derivatives_v,interpolation)
end
function DiscontinuousTriangularLinear(SF::Triangular,beta)
    weights = SF.weights
    gauss_u = SF.gauss_u
    gauss_v = SF.gauss_v
    interp  = discontinuousTriangularLinear(gauss_u',gauss_v',beta)
    dX,dY   = discontinuousTriangularLinearDerivatives(gauss_u',gauss_v',beta)
    return DiscontinuousTriangularLinear(weights,gauss_u,gauss_v,dX,dY,interp,beta)
end
function DiscontinuousTriangularQuadratic(SF::Triangular,beta)
    weights = SF.weights
    gauss_u = SF.gauss_u
    gauss_v = SF.gauss_v
    interp  = discontinuousTriangularQuadratic(gauss_u',gauss_v',beta)
    dX,dY   = discontinuousTriangularQuadraticDerivatives(gauss_u',gauss_v',beta)
    return DiscontinuousTriangularQuadratic(weights,gauss_u,gauss_v,dX,dY,interp,beta)
end

#==========================================================================================
                            Quadrilateral Constructors                                 
==========================================================================================#
function QuadrilateralLinearSerendipity(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralLinearSerendipity(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralLinearSerendipity(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralLinear4(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralLinear4(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralLinear4(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralQuadratic(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralQuadratic(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralQuadratic(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralQuadraticLagrange(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralQuadraticLagrange(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralQuadraticLagrange(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function DiscontinuousQuadrilateralConstant(n::Real,m::Real,alpha_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    alpha = get_beta_quad_linear(alpha_type)
    TMP = DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,dX,dY,interpolation)
end
function DiscontinuousQuadrilateralLinear4(n::Real,m::Real,alpha_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    alpha = get_beta_quad_linear(alpha_type)
    TMP = DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),alpha)
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,dX,dY,interpolation,alpha)
end
function DiscontinuousQuadrilateralQuadraticLagrange(n::Real,m::Real,alpha_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    alpha = get_beta_quad_quadratic(alpha_type)
    TMP = DiscontinuousQuadrilateralQuadraticLagrange(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),alpha)
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralQuadraticLagrange(weights,nodes_u, nodes_v,dX,dY,interpolation,alpha)
end
function DiscontinuousQuadrilateralConstant(SF::Quadrilateral)
    nodes_u = SF.gauss_u
    nodes_v = SF.gauss_v
    weights = SF.weights
    TMP = DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation = TMP(nodes_u',nodes_v')
    dX, dY        = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,dX,dY,interpolation)
end
function DiscontinuousQuadrilateralLinear4(SF::Quadrilateral,alpha)
    nodes_u = SF.gauss_u
    nodes_v = SF.gauss_v
    weights = SF.weights
    TMP = DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),alpha)
    interpolation = TMP(nodes_u',nodes_v')
    dX, dY        = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,dX,dY,interpolation,alpha)
end
function DiscontinuousQuadrilateralQuadraticLagrange(SF::Quadrilateral,alpha)
    nodes_u = SF.gauss_u
    nodes_v = SF.gauss_v
    weights = SF.weights
    TMP = DiscontinuousQuadrilateralQuadraticLagrange(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),alpha)
    interpolation = TMP(nodes_u',nodes_v')
    dX, dY        = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralQuadraticLagrange(weights,nodes_u, nodes_v,dX,dY,interpolation,alpha)
end

