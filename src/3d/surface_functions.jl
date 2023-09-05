#==========================================================================================
                                QuadrilateralLinear Surface
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        4 -- 3
                                        |    |
                                        1 -- 2
==========================================================================================#
quadrilateralLinearSerendipity(u,v) =  [(1 .- u).*(1 .- v)/4;
                                        (1 .+ u).*(1 .- v)/4;
                                        (1 .+ u).*(1 .+ v)/4;
                                        (1 .- u).*(1 .+ v)/4]
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
quadrilateralQuadratic(u,v) = [(1 .- u).*(v .- 1).*(u + v .+ 1)/4;
                               (1 .+ u).*(v .- 1).*(v - u .+ 1)/4;
                               (1 .+ u).*(v .+ 1).*(u + v .- 1)/4;
                               (u .- 1).*(v .+ 1).*(u - v .+ 1)/4;
                               (1 .- v).*(1 .- u.^2)/2;
                               (1 .+ u).*(1 .- v.^2)/2;
                               (1 .+ v).*(1 .- u.^2)/2;
                               (1 .- u).*(1 .- v.^2)/2]
function basisFunction(surfaceFunction::QuadrilateralQuadratic,u,v)
    return quadrilateralQuadratic(u,v)
end
function quadrilateralQuadraticDerivative(u,v)
    dNu = [ -(v .- 1).*(u + v .+ 1)/4 .+ (1 .- u).*(v .- 1)/4;
             (v .- 1).*(v - u .+ 1)/4 .- (1 .+ u).*(v .- 1)/4;
             (v .+ 1).*(u + v .- 1)/4 .+ (1 .+ u).*(v .+ 1)/4;
             (v .+ 1).*(u - v .+ 1)/4 .- (1 .- u).*(v .+ 1)/4;
            -(1 .- v).* u;
             (1 .- v.^2)/2;
            -(1 .+ v).*u;
            -(1 .- v.^2)/2]
    dNv = [ (1 .- u).*(u + v .+ 1)/4 + (1 .- u).*(v .- 1)/4;
            (1 .+ u).*(v - u .+ 1)/4 + (1 .+ u).*(v .- 1)/4;
            (1 .+ u).*(u + v .- 1)/4 + (1 .+ u).*(v .+ 1)/4;
            (u .- 1).*(u - v .+ 1)/4 + (1 .- u).*(v .+ 1)/4;
            (u.^2 .- 1)/2;
            -(1 .+ u).*v;
             (1 .- u.^2)/2;
            -(1 .- u).*v]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::QuadrilateralQuadratic,u,v)
    return quadrilateralQuadraticDerivative(u,v)
end
#==========================================================================================
                            QuadrilateralQuadraticLagrange (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3  9  4
                                        6  7  8
                                        1  5  2
==========================================================================================#
quadrilateralQuadraticLagrange(u,v) = [u.*(1 .- u).*v.*(1 .- v)/4;
                                      -u.*(1 .+ u).*v.*(1 .- v)/4;
                                      -u.*(1 .- u).*v.*(1 .+ v)/4;
                                       u.*(1 .+ u).*v.*(1 .+ v)/4;
                                      -(1 .+ u).*(1 .- u).*v.*(1 .- v)/2;
                                      -u.*(1 .- u).*(1 .+ v).*(1 .- v)/2;
                                        (1 .- u.^2).*(1 .- v.^2);
                                       u.*(1 .+ u).*(1 .+ v).*(1 .- v)/2;
                                       (1 .+ u).*(1 .- u).*(1 .+ v).*v/2]
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
quadrilateralLinear4(u,v) = [(1 .- u).*(1 .- v)/4;
                             (1 .+ u).*(1 .- v)/4;
                             (1 .- u).*(1 .+ v)/4;
                             (1 .+ u).*(1 .+ v)/4]
function basisFunction(surfaceFunction::QuadrilateralLinear4,u,v)
    return quadrilateralLinear4(u,v)
end
function quadrilateralLinear4Derivatives(u,v)
    dNu = [-(1 .- v)/4;  (1 .- v)/4; -(1 .+ v)/4; (1 .+ v)/4]
    dNv = [-(1 .- u)/4; -(1 .+ u)/4;  (1 .- u)/4; (1 .+ u)/4]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::QuadrilateralLinear4,u,v)
    return quadrilateralLinear4Derivatives(u,v)
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
    dNu = zeros(eltype(u),1,length(u))
    dNv = zeros(eltype(v),1,length(v))
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
discontinuousQuadrilateralLinear4(u,v,beta) = quadrilateralLinear4(u./(1-beta),v./(1-beta))
function basisFunction(surfaceFunction::DiscontinuousQuadrilateralLinear4,u,v)
    return discontinuousQuadrilateralLinear4(u,v,surfaceFunction.beta)
end
function discontinuousQuadrilateralLinear4Derivatives(u,v,beta)
    return quadrilateralLinear4Derivatives(u./(1 - beta),v./(1 - beta))./(1 - beta)
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousQuadrilateralLinear4,u,v)
    return discontinuousQuadrilateralLinear4Derivatives(u,v,surfaceFunction.beta)
end
#==========================================================================================
                DiscontinuousQuadrilateralQuadraticLagrange (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3  9  4
                                        6  7  8
                                        1  5  2
==========================================================================================#
discontinuousQuadrilateralQuadraticLagrange(u,v,beta) = quadrilateralQuadraticLagrange(u./(1-beta),v./(1-beta))
function basisFunction(surfaceFunction::DiscontinuousQuadrilateralQuadraticLagrange,u,v)
    return discontinuousQuadrilateralQuadraticLagrange(u,v,surfaceFunction.beta)
end
# function discontinuousQuadrilateralQuadraticLagrangeDerivative(u,v,beta)
#     return quadrilateralQuadraticLagrangeDerivative(u./(1.0-beta),v./(1.0-beta))./(1.0-beta)
# end
# function basisFunctionDerivative(surfaceFunction::DiscontinuousQuadrilateralQuadraticLagrange,u,v)
#     return discontinuousQuadrilateralQuadraticLagrangeDerivative(u,v,surfaceFunction.beta)
# end
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
function basisFunction(surfaceFunction::QuadrilateralLegendre,u,v)
    return quadrilateralLegendre(surfaceFunction.N,surfaceFunction.M,u,v)
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
triangularLinear(u,v) = [1 .- u .- v; u; v];
function basisFunction(surfaceFunction::TriangularLinear,u,v)
    return triangularLinear(u,v)
end
triangularLinearDerivatives(u,v) = [-1;1;0].*ones(eltype(u),1,length(u)),
                                   [-1;0;1].*ones(eltype(v),1,length(v));
function basisFunctionDerivative(surfaceFunction::TriangularLinear,u,v)
    return triangularLinearDerivatives(u,v)
end
#==========================================================================================
                                TriangularQuadratic Surface
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3
                                          6  5
                                          1  4  2
==========================================================================================#
triangularQuadratic(u,v) = [(1 .- v .- u).*(1 .- 2*v .- 2*u);
                            u.*(2*u .- 1);
                            v.*(2*v .- 1);
                            4*u.*(1 .- v .- u);
                            4*u.*v;
                            4*v.*(1 .- v .- u)]
function basisFunction(surfaceFunction::TriangularQuadratic,u,v)
    return triangularQuadratic(u,v)
end
function triangularQuadraticDerivatives(u,v)
    dNu = [-(2*(1 .- v .- u) .- 1) - 2*(1 .- v .- u);
            4*u .- 1;
            zeros(eltype(u),1,length(u));
            4*(1 .- v .- u) - 4*u;
            4*v;
           -4*v]
    dNv = [-(2*(1 .- v .- u) .- 1) - 2*(1 .- v .- u);
            zeros(eltype(u),1,length(u));
            4*v .- 1;
           -4*u;
            4*u;
            4*(1 .- v .- u) - 4*v]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::TriangularQuadratic,u,v)
    return triangularQuadraticDerivatives(u,v)
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
    dNu = zeros(eltype(u),1,length(u))
    dNv = zeros(eltype(v),1,length(v))
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
gamma3(u,v)  = 1 .- u .- v
psik(s,β)    = (s .- β)/(1 - 3*β)
psi3(u,v,β)  =  psik(gamma3(u,v),β)
psikm(s,β)   =  ones(1,length(s))/(1 - 3*β)
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
discontinuousTriangularQuadratic(u,v,β) = [psi3(u,v,β).*(2*psi3(u,v,β) .- 1);
                                            psik(u,β).*(2*psik(u,β) .- 1);
                                            psik(v,β).*(2*psik(v,β) .- 1);
                                          4*psik(u,β).*psi3(u,v,β);
                                          4*psik(u,β).*psik(v,β);
                                          4*psik(v,β).*psi3(u,v,β)]
function basisFunction(surfaceFunction::DiscontinuousTriangularQuadratic,u,v)
    return discontinuousTriangularQuadratic(u,v,surfaceFunction.beta)
end
function discontinuousTriangularQuadraticDerivatives(u,v,β)
    dNu = [ (4*psi3(u,v,β) .- 1) .* psi3m(u,v,β);
            (4*psik(u,β)   .- 1) .* psikm(u,β);
            zeros(eltype(u),1,length(u));
            4*psikm(u,β).*psi3(u,v,β) + 4*psik(u,β).*psi3m(u,v,β);
            4*psikm(u,β).*psik(v,β);
            4*psik(v,β).*psi3m(u,v,β)]
    dNv = [(4*psi3(u,v,β) .- 1) .* psi3m(u,v,β);
            zeros(eltype(v),1,length(v));
           (4*psik(v,β)   .- 1) .* psikm(v,β);
            4*psik(u,β).*psi3m(u,v,β);
            4*psik(u,β).*psikm(v,β);
            4*psikm(v,β).*psi3(u,v,β) + 4*psik(v,β).*psi3m(u,v,β)]
    return dNu, dNv
end
function basisFunctionDerivative(surfaceFunction::DiscontinuousTriangularQuadratic,u,v)
    return discontinuousTriangularQuadraticDerivatives(u,v,surfaceFunction.beta)
end
#==========================================================================================
                                Triangular Constructors
==========================================================================================#
function TriangularLinear(n::Real)
    nodes_u, nodes_v, weights = gauss_points_triangle(n)
    TMP = TriangularLinear(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularLinear(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
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
function TriangularQuadratic(n::Real)
    nodes_u, nodes_v, weights = gauss_points_triangle(n)
    TMP = TriangularQuadratic(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularQuadratic(weights,nodes_u,nodes_v,dX,dY,interpolation)
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
function DiscontinuousQuadrilateralConstant(n::Real,m::Real,beta_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    beta = get_beta_quad_linear(beta_type)
    TMP = DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralConstant(weights, nodes_u, nodes_v,dX,dY,interpolation)
end
function DiscontinuousQuadrilateralLinear4(n::Real,m::Real,beta_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    beta = get_beta_quad_linear(beta_type)
    TMP = DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),beta)
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,dX,dY,interpolation,beta)
end
function DiscontinuousQuadrilateralQuadraticLagrange(n::Real,m::Real,beta_type=:legendre)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    beta = get_beta_quad_quadratic(beta_type)
    TMP = DiscontinuousQuadrilateralQuadraticLagrange(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),beta)
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralQuadraticLagrange(weights,nodes_u, nodes_v,dX,dY,interpolation,beta)
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
function DiscontinuousQuadrilateralLinear4(SF::Quadrilateral,beta)
    nodes_u = SF.gauss_u
    nodes_v = SF.gauss_v
    weights = SF.weights
    TMP = DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),beta)
    interpolation = TMP(nodes_u',nodes_v')
    dX, dY        = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralLinear4(weights,nodes_u, nodes_v,dX,dY,interpolation,beta)
end
function DiscontinuousQuadrilateralQuadraticLagrange(SF::Quadrilateral,beta)
    nodes_u = SF.gauss_u
    nodes_v = SF.gauss_v
    weights = SF.weights
    TMP = DiscontinuousQuadrilateralQuadraticLagrange(weights, nodes_u, nodes_v,rand(3,2),rand(3,2),rand(3,2),beta)
    interpolation = TMP(nodes_u',nodes_v')
    dX, dY        = TMP'(nodes_u',nodes_v')
    return DiscontinuousQuadrilateralQuadraticLagrange(weights,nodes_u, nodes_v,dX,dY,interpolation,beta)
end
