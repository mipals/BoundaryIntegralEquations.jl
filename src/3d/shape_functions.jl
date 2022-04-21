#==========================================================================================
                                QuadrilateralLinear Surface                            
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        4 -- 3                                      
                                        |    |                                      
                                        1 -- 2                                       
==========================================================================================#
quadrilateralLinear(u,v) = [0.25*(1.0 .- u).*(1.0 .- v);
                            0.25*(1.0 .+ u).*(1.0 .- v);
                            0.25*(1.0 .+ u).*(1.0 .+ v);
                            0.25*(1.0 .- u).*(1.0 .+ v)];
function basisFunction(surfaceFunction::QuadrilateralLinear,u,v)
    return quadrilateralLinear(u,v)
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
                            QuadrilateralQuadratic9 (COMSOL Layout)
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                        3  9  4                                      
                                        6  7  8                                      
                                        1  5  2                                      
==========================================================================================#
quadrilateralQuadratic9(u,v) = [0.25*u.*(1.0 .- u).*v.*(1.0 .- v);
                               -0.25*u.*(1.0 .+ u).*v.*(1.0 .- v);
                               -0.25*u.*(1.0 .- u).*v.*(1.0 .+ v);
                                0.25*u.*(1.0 .+ u).*v.*(1.0 .+ v);
                               -0.50*(1.0 .+ u).*(1.0 .- u).*v.*(1.0 .- v);
                               -0.50*u.*(1.0 .- u).*(1.0 .+ v).*(1.0 .- v);
                               (1.0 .- u.^2).*(1.0 .- v.^2);
                                0.50*u.*(1.0 .+ u).*(1.0 .+ v).*(1.0 .- v);
                                0.50*(1.0 .+ u).*(1.0 .- u).*(1.0 .+ v).*v]
function basisFunction(surfaceFunction::QuadrilateralQuadratic9,u,v)
    return quadrilateralQuadratic9(u,v)
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
                             0.25*(1.0 .+ u).*(1.0 .+ v)];
function basisFunction(surfaceFunction::QuadrilateralLinear4,u,v)
    return quadrilateralLinear4(u,v)
end

#==========================================================================================
                                   TriangularLinear Surface                             
 ——————————————————————————————————————  Grid  ———————————————————————————————————————————
                                          3                                             
                                          | \                                           
                                          1 - 2                                         
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
                                           1                                         
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
γ3(u,v) = 1.0 .- u .- v
psik(s,β) = (s .- β)/(1.0 - 3.0*β)
psi3(u,v,β) = psik(γ3(u,v),β)
psikm(s,β) =  ones(1,length(s))/(1.0 - 3.0*β)
psi3m(u,v,β) = -psikm(γ3(u,v),β)

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
                                Constructors                                 
==========================================================================================#
function TriangularLinear(n::Real,m::Real)
    nodes_u, nodes_v, weights = triangularQuadpoints(n,m)
    TMP = TriangularLinear(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularLinear(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function TriangularQuadratic(n::Real,m::Real)
    nodes_u, nodes_v, weights = triangularQuadpoints(n,m)
    TMP = TriangularQuadratic(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return TriangularQuadratic(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralLinear(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralLinear(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralLinear(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralLinear4(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralLinear4(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralLinear4(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralQuadratic(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralQuadratic(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralQuadratic(weights,nodes_u,nodes_v,dX,dY,interpolation)
end
function QuadrilateralQuadratic9(n::Real,m::Real)
    nodes_u, nodes_v, weights = quadrilateralQuadpoints(n,m)
    TMP = QuadrilateralQuadratic9(nodes_u, nodes_v, weights,rand(3,2),rand(3,2),rand(3,2))
    interpolation           = TMP(nodes_u',nodes_v')
    dX, dY                  = TMP'(nodes_u',nodes_v')
    return QuadrilateralQuadratic9(weights,nodes_u,nodes_v,dX,dY,interpolation)
end