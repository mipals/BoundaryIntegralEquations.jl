#==========================================================================================
This file contains the adaptive integration from OpenBEM with only very minor changes.
==========================================================================================#
"""
    gauss_points_triangle(n)

Returns u,v,w where (u,v) are local coordinates on a triangle with corners
                |\
                | \
                |__\
and w are the weights. The sum of weights are equal to the area of the triangle: 0.5
"""
function gauss_points_triangle(n)
    if n == 1
        IP=[1/3 1/3 0.5]
    elseif n == 3
        IP=[1/6 1/6 1/6;
            2/3 1/6 1/6;
            1/6 2/3 1/6]
    elseif n == 4
        IP=[1/3 1/3 -9/32;
            0.2 0.6 25/96;
            0.2 0.2 25/96;
            0.6 0.2 25/96]
    elseif n == 7
        IP=[1/3 1/3 0.1125;
            0.797426985353087 0.101286507323456 0.0629695902724135;
            0.101286507323456 0.797426985353087 0.0629695902724135;
            0.101286507323456 0.101286507323456 0.0629695902724135;
            0.059715871789770 0.470142064105115 0.066197076394253;
            0.470142064105115 0.059715871789770 0.066197076394253;
            0.470142064105115 0.470142064105115 0.066197076394253]
    elseif n == 12
        IP=[0.873821971016996 0.063089014491502 0.0254224531851035;
            0.063089014491502 0.873821971016996 0.0254224531851035;
            0.063089014491502 0.063089014491502 0.0254224531851035;
            0.501426509658279 0.249286745170910 0.0583931378631895;
            0.249286745170910 0.501426509658279 0.0583931378631895;
            0.249286745170910 0.249286745170910 0.0583931378631895;
            0.636502499121399 0.310352451033785 0.041425537809187;
            0.636502499121399 0.053145049844816 0.041425537809187;
            0.310352451033785 0.636502499121399 0.041425537809187;
            0.310352451033785 0.053145049844816 0.041425537809187;
            0.053145049844816 0.636502499121399 0.041425537809187;
            0.053145049844816 0.310352451033785 0.041425537809187]
    elseif n == 13
        IP=[1/3                 1/3            -0.074785022233835;
            0.479308067841923 0.260345966079038 0.087807628716602;
            0.260345966079038 0.479308067841923 0.087807628716602;
            0.260345966079038 0.260345966079038 0.087807628716602;
            0.869739794195568 0.065130102902216 0.0266736178044195;
            0.065130102902216 0.869739794195568 0.0266736178044195;
            0.065130102902216 0.065130102902216 0.0266736178044195;
            0.638444188569809 0.312865496004875 0.0385568804451285;
            0.638444188569809 0.048690315425316 0.0385568804451285;
            0.312865496004875 0.638444188569809 0.0385568804451285;
            0.312865496004875 0.048690315425316 0.0385568804451285;
            0.048690315425316 0.638444188569809 0.0385568804451285;
            0.048690315425316 0.312865496004875 0.0385568804451285]
    else
        error("Triangular quadrature not implemented for n=$(n)")
    end
    return IP[:,1],IP[:,2],IP[:,3] # x,y,weights
end
#==========================================================================================
                                Utility Functions
==========================================================================================#
"""
    maximum_side_length(element_coordinates)

Computes the maximum norm of distances between nodes in cycle.
"""
function maximum_side_length(element_coordinates)
    return maximum(sqrt.(sum(diff([element_coordinates element_coordinates[:,1]],dims=2).^2,dims=1)))
end

"""
    minimum_side_length(element_coordinates)

Computes the minimum norm of distances between nodes in cycle.
"""
function minimum_side_length(element_coordinates)
    return minimum(sqrt.(sum(diff([element_coordinates element_coordinates[:,1]],dims=2).^2,dims=1)))
end

"""
    minimum_distance(element_coordinates,source)

Computes the minimum distance from source to each node in element.
"""
function minimum_distance(element_coordinates,source)
    return minimum(sqrt.(sum((element_coordinates .- source).^2,dims=1)))
end

"""
    compute_center(element_coordinates)

Computes the center of element_coordinates.
"""
function compute_center(element_coordinates)
    return sum(element_coordinates,dims=2)/3;
end

"""
    singularity_check(p,coords,sing_check,factor=2.0)

Check distance from `p` to `coords` and compares it to element size.
If the distance if larger than a certain `factor` we need to subdivide element.
"""
function singularity_check(p,coords,sing_check,factor=2.0)
    if sing_check
        max_node_dist = max_node_distance(coords)
        max_dist      = max_distance(coords,p)
        if max_dist > max_node_dist*factor
            return false
        else
            return true
        end
    else
        return false
    end
end

"""
    max_distance(coords,p)

Computes maximum distance from `p` to the coordinates in `coords`.
"""
function max_distance(coords,p)
    max_dist = zero(eltype(p))
    @inbounds for i = 1:size(coords,2)
        tmp = hypot(coords[1,i] - p[1], coords[2,i] - p[2], coords[3,i] - p[3])
        max_dist = tmp > max_dist ? tmp : max_dist
    end
    return max_dist
end

"""
    min_distance(coords,p)

Computes minimum distance from `p` to the coordinates in `coords`.
"""
function min_distance(coords,p)
    min_dist = Inf
    @inbounds for i = 1:size(coords,2)
        tmp = hypot(coords[1,i] - p[1], coords[2,i] - p[2], coords[3,i] - p[3])
        min_dist = min_dist > tmp ? tmp : min_dist
    end
    return min_dist
end

"""
    max_node_distance(coords)

Computes maximum distance between nodes in `coords`.
"""
function max_node_distance(coords)
    max_node_dist = zero(eltype(coords))
    @inbounds for i = 2:size(coords,2)
        tmp = hypot(coords[1,i] - coords[1,i-1],
                    coords[2,i] - coords[2,i-1],
                    coords[3,i] - coords[3,i-1])
        max_node_dist = tmp > max_node_dist ? tmp : max_node_dist
    end
    return max_node_dist
end
#==========================================================================================
                            Functions for Adaptive Gaussian integration
==========================================================================================#
function subdivide_triangle(shape_function::Triangular,
                            DivsIN,
                            source,
                            element_coordinates,
                            element_size,
                            TolE)
    DivsOUT = zeros(6,0)
    for i = 1:size(DivsIN,2)
        P1 = DivsIN[1:2,i]
        P2 = DivsIN[3:4,i]
        P3 = DivsIN[5:6,i]

        P12 = (P1 + P2)/2
        P23 = (P2 + P3)/2
        P31 = (P3 + P1)/2

        center = (P1 + P2 + P3)/3.0

        subelement = element_coordinates * shape_function([center[1] P1[1] P2[1] P3[1]],
                                                          [center[2] P1[2] P2[2] P3[2]])

        distSub = minimum_distance(subelement,source)
        size_of_sub_element = maximum_side_length(subelement[:,2:end])

        if size_of_sub_element < element_size*TolE
            size_of_sub_element = 0
        end

        if size_of_sub_element > distSub/4
            NewDivs = [P1  P12 P31 P12;
                       P12 P2  P23 P23;
                       P31 P23 P3  P31]
            NewDivs = subdivide_triangle(shape_function,
                                        NewDivs,
                                        source,
                                        element_coordinates,
                                        element_size,TolE)
            DivsOUT = [DivsOUT NewDivs]
        else
            DivsOUT = [DivsOUT DivsIN[:,i]]
        end
    end
    return DivsOUT
end

function singular_triangle_integration(shape_function::SurfaceFunction,
                                     n,source,element_coordinates,TolE=1e-6)
    element_size = maximum_side_length(element_coordinates)

    aThird = 1.0/3.0
    divsIN = [0.0; 0.0; 1.0; 0.0; 0.0; 1.0] # Standard
    divsOUT = subdivide_triangle(shape_function,divsIN,source,element_coordinates,element_size,TolE)

    nodesx,nodesy,w = gauss_points_triangle(n)

    nGauss = length(nodesx)
    nSubdivisons = size(divsOUT,2)

    nodesX  = zeros(nSubdivisons * nGauss)
    nodesY  = zeros(nSubdivisons * nGauss)
    weights = zeros(nSubdivisons * nGauss)
    for i = 1:nSubdivisons
        coords  = reshape(divsOUT[:,i],2,3)
        Dim     = minimum_side_length(coords)
        center  = compute_center(coords)

        idx = ((i - 1)*nGauss + 1):(i * nGauss)
        nodesX[idx]  = (nodesx .- aThird)*Dim .+ center[1]
        nodesY[idx]  = (nodesy .- aThird)*Dim .+ center[2]
        weights[idx] =  w * Dim^2
    end

    return nodesX,nodesY,weights
end
