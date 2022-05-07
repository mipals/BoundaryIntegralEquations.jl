function create_pointmat(segments)
    nnodes = 2*Int(sum(segments[:,5]));
    pointmat = fill(1.0,nnodes, 3)
    counter = 1;
    for segment in eachrow(segments)
        nel  = Int(segment[5]);
        xtmp = range(segment[1],segment[3],length=2*nel+1)
        ytmp = range(segment[2],segment[4],length=2*nel+1)
        pointmat[counter:counter+2*nel-1,1] = xtmp[1:end-1];
        pointmat[counter:counter+2*nel-1,2] = ytmp[1:end-1];
        counter += 2*nel;
    end
    return pointmat
end

function create_topology(segments)
    n_elements = Int(sum(segments[:,5]));
    topology = fill(1, n_elements, 4)
    idx = 1
    for i = 1:n_elements
        for j = 1:2
            topology[i,j] = idx;
            idx += 1;
        end
        topology[i,3] = idx;
    end
    topology[end,3] = 1
    return topology
end

function create_circle(shape_function::CurveLinear,n_elements,radius=1.0, angle_output=0)
    θ = reverse(collect(range(-pi,pi,length=n_elements+1)));
    θ = θ[1:end-1]'

    if angle_output == 0
        return radius*[cos.(θ); sin.(θ)];
    elseif angle_output == 1
        return radius*[cos.(θ); sin.(θ)], θ' 
    end
end

function create_circle(shape_function::CurveQuadratic,n_elements,radius=1.0, angle_output=0)
    θ = reverse(collect(range(-pi,pi,length=2*n_elements+1)));
    θ = θ[1:end-1]'

    if angle_output == 0
        return radius*[cos.(θ); sin.(θ)];
    elseif angle_output == 1
        return radius*[cos.(θ); sin.(θ)], θ' 
    end
    
end

function create_top(shape_function::CurveLinear,n_elements)
    topology = ones(Int64, 2, n_elements)
    topology[1,:] = 1:1:n_elements
    topology[2,1:end-1] = 2:1:n_elements
    return topology
end

function create_top(shape_function::CurveQuadratic,n_elements)
    topology = zeros(Int64, 3, n_elements)
    idx = 1
    for i = 1:n_elements
        for j = 1:2
            topology[j,i] = idx;
            idx += 1;
        end
        topology[3,i] = idx;
    end
    topology[3,end] = 1
    return topology
end

function circle_nodegen(shape_function,n_elements,radius=1.0)
    coordinates = create_circle(shape_function,n_elements, radius);
    topology    = create_top(shape_function,n_elements);
    return coordinates, topology
end

function meshCircle(shape_function::CurveLinear, n_elements,radius=1.0)
    coordinates, topology = circle_nodegen(shape_function,n_elements,radius)
    return Mesh2d(coordinates,topology,CurveLinear(1))
end
function meshCircle(shape_function::CurveQuadratic,n_elements,radius=1.0)
    coordinates, topology = circle_nodegen(shape_function,n_elements,radius)
    return Mesh2d(coordinates,topology,CurveQuadratic(1))
end
