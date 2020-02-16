using NamedDims

function buildLocalNodeCoordinates_1D(deg)
    num_nodes = deg+1
    ξ = NamedDimsArray{(:local_node_id, :ℝᴺ,)}(zeros(Float64,num_nodes,1))
    n = 0
    for i = 1:deg+1
        n += 1
        ξ[n] = range(-1,stop=1,length=deg+1)[i]
    end
    return ξ
end

function buildLocalNodeCoordinates_2D(deg)
    num_nodes = (deg+1)^2
    ξ = NamedDimsArray{(:local_node_id, :ℝᴺ,)}(zeros(Float64,num_nodes,2))
    n = 0
    for j = 1:deg+1
        for i = 1:deg+1
            n += 1
            ξ[n,:] = [range(-1,stop=1,length=deg+1)[i], range(-1,stop=1,length=deg+1)[j]]
        end
    end
    return ξ
end

function buildLocalNodeCoordinates_3D(deg)
    num_nodes = (deg+1)^3
    ξ = NamedDimsArray{(:local_node_id, :ℝᴺ,)}(zeros(Float64,num_nodes,3))
    n = 0
    for k = 1:deg+1
        for j = 1:deg+1
            for i = 1:deg+1
                n += 1
                ξ[n,:] = [range(-1,stop=1,length=deg+1)[i], range(-1,stop=1,length=deg+1)[j], range(-1,stop=1,length=deg+1)[k]]
            end
        end
    end
    return ξ
end