module feMesh

include("feDatastruct.jl")
import .feDatastruct

function buildGlobalNodeCoordinateArray(G)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    # Load Nodal Positions
    if num_dim == 1
        coordx = G["coordx"].var[:]
        NodeCoords = coordx
    elseif num_dim == 2
        coordx = G["coordx"].var[:]
        coordy = G["coordy"].var[:]
        NodeCoords = [coordx coordy]
    elseif num_dim == 3
        coordx = G["coordx"].var[:]
        coordy = G["coordy"].var[:]
        coordz = G["coordz"].var[:]
        NodeCoords = [coordx coordy coordz]
    end
    
    return NodeCoords
end

function buildGlobalElementNodalConnectivityArray(ELEMS)
    num_elems = length(ELEMS)
    # Find max number of nodes in the elements
    max_loc_nodes = 0
    for e = 1:num_elems
        num_loc_nodes = ELEMS[e].NumNodes
        if num_loc_nodes > max_loc_nodes
            max_loc_nodes = ELEMS[e].NumNodes
        end
    end

    ElemConnect = fill(-1,max_loc_nodes,num_elems)
    for e = 1:num_elems
        num_loc_nodes = ELEMS[e].NumNodes
        ElemConnect[1:num_loc_nodes,e] = ELEMS[e].ChildNodes
    end

    return ElemConnect
end

function buildNodeGlobalDOFConnectivityArray(NODES)
    num_nodes = length(NODES)
    # Find max number of dofs in the nodes
    max_num_dofs = 0
    for n = 1:num_nodes
        nodeLocalDOFS = NODES[n].ChildDOFS
        num_local_dofs = length(nodeLocalDOFS)
        if num_local_dofs > max_num_dofs
            max_num_dofs = num_local_dofs
        end
    end

    # Populate the array
    NodeConnect = fill(-1,max_num_dofs,num_nodes)
    for n = 1:num_nodes
        nodeLocalDOFS = NODES[n].ChildDOFS
        num_local_dofs = length(nodeLocalDOFS)
        for ldof = 1:num_local_dofs
            NodeConnect[ldof,n] = nodeLocalDOFS[ldof]
        end
    end
    return NodeConnect
end

function buildConstrainedDOFList(NODES, NS::Array{feDatastruct.feNodeSet,1})
    NodeConnect = buildNodeGlobalDOFConnectivityArray(NODES)
    num_gdof = max(NodeConnect...)
    isConstrainedGDOF = fill(false,num_gdof)
    num_nodesets = length(NS)
    for i = 1:num_nodesets
        constrained_ldof = NS[i].ConstrainedDOF
        num_ns_nodes = length(NS[i].ChildNodes)
        for n = 1:num_ns_nodes
            for ldof = 1:length(constrained_ldof)
                gdofID = NodeConnect[constrained_ldof,NS[i].ChildNodes[n]]
                isConstrainedGDOF[gdofID] = true
            end
        end
    end

    return isConstrainedGDOF
end

end