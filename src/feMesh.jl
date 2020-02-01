module feMesh

include("feDatastruct.jl")

export buildGlobalNodeCoordinateArray
export buildGlobalElementNodalConnectivityArray
export buildNodeGlobalDOFConnectivityArray
export buildConstrainedDOFList



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

function buildElementQuadrature(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        eDegree = ELEMS[e].Degree
        num_dims = length(eDegree)
        ELEMS[e].Quadrature = feDatastruct.Quadrature()
        ELEMS[e].Quadrature.Type = "Gauss-Legendre"
        ELEMS[e].Quadrature.Points = Array{Array{Float64,1},1}(undef,num_dims) #Array{Any,1}(undef,num_dims)
        ELEMS[e].Quadrature.Weights = Array{Array{Float64,1},1}(undef,num_dims) #Array{Any,1}(undef,num_dims)
        for dim = 1:num_dims
            nPts = Int(ceil((eDegree[dim] + 1)/2))
            ξ, W = feQuadrature.computeGaussQuadrature(nPts)
            ELEMS[e].Quadrature.Points[dim] = ξ
            ELEMS[e].Quadrature.Weights[dim] = W
        end
        
    end
    return ELEMS
end

function getDimension(elem_type::String)
    ETYPE_1D = ["BAR2"]
    ETYPE_2D = ["TRI3","QUAD4"]
    ETYPE_3D = ["TETRA4","PYRAMID5","WEDGE6","HEX8"]
    if      elem_type in ETYPE_1D
        Dimension = 1
    elseif  elem_type in ETYPE_2D
        Dimension = 2
    elseif  elem_type in ETYPE_3D
        Dimension = 3
    end
    return Dimension
end

function assignNodeDOFS(NODES, num_dim)
    dofID = 0
    for n = 1:length(NODES)
        NODES[n].ChildDOFS = zeros(Int64,num_dim)
        for d = 1:num_dim
            dofID+=1
            NODES[n].ChildDOFS[d] = dofID
        end
    end
    return NODES
end

function setParentElement(ELEMS,NODES)
    num_nodes = length(NODES)
    for n = 1:num_nodes
        NODES[n].ParentElements = Int64[]
    end

    num_elem = length(ELEMS)
    for e = 1:num_elem
        num_nodes = length(ELEMS[e].ChildNodes)
        for n = 1:num_nodes
            gnID = ELEMS[e].ChildNodes[n]
            if (e in NODES[gnID].ParentElements) == false
                append!(NODES[gnID].ParentElements,e)
            end
        end
    end
    return NODES
end

function setElementNodeTypes(ELEMS,NODES)
    num_elem = length(ELEMS)
    for e = 1:num_elem
        num_nodes = length(ELEMS[e].ChildNodes)
        ELEMS[e].BoundaryNodes = zeros(Int64,num_nodes)
        ELEMS[e].CornerNodes = zeros(Int64,num_nodes)
        ELEMS[e].FaceNodes = zeros(Int64,num_nodes)
        ELEMS[e].InternalNodes = zeros(Int64,num_nodes)
        for n = 1:num_nodes
            gnID = ELEMS[e].ChildNodes[n]
            if NODES[gnID].isElementBoundaryNode
                ELEMS[e].BoundaryNodes[n] = gnID
            end
            if NODES[gnID].isElementCornerNode
                ELEMS[e].CornerNodes[n] = gnID
            end
            if NODES[gnID].isElementFaceNode
                ELEMS[e].FaceNodes[n] = gnID
            end
            if NODES[gnID].isElementInternalNode
                ELEMS[e].InternalNodes[n] = gnID
            end
        end
    end
    
    return ELEMS
end

end