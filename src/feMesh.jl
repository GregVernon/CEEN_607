using NamedDims

function buildGlobalNodeCoordinateArray(NODES)
    # Get global information about the Genesis file
    num_nodes = size(NODES,:global_node_id)
    num_dim = size(NODES[1].Coordinates,:ℝᴺ)
    
    NodeCoords = NamedDimsArray{(:global_node_id, :ℝᴺ,)}(zeros(Float64,num_nodes,num_dim))
    for n = 1:num_nodes
            NodeCoords[n,:] = NODES[n].Coordinates
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

    ElemConnect = NamedDimsArray{(:local_node_id,:global_elem_id,)}(fill(-1,max_loc_nodes,num_elems))
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
    NodeConnect = NamedDimsArray{(:local_dof_id,:global_node_id,)}(fill(-1,max_num_dofs,num_nodes))
    for n = 1:num_nodes
        nodeLocalDOFS = NODES[n].ChildDOFS
        num_local_dofs = length(nodeLocalDOFS)
        for ldof = 1:num_local_dofs
            NodeConnect[ldof,n] = nodeLocalDOFS[ldof]
        end
    end
    return NodeConnect
end

function buildConstrainedDOFList(NODES, NS)
    NodeConnect = buildNodeGlobalDOFConnectivityArray(NODES)
    num_gdof = max(NodeConnect...)
    isConstrainedGDOF = NamedDimsArray{(:global_dof_id,)}(fill(false,num_gdof))
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

function buildDOFConstraintTypeArray(NODES,NODESETS)
    num_nodesets = length(NODESETS)
    num_nodes = length(NODES)
    num_loc_dof = length(NODES[1].ChildDOFS)

    DOFConstraintType = NamedDimsArray{(:global_node_id,:local_dof_id,)}([feEnumerate.unconstrained for i = 1:num_loc_dof, j=1:num_nodes])
    for ns = 1:num_nodesets
        if isempty(NODESETS[ns].BC_Type) == false
            for n = 1:length(NODESETS[ns].ChildNodes)
                for i = 1:length(NODESETS[ns].BC_DOF)
                    DOFConstraintType[NODESETS[ns].BC_DOF[i], NODESETS[ns].ChildNodes[n]] = NODESETS[ns].BC_Type[i]
                end
            end
        end
    end
    return DOFConstraintType
end

function buildDOFConstraintValueArray(NODES,NODESETS)
    num_nodesets = length(NODESETS)
    num_nodes = length(NODES)
    num_loc_dof = length(NODES[1].ChildDOFS)
    DOFConstraintValue = NamedDimsArray{(:global_node_id,:local_dof_id,)}(Array{Any,2}(fill(nothing,num_loc_dof,num_nodes)))
    for ns = 1:num_nodesets
        if isempty(NODESETS[ns].BC_Type) == false
            for n = 1:length(NODESETS[ns].ChildNodes)
                for i = 1:length(NODESETS[ns].BC_DOF)
                    DOFConstraintValue[NODESETS[ns].BC_DOF[i], NODESETS[ns].ChildNodes[n]] = NODESETS[ns].BC_Value[i]
                end
            end
        end
    end
    return DOFConstraintValue
end

function buildLoadTypeArray(ELEMS,SURFACESETS)
    num_surfacesets = length(SURFACESETS)
    num_elems = length(ELEMS)
    num_loc_sides = length(ELEMS[1].SideNodes)
    LoadTypeArray = NamedDimsArray{(:local_side_id, :global_elem_id,)}(Array{Any,2}(fill(nothing,num_loc_sides,num_elems)))
    for ss = 1:num_surfacesets
        if isempty(SURFACESETS[ss].LC_Type) == false
            for e = 1:length(SURFACESETS[ss].ChildElements)
                if SURFACESETS[ss].LC_Type == feEnumerate.body
                    LoadTypeArray[1,SURFACESETS[ss].ChildElements[e]] = feEnumerate.body
                else
                    for i = 1:length(SURFACESETS[ss].LC_Type)
                        LoadTypeArray[SURFACESETS[ss].ChildElements_LocalFace[i] + 1, SURFACESETS[ss].ChildElements[e]] = SURFACESETS[ss].LC_Type[i]
                    end
                end
            end
        end
    end
    return LoadTypeArray
end

function buildElementQuadrature(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        eDegree = ELEMS[e].Degree
        num_dims = ELEMS[e].Dimension
        # Preallocate Quadrature Array for Element Sides
        num_loc_sides = size(ELEMS[1].SideNodes, :local_side_id)
        for i = 1:num_loc_sides
            if i == 1
                ELEMS[e].Quadrature = NamedDimsArray{(:local_side_id,)}(Array{feQuadrature,1}())
                push!(ELEMS[e].Quadrature, feQuadrature())
            else
                push!(ELEMS[e].Quadrature, feQuadrature())
            end
        end

        # Preallocate Quadrature Points for each Element Side
        for side_id = 1:num_loc_sides
            if side_id == 1
                nPts = (Int(ceil((max(eDegree...))/2)) + 1)^2
                for qp_id = 1:nPts
                    if qp_id == 1
                        ELEMS[e].Quadrature[side_id].QuadraturePoints = NamedDimsArray{(:local_qp_id,)}(Array{feQuadraturePoint,1}())
                        push!(ELEMS[e].Quadrature[side_id].QuadraturePoints, feQuadraturePoint())
                    else
                        push!(ELEMS[e].Quadrature[side_id].QuadraturePoints, feQuadraturePoint())
                    end
                end
            else
                nPts = Int(ceil((max(eDegree...))/2)) + 1
                for qp_id = 1:nPts
                    if qp_id == 1
                        ELEMS[e].Quadrature[side_id].QuadraturePoints = NamedDimsArray{(:local_qp_id,)}(Array{feQuadraturePoint,1}())
                        push!(ELEMS[e].Quadrature[side_id].QuadraturePoints, feQuadraturePoint())
                    else
                        push!(ELEMS[e].Quadrature[side_id].QuadraturePoints, feQuadraturePoint())
                    end
                end
            end
        end

        # Define Quadrature on the element domain
        ELEMS[e].Quadrature[1].Type = "Gauss-Legendre"
        nPts = Int(ceil((max(eDegree...))/2)) + 1
        if num_dims == 2
            side_id = 1
            ξ,W = GaussQuadratureRule_2D(nPts)
            xₐ = buildLocalNodeCoordinates_2D(eDegree[1])
            Nₐ = ξ->LagrangeBasis_2D(eDegree[1],ξ)
            ∇Nₐ = ξ->∇LagrangeBasis_2D(eDegree[1],ξ)
            for qp = 1:size(ξ, :local_qp_id)
                Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
                
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Coordinates = ξ[qp,:]
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Weights = W[qp]
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Nₐ = Nₐ(ξ[qp,:])
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].∇Nₐ = ∇Nₐ(ξ[qp,:])
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Jᵢⱼ = Jᵢⱼ(ξ[qp,:])
                ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ[qp,:]), Jᵢⱼ(ξ[qp,:]))
            end
            
            
            for side_id = 2:num_loc_sides
                ξ,W = GaussQuadratureRule_1D(nPts)
                xₐ = buildLocalNodeCoordinates_1D(eDegree[1])
                Nₐ = ξ->LagrangeBasis_1D(eDegree[1],ξ)
                ∇Nₐ = ξ->∇LagrangeBasis_1D(eDegree[1],ξ)
                for qp = 1:size(ξ, :local_qp_id)
                    Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)

                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Coordinates = ξ[qp]
                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Weights = W[qp]
                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Nₐ = Nₐ([ξ[qp]])
                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].∇Nₐ = ∇Nₐ([ξ[qp]])
                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].Jᵢⱼ = Jᵢⱼ([ξ[qp]])
                    ELEMS[e].Quadrature[side_id].QuadraturePoints[qp].∇ₓNₐ = compute∇ₓNₐ(∇Nₐ([ξ[qp]]), Jᵢⱼ([ξ[qp]]))
                end
            end

        elseif num_dims == 3
            ξ,W = GaussQuadratureRule_3D(nPts)
            ELEMS[e].Quadrature[1].Points = ξ
            ELEMS[e].Quadrature[1].Weights = W
            for side_id = 2:num_loc_sides+1
                ξ, W = GaussQuadratureRule_2D(nPts)
                ELEMS[e].Quadrature[side_id].Points = ξ
                ELEMS[e].Quadrature[side_id].Weights = W
            end
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
        NODES[n].ChildDOFS = NamedDimsArray{(:local_dof_id,)}(zeros(Int64,num_dim))
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
        ELEMS[e].BoundaryNodes = NamedDimsArray{(:local_node_id,)}(zeros(Int64,num_nodes))
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

function makeExodusElement(elem_type)
    ElementType = elem_type
    if elem_type == "BAR2"
        ElementNodeOrder = [1 2]
        ElementFaceOrder = [1]
        FaceNodeOrder = [1 2]
        isBoundaryNode = fill(true,2)
        isCornerNode = fill(true,2)
        isFaceNode = fill(true,2)
        isInternalNode = fill(false,2)
    end

    if elem_type == "QUAD4"
        ElementNodeOrder = [1, 2, 4, 3]
        ElementFaceOrder = [1]
        FaceNodeOrder = [3 1; 2 4; 1 2; 4 3]
        isBoundaryNode = fill(true,4)
        isCornerNode = fill(true,4)
        isFaceNode = fill(true,4)
        isInternalNode = fill(false,4)
    end

    if elem_type == "HEX8"
        ElementNodeOrder = [1, 2, 4, 3, 5, 6, 8, 7]
        ElementFaceOrder = [4, 2, 1, 3, 5, 6]
        FaceNodeOrder = [1 2 6 5; 2 4 8 6; 4 3 7 8; 1 5 7 3; 1 3 4 2; 5 6 8 7]
        isBoundaryNode = fill(true,8)
        isCornerNode = fill(true,8)
        isFaceNode = fill(true,8)
        isInternalNode = fill(false,8)
    end
    ElementFaceOrder = NamedDimsArray{(:local_face_id,)}(ElementFaceOrder)
    ElementNodeOrder = NamedDimsArray{(:local_node_id,)}(ElementNodeOrder)
    FaceNodeOrder = NamedDimsArray{(:local_face_id,:local_node_id,)}(FaceNodeOrder)
    isBoundaryNode = NamedDimsArray{(:local_node_id,)}(isBoundaryNode)
    isCornerNode = NamedDimsArray{(:local_node_id,)}(isCornerNode)
    isFaceNode = NamedDimsArray{(:local_node_id,)}(isFaceNode)
    isInternalNode = NamedDimsArray{(:local_node_id,)}(isInternalNode)
    GE = ExodusElement(ElementType,ElementFaceOrder,ElementNodeOrder,FaceNodeOrder,isBoundaryNode,isCornerNode,isFaceNode,isInternalNode)
    return GE
end