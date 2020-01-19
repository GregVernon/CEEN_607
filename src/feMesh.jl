module feMesh

import NCDatasets
include("feDatastruct.jl")
import .feDatastruct

export importGenomat

function importGenomat(filename)
    G = NCDatasets.Dataset(filename,"r")

    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    ELEMS = initElements(G)
    ELEMS,NODES = initNodes(G,ELEMS)
    
    NS = initNodeSets(G)
    SS = initSurfaceSets(G,ELEMS)
    
    return ELEMS, NODES,  NS, SS
end

mutable struct Genomat
    num_dim
    num_nodes
    num_elem
    num_el_blk
    num_node_sets
    num_side_sets
end

function initSurfaceSets(G,ELEMS)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    SS = [feDatastruct.feSurfaceSet() for s = 1:num_side_sets]
    for s = 1:num_side_sets
        es_name = join(["elem_ss" string(s)])
        ss_name = join(["side_ss" string(s)])
        SS[s].ChildElements = G[es_name].var[:]
        SS[s].ChildElements_LocalFace = G[ss_name].var[:]
        SS[s].ChildNodes = Int64[]
        for e = 1:length(SS[s].ChildElements)
            geID = SS[s].ChildElements[e]
            elem_type = ELEMS[geID].ElementFamily
            locFaceID = SS[s].ChildElements_LocalFace[e]
            R = feDatastruct.makeExodusElement(elem_type)
            locNodeID = R.FaceNodeOrder[R.ElementFaceOrder[locFaceID]] # FIXME
            append!(SS[s].ChildNodes, ELEMS[geID].ChildNodes[locNodeID])
        end
        unique!(SS[s].ChildNodes)
    end

    return SS
end

function initNodeSets(G)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
     
    NS = [feDatastruct.feNodeSet() for n = 1:num_node_sets]
    for n = 1:num_node_sets
        ns_name = join(["node_ns" string(n)])
        NS[n].ChildNodes = G[ns_name].var[:]
    end

    return NS
end

function initElements(G)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]

    # Create the Elements
    ELEMS = [feDatastruct.feElement() for e = 1:num_elem]
    e = 0
    for b = 1:num_el_blk
        blk_name = join(["connect" string(b)])
        # Determine Dimensionality of the Elements
        elem_type = G[blk_name].attrib["elem_type"]
        referenceElement = feDatastruct.makeExodusElement(elem_type)
        E2G_nodeOrder = referenceElement.ElementNodeOrder
        Dimension = getDimension(elem_type)
        for blk_e = 1:num_elem
            e+=1
            # ELEMS[e] = feDatastruct.feElement()
            ELEMS[e].Dimension = Dimension
            ELEMS[e].Degree = ones(Int8, Dimension)
            ELEMS[e].ElementFamily = elem_type
            ELEMS[e].ChildNodes = G[blk_name].var[E2G_nodeOrder,e]
            ELEMS[e].GlobalID = e
            ELEMS[e].ParentBlocks = b
            ELEMS[e].NumNodes = length(ELEMS[e].ChildNodes)
        end
    end
    return ELEMS
end

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

function initNodes(G,ELEMS)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    # Create the Nodes
    NODES = [feDatastruct.feNode() for n = 1:num_nodes]
    for n = 1:num_nodes
        # A Genomat mesh is linear only
        NODES[n].isElementBoundaryNode = true
        NODES[n].isElementCornerNode = true
        NODES[n].isElementFaceNode = true
        NODES[n].isElementInternalNode = false
        # Load Nodal Positions
        if num_dim == 1
            coordx = G["coordx"].var[n]
            NODES[n].Coordinates = coordx
        elseif num_dim == 2
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            NODES[n].Coordinates = [coordx coordy]
        elseif num_dim == 3
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            coordz = G["coordz"].var[n]
            NODES[n].Coordinates = [coordx coordy coordz]
        end
    end

    # Assign global DOF ID to each node
    NODES = assignNodeDOFS(NODES,num_dim)

    # Get parent element ID for each node
    NODES = setParentElement(ELEMS,NODES)

    # Backfill node boundary information for each element
    ELEMS = setElementNodeTypes(ELEMS,NODES)

    return ELEMS,NODES
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