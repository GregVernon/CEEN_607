import NCDatasets

function importGenesis(filename)
    G = NCDatasets.Dataset(filename,"r")

    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    ELEM = initElements(G)
    ELEM,NODES = initNodes(G,ELEM)
    
    NS = initNodeSets(G)
    SS = initSurfaceSets(G,ELEM)
    ES = initElementSets(G)

    GEOM = MESH()
    GEOM.Elements = ELEM
    GEOM.Nodes = NODES
    GEOM.NodeSets = NS
    GEOM.SurfaceSets = SS
    GEOM.ElementSets = ES

    NCDatasets.close(G)
    return GEOM
end

function initElementSets(G)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]

    ES = NamedDimsArray{(:global_elementset_id,)}([feElementSet() for b = 1:num_el_blk])
    num_last_elem = 0
    for b = 1:num_el_blk
        eb_name = join(["connect" string(b)])
        num_elem_in_block = size(G[eb_name].var[:],2)
        ES[b].ChildElements = collect([(num_last_elem+1) : (num_last_elem+num_elem_in_block)]...)
        num_last_elem += num_elem_in_block
        ES[b].ElementFamily = G[eb_name].attrib["elem_type"]
    end

    if "eb_names" in keys(G)
        # Get the names of the surfacesets
        eb_names = G["eb_names"].var[:]
        for n = 1:num_el_blk
            eb_name = eb_names[:,n]
            eb_name = join(eb_name[eb_name .!= '\0'])
            ES[n].Name = eb_name
        end
    end

    return ES
end

function initSurfaceSets(G,ELEM)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    SS = NamedDimsArray{(:global_surfaceset_id,)}([feSurfaceSet() for s = 1:num_side_sets])
    for s = 1:num_side_sets
        es_name = join(["elem_ss" string(s)])
        ss_name = join(["side_ss" string(s)])
        SS[s].ChildElements = NamedDimsArray{(:child_elem_id,)}(G[es_name].var[:])
        SS[s].ChildElements_LocalFace = NamedDimsArray{(:child_elem_id,)}(G[ss_name].var[:])
        SS[s].ChildNodes = Int64[]
        for e = 1:length(SS[s].ChildElements)
            geID = SS[s].ChildElements[e]
            elem_type = ELEM[geID].ElementFamily
            locFaceID = SS[s].ChildElements_LocalFace[e]
            R = makeExodusElement(elem_type)
            locNodeID = R.FaceNodeOrder[locFaceID]
            append!(SS[s].ChildNodes, ELEM[geID].ChildNodes[locNodeID])
        end
        SS[s].ChildNodes = NamedDimsArray{(:child_node_id,)}(unique(SS[s].ChildNodes))
    end

    if "ss_names" in keys(G)
        # Get the names of the surfacesets
        ss_names = G["ss_names"].var[:]
        for n = 1:num_side_sets
            ss_name = ss_names[:,n]
            ss_name = join(ss_name[ss_name .!= '\0'])
            SS[n].Name = ss_name
        end
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
     
    NS = NamedDimsArray{(:global_nodeset_id,)}([feNodeSet() for n = 1:num_node_sets])
    for n = 1:num_node_sets
        ns_name = join(["node_ns" string(n)])
        NS[n].ChildNodes = NamedDimsArray{(:child_node_id,)}(G[ns_name].var[:])
    end

    if "ns_names" in keys(G)
        # Get the names of the nodesets
        ns_names = G["ns_names"].var[:]
        for n = 1:num_node_sets
            ns_name = ns_names[:,n]
            ns_name = join(ns_name[ns_name .!= '\0'])
            NS[n].Name = ns_name
        end
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
    ELEM = NamedDimsArray{(:global_element_id,)}([feElement() for e = 1:num_elem])
    e = 0
    for b = 1:num_el_blk
        blk_name = join(["connect" string(b)])
        # Determine Dimensionality of the Elements
        elem_type = G[blk_name].attrib["elem_type"]
        referenceElement = makeExodusElement(elem_type)
        E2G_nodeOrder = referenceElement.ElementNodeOrder
        Dimension = getDimension(elem_type)
        num_variates = Dimension
        num_elem_in_block = size(G[blk_name].var[:],2)
        for blk_e = 1:num_elem_in_block
            e+=1
            # ELEM[e] = feElement()
            ELEM[e].Dimension = Dimension
            ELEM[e].NumVariates = num_variates
            ELEM[e].Degree = ones(Int8, Dimension)
            ELEM[e].ElementFamily = elem_type
            ELEM[e].ChildNodes = G[blk_name].var[E2G_nodeOrder,blk_e]
            ELEM[e].GlobalID = e
            ELEM[e].ParentBlock = b
            ELEM[e].NumNodes = length(ELEM[e].ChildNodes)
        end
    end
    return ELEM
end

function initNodes(G,ELEM)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    num_elem      = G.dim["num_elem"]
    num_el_blk    = G.dim["num_el_blk"]
    num_node_sets = G.dim["num_node_sets"]
    num_side_sets = G.dim["num_side_sets"]
    
    # Create the Nodes
    NODES = NamedDimsArray{(:global_node_id,)}([feNode() for n = 1:num_nodes])
    for n = 1:num_nodes
        # Mesh input is linear only
        NODES[n].isElementBoundaryNode = true
        NODES[n].isElementCornerNode = true
        NODES[n].isElementFaceNode = true
        NODES[n].isElementInternalNode = false
        # Load Nodal Positions
        if num_dim == 1
            coordx = G["coordx"].var[n]
            NODES[n].Coordinates = NamedDimsArray{(:ℝᴺ,)}(coordx)
        elseif num_dim == 2
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            NODES[n].Coordinates = NamedDimsArray{(:ℝᴺ,)}([coordx, coordy])
        elseif num_dim == 3
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            coordz = G["coordz"].var[n]
            NODES[n].Coordinates = NamedDimsArray{(:ℝᴺ,)}([coordx, coordy, coordz])
        end
    end

    # Assign global DOF ID to each node
    NODES = assignNodeDOFS(NODES,num_dim)

    # Get parent element ID for each node
    NODES = setParentElement(ELEM,NODES)

    # Backfill node boundary information for each element
    ELEM = setElementNodeTypes(ELEM,NODES)

    return ELEM,NODES
end