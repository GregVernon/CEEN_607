module feMeshImport

import NCDatasets
include("feDatastruct.jl")
include("feMesh.jl")
import .feDatastruct
import .feMesh

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
    
    GEOM = MESH()
    GEOM.Elements = ELEMS
    GEOM.Nodes = NODES
    GEOM.NodeSets = NS
    GEOM.SurfaceSets = SS
    return GEOM
end

mutable struct MESH  
    Elements
    Nodes
    NodeSets
    SurfaceSets
    MESH() = new()
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
            locNodeID = R.FaceNodeOrder[locFaceID]
            append!(SS[s].ChildNodes, ELEMS[geID].ChildNodes[locNodeID])
        end
        unique!(SS[s].ChildNodes)
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

function initNodeSets(G,PARAMS)
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

    if "ns_names" in keys(G)
        # Get the names of the nodesets
        ns_names = G["ns_names"].var[:]
        for n = 1:num_node_sets
            ns_name = ns_names[:,n]
            ns_name = join(ns_name[ns_name .!= '\0'])
            NS[n].Name = ns_name
        end

        # Set the constrained DOF for the nodesets
        for n = 1:num_node_sets
            for m = 1:num_node_sets
                if lowercase(NS[n].Name) == lowercase(PARAMS.NodeSets[m].NodeSetName)
                    NS[n].ConstrainedDOF = PARAMS.NodeSets[m].DOF
                end
            end
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
    ELEMS = [feDatastruct.feElement() for e = 1:num_elem]
    e = 0
    for b = 1:num_el_blk
        blk_name = join(["connect" string(b)])
        # Determine Dimensionality of the Elements
        elem_type = G[blk_name].attrib["elem_type"]
        referenceElement = feDatastruct.makeExodusElement(elem_type)
        E2G_nodeOrder = referenceElement.ElementNodeOrder
        Dimension = feMesh.getDimension(elem_type)
        num_elem_in_block = size(G[blk_name].var[:],2)
        for blk_e = 1:num_elem_in_block
            e+=1
            # ELEMS[e] = feDatastruct.feElement()
            ELEMS[e].Dimension = Dimension
            ELEMS[e].Degree = ones(Int8, Dimension)
            ELEMS[e].ElementFamily = elem_type
            ELEMS[e].ChildNodes = G[blk_name].var[E2G_nodeOrder,blk_e]
            ELEMS[e].GlobalID = e
            ELEMS[e].ParentBlocks = b
            ELEMS[e].NumNodes = length(ELEMS[e].ChildNodes)
        end
    end
    return ELEMS
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
    NODES = feMesh.assignNodeDOFS(NODES,num_dim)

    # Get parent element ID for each node
    NODES = feMesh.setParentElement(ELEMS,NODES)

    # Backfill node boundary information for each element
    ELEMS = feMesh.setElementNodeTypes(ELEMS,NODES)

    return ELEMS,NODES
end

end