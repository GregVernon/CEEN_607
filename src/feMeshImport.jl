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
        
        E2G_ElementFaceOrder = makeExodusElement(ELEM[1].ElementFamily).ElementFaceOrder
        SS[s].ChildElements = NamedDimsArray{(:child_elem_id,)}(G[es_name].var[:])
        SS[s].ChildElements_LocalFace = NamedDimsArray{(:child_elem_id,)}([findfirst(G[ss_name].var[i] .== collect(E2G_ElementFaceOrder)) for i = 1:length(G[ss_name].var[:])])
        SS[s].ChildNodes = Int64[]
        for e = 1:length(SS[s].ChildElements)
            geID = SS[s].ChildElements[e]
            elem_type = ELEM[geID].ElementFamily
            locFaceID = SS[s].ChildElements_LocalFace[e]
            R = makeExodusElement(elem_type)
            locNodeID = unname(R.FaceNodeOrder[locFaceID,:])
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
    node_num_map  = G["node_num_map"]

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
    node_num_map  = G["node_num_map"]
    
    # Create the Nodes
    NODES = NamedDimsArray{(:global_node_id,)}([feNode() for n = 1:num_nodes])
    for n = 1:num_nodes
        # Load Nodal Positions
        if num_dim == 1
            coordx = G["coordx"].var[n]
            NODES[node_num_map[n]].Coordinates = NamedDimsArray{(:ℝᴺ,)}(coordx)
        elseif num_dim == 2
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            NODES[node_num_map[n]].Coordinates = NamedDimsArray{(:ℝᴺ,)}([coordx, coordy])
        elseif num_dim == 3
            coordx = G["coordx"].var[n]
            coordy = G["coordy"].var[n]
            coordz = G["coordz"].var[n]
            NODES[node_num_map[n]].Coordinates = NamedDimsArray{(:ℝᴺ,)}([coordx, coordy, coordz])
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

function assignLoadConditions(GEOM, LC)
    SS = GEOM.SurfaceSets
    num_element_sets = length(GEOM.ElementSets)
    num_surface_sets = length(GEOM.SurfaceSets)
    num_load_conditions = length(LC)
    # Set the applied loads for the surface sets
    for n = 1:num_surface_sets
        if isdefined(GEOM.SurfaceSets[n], :LC_Type) == false
            GEOM.SurfaceSets[n].LC_Type = NamedDimsArray{(:lc_id,)}(Array{feEnumerate.enum_LoadCondition,1}())
            GEOM.SurfaceSets[n].LC_Direction = NamedDimsArray{(:lc_id,)}(Array{Float64,1}[])
            GEOM.SurfaceSets[n].LC_Magnitude = NamedDimsArray{(:lc_id,)}(Float64[])
        end
        for m = 1:num_load_conditions
            if LC[m].Type in [feEnumerate.pressure, feEnumerate.traction, feEnumerate.force]
                if lowercase(GEOM.SurfaceSets[n].Name) == lowercase(LC[m].SurfaceSetName)
                    push!(GEOM.SurfaceSets[n].LC_Type, LC[m].Type)
                    push!(GEOM.SurfaceSets[n].LC_Direction, LC[m].Direction)
                    push!(GEOM.SurfaceSets[n].LC_Magnitude, LC[m].Magnitude)
                end
            end
        end
    end

    # Set the applied loads for the element sets
    for n = 1:num_element_sets
        if isdefined(GEOM.ElementSets[n], :LC_Type) == false
            GEOM.ElementSets[n].LC_Type = NamedDimsArray{(:lc_id,)}(Array{feEnumerate.enum_LoadCondition,1}())
            GEOM.ElementSets[n].LC_Direction = NamedDimsArray{(:lc_id,)}(Array{Float64,1}[])
            GEOM.ElementSets[n].LC_Magnitude = NamedDimsArray{(:lc_id,)}(Float64[])
        end
        for m = 1:num_load_conditions
            if LC[m].Type in [feEnumerate.body]
                if lowercase(GEOM.ElementSets[n].Name) == lowercase(LC[m].ElementSetName)
                    push!(GEOM.ElementSets[n].LC_Type, LC[m].Type)
                    push!(GEOM.ElementSets[n].LC_Direction, LC[m].Direction)
                    push!(GEOM.ElementSets[n].LC_Magnitude, LC[m].Magnitude)
                end
            end
        end
    end

    return GEOM
end

function assignBoundaryConditions(GEOM, BC)
    num_node_sets = length(GEOM.NodeSets)
    num_boundary_conditions = length(BC)
    # Set the constrained DOF for the nodesets
    for n = 1:num_node_sets
        if isdefined(GEOM.NodeSets[n], :BC_Type) == false
            GEOM.NodeSets[n].BC_Type = NamedDimsArray{(:bc_id,)}(Array{feEnumerate.enum_BoundaryCondition,1}())
            GEOM.NodeSets[n].BC_DOF = NamedDimsArray{(:bc_id,)}(Int[])
            GEOM.NodeSets[n].BC_Value = NamedDimsArray{(:bc_id,)}(Float64[])
        end
        for m = 1:num_boundary_conditions
            if lowercase(GEOM.NodeSets[n].Name) == lowercase(BC[m].NodeSetName)
                push!(GEOM.NodeSets[n].BC_Type, BC[m].Type)
                push!(GEOM.NodeSets[n].BC_DOF, BC[m].DOF)
                push!(GEOM.NodeSets[n].BC_Value, BC[m].Value)
            end
        end
    end
    return GEOM
end