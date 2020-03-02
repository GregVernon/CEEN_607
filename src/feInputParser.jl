function importInput(filename)
    fLines = readlines(filename)
    fLines = strip.(fLines)
    fLines = lowercase.(fLines)

    GC, BC, LC = findInputCards(fLines)
    GEOM = importGenesis(GC.Filename)
    GEOM = assignLoadConditions(GEOM, LC)
    GEOM = assignBoundaryConditions(GEOM, BC)

    PARAMS = InputParams()
    PARAMS.BoundaryConditions = BC
    PARAMS.LoadConditions = LC

    return GEOM, PARAMS
end

function findInputCards(fLines)
    # First get the geometry
    geomLine = findall(occursin.("geometry",fLines))[1]
    geomCard = fLines[geomLine]
    GC = parseGeometryCard(geomCard)
    # G = importGenesis(geomFile)

    # Next get parameters
    beginLines = findall(occursin.("begin",fLines))
    endLines = findall(occursin.("end",fLines))

    # Find Boundary Condition Cards
    bcStartLines = findall(occursin.("begin boundary condition",fLines))
    num_bc_cards = length(bcStartLines)
    bcEndLines = zeros(Int,num_bc_cards)
    for bc = 1:num_bc_cards
        bc_endline_index = findall(endLines .> bcStartLines[bc])[1]
        bcEndLines[bc] = endLines[bc_endline_index]
    end
    
    BC = [BoundaryCondition() for bc = 1:num_bc_cards]
    for bc = 1:num_bc_cards
        bcCard = fLines[(bcStartLines[bc]+1) : (bcEndLines[bc]-1)]
        BC[bc] = parseBoundaryCondition(bcCard)
    end

    # Find Load Condition Cards
    lcStartLines = findall(occursin.("begin load condition",fLines))
    num_lc_cards = length(lcStartLines)
    lcEndLines = zeros(Int,num_lc_cards)
    for lc = 1:num_lc_cards
        lc_endline_index = findall(endLines .> lcStartLines[lc])[1]
        lcEndLines[lc] = endLines[lc_endline_index]
    end

    LC = [LoadCondition() for lc = 1:num_lc_cards]
    for lc = 1:num_lc_cards
        lcCard = fLines[(lcStartLines[lc]+1) : (lcEndLines[lc]-1)]
        LC[lc] = parseLoadCondition(lcCard)
    end
   
    return GC, BC, LC
end

function parseGeometryCard(Card)
    geomFile = strip(split(Card)[end])

    GC = GeometryCard()
    GC.Filename = geomFile
    return GC
end

function parseBoundaryCondition(Card)
    # Get BC Type
    bc_type_index = findall(occursin.("type",Card))
    bc_type_line = Card[bc_type_index][1]
    bc_type = strip(split(bc_type_line,"=")[end])
    bc_type = getproperty(feEnumerate,Symbol(bc_type))
    
    # Get BC Nodeset
    bc_ns_index = findall(occursin.("node set",Card))
    bc_ns_line = Card[bc_ns_index][1]
    bc_ns = strip(split(bc_ns_line,"=")[end])
    
    # Get BC DOF ID
    bc_dofid_index = findall(occursin.("dof id",Card))
    bc_dofid_line = Card[bc_dofid_index][1]
    bc_dofid = strip(split(bc_dofid_line,"=")[end])
    bc_dofid = parse(Int,bc_dofid)
    
    # Get BC Value
    bc_value_index = findall(occursin.("value",Card))
    bc_value_line = Card[bc_value_index][1]
    bc_value = strip(split(bc_value_line,"=")[end])
    bc_value = parse(Float64,bc_value)
    
    # Create a BoundaryCondition instance and set the values
    BC = BoundaryCondition()
    BC.Type = bc_type
    BC.NodeSetName = bc_ns
    BC.DOF = bc_dofid
    BC.Value = bc_value

    return BC
end

function parseLoadCondition(Card)
    # Get LC Type
    lc_type_index = findall(occursin.("type",Card))
    lc_type_line = Card[lc_type_index][1]
    lc_type = strip(split(lc_type_line,"=")[end])
    lc_type = getproperty(feEnumerate,Symbol(lc_type))
    # Get LC element set or surface set
    if lc_type == feEnumerate.body
        lc_es_index = findall(occursin.("element set",Card))
        lc_es_line = Card[lc_es_index][1]
        lc_es = strip(split(lc_es_line,"=")[end])
    else
        lc_ss_index = findall(occursin.("surface set",Card))
        lc_ss_line = Card[lc_ss_index][1]
        lc_ss = strip(split(lc_ss_line,"=")[end])
    end

    # Get LC Direction (if applicable)
    lc_direction_index = findall(occursin.("direction",Card))
    if isempty(lc_direction_index) == false
        lc_direction_line = Card[lc_direction_index][1]
        lc_direction = strip(split(lc_direction_line,"=")[end])
        lc_direction = reduce(replace, ["["=>"", "]"=>""], init=lc_direction)
        lc_direction = [parse(Float64,val) for val in split(lc_direction,",")]
    else
        lc_direction = []
    end
    
    # Get LC Magnitude (if applicable)
    lc_magnitude_index = findall(occursin.("magnitude",Card))
    if isempty(lc_magnitude_index) == false
        lc_magnitude_line = Card[lc_magnitude_index][1]
        lc_magnitude = strip(split(lc_magnitude_line,"=")[end])
        lc_magnitude = parse(Float64,lc_magnitude)
    else
        lc_magnitude = []
    end
    
    # Create a LoadCondition instance and set the values
    LC = LoadCondition()
    LC.Type = lc_type
    if Int(lc_type) == Int(feEnumerate.body)
        LC.ElementSetName = lc_es
    else
        LC.SurfaceSetName = lc_ss
    end
    LC.Direction = lc_direction
    LC.Magnitude = lc_magnitude

    return LC
end
