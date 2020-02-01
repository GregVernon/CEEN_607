module feInputParser

include("feMeshImport.jl")
include("feEnumerations.jl")

export importInput
export findInputCards
export parseBoundaryCondition
export parseLoadCondition
export InputParams
export BoundaryCondition
export LoadCondition
export BodyCondition

function importInput(filename)
    fLines = readlines(filename)
    fLines = strip.(fLines)
    fLines = lowercase.(fLines)

    GEOM, BC, LC = findInputCards(fLines)
    
    PARAMS = InputParams()
    PARAMS.BoundaryConditions = BC
    PARAMS.LoadConditions = LC

    return GEOM, PARAMS
end

function findInputCards(fLines)
    # First get the geometry
    geomLine = findall(occursin.("geometry",fLines))[1]
    geomFile = strip(split(fLines[geomLine])[end])
    G = feMeshImport.importGenomat(geomFile)

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

    LC = [LoadCondition() for bc = 1:num_lc_cards]
    for lc = 1:num_lc_cards
        lcCard = fLines[(lcStartLines[lc]+1) : (lcEndLines[lc]-1)]
        LC[lc] = parseLoadCondition(lcCard)
    end
   
    return G, BC, LC
end

function parseBoundaryCondition(Card)
    # Get BC Type
    bc_type_index = findall(occursin.("type",Card))
    bc_type_line = Card[bc_type_index][1]
    bc_type = strip(split(bc_type_line,"=")[end])
    bc_type = getproperty(feEnumerations,Symbol(bc_type))
    
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
    lc_type = getproperty(feEnumerations,Symbol(lc_type))
    
    # Get LC Nodeset
    lc_ns_index = findall(occursin.("surface set",Card))
    lc_ns_line = Card[lc_ns_index][1]
    lc_ns = strip(split(lc_ns_line,"=")[end])
      
    # Get LC Value
    lc_value_index = findall(occursin.("value",Card))
    lc_value_line = Card[lc_value_index][1]
    lc_value = strip(split(lc_value_line,"=")[end])
    lc_value = parse(Float64,lc_value)
    
    # Create a LoadCondition instance and set the values
    LC = LoadCondition()
    LC.Type = lc_type
    LC.SurfaceSetName = lc_ns
    LC.Value = lc_value

    return LC
end

mutable struct InputParams
    BoundaryConditions
    LoadConditions
    InputParams() = new()
end

mutable struct BoundaryCondition
    Type
    NodeSetName
    DOF
    Value
    BoundaryCondition() = new()
end

mutable struct LoadCondition
    Type
    SurfaceSetName
    Value
    LoadCondition() = new()
end

mutable struct BodyCondition
    Type
    ElementSetName
    Value
    BodyCondition() = new()
end

end