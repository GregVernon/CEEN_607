function [MESH, PARAMS] = importInput(filename)

    fLines = fileread(filename);
    fLines = strsplit(fLines,["\n","\r"])';
    fLines = strip(fLines);
    fLines = lower(fLines);

    [GC, BC, LC] = findInputCards(fLines);
    MESH = importMesh(GC.Filename);
    MESH = MESH.assignLoadConditions(LC);
    MESH = MESH.assignBoundaryConditions(BC);

    PARAMS = struct("BoundaryConditions",BC,"LoadConditions",LC);   
end


function [GC, BC, LC] = findInputCards(fLines)
    % First get the geometry
    geomLine = contains(fLines,"geometry");
    geomCard = fLines(geomLine);
    GC = parseGeometryCard(geomCard);

    % Next get parameters
    beginLines = find(contains(fLines, "begin"));
    endLines = find(contains(fLines, "end"));

    % Find Boundary Condition Cards
    bcStartLines = find(contains(fLines,"begin boundary condition"));
    num_bc_cards = length(bcStartLines);
    bcEndLines = zeros(num_bc_cards,1);
    for bc = 1:num_bc_cards
        bc_endline_index = find(endLines > bcStartLines(bc),1,'first');
        bcEndLines(bc) = endLines(bc_endline_index);
    end
    
    BC = repmat(feBoundaryCondition(), num_bc_cards, 1);
    for bc = 1:num_bc_cards
        bcCard = fLines((bcStartLines(bc)+1) : (bcEndLines(bc)-1));
        BC(bc) = parseBoundaryCondition(bcCard);
    end

    % Find Load Condition Cards
    lcStartLines = find(contains(fLines,"begin load condition"));
    num_lc_cards = length(lcStartLines);
    lcEndLines = zeros(num_lc_cards,1);
    for lc = 1:num_lc_cards
        lc_endline_index = find(endLines > lcStartLines(lc),1,'first');
        lcEndLines(lc) = endLines(lc_endline_index);
    end

    LC = repmat(feLoadCondition(), num_lc_cards, 1);
    for lc = 1:num_lc_cards
        lcCard = fLines((lcStartLines(lc)+1) : (lcEndLines(lc)-1));
        LC(lc) = parseLoadCondition(lcCard);
    end
end

function GC = parseGeometryCard(Card)
    geomFile = strip(split(Card));
    geomFile = geomFile{end};
    
    GC.Filename = geomFile;
end

function BC = parseBoundaryCondition(Card)
    % Get BC Type
    bc_type_index = contains(Card,"type");
    bc_type_line = Card(bc_type_index);
    bc_type = strip(split(bc_type_line,"="));
    bc_type = bc_type(end);
    
    % Get BC Nodeset
    bc_ns_index = contains(Card,"node set");
    bc_ns_line = Card(bc_ns_index);
    bc_ns = strip(split(bc_ns_line,"="));
    bc_ns = bc_ns(end);
    
    % Get BC DOF ID
    bc_dofid_index = contains(Card,"dof id");
    bc_dofid_line = Card(bc_dofid_index);
    bc_dofid = strip(split(bc_dofid_line,"="));
    bc_dofid = str2double(bc_dofid(end));
    
    % Get BC Value
    bc_value_index = contains(Card,"value");
    bc_value_line = Card(bc_value_index);
    bc_value = strip(split(bc_value_line,"="));
    bc_value = str2double(bc_value(end));
    
    % Create a BoundaryCondition instance and set the values
    BC = feBoundaryCondition();
    BC.Type = bc_type;
    BC.NodeSetName = bc_ns;
    BC.DOF = bc_dofid;
    BC.Value = bc_value;
end

function LC = parseLoadCondition(Card)
    % Get LC Type
    lc_type_index = contains(Card,"type");
    lc_type_line = Card(lc_type_index);
    lc_type = strip(split(lc_type_line,"="));
    lc_type = lc_type(end);
    
    % Get LC element set or surface set
    if strcmpi(lc_type, "body")
        lc_es_index = contains(Card,"element set");
        lc_es_line = Card(lc_es_index);
        lc_es = strip(split(lc_es_line,"="));
        lc_es = lc_es(end);
    else
        lc_ss_index = contains(Card,"surface set");
        lc_ss_line = Card(lc_ss_index);
        lc_ss = strip(split(lc_ss_line,"="));
        lc_ss = lc_ss(end);
    end

    % Get LC Direction (if applicable)
    lc_direction_index = find(contains(Card,"direction"));
    if isempty(lc_direction_index) == false
        lc_direction_line = Card(lc_direction_index);
        lc_direction = strip(split(lc_direction_line,"="));
        lc_direction = str2num(lc_direction{end});  %#ok<ST2NM>
    else
        lc_direction = [];
    end
    
    % Get LC Magnitude (if applicable)
    lc_magnitude_index = find(contains(Card,"magnitude"));
    if isempty(lc_magnitude_index) == false
        lc_magnitude_line = Card(lc_magnitude_index);
        lc_magnitude = strip(split(lc_magnitude_line,"="));
        lc_magnitude = str2double(lc_magnitude(end));
    else
        lc_magnitude = [];
    end
    
    % Create a LoadCondition instance and set the values
    LC = feLoadCondition();
    LC.Type = lc_type;
    if strcmpi(lc_type, "body")
        LC.ElementSetName = lc_es;
    else
        LC.SurfaceSetName = lc_ss;
    end
    LC.Direction = lc_direction;
    LC.Magnitude = lc_magnitude;
end