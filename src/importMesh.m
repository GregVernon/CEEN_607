function MESH = importMesh(filename)
%%%%% Import the mesh file
EXO = load(filename);

%%%%% Initialize Mesh class
MESH = feMesh(EXO.nelems);

%%%%% Initialize Elements
eID = 0;
nBlk = size(EXO.element_blocks,2);
for blk = 1:nBlk
    nBlkElems = size(EXO.element_blocks{4,blk},2);
    for e = 1:nBlkElems
        eID = eID + 1;
        %%%% Element Properties
        MESH.Elements(eID).Degree = [1 1];
        MESH.Elements(eID).Dimension = 2;
        MESH.Elements(eID).GlobalID = eID;
        MESH.Elements(eID).Type = EXO.element_blocks{3,blk};
        if     strcmpi(MESH.Elements(eID).Type, "BAR") == true
            MESH.Elements(eID).NodeConnectivity = EXO.element_blocks{4,blk}(:,e);
        elseif strcmpi(MESH.Elements(eID).Type, "QUAD4") == true
            MESH.Elements(eID).NodeConnectivity = double(EXO.element_blocks{4,blk}([1 2 4 3],e));
            MESH.Elements(eID).DOFConnectivity = 2*(MESH.Elements(eID).NodeConnectivity' - 1) + [1; 2];
        elseif strcmpi(MESH.Elements(eID).Type, "HEX8") == true
            MESH.Elements(eID).NodeConnectivity = EXO.element_blocks{4,blk}([1 2 4 3 5 6 8 7],e);
        else
            error("Element Type not Supported for Import");
        end
        MESH.NodeConnectivity = MESH.gather_node_connectivity();
        MESH.DOFConnectivity = MESH.gather_dof_connectivity();
        
        MESH.Elements(eID).Parametric = feElementParametric();
        MESH.Elements(eID).Reference = feElementReference();
               
        %%%% Initialize Parametric Configuration
        MESH.Elements(eID).Parametric.Degree = MESH.Elements(eID).Degree;
        MESH.Elements(eID).Parametric.Dimension = MESH.Elements(eID).Dimension;
        MESH.Elements(eID).Parametric.GlobalID = MESH.Elements(eID).GlobalID;
        MESH.Elements(eID).Parametric.Type = MESH.Elements(eID).Type;
        
        %%% Parametric Nodes
        num_loc_nodes = length(MESH.Elements(eID).NodeConnectivity);
        num_loc_dof = size(MESH.Elements(eID).DOFConnectivity,1);
        MESH.Elements(eID).Parametric.Nodes = repmat(feNode(),num_loc_nodes,1);
        for n = 1:num_loc_nodes
            MESH.Elements(eID).Parametric.Nodes(n).ConfigurationType = "Parametric";
            MESH.Elements(eID).Parametric.Nodes(n).Coordinates = feElementParametric.computeNodeCoordinates(MESH.Elements(eID).Degree, n);
            MESH.Elements(eID).Parametric.Nodes(n).Dimension = MESH.Elements(eID).Dimension;
            MESH.Elements(eID).Parametric.Nodes(n).GlobalID = MESH.Elements(eID).NodeConnectivity(n);
            MESH.Elements(eID).Parametric.Nodes(n).ChildDOF_ID = (MESH.Elements(eID).Parametric.Nodes(n).GlobalID - 1) * MESH.Elements(eID).Parametric.Nodes(n).Dimension + n;
        end
        
        %%% Parametric Quadrature
        if     MESH.Elements(eID).Dimension == 1
            Quadrature = feQuadrature("Solid");
        elseif MESH.Elements(eID).Dimension == 2
            Quadrature = [feQuadrature("Solid","QP1"); repmat(feQuadrature("Boundary","QP1"),4,1)];
        elseif MESH.Elements(eID).Dimension == 3
            Quadrature = [feQuadrature("Solid","QP1"); repmat(feQuadrature("Boundary","QP1"),6,1)];
        end
        MESH.Elements(eID).Parametric.Quadrature = Quadrature;
        
        %%%% Initialize Reference Configuration
        MESH.Elements(eID).Reference.Degree = MESH.Elements(eID).Degree;
        MESH.Elements(eID).Reference.Dimension = MESH.Elements(eID).Dimension;
        MESH.Elements(eID).Reference.GlobalID = MESH.Elements(eID).GlobalID;
        MESH.Elements(eID).Reference.Type = MESH.Elements(eID).Type;
        
        MESH.Elements(eID).Reference.DirichletConditions = repmat({nan(num_loc_dof,1)},num_loc_nodes,1);
        MESH.Elements(eID).Reference.BodyForce = zeros(num_loc_dof,1);
        MESH.Elements(eID).Reference.NodeForce = repmat({zeros(num_loc_dof,1)},num_loc_nodes,1);
        MESH.Elements(eID).Reference.SurfacePressure = repmat({0},num_loc_nodes,1);
        
        %%% Reference Nodes
        MESH.Elements(eID).Reference.Nodes = repmat(feNode(),num_loc_nodes,1);
        for n = 1:num_loc_nodes
            MESH.Elements(eID).Reference.Nodes(n).ConfigurationType = "Reference";
            MESH.Elements(eID).Reference.Nodes(n).Coordinates = getExodusNodeCoordinates(EXO, MESH.Elements(eID).Parametric.Nodes(n).GlobalID);
            MESH.Elements(eID).Reference.Nodes(n).Dimension = MESH.Elements(eID).Dimension;
            MESH.Elements(eID).Reference.Nodes(n).GlobalID = MESH.Elements(eID).NodeConnectivity(n);
            MESH.Elements(eID).Reference.Nodes(n).ChildDOF_ID = (MESH.Elements(eID).Reference.Nodes(n).GlobalID - 1) * MESH.Elements(eID).Reference.Nodes(n).Dimension + n;
        end
        
        %%% Parametric Quadrature
        if     MESH.Elements(eID).Dimension == 1
            Quadrature = feQuadrature("Solid");
        elseif MESH.Elements(eID).Dimension == 2
            Quadrature = [feQuadrature("Solid","QP1"); repmat(feQuadrature("Boundary","QP1"),4,1)];
        elseif MESH.Elements(eID).Dimension == 3
            Quadrature = [feQuadrature("Solid","QP1"); repmat(feQuadrature("Boundary","QP1"),6,1)];
        end
        MESH.Elements(eID).Reference.Quadrature = Quadrature;
    end
end

%%%%% Initialize Element Sets
num_elem_sets = EXO.nblks;
MESH.ElementSets = repmat(feElementSet(), num_elem_sets, 1);
eCount = 0;
for es = 1:num_elem_sets
    num_blk_elems = size(EXO.element_blocks{4,es},2);
    MESH.ElementSets(es).Name = EXO.element_blocks{1,es};
    MESH.ElementSets(es).ElementID = (eCount+1) : (eCount+num_blk_elems);
end

%%%% Initialize Surface Sets
num_surf_sets = EXO.nssets;
MESH.SurfaceSets = repmat(feSurfaceSet(), num_surf_sets, 1);
for ss = 1:num_surf_sets
    MESH.SurfaceSets(ss).Name = EXO.side_sets{1,ss};
    MESH.SurfaceSets(ss).ElementID = EXO.side_sets{3,ss};
    MESH.SurfaceSets(ss).LocalSideID = feSurfaceSet.ExodusSideID_to_TensorSideID(EXO.element_blocks{3,1}, EXO.side_sets{4,ss});
    MESH.SurfaceSets(ss).GlobalNodeID = EXO.side_sets{6,ss};
end

%%%% Initialize Node Sets
num_node_sets = EXO.nnsets;
MESH.NodeSets = repmat(feNodeSet(), num_node_sets, 1);
for ns = 1:num_node_sets
    MESH.NodeSets(ns).Name = EXO.node_sets{1,ns};
    MESH.NodeSets(ns).GlobalNodeID = EXO.node_sets{3,ns};
    % Determine which elements contain each node in nodeset
    num_node_in_set = length(MESH.NodeSets(ns).GlobalNodeID);
    MESH.NodeSets(ns).ElementID = cell(num_node_in_set,1);
    for n = 1:num_node_in_set
        global_node_id =  MESH.NodeSets(ns).GlobalNodeID(n);
        MESH.NodeSets(ns).ElementID{n} = find(any(MESH.NodeConnectivity == global_node_id,1));
    end
end

%%%% After initialization, precompute values and cache them
MESH = MESH.precomputeMesh();
end


function NodeCoords = getExodusNodeCoordinates(EXO, node_id)
nDim = EXO.naxes;
if     nDim == 1
    NodeCoords =  EXO.x0(node_id);
elseif nDim == 2
    NodeCoords = [EXO.x0(node_id); EXO.y0(node_id)];
elseif nDim == 3
    NodeCoords = [EXO.x0(node_id); EXO.y0(node_id); EXO.z0(node_id)];
end
end