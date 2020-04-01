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
    for e = 1:nBlk
        eID = eID + 1;
        %%%% Element Properties
        MESH.Elements(eID).Degree = [1 1];
        MESH.Elements(eID).Dimension = 2;
        MESH.Elements(eID).GlobalID = eID;
        MESH.Elements(eID).Type = EXO.element_blocks{3,blk};
        MESH.Elements(eID).NodeConnectivity = EXO.element_blocks{4,blk}(:,e);
        MESH.Elements(eID).Parametric = feElementConfiguration("Parametric");
        MESH.Elements(eID).Reference = feElementConfiguration("Reference");
               
        %%%% Initialize Parametric Configuration
        MESH.Elements(eID).Parametric.Degree = MESH.Elements(eID).Degree;
        MESH.Elements(eID).Parametric.Dimension = MESH.Elements(eID).Dimension;
        MESH.Elements(eID).Parametric.GlobalID = MESH.Elements(eID).GlobalID;
        MESH.Elements(eID).Parametric.Type = MESH.Elements(eID).Type;
        
        %%% Parametric Nodes
        num_loc_nodes = length(MESH.Elements(eID).NodeConnectivity);
        MESH.Elements(eID).Parametric.Nodes = repmat(feNode(),num_loc_nodes,1);
        for n = 1:num_loc_nodes
            MESH.Elements(eID).Parametric.Nodes(n).ConfigurationType = "Parametric";
            MESH.Elements(eID).Parametric.Nodes(n).Coordinates = feElementParametric(MESH.Elements(eID).Degree, n);
            MESH.Elements(eID).Parametric.Nodes(n).Dimension = MESH.Elements(eID).Dimension;
            MESH.Elements(eID).Parametric.Nodes(n).GlobalID = MESH.Elements(eID).NodeConnectivity(n);
            MESH.Elements(eID).Parametric.Nodes(n).ChildDOF_ID = (MESH.Elements.Parametric.Nodes(n).GlobalID - 1) * MESH.Elements.Parametric.Nodes(n).Dimension + n;
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
               
        %%% Reference Nodes
        MESH.Elements.Reference.Nodes = repmat(feNode(),num_loc_nodes,1);
        for n = 1:num_loc_nodes
            MESH.Elements(eID).Reference.Nodes(n).ConfigurationType = "Reference";
            MESH.Elements(eID).Reference.Nodes(n).Coordinates = getExodusNodeCoordinates(EXO, MESH.Elements.Reference.Nodes(n).GlobalID);
            MESH.Elements(eID).Reference.Nodes(n).Dimension = MESH.Elements(eID).Dimension;
            MESH.Elements(eID).Reference.Nodes(n).GlobalID = MESH.Elements(eID).NodeConnectivity(n);
            MESH.Elements(eID).Reference.Nodes(n).ChildDOF_ID = (MESH.Elements.Parametric.Nodes(n).GlobalID - 1) * MESH.Elements.Parametric.Nodes(n).Dimension + n;
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
end


function NodeCoords = getExodusNodeCoordinates(EXO, node_id)
nDim = EXO.naxes;
if     nDim == 1
    NodeCoords =  x0(node_id);
elseif nDim == 2
    NodeCoords = [x0(node_id); y0(node_id)];
elseif nDim == 3
    NodeCoords = [x0(node_id); y0(node_id); z0(node_id)];
end

end